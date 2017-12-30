/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <cmath>  // sqrt
#include <string>
#include <memory>
#include "viscousFlux.hpp"
#include "physicsModels.hpp"
#include "primitive.hpp"       // primitive
#include "matrix.hpp"         // squareMatrix
#include "utility.hpp"        // TauNormal
#include "wallData.hpp"       // wallVars

using std::cout;
using std::endl;
using std::cerr;
using std::max;
using std::unique_ptr;

// constructor -- initialize flux from velocity gradient
/*
Viscous flux normal to face:
F = [ 0,
      taux,
      tauy,
      tauz,
      tau (dot) vel + K * tGrad (dot) area
      (mu + sk *mut) * tkeGrad (dot) area
      (mu + sw *mut) * omegaGrad (dot) area ]

In the above equation tau is the wall shear stress. Taux, tauy, and tauz are the
rows of the wall shear stress tensor i.e. taux = tauxx + tauxy + tauxz. K is the
thermal conductivity, tGrad is the temperature gradient, and area is the
normalized face area. Sk and sw are turbulence model coefficients

Wall shear stress:
tau = lambda * velGradTrace * area + mu * ( velGrad * area + velGrad' * area)

In the above equation lambda is the bulk viscosity, velGradTrace is the trace of
the velocity gradient, area is the normalized face area, mu is the dynamic
viscosity, and velGrad is the velocity gradient tensor.
*/
void viscousFlux::CalcFlux(const tensor<double> &velGrad, const physics &phys,
                           const vector3d<double> &tGrad,
                           const vector3d<double> &normArea,
                           const vector3d<double> &tkeGrad,
                           const vector3d<double> &omegaGrad,
                           const primitive &state, const double &lamVisc,
                           const double &turbVisc, const double &f1) {
  // velGrad -- velocity gradient tensor
  // phys -- physics models
  // tGrad -- temperature gradient
  // normArea -- unit area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // state -- primitive variables at face
  // lamVisc -- laminar viscosity
  // turbVisc -- turbulent viscosity
  // f1 -- first blending coefficient

  // get viscosity with nondimensional normalization
  const auto mu = phys.Transport()->NondimScaling() * lamVisc;
  const auto mut = phys.Transport()->NondimScaling() * turbVisc;

  // wall shear stress
  const auto tau = TauNormal(velGrad, normArea, mu, mut, phys.Transport());

  const auto t = state.Temperature(phys.EoS());

  (*this)[this->MomentumXIndex()] = tau.X();
  (*this)[this->MomentumYIndex()] = tau.Y();
  (*this)[this->MomentumZIndex()] = tau.Z();
  (*this)[this->EnergyIndex()] =
      tau.DotProd(state.Velocity()) +
      (phys.Transport()->EffectiveConductivity(t, state.MassFractions()) +
       phys.Transport()->TurbConductivity(
           mut, phys.Turbulence()->TurbPrandtlNumber(), t, phys.Thermodynamic(),
           state.MassFractions())) *
          tGrad.DotProd(normArea);

  // turbulence viscous flux
  if (this->HasTurbulenceData()) {
    // get molecular diffusion coefficients for turbulence equations
    const auto tkeCoeff = phys.Turbulence()->SigmaK(f1);
    const auto omgCoeff = phys.Turbulence()->SigmaW(f1);

    // some turbulence models use the unlimited eddy viscosity for the
    // turbulence viscous flux instead of the limited eddy viscosity
    const auto mutt = phys.Turbulence()->UseUnlimitedEddyVisc()
                          ? phys.Transport()->NondimScaling() *
                                phys.Turbulence()->EddyViscNoLim(state)
                          : mut;
    (*this)[this->TurbulenceIndex()] =
        (mu + tkeCoeff * mutt) * tkeGrad.DotProd(normArea);
    (*this)[this->TurbulenceIndex() + 1] =
        (mu + omgCoeff * mutt) * omegaGrad.DotProd(normArea);
  }
}

wallVars viscousFlux::CalcWallFlux(
    const tensor<double> &velGrad, const physics &phys,
    const vector3d<double> &tGrad, const vector3d<double> &normArea,
    const vector3d<double> &tkeGrad, const vector3d<double> &omegaGrad,
    const primitive &state, const double &lamVisc, const double &turbVisc,
    const double &f1) {
  // velGrad -- velocity gradient tensor
  // phys -- physics models
  // tGrad -- temperature gradient
  // normArea -- unit area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // state -- primitive variables at face
  // lamVisc -- laminar viscosity
  // turbVisc -- turbulent viscosity
  // f1 -- first blending coefficient

  wallVars wVars(state.NumSpecies());

  // get viscosity with nondimensional normalization
  wVars.viscosity_ = phys.Transport()->NondimScaling() * lamVisc;
  wVars.turbEddyVisc_ = phys.Transport()->NondimScaling() * turbVisc;

  // wall shear stress
  wVars.shearStress_ = TauNormal(velGrad, normArea, wVars.viscosity_,
                                 wVars.turbEddyVisc_, phys.Transport());

  const auto t = state.Temperature(phys.EoS());

  // wall heat flux
  wVars.heatFlux_ =
      (phys.Transport()->EffectiveConductivity(t, state.MassFractions()) +
       phys.Transport()->TurbConductivity(
           wVars.turbEddyVisc_, phys.Turbulence()->TurbPrandtlNumber(), t,
           phys.Thermodynamic(), state.MassFractions())) *
      tGrad.DotProd(normArea);

  (*this)[this->MomentumXIndex()] = wVars.shearStress_.X();
  (*this)[this->MomentumYIndex()] = wVars.shearStress_.Y();
  (*this)[this->MomentumZIndex()] = wVars.shearStress_.Z();
  (*this)[this->EnergyIndex()] =
      wVars.shearStress_.DotProd(state.Velocity()) + wVars.heatFlux_;

  // calculate other wall data
  wVars.density_ = state.Rho();
  wVars.temperature_ = t;
  wVars.frictionVelocity_ = sqrt(wVars.shearStress_.Mag() / wVars.density_);

  // turbulence viscous flux
  if (this->HasTurbulenceData()) {
    wVars.tke_ = state.Tke();
    wVars.sdr_ = state.Omega();

    // get molecular diffusion coefficients for turbulence equations
    const auto tkeCoeff = phys.Turbulence()->SigmaK(f1);
    const auto omgCoeff = phys.Turbulence()->SigmaW(f1);

    // some turbulence models use the unlimited eddy viscosity for the
    // turbulence viscous flux instead of the limited eddy viscosity
    const auto mutt = phys.Turbulence()->UseUnlimitedEddyVisc()
                          ? phys.Transport()->NondimScaling() *
                                phys.Turbulence()->EddyViscNoLim(state)
                          : wVars.turbEddyVisc_;
    (*this)[this->TurbulenceIndex()] =
        (wVars.viscosity_ + tkeCoeff * mutt) * tkeGrad.DotProd(normArea);
    (*this)[this->TurbulenceIndex() + 1] =
        (wVars.viscosity_ + omgCoeff * mutt) * omegaGrad.DotProd(normArea);
  }

  return wVars;
}

void viscousFlux::CalcWallLawFlux(
    const vector3d<double> &tauWall, const double &qWall, const double &muWall,
    const double &mutWall, const vector3d<double> &velWall,
    const vector3d<double> &normArea, const vector3d<double> &tkeGrad,
    const vector3d<double> &omegaGrad, const unique_ptr<turbModel> &turb) {
  // tauWall -- wall shear stress
  // qWall -- wall heat flux
  // muWall -- wall viscosity
  // mutWall -- wall eddy viscosity
  // normArea -- unit area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // turb -- turbulence model

  (*this)[this->MomentumXIndex()] = tauWall.X();
  (*this)[this->MomentumYIndex()] = tauWall.Y();
  (*this)[this->MomentumZIndex()] = tauWall.Z();
  (*this)[this->EnergyIndex()] = tauWall.DotProd(velWall) + qWall;

  // turbulence viscous flux
  if (this->HasTurbulenceData()) {
    // get molecular diffusion coefficients for turbulence equations
    const auto tkeCoeff = turb->WallSigmaK();
    const auto omgCoeff = turb->WallSigmaW();

    // some turbulence models use the unlimited eddy viscosity for the
    // turbulence viscous flux instead of the limited eddy viscosity
    // for wall laws, eddy viscosity is prescribed
    (*this)[this->TurbulenceIndex()] =
        (muWall + tkeCoeff * mutWall) * tkeGrad.DotProd(normArea);
    (*this)[this->TurbulenceIndex() + 1] =
        (muWall + omgCoeff * mutWall) * omegaGrad.DotProd(normArea);
  }
}

// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, viscousFlux &flux) {
  for (auto rr = 0; rr < flux.Size(); rr++) {
    os << flux[rr] << endl;
  }
  return os;
}

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
#include "eos.hpp"            // equation of state
#include "transport.hpp"      // transport model
#include "thermodynamic.hpp"  // thermodynamic model
#include "primVars.hpp"       // primVars
#include "turbulence.hpp"     // turbModel
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
void viscousFlux::CalcFlux(
    const tensor<double> &velGrad, const unique_ptr<transport> &trans,
    const unique_ptr<thermodynamic> &thermo, const vector3d<double> &tGrad,
    const vector3d<double> &normArea, const vector3d<double> &tkeGrad,
    const vector3d<double> &omegaGrad, const unique_ptr<turbModel> &turb,
    const primVars &state, const double &lamVisc, const double &turbVisc,
    const double &f1) {
  // velGrad -- velocity gradient tensor
  // trans -- viscous transport model
  // thermo -- thermodynamic model
  // tGrad -- temperature gradient
  // normArea -- unit area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // turb -- turbulence model
  // state -- primative variables at face
  // lamVisc -- laminar viscosity
  // turbVisc -- turbulent viscosity
  // f1 -- first blending coefficient

  // get viscosity with nondimensional normalization
  const auto mu = trans->NondimScaling() * lamVisc;
  const auto mut = trans->NondimScaling() * turbVisc;

  // wall shear stress
  const auto tau = TauNormal(velGrad, normArea, mu, mut, trans);

  data_[0] = tau.X();
  data_[1] = tau.Y();
  data_[2] = tau.Z();
  data_[3] = tau.DotProd(state.Velocity()) +
      (trans->Conductivity(mu, thermo) +
       trans->TurbConductivity(mut, turb->TurbPrandtlNumber(), thermo)) *
      tGrad.DotProd(normArea);

  // turbulence viscous flux
  // get molecular diffusion coefficients for turbulence equations
  const auto tkeCoeff = turb->SigmaK(f1);
  const auto omgCoeff = turb->SigmaW(f1);

  // some turbulence models use the unlimited eddy viscosity for the
  // turbulence viscous flux instead of the limited eddy viscosity
  const auto mutt = turb->UseUnlimitedEddyVisc() ?
      trans->NondimScaling() * turb->EddyViscNoLim(state) : mut;
  data_[4] = (mu + tkeCoeff * mutt) * tkeGrad.DotProd(normArea);
  data_[5] = (mu + omgCoeff * mutt) * omegaGrad.DotProd(normArea);
}

wallVars viscousFlux::CalcWallFlux(
    const tensor<double> &velGrad, const unique_ptr<transport> &trans,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<eos> &eqnState,
    const vector3d<double> &tGrad, const vector3d<double> &normArea,
    const vector3d<double> &tkeGrad, const vector3d<double> &omegaGrad,
    const unique_ptr<turbModel> &turb, const primVars &state,
    const double &lamVisc, const double &turbVisc, const double &f1) {
  // velGrad -- velocity gradient tensor
  // trans -- viscous transport model
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // tGrad -- temperature gradient
  // normArea -- unit area vector of face
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // turb -- turbulence model
  // state -- primative variables at face
  // lamVisc -- laminar viscosity
  // turbVisc -- turbulent viscosity
  // f1 -- first blending coefficient

  wallVars wVars;

  // get viscosity with nondimensional normalization
  wVars.viscosity_ = trans->NondimScaling() * lamVisc;
  wVars.turbEddyVisc_ = trans->NondimScaling() * turbVisc;

  // wall shear stress
  wVars.shearStress_ =
      TauNormal(velGrad, normArea, wVars.viscosity_, wVars.turbEddyVisc_, trans);

  // wall heat flux
  wVars.heatFlux_ =
      (trans->Conductivity(wVars.viscosity_, thermo) +
       trans->TurbConductivity(wVars.turbEddyVisc_, turb->TurbPrandtlNumber(),
                               thermo)) *
      tGrad.DotProd(normArea);

  data_[0] = wVars.shearStress_.X();
  data_[1] = wVars.shearStress_.Y();
  data_[2] = wVars.shearStress_.Z();
  data_[3] = wVars.shearStress_.DotProd(state.Velocity()) + wVars.heatFlux_;

  // calculate other wall data
  wVars.density_ = state.Rho();
  wVars.temperature_ = state.Temperature(eqnState);
  wVars.tke_ = state.Tke();
  wVars.sdr_ = state.Omega();
  wVars.frictionVelocity_ = sqrt(wVars.shearStress_.Mag() / wVars.density_);

  // turbulence viscous flux
  // get molecular diffusion coefficients for turbulence equations
  const auto tkeCoeff = turb->SigmaK(f1);
  const auto omgCoeff = turb->SigmaW(f1);

  // some turbulence models use the unlimited eddy viscosity for the
  // turbulence viscous flux instead of the limited eddy viscosity
  const auto mutt = turb->UseUnlimitedEddyVisc() ?
      trans->NondimScaling() * turb->EddyViscNoLim(state) : wVars.turbEddyVisc_;
  data_[4] = (wVars.viscosity_ + tkeCoeff * mutt) * tkeGrad.DotProd(normArea);
  data_[5] = (wVars.viscosity_ + omgCoeff * mutt) * omegaGrad.DotProd(normArea);

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

  data_[0] = tauWall.X();
  data_[1] = tauWall.Y();
  data_[2] = tauWall.Z();
  data_[3] = tauWall.DotProd(velWall) + qWall;

  // turbulence viscous flux
  // get molecular diffusion coefficients for turbulence equations
  const auto tkeCoeff = turb->WallSigmaK();
  const auto omgCoeff = turb->WallSigmaW();

  // some turbulence models use the unlimited eddy viscosity for the
  // turbulence viscous flux instead of the limited eddy viscosity
  // for wall laws, eddy viscosity is prescribed
  data_[4] = (muWall + tkeCoeff * mutWall) * tkeGrad.DotProd(normArea);
  data_[5] = (muWall + omgCoeff * mutWall) * omegaGrad.DotProd(normArea);
}


// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, viscousFlux &flux) {
  os << "0.0" << endl;
  os << flux.MomX() << endl;
  os << flux.MomY() << endl;
  os << flux.MomZ() << endl;
  os << flux.Engy() << endl;
  os << flux.MomK() << endl;
  os << flux.MomO() << endl;
  return os;
}

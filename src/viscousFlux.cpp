/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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
#include "eos.hpp"         // idealGas
#include "primVars.hpp"    // primVars
#include "turbulence.hpp"  // turbModel
#include "matrix.hpp"      // squareMatrix
#include "utility.hpp"     // TauNormal

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
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
viscousFlux::viscousFlux(
    const tensor<double> &velGrad, const sutherland &suth,
    const idealGas &eqnState, const vector3d<double> &tGrad,
    const vector3d<double> &normArea, const vector3d<double> &tkeGrad,
    const vector3d<double> &omegaGrad, const unique_ptr<turbModel> &turb,
    const primVars &state, const double &lamVisc, const double &turbVisc,
    const double &f1) {
  // velGrad -- velocity gradient tensor
  // suth -- method to get viscosity (Sutherland's law)
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

  // get viscosity with nondimensional normalization
  const auto mu = suth.NondimScaling() * lamVisc;
  const auto mut = suth.NondimScaling() * turbVisc;

  // get molecular diffusion coefficients for turbulence equations
  const auto tkeCoeff = turb->SigmaK(f1);
  const auto omgCoeff = turb->SigmaW(f1);

  // wall shear stress
  const auto tau = TauNormal(velGrad, normArea, mu, mut, suth);

  data_[0] = tau.X();
  data_[1] = tau.Y();
  data_[2] = tau.Z();
  data_[3] = tau.DotProd(state.Velocity()) +
      (eqnState.Conductivity(mu) +
       eqnState.TurbConductivity(mut, turb->TurbPrandtlNumber())) *
      tGrad.DotProd(normArea);

  // turbulence viscous flux
  // some turbulence models use the unlimited eddy viscosity for the
  // turbulence viscous flux instead of the limited eddy viscosity
  const auto mutt = turb->UseUnlimitedEddyVisc() ?
      suth.NondimScaling() * turb->EddyViscNoLim(state) : mut;
  data_[4] = (mu + tkeCoeff * mutt) * tkeGrad.DotProd(normArea);
  data_[5] = (mu + omgCoeff * mutt) * omegaGrad.DotProd(normArea);
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


// function to calculate the velocity gradients at a cell face using the Thin
// Shear Layer approximation
tensor<double> CalcVelGradTSL(const primVars &left, const primVars &right,
                              const vector3d<double> &normArea,
                              const double &dist) {
  // left -- left state (primative)
  // right -- right state (primative)
  // normArea -- unit area vector of face
  // dist -- distance between centroid of left cell and right cell

  // calculate velocity derivatives
  const auto velDeriv = (right.Velocity() - left.Velocity()) / dist;

  // populate velocity gradient tensor
  tensor<double> velGrad(
      velDeriv.X() * normArea.X(),
      velDeriv.Y() * normArea.X(),
      velDeriv.Z() * normArea.X(),
      velDeriv.X() * normArea.Y(),
      velDeriv.Y() * normArea.Y(),
      velDeriv.Z() * normArea.Y(),
      velDeriv.X() * normArea.Z(),
      velDeriv.Y() * normArea.Z(),
      velDeriv.Z() * normArea.Z());

  return velGrad;
}

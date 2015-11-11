/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <cmath>  // sqrt
#include <string>
#include <memory>
#include "inviscidFlux.hpp"
#include "eos.hpp"
#include "primVars.hpp"
#include "genArray.hpp"
#include "matrix.hpp"
#include "turbulence.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::copysign;
using std::unique_ptr;

// constructor -- initialize flux from state vector
// flux is a 3D flux in the normal direction of the given face
/*

F = [rho * vel (dot) area
     rho * vel (dot) area * velx + P * areax
     rho * vel (dot) area * vely + P * areay
     rho * vel (dot) area * velz + P * areaz
     rho * vel (dot) area * H]

rho -- density
vel -- velocity vector (3D)
area -- area vector (3D)
P -- pressure
H -- enthalpy
velx, vely, velz -- velocity components
areax, areay, areaz -- area components

Constructor is put in a private member function because identical code is
used for constructing from primative variables and conservative variables
once the conservative variables have been changed to primative variables.
The C++11 way of delegating constructors is not used because the primVars
class is not fully defined in the inviscidFlux.hpp header. This way both
constructors (primVars version and genArray version) can call this function
to avoid code duplication.
*/
void inviscidFlux::ConstructFromPrim(const primVars &state,
                                     const idealGas &eqnState,
                                     const vector3d<double> &normArea) {
  // state -- primative variables
  // eqnState -- equation of state
  // normArea -- unit area vector of face

  const auto vel = state.Velocity();
  const auto velNorm = vel.DotProd(normArea);

  data_[0] = state.Rho() * velNorm;
  data_[1] = state.Rho() * velNorm * vel.X() + state.P() * normArea.X();
  data_[2] = state.Rho() * velNorm * vel.Y() + state.P() * normArea.Y();
  data_[3] = state.Rho() * velNorm * vel.Z() + state.P() * normArea.Z();
  data_[4] = state.Rho() * velNorm * state.Enthalpy(eqnState);

  data_[5] = state.Rho() * velNorm * state.Tke();
  data_[6] = state.Rho() * velNorm * state.Omega();
}

inviscidFlux::inviscidFlux(const primVars &state, const idealGas &eqnState,
                           const vector3d<double> &normArea) {
  // state -- primative variables
  // eqnState -- equation of state
  // normArea -- unit area vector of face

  this->ConstructFromPrim(state, eqnState, normArea);
}

// constructor -- initialize flux from state vector using conservative variables
// flux is a 3D flux in the normal direction of the given face
inviscidFlux::inviscidFlux(const genArray &cons, const idealGas &eqnState,
                           const unique_ptr<turbModel> &turb,
                           const vector3d<double> &normArea) {
  // cons -- genArray of conserved variables
  // eqnState -- equation of state
  // turb -- turbulence model
  // normArea -- unit area vector of face

  // convert conserved variables to primative variables
  const primVars state(cons, false, eqnState, turb);

  this->ConstructFromPrim(state, eqnState, normArea);
}

/* Function to calculate inviscid flux using Roe's approximate Riemann solver.
The function takes in the primative varibles constructed
from the left and right states, an equation of state, a face area vector, and
outputs the inviscid flux as well as the maximum wave speed.
The function uses Harten's entropy fix to correct wave speeds near 0 and near
sonic.
__________________________________________
|         |       Ul|Ur        |         |
|         |         |          |         |
|   Ui-1  |   Ui    Ua  Ui+1   |   Ui+2  |
|         |         |          |         |
|         |         |          |         |
|_________|_________|__________|_________|

In the diagram above, Ul and Ur represent the reconstructed states, which for
the MUSCL scheme (1st and 2nd order) come from the stencil shown
above. Ua represents the average state at the face at which the flux will be
calculated. In this case it is a Roe average. Once the average
state has been calculated, the Roe flux can be calculated using the following
equation.

F = 1/2 * (Fl + Fr - D)

F represents the calculated Roe flux at the given face. Fl is the inviscid flux
calculated from the reconstructed state Ul. Fr is the inviscid
flux calculated from the reconstructed state Ur. D is the dissipation term which
is calculated using the Roe averaged state, as well as the
eigen values and eigen vectors resulting from the Roe linearization.

D = A * (Ur - Ul) = T * L * T^-1 * (Ur - Ul) = T * L * (Cr - Cl)

A is the linearized Roe matrix. It is equal to the convective flux jacobian
(dF/dU) calculated with the Roe averaged state. The linearization
essentially states that the flux jacobian (which is the change in flux over
change in state) mulitplied by the change in state is equal to the
change in flux. The change in flux is then added to the average of the physical
right and left fluxes (central difference). The Roe matrix, A,
can be diagonalized where T and T^-1 are the right and left eigenvectors
respectively and L is the eigenvalues. T^-1 multiplied by the change in
state results in the change in characteristic wave amplitude (Cr - Cl), or wave
strength. In its final form (right most) T represents the
characteristic waves, L represents the wave speeds, and (Cr - Cl) represents the
wave strength across the face.

*/
inviscidFlux RoeFlux(const primVars &left, const primVars &right,
                     const idealGas &eqnState,
                     const vector3d<double> &areaNorm) {
  // left -- primative variables from left
  // right -- primative variables from right
  // eqnState -- equation of state
  // areaNorm -- norm area vector of face

  // compute Rho averaged quantities
  // density ratio
  const auto denRatio = sqrt(right.Rho() / left.Rho());
  // Roe averaged density
  const auto rhoR = left.Rho() * denRatio;
  // Roe averaged velocities - u, v, w
  const auto uR = (left.U() + denRatio * right.U()) / (1.0 + denRatio);
  const auto vR = (left.V() + denRatio * right.V()) / (1.0 + denRatio);
  const auto wR = (left.W() + denRatio * right.W()) / (1.0 + denRatio);
  // Roe averaged total enthalpy
  const auto hR = (left.Enthalpy(eqnState) + denRatio *
                   right.Enthalpy(eqnState)) / (1.0 + denRatio);
  // Roe averaged speed of sound
  const auto aR = sqrt((eqnState.Gamma() - 1.0) *
                       (hR - 0.5 * (uR * uR + vR * vR + wR * wR)));
  // Roe averaged tke
  const auto kR = (left.Tke() + denRatio * right.Tke()) / (1.0 + denRatio);
  // Roe averaged specific dissipation (omega)
  const auto omR = (left.Omega() + denRatio * right.Omega()) / (1.0 + denRatio);

  // Roe averaged face normal velocity
  const vector3d<double> velR(uR, vR, wR);

  // Roe velocity dotted with normalized area vector
  const auto velRSum = velR.DotProd(areaNorm);

  // normal velocity difference between left and right states
  const auto normVelDiff = right.Velocity().DotProd(areaNorm) -
      left.Velocity().DotProd(areaNorm);

  // calculate wave strengths (Cr - Cl)
  double waveStrength[NUMVARS - 1] = {
    ((right.P() - left.P()) - rhoR * aR * normVelDiff) / (2.0 * aR * aR),
    (right.Rho() - left.Rho()) - (right.P() - left.P()) / (aR * aR),
    ((right.P() - left.P()) + rhoR * aR * normVelDiff) / (2.0 * aR * aR),
    rhoR,
    rhoR * (right.Tke() - left.Tke()) + kR * (right.Rho() - left.Rho()) -
    (right.P() - left.P()) * kR / (aR * aR),
    rhoR * (right.Omega() - left.Omega()) + omR * (right.Rho() - left.Rho()) -
    (right.P() - left.P()) * omR / (aR * aR)};

  // calculate absolute value of wave speeds (L)
  double waveSpeed[NUMVARS - 1] = {
    fabs(velRSum - aR),  // left moving acoustic wave speed
    fabs(velRSum),  // entropy wave speed
    fabs(velRSum + aR),  // right moving acoustic wave speed
    fabs(velRSum),  // shear wave speed
    fabs(velRSum),  // turbulent eqn 1 (k) wave speed
    fabs(velRSum)};  // turbulent eqn 2 (omega) wave speed

  // calculate entropy fix (Harten) and adjust wave speeds if necessary
  // default setting for entropy fix to kick in
  constexpr auto entropyFix = 0.1;

  if (waveSpeed[0] < entropyFix) {
    waveSpeed[0] = 0.5 * (waveSpeed[0] * waveSpeed[0] /
                          entropyFix + entropyFix);
  }
  if (waveSpeed[2] < entropyFix) {
    waveSpeed[2] = 0.5 * (waveSpeed[2] * waveSpeed[2] /
                          entropyFix + entropyFix);
  }

  // calculate right eigenvectors (T)
  // calculate eigenvector due to left acoustic wave
  const genArray lAcousticEigV(1.0, uR - aR * areaNorm.X(),
                               vR - aR * areaNorm.Y(), wR - aR * areaNorm.Z(),
                               hR - aR * velRSum, kR, omR);

  // calculate eigenvector due to entropy wave
  const genArray entropyEigV(1.0, uR, vR, wR,
                             0.5 * (uR * uR + vR * vR + wR * wR),
                             0.0, 0.0);

  // calculate eigenvector due to right acoustic wave
  const genArray rAcousticEigV(1.0, uR + aR * areaNorm.X(),
                               vR + aR * areaNorm.Y(), wR + aR * areaNorm.Z(),
                               hR + aR * velRSum, kR, omR);

  // calculate eigenvector due to shear wave
  const genArray shearEigV(0.0,
                           (right.U() - left.U()) - normVelDiff * areaNorm.X(),
                           (right.V() - left.V()) - normVelDiff * areaNorm.Y(),
                           (right.W() - left.W()) - normVelDiff * areaNorm.Z(),
                           uR * (right.U() - left.U()) +
                           vR * (right.V() - left.V()) +
                           wR * (right.W() - left.W()) - velRSum * normVelDiff,
                           0.0, 0.0);

  // calculate eigenvector due to turbulent equation 1
  const genArray tkeEigV(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

  // calculate eigenvector due to turbulent equation 2
  const genArray omgEigV(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

  // calculate dissipation term ( eigenvector * wave speed * wave strength)
  genArray dissipation(0.0);
  for (auto ii = 0; ii < NUMVARS; ii++) {
    dissipation[ii] =
        waveSpeed[0] * waveStrength[0] *
            lAcousticEigV[ii]  // contribution from left acoustic wave
        +
        waveSpeed[1] * waveStrength[1] *
            entropyEigV[ii]  // contribution from entropy wave
        +
        waveSpeed[2] * waveStrength[2] *
            rAcousticEigV[ii]  // contribution from right acoustic wave
        +
        waveSpeed[3] * waveStrength[3] *
            shearEigV[ii]  // contribution from shear wave
        +
        waveSpeed[4] * waveStrength[4] *
            tkeEigV[ii]  // contribution from turbulent wave 1
        +
        waveSpeed[5] * waveStrength[5] *
            omgEigV[ii];  // contribution from turbulent wave 2
  }

  // calculate left/right physical flux
  inviscidFlux leftFlux(left, eqnState, areaNorm);
  inviscidFlux rightFlux(right, eqnState, areaNorm);

  // calculate numerical Roe flux
  leftFlux.RoeFlux(rightFlux, dissipation);

  return leftFlux;
}

/* Member function to calculate the Roe flux, given the left and right
 * convective fluxes as well as the dissipation term.
*/
void inviscidFlux::RoeFlux(const inviscidFlux &right, const genArray &diss) {
  // right -- right convective flux
  // diss -- dissipation term

  for (auto ii = 0; ii < NUMVARS; ii++) {
    data_[ii] = 0.5 * (data_[ii] + right.data_[ii] - diss[ii]);
  }
}

// function to calculate approximate Roe flux jacobians -- NOT USED in LUSGS
// method
void ApproxRoeFluxJacobian(const primVars &left, const primVars &right,
                           const idealGas &eqnState,
                           const vector3d<double> &areaNorm, double &maxWS,
                           squareMatrix &dF_dUl, squareMatrix &dF_dUr) {
  // left --> primative variables from left side
  // right --> primative variables from right side
  // eqnStat --> ideal gas equation of state
  // areaNorm --> face unit area vector
  // maxWS --> maximum wave speed
  // dF_dUl --> dF/dUl, derivative of the Roe flux wrt the left state
  // (conservative variables)
  // dF_dUr --> dF/dUlr, derivative of the Roe flux wrt the right state
  // (conservative variables)

  // check to see that output matricies are correct size
  if ((dF_dUl.Size() != 5) || (dF_dUr.Size() != 5)) {
    cerr << "ERROR: Input matricies to RoeFLuxJacobian function are not the "
            "correct size!" << endl;
  }

  // compute Rho averaged quantities
  const auto denRatio = sqrt(right.Rho() / left.Rho());
  const auto uR = (left.U() + denRatio * right.U()) /
      (1.0 + denRatio);  // Roe averaged u-velocity
  const auto vR = (left.V() + denRatio * right.V()) /
      (1.0 + denRatio);  // Roe averaged v-velocity
  const auto wR = (left.W() + denRatio * right.W()) /
      (1.0 + denRatio);  // Roe averaged w-velocity
  const auto hR = (left.Enthalpy(eqnState) + denRatio *
                   right.Enthalpy(eqnState)) /
      (1.0 + denRatio);  // Roe averaged total enthalpy
  // Roe averaged speed of sound
  const auto aR = sqrt((eqnState.Gamma() - 1.0) *
                       (hR - 0.5 * (uR * uR + vR * vR + wR * wR)));
  // Roe averaged face normal velocity
  const vector3d<double> velR(uR, vR, wR);

  // dot product of velocities (Roe, left, right) with unit area vector
  const auto velRNorm = velR.DotProd(areaNorm);
  const auto velLeftNorm = left.Velocity().DotProd(areaNorm);
  const auto velRightNorm = right.Velocity().DotProd(areaNorm);

  // calculate diagonal eigenvalue matrix |lambda|
  squareMatrix lambda(5);
  lambda.Zero();
  lambda.SetData(0, 0, fabs(velRNorm - aR));
  lambda.SetData(1, 1, fabs(velRNorm));
  lambda.SetData(2, 2, fabs(velRNorm + aR));
  lambda.SetData(3, 3, fabs(velRNorm));
  lambda.SetData(4, 4, fabs(velRNorm));

  // calculate Roe jacobian matrix A
  // contribution due to normal velocity eigenvalues
  squareMatrix A(5);

  // column zero
  A.SetData(0, 0,
            1.0 - 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR));
  A.SetData(
      1, 0,
      -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * uR +
          velRNorm * areaNorm.X());
  A.SetData(
      2, 0,
      -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * vR +
          velRNorm * areaNorm.Y());
  A.SetData(
      3, 0,
      -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * wR +
          velRNorm * areaNorm.Z());
  A.SetData(4, 0,
            velRNorm * velRNorm -
                0.5 * velR.MagSq() * (1.0 + 0.5 * (eqnState.Gamma() - 1.0) *
                                                velR.MagSq() / (aR * aR)));

  // column one
  A.SetData(0, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * uR);
  A.SetData(1, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * uR + 1.0 -
                      areaNorm.X() * areaNorm.X());
  A.SetData(2, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * uR -
                      areaNorm.Y() * areaNorm.X());
  A.SetData(3, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * uR -
                      areaNorm.Z() * areaNorm.X());
  A.SetData(
      4, 1,
      (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * uR -
          velRNorm * areaNorm.X());

  // column two
  A.SetData(0, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * vR);
  A.SetData(1, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * vR -
                      areaNorm.X() * areaNorm.Y());
  A.SetData(2, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * vR + 1.0 -
                      areaNorm.Y() * areaNorm.Y());
  A.SetData(3, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * vR -
                      areaNorm.Z() * areaNorm.Y());
  A.SetData(
      4, 2,
      (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * vR -
          velRNorm * areaNorm.Y());

  // column three
  A.SetData(0, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * wR);
  A.SetData(1, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * wR -
                      areaNorm.X() * areaNorm.Z());
  A.SetData(2, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * wR -
                      areaNorm.Y() * areaNorm.Z());
  A.SetData(3, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * wR + 1.0 -
                      areaNorm.Z() * areaNorm.Z());
  A.SetData(
      4, 3,
      (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * wR -
          velRNorm * areaNorm.Z());

  // column four
  A.SetData(0, 4, -(eqnState.Gamma() - 1.0) / (aR * aR));
  A.SetData(1, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * uR);
  A.SetData(2, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * vR);
  A.SetData(3, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * wR);
  A.SetData(4, 4,
            -(eqnState.Gamma() - 1.0) / (aR * aR) * velR.MagSq() / (aR * aR));

  A = fabs(velRNorm) * A;

  // contribution due to u - c wave
  squareMatrix temp(5);

  // column zero
  temp.SetData(0, 0,
               0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) +
                   0.5 * velRNorm / aR);
  temp.SetData(1, 0,
               (uR - aR * areaNorm.X()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) +
                    0.5 * velRNorm / aR));
  temp.SetData(2, 0,
               (vR - aR * areaNorm.Y()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) +
                    0.5 * velRNorm / aR));
  temp.SetData(3, 0,
               (wR - aR * areaNorm.Z()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) +
                    0.5 * velRNorm / aR));
  temp.SetData(4, 0,
               (hR - velRNorm * aR * areaNorm.X()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) +
                    0.5 * velRNorm / aR));

  // column one
  temp.SetData(0, 1, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR -
                         areaNorm.X() / (2.0 * aR));
  temp.SetData(1, 1, (uR - aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR -
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(2, 1, (vR - aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR -
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(3, 1, (wR - aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR -
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(4, 1, (hR - aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR -
                          areaNorm.X() / (2.0 * aR)));

  // column two
  temp.SetData(0, 2, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR -
                         areaNorm.Y() / (2.0 * aR));
  temp.SetData(1, 2, (uR - aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR -
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(2, 2, (vR - aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR -
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(3, 2, (wR - aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR -
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(4, 2, (hR - aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR -
                          areaNorm.Y() / (2.0 * aR)));

  // column three
  temp.SetData(0, 3, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR -
                         areaNorm.Z() / (2.0 * aR));
  temp.SetData(1, 3, (uR - aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR -
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(2, 3, (vR - aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR -
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(3, 3, (wR - aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR -
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(4, 3, (hR - aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR -
                          areaNorm.Z() / (2.0 * aR)));

  // column four
  temp.SetData(0, 4, (eqnState.Gamma() - 1.0) / (2.0 * aR * aR));
  temp.SetData(1, 4, (uR - aR * areaNorm.X()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(2, 4, (vR - aR * areaNorm.Y()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(3, 4, (wR - aR * areaNorm.Z()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(
      4, 4, (hR - aR * velRNorm) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR));

  A = A + (fabs(velRNorm - aR) * temp);

  // contribution due to u + c wave

  // column zero
  temp.SetData(0, 0,
               0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) -
                   0.5 * velRNorm / aR);
  temp.SetData(1, 0,
               (uR + aR * areaNorm.X()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) -
                    0.5 * velRNorm / aR));
  temp.SetData(2, 0,
               (vR + aR * areaNorm.Y()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) -
                    0.5 * velRNorm / aR));
  temp.SetData(3, 0,
               (wR + aR * areaNorm.Z()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) -
                    0.5 * velRNorm / aR));
  temp.SetData(4, 0,
               (hR + velRNorm * aR * areaNorm.X()) *
                   (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) -
                    0.5 * velRNorm / aR));

  // column one
  temp.SetData(0, 1, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR +
                         areaNorm.X() / (2.0 * aR));
  temp.SetData(1, 1, (uR + aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR +
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(2, 1, (vR + aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR +
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(3, 1, (wR + aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR +
                          areaNorm.X() / (2.0 * aR)));
  temp.SetData(4, 1, (hR + aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR +
                          areaNorm.X() / (2.0 * aR)));

  // column two
  temp.SetData(0, 2, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR +
                         areaNorm.Y() / (2.0 * aR));
  temp.SetData(1, 2, (uR + aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR +
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(2, 2, (vR + aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR +
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(3, 2, (wR + aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR +
                          areaNorm.Y() / (2.0 * aR)));
  temp.SetData(4, 2, (hR + aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR +
                          areaNorm.Y() / (2.0 * aR)));

  // column three
  temp.SetData(0, 3, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR +
                         areaNorm.Z() / (2.0 * aR));
  temp.SetData(1, 3, (uR + aR * areaNorm.X()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR +
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(2, 3, (vR + aR * areaNorm.Y()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR +
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(3, 3, (wR + aR * areaNorm.Z()) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR +
                          areaNorm.Z() / (2.0 * aR)));
  temp.SetData(4, 3, (hR + aR * velRNorm) *
                         (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR +
                          areaNorm.Z() / (2.0 * aR)));

  // column four
  temp.SetData(0, 4, (eqnState.Gamma() - 1.0) / (2.0 * aR * aR));
  temp.SetData(1, 4, (uR + aR * areaNorm.X()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(2, 4, (vR + aR * areaNorm.Y()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(3, 4, (wR + aR * areaNorm.Z()) * (eqnState.Gamma() - 1.0) /
                         (2.0 * aR * aR));
  temp.SetData(
      4, 4, (hR + aR * velRNorm) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR));

  A = A + (fabs(velRNorm + aR) * temp);

  // begin jacobian calculation

  // derivative of Roe flux wrt left conservative variables
  dF_dUl.Zero();

  // column zero
  dF_dUl.SetData(0, 0, 0.0);
  dF_dUl.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               left.Velocity().MagSq() * areaNorm.X() -
                           left.U() * velLeftNorm);
  dF_dUl.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               left.Velocity().MagSq() * areaNorm.Y() -
                           left.V() * velLeftNorm);
  dF_dUl.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               left.Velocity().MagSq() * areaNorm.Z() -
                           left.W() * velLeftNorm);
  dF_dUl.SetData(4, 0,
                 (0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() -
                  left.Enthalpy(eqnState)) *
                     velLeftNorm);

  // column one
  dF_dUl.SetData(1, 0, areaNorm.X());
  dF_dUl.SetData(1, 1, left.U() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * left.U() * areaNorm.X() +
                           velLeftNorm);
  dF_dUl.SetData(1, 2, left.V() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Y());
  dF_dUl.SetData(1, 3, left.W() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Z());
  dF_dUl.SetData(1, 4, left.Enthalpy(eqnState) * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * left.U() * velLeftNorm);

  // column two
  dF_dUl.SetData(2, 0, areaNorm.Y());
  dF_dUl.SetData(2, 1, left.U() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * left.V() * areaNorm.X());
  dF_dUl.SetData(2, 2, left.V() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Y() +
                           velLeftNorm);
  dF_dUl.SetData(2, 3, left.W() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Z());
  dF_dUl.SetData(2, 4, left.Enthalpy(eqnState) * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * left.V() * velLeftNorm);

  // column three
  dF_dUl.SetData(3, 0, areaNorm.Z());
  dF_dUl.SetData(3, 1, left.U() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * left.W() * areaNorm.X());
  dF_dUl.SetData(3, 2, left.V() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Y());
  dF_dUl.SetData(3, 3, left.W() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Z() +
                           velLeftNorm);
  dF_dUl.SetData(3, 4, left.Enthalpy(eqnState) * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * left.W() * velLeftNorm);

  // column four
  dF_dUl.SetData(4, 0, 0.0);
  dF_dUl.SetData(4, 1, (eqnState.Gamma() - 1.0) * areaNorm.X());
  dF_dUl.SetData(4, 2, (eqnState.Gamma() - 1.0) * areaNorm.Y());
  dF_dUl.SetData(4, 3, (eqnState.Gamma() - 1.0) * areaNorm.Z());
  dF_dUl.SetData(4, 4, eqnState.Gamma() * velLeftNorm);

  dF_dUl = 0.5 * (dF_dUl + A);

  // Compute derivative of flux wrt right conservative variables

  dF_dUr.Zero();

  // calculate flux derivatives
  // column zero
  dF_dUr.SetData(0, 0, 0.0);
  dF_dUr.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               right.Velocity().MagSq() * areaNorm.X() -
                           right.U() * velRightNorm);
  dF_dUr.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               right.Velocity().MagSq() * areaNorm.Y() -
                           right.V() * velRightNorm);
  dF_dUr.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) *
                               right.Velocity().MagSq() * areaNorm.Z() -
                           right.W() * velRightNorm);
  dF_dUr.SetData(4, 0,
                 (0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() -
                  right.Enthalpy(eqnState)) *
                     velRightNorm);

  // column one
  dF_dUr.SetData(1, 0, areaNorm.X());
  dF_dUr.SetData(1, 1, right.U() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * right.U() * areaNorm.X() +
                           velRightNorm);
  dF_dUr.SetData(1, 2, right.V() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Y());
  dF_dUr.SetData(1, 3, right.W() * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Z());
  dF_dUr.SetData(1, 4, right.Enthalpy(eqnState) * areaNorm.X() -
                           (eqnState.Gamma() - 1.0) * right.U() * velRightNorm);

  // column two
  dF_dUr.SetData(2, 0, areaNorm.Y());
  dF_dUr.SetData(2, 1, right.U() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * right.V() * areaNorm.X());
  dF_dUr.SetData(2, 2, right.V() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Y() +
                           velRightNorm);
  dF_dUr.SetData(2, 3, right.W() * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Z());
  dF_dUr.SetData(2, 4, right.Enthalpy(eqnState) * areaNorm.Y() -
                           (eqnState.Gamma() - 1.0) * right.V() * velRightNorm);

  // column three
  dF_dUr.SetData(3, 0, areaNorm.Z());
  dF_dUr.SetData(3, 1, right.U() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * right.W() * areaNorm.X());
  dF_dUr.SetData(3, 2, right.V() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Y());
  dF_dUr.SetData(3, 3, right.W() * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Z() +
                           velRightNorm);
  dF_dUr.SetData(3, 4, right.Enthalpy(eqnState) * areaNorm.Z() -
                           (eqnState.Gamma() - 1.0) * right.W() * velRightNorm);

  // column four
  dF_dUr.SetData(4, 0, 0.0);
  dF_dUr.SetData(4, 1, (eqnState.Gamma() - 1.0) * areaNorm.X());
  dF_dUr.SetData(4, 2, (eqnState.Gamma() - 1.0) * areaNorm.Y());
  dF_dUr.SetData(4, 3, (eqnState.Gamma() - 1.0) * areaNorm.Z());
  dF_dUr.SetData(4, 4, eqnState.Gamma() * velRightNorm);

  dF_dUr = 0.5 * (dF_dUr - A);
}

// //function to calculate Lax-Friedrichs flux jacobians
// void LaxFriedrichsFluxJacobian( const primVars &left, const primVars &right,
// const idealGas &eqnState, const vector3d<double>& areaVec,
// double &specRadL, double &specRadR, squareMatrix &dF_dUl, squareMatrix
// &dF_dUr){

//   //left --> primative variables from left side
//   //right --> primative variables from right side
//   //eqnStat --> ideal gas equation of state
//   //areaVec --> face area vector
//   //maxWS --> maximum wave speed
//   //dF_dUl --> dF/dUl, derivative of the Roe flux wrt the left state
// (conservative variables)
//   //dF_dUr --> dF/dUlr, derivative of the Roe flux wrt the right state
// (conservative variables)

//   //check to see that output matricies are correct size
//   if( (dF_dUl.Size() != 5) || (dF_dUr.Size() != 5)){
//     cerr << "ERROR: Input matricies to LaxFreidrichsFLuxJacobian function are
// not the correct size!" << endl;
//   }

//   vector3d<double> areaNorm = areaVec / areaVec.Mag();  //normalize area
// vector to unit vector

//   //dot product of velocities with unit area vector
//   double maxWS = ConvSpecRad(areaNorm, left, right, eqnState);

//   double velLeftNorm = left.Velocity().DotProd(areaNorm);
//   double velRightNorm = right.Velocity().DotProd(areaNorm);

//   //calculate spectral radii
//   specRadL = maxWS;
//   specRadR = maxWS;

//   //form spectral radii identity matrices
//   squareMatrix dissLeft(5);
//   dissLeft.Identity();
//   dissLeft = specRadL * dissLeft;

//   squareMatrix dissRight(5);
//   dissRight.Identity();
//   dissRight = specRadR * dissRight;

//   //begin jacobian calculation

//   dF_dUl.Zero();
//   dF_dUr.Zero();

//   //calculate flux derivatives wrt left state
//   //column zero
//   dF_dUl.SetData(0, 0, 0.0);
//   dF_dUl.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// left.Velocity().MagSq() * areaNorm.X() - left.U() * velLeftNorm);
//   dF_dUl.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// left.Velocity().MagSq() * areaNorm.Y() - left.V() * velLeftNorm);
//   dF_dUl.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// left.Velocity().MagSq() * areaNorm.Z() - left.W() * velLeftNorm);
//   dF_dUl.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) *
// left.Velocity().MagSq() - left.Enthalpy(eqnState)) * velLeftNorm);

//   //column one
//   dF_dUl.SetData(0, 1, areaNorm.X());
//   dF_dUl.SetData(1, 1, left.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// left.U() * areaNorm.X() + velLeftNorm);
//   dF_dUl.SetData(2, 1, left.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// left.U() * areaNorm.Y());
//   dF_dUl.SetData(3, 1, left.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// left.U() * areaNorm.Z());
//   dF_dUl.SetData(4, 1, left.Enthalpy(eqnState) * areaNorm.X() -
// (eqnState.Gamma() - 1.0) * left.U() * velLeftNorm);

//   //column two
//   dF_dUl.SetData(0, 2, areaNorm.Y());
//   dF_dUl.SetData(1, 2, left.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// left.V() * areaNorm.X());
//   dF_dUl.SetData(2, 2, left.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// left.V() * areaNorm.Y() + velLeftNorm);
//   dF_dUl.SetData(3, 2, left.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// left.V() * areaNorm.Z());
//   dF_dUl.SetData(4, 2, left.Enthalpy(eqnState) * areaNorm.Y() -
// (eqnState.Gamma() - 1.0) * left.V() * velLeftNorm);

//   //column three
//   dF_dUl.SetData(0, 3, areaNorm.Z());
//   dF_dUl.SetData(1, 3, left.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// left.W() * areaNorm.X());
//   dF_dUl.SetData(2, 3, left.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// left.W() * areaNorm.Y());
//   dF_dUl.SetData(3, 3, left.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// left.W() * areaNorm.Z() + velLeftNorm);
//   dF_dUl.SetData(4, 3, left.Enthalpy(eqnState) * areaNorm.Z() -
// (eqnState.Gamma() - 1.0) * left.W() * velLeftNorm);

//   //column four
//   dF_dUl.SetData(0, 4, 0.0);
//   dF_dUl.SetData(1, 4, (eqnState.Gamma() - 1.0) * areaNorm.X());
//   dF_dUl.SetData(2, 4, (eqnState.Gamma() - 1.0) * areaNorm.Y());
//   dF_dUl.SetData(3, 4, (eqnState.Gamma() - 1.0) * areaNorm.Z());
//   dF_dUl.SetData(4, 4, eqnState.Gamma() * velLeftNorm);

//   dF_dUl = 0.5 * (dF_dUl + dissLeft);

//   //calculate flux derivatives wrt right state
//   //column zero
//   dF_dUr.SetData(0, 0, 0.0);
//   dF_dUr.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// right.Velocity().MagSq() * areaNorm.X() - right.U() * velRightNorm);
//   dF_dUr.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// right.Velocity().MagSq() * areaNorm.Y() - right.V() * velRightNorm);
//   dF_dUr.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) *
// right.Velocity().MagSq() * areaNorm.Z() - right.W() * velRightNorm);
//   dF_dUr.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) *
// right.Velocity().MagSq() - right.Enthalpy(eqnState)) * velRightNorm);

//   //column one
//   dF_dUr.SetData(0, 1, areaNorm.X());
//   dF_dUr.SetData(1, 1, right.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// right.U() * areaNorm.X() + velRightNorm);
//   dF_dUr.SetData(2, 1, right.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// right.U() * areaNorm.Y());
//   dF_dUr.SetData(3, 1, right.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) *
// right.U() * areaNorm.Z());
//   dF_dUr.SetData(4, 1, right.Enthalpy(eqnState) * areaNorm.X() -
// (eqnState.Gamma() - 1.0) * right.U() * velRightNorm);

//   //column two
//   dF_dUr.SetData(0, 2, areaNorm.Y());
//   dF_dUr.SetData(1, 2, right.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// right.V() * areaNorm.X());
//   dF_dUr.SetData(2, 2, right.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// right.V() * areaNorm.Y() + velRightNorm);
//   dF_dUr.SetData(3, 2, right.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) *
// right.V() * areaNorm.Z());
//   dF_dUr.SetData(4, 2, right.Enthalpy(eqnState) * areaNorm.Y() -
// (eqnState.Gamma() - 1.0) * right.V() * velRightNorm);

//   //column three
//   dF_dUr.SetData(0, 3, areaNorm.Z());
//   dF_dUr.SetData(1, 3, right.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// right.W() * areaNorm.X());
//   dF_dUr.SetData(2, 3, right.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// right.W() * areaNorm.Y());
//   dF_dUr.SetData(3, 3, right.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) *
// right.W() * areaNorm.Z() + velRightNorm);
//   dF_dUr.SetData(4, 3, right.Enthalpy(eqnState) * areaNorm.Z() -
// (eqnState.Gamma() - 1.0) * right.W() * velRightNorm);

//   //column four
//   dF_dUr.SetData(0, 4, 0.0);
//   dF_dUr.SetData(1, 4, (eqnState.Gamma() - 1.0) * areaNorm.X());
//   dF_dUr.SetData(2, 4, (eqnState.Gamma() - 1.0) * areaNorm.Y());
//   dF_dUr.SetData(3, 4, (eqnState.Gamma() - 1.0) * areaNorm.Z());
//   dF_dUr.SetData(4, 4, eqnState.Gamma() * velRightNorm);

//   dF_dUr = 0.5 * (dF_dUr - dissRight);

// }

// non-member functions
// -----------------------------------------------------------

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, inviscidFlux &flux) {
  for (auto ii = 0; ii < NUMVARS; ii++) {
    os << flux.data_[ii];
    if (ii != NUMVARS - 1) {
      os << ", ";
    }
  }
  return os;
}

// member function for scalar multiplication
inviscidFlux inviscidFlux::operator*(const double &scalar) const {
  inviscidFlux temp = *this;
  for (auto &val : temp.data_) {
    val *= scalar;
  }
  return temp;
}

// friend function to allow multiplication (elementwise) from either direction
inviscidFlux operator*(const double &scalar, const inviscidFlux &flux) {
  inviscidFlux temp;
  for (auto ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] = flux.data_[ii] * scalar;
  }
  return temp;
}

// operator overload for addition
inviscidFlux inviscidFlux::operator+(const inviscidFlux &invf2) const {
  inviscidFlux invf1 = *this;
  for (auto ii = 0; ii < NUMVARS; ii++) {
    invf1.data_[ii] += invf2.data_[ii];
  }
  return invf1;
}

// operator overload for subtraction
inviscidFlux inviscidFlux::operator-(const inviscidFlux &invf2) const {
  inviscidFlux invf1 = *this;
  for (auto ii = 0; ii < NUMVARS; ii++) {
    invf1.data_[ii] -= invf2.data_[ii];
  }
  return invf1;
}

// member function for scalar division
inviscidFlux inviscidFlux::operator/(const double &scalar) const {
  inviscidFlux temp = *this;
  for (auto &val : temp.data_) {
    val /= scalar;
  }
  return temp;
}

// friend function to allow division (elementwise) from either direction
inviscidFlux operator/(const double &scalar, const inviscidFlux &flux) {
  inviscidFlux temp;
  for (auto ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] = scalar / flux.data_[ii];
  }
  return temp;
}

// convert the inviscid flux to a genArray
genArray inviscidFlux::ConvertToGenArray() const {
  return genArray(this->RhoVel(), this->RhoVelU(), this->RhoVelV(),
                  this->RhoVelW(), this->RhoVelH(), this->RhoVelK(),
                  this->RhoVelO());
}

// function to take in the primative variables, equation of state, face area
// vector, and conservative variable update and calculate the change in the
// convective flux
genArray ConvectiveFluxUpdate(const primVars &state, const idealGas &eqnState,
                              const unique_ptr<turbModel> &turb,
                              const vector3d<double> &normArea,
                              const genArray &du) {
  // get inviscid flux of old state
  const inviscidFlux oldFlux(state, eqnState, normArea);

  // get updated state in primative variables
  const primVars stateUpdate = state.UpdateWithConsVars(eqnState, du, turb);

  // get updated inviscid flux
  const inviscidFlux newFlux(stateUpdate, eqnState, normArea);

  // calculate difference in flux
  inviscidFlux dFlux = newFlux - oldFlux;

  return dFlux.ConvertToGenArray();
}

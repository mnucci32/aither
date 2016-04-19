/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

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
#include <algorithm>  // max
#include "fluxJacobian.hpp"
#include "turbulence.hpp"    // turbModel
#include "input.hpp"         // input
#include "primVars.hpp"      // primVars
#include "genArray.hpp"      // genArray
#include "matrix.hpp"        // squareMatrix
#include "inviscidFlux.hpp"  // ConvectiveFluxUpdate

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std:: unique_ptr;

// constructor
fluxJacobian::fluxJacobian(const primVars &state,
                           const unitVec3dMag<double> &fAreaL,
                           const unitVec3dMag<double> &fAreaR,
                           const idealGas &eos, const sutherland &suth,
                           const double &vol,
                           const unique_ptr<turbModel> &turb,
                           const input &inp, const bool &mainDiagonal) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model
  // inp -- input variables
  // mainDiagonal -- flag to determine if jacobian is for main diagonal of
  //                 implicit matrix

  // initialize jacobians
  flowJacobian_ = 0.0;
  turbJacobian_ = 0.0;

  const auto isTurbulent = inp.IsTurbulent();

  this->AddInviscidJacobian(state, fAreaL, fAreaR, eos, turb, isTurbulent);

  // if viscous, add viscous jacobians
  if (inp.IsViscous()) {
    this->AddViscousJacobian(state, fAreaL, fAreaR, eos, suth, vol,
                             turb, isTurbulent);
  }

  // source term is only added to main diagonal,
  // and only for turbulence equations
  if (mainDiagonal && isTurbulent) {
    this->AddTurbSourceJacobian(state, suth, vol, turb);
  }
}

// member functions
// member function to add inviscid jacobians
void fluxJacobian::AddInviscidJacobian(const primVars &state,
                                       const unitVec3dMag<double> &fAreaL,
                                       const unitVec3dMag<double> &fAreaR,
                                       const idealGas &eos,
                                       const unique_ptr<turbModel> &turb,
                                       const bool &isTurbulent) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // turb -- turbulence model
  // isTurbulent -- flag to determine if simulation is turbulent

  flowJacobian_ += state.InvCellSpectralRadius(fAreaL, fAreaR, eos);

  if (isTurbulent) {
    turbJacobian_ += turb->InviscidSpecRad(state, fAreaL, fAreaR);
  }
}

// member function to add viscous jacobians
void fluxJacobian::AddViscousJacobian(const primVars &state,
                                      const unitVec3dMag<double> &fAreaL,
                                      const unitVec3dMag<double> &fAreaR,
                                      const idealGas &eos,
                                      const sutherland &suth,
                                      const double &vol,
                                      const unique_ptr<turbModel> &turb,
                                      const bool &isTurbulent) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model
  // isTurbulent -- flag to determine if simulation is turbulent

  // factor of 2 because viscous spectral radius is not halved (Blazek 6.53)
  flowJacobian_ += 2.0 * state.ViscCellSpectralRadius(fAreaL, fAreaR, eos, suth,
                                                      vol, turb);

  if (isTurbulent) {
    // factor of 2 because viscous spectral radius is not halved (Blazek 6.53)
    turbJacobian_ += 2.0 * turb->ViscSpecRad(state, fAreaL, fAreaR, eos, suth,
                                             vol);
  }
}

// member function to add source jacobians
// this should only be used for the main diagonal
void fluxJacobian::AddTurbSourceJacobian(const primVars &state,
                                         const sutherland &suth,
                                         const double &vol,
                                         const unique_ptr<turbModel> &turb) {
  // state -- primative variables
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model

  turbJacobian_ -= turb->SrcSpecRad(state, suth, vol);
}

// member function to multiply the flux jacobians with a genArray
genArray fluxJacobian::ArrayMult(genArray arr) const {
  arr[0] *= flowJacobian_;
  arr[1] *= flowJacobian_;
  arr[2] *= flowJacobian_;
  arr[3] *= flowJacobian_;
  arr[4] *= flowJacobian_;

  arr[5] *= turbJacobian_;
  arr[6] *= turbJacobian_;

  return arr;
}

// function to take the inverse of a flux jacobian
fluxJacobian fluxJacobian::Inverse(const bool &isTurbulent) const {
  auto inv = *this;
  inv.flowJacobian_ = 1.0 / inv.flowJacobian_;

  if (isTurbulent) {
    inv.turbJacobian_ = 1.0 / inv.turbJacobian_;
  }
  return inv;
}

// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}

/* Function to calculate Rusanov flux jacobian. The Rusanov flux is defined as
shown below.

  F = 0.5 * (F(Ul) + F(Ur) - L(Ul, Ur) * (Ur - Ul)

Differentiating by the left and right states gives the left and right flux
jacobians.

  dF_Ul = 0.5 * (A(Ul) + L(Ul, Ur))
  dF_Ur = 0.5 * (A(Ur) - L(Ul, Ur))

In the above equations the dissipation term L is held constant during
differentiation. A represents the convective flux jacobian matrix.
 */
squareMatrix RusanovFluxJacobian(const primVars &left, const primVars &right,
                                 const idealGas &eos,
                                 const vector3d<double> &areaNorm,
                                 const bool &positive) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eos -- ideal gas equation of state
  // areaNorm -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation

  // dot product of velocities with unit area vector
  const auto specRad = std::max(fabs(left.Velocity().DotProd(areaNorm)) +
                                left.SoS(eos),
                                fabs(right.Velocity().DotProd(areaNorm)) +
                                right.SoS(eos));

  // form dissipation matrix based on spectral radius
  squareMatrix dissipation(5);
  dissipation.Identity();
  dissipation *= specRad;

  // begin jacobian calculation
  const auto fluxJac = positive ?
      InvFluxJacobian(left, eos, areaNorm) :
      InvFluxJacobian(right, eos, areaNorm);

  return positive ? 0.5 * (fluxJac + dissipation) :
      0.5 * (fluxJac - dissipation);
}

// function to calculate inviscid flux jacobian
squareMatrix InvFluxJacobian(const primVars &state,
                             const idealGas &eqnState,
                             const vector3d<double> &areaNorm) {
  // state -- primative variables from left side
  // eqnState -- ideal gas equation of state
  // areaNorm -- face area vector

  const auto velNorm = state.Velocity().DotProd(areaNorm);
  const auto gammaMinusOne = eqnState.Gamma() - 1.0;
  const auto phi = 0.5 * gammaMinusOne * state.Velocity().MagSq();
  const auto a1 = eqnState.Gamma() * state.Energy(eqnState) - phi;
  const auto a3 = eqnState.Gamma() - 2.0;

  // begin jacobian calculation
  squareMatrix A(5);

  // calculate flux derivatives wrt left state
  // column zero
  A(0, 0) = 0.0;
  A(1, 0) = phi * areaNorm.X() - state.U() * velNorm;
  A(2, 0) = phi * areaNorm.Y() - state.V() * velNorm;
  A(3, 0) = phi * areaNorm.Z() - state.W() * velNorm;
  A(4, 0) = velNorm * (phi - a1);

  // column one
  A(0, 1) = areaNorm.X();
  A(1, 1) = velNorm - a3 * areaNorm.X() * state.U();
  A(2, 1) = state.V() * areaNorm.X() - gammaMinusOne * state.U() * areaNorm.Y();
  A(3, 1) = state.W() * areaNorm.X() - gammaMinusOne * state.U() * areaNorm.Z();
  A(4, 1) = a1 * areaNorm.X() - gammaMinusOne * state.U() * velNorm;

  // column two
  A(0, 2) = areaNorm.Y();
  A(1, 2) = state.U() * areaNorm.Y() - gammaMinusOne * state.V() * areaNorm.X();
  A(2, 2) = velNorm - a3 * areaNorm.Y() * state.V();
  A(3, 2) = state.W() * areaNorm.Y() - gammaMinusOne * state.V() * areaNorm.Z();
  A(4, 2) = a1 * areaNorm.Y() - gammaMinusOne * state.V() * velNorm;

  // column three
  A(0, 3) = areaNorm.Z();
  A(1, 3) = state.U() * areaNorm.Z() - gammaMinusOne * state.W() * areaNorm.X();
  A(2, 3) = state.V() * areaNorm.Z() - gammaMinusOne * state.W() * areaNorm.Y();
  A(3, 3) = velNorm - a3 * areaNorm.Z() * state.W();
  A(4, 3) = a1 * areaNorm.Z() - gammaMinusOne * state.W() * velNorm;

  // column four
  A(0, 4) = 0.0;
  A(1, 4) = gammaMinusOne * areaNorm.X();
  A(2, 4) = gammaMinusOne * areaNorm.Y();
  A(3, 4) = gammaMinusOne * areaNorm.Z();
  A(4, 4) = eqnState.Gamma() * velNorm;

  return A;
}

/* Function to calculate approximate Roe flux jacobian. The Roe flux is
defined as shown below.

  F = 0.5 * (F(Ul) + F(Ur) - Aroe(Ul, Ur) * (Ur - Ul)

Differentiating by the left and right states gives the left and right flux
jacobians.

  dF_Ul = 0.5 * (A(Ul) + Aroe(Ul, Ur))
  dF_Ur = 0.5 * (A(Ur) - Aroe(Ul, Ur))

In the above equations the Roe matrix Aroe is held constant during
differentiation. A represents the convective flux jacobian matrix.
 */
squareMatrix ApproxRoeFluxJacobian(const primVars &left, const primVars &right,
                                   const idealGas &eos,
                                   const vector3d<double> &areaNorm,
                                   const bool &positive) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eos -- ideal gas equation of state
  // areaNorm -- face unit area vector
  // positive -- flag to determine whether to add or subtract dissipation

  // compute Roe averaged state
  const auto roeAvg = RoeAveragedState(left, right, eos);

  // compute Roe matrix
  const auto Aroe = InvFluxJacobian(roeAvg, eos, areaNorm);

  // compute convective flux jacobian
  const auto fluxJac = positive ? InvFluxJacobian(left, eos, areaNorm) :
      InvFluxJacobian(right, eos, areaNorm);

  return positive ? 0.5 * (fluxJac + Aroe) : 0.5 * (fluxJac - Aroe);
}

genArray RusanovOffDiagonal(const primVars &state, const genArray &update,
                            const unitVec3dMag<double> &fAreaL,
                            const unitVec3dMag<double> &fAreaR,
                            const double &vol, const idealGas &eos,
                            const sutherland &suth,
                            const unique_ptr<turbModel> &turb,
                            const input &inp, const bool &positive) {
  // state -- primative variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fAreaL -- face area vector on off diagonal boundary
  // fAreaR -- face area vector opposite off diagonal boundary
  // vol -- cell volume
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate = state.UpdateWithConsVars(eos, update, turb);

  // calculate updated convective flux
  const auto fluxChange = ConvectiveFluxUpdate(state, stateUpdate, eos,
                                               fAreaL.UnitVector());

  // can't use stored cell spectral radius b/c it has contribuitons from i, j, k
  const uncoupledScalar specRad(state.CellSpectralRadius(fAreaL, fAreaR, eos,
                                                         suth, vol, turb,
                                                         inp.IsViscous()),
                                turb->SpectralRadius(state, fAreaL, fAreaR, eos,
                                                     suth, vol, false));

  return positive ?
    0.5 * (fAreaL.Mag() * fluxChange + specRad.ArrayMult(update)) :
    0.5 * (fAreaL.Mag() * fluxChange - specRad.ArrayMult(update));
}


genArray RoeOffDiagonal(const primVars &left, const primVars &right,
                        const genArray &update,
                        const unitVec3dMag<double> &fAreaL,
                        const unitVec3dMag<double> &fAreaR,
                        const double &vol, const idealGas &eos,
                        const sutherland &suth,
                        const unique_ptr<turbModel> &turb,
                        const input &inp, const bool &positive) {
  // left -- primative variables at left side
  // right -- primative variables at right side
  // update -- conserved variable update at off diagonal
  // fAreaL -- face area vector on off diagonal boundary
  // fAreaR -- face area vector opposite off diagonal boundary
  // vol -- cell volume
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  const auto areaNorm = fAreaL.UnitVector();

  // calculate Roe flux with old variables
  const auto oldFlux = RoeFlux(left, right, eos, areaNorm);

  // calculate updated Roe flux
  const auto stateUpdate = positive ?
      left.UpdateWithConsVars(eos, update, turb) :
      right.UpdateWithConsVars(eos, update, turb);

  const auto newFlux = positive ?
    RoeFlux(stateUpdate, right, eos, areaNorm) :
    RoeFlux(left, stateUpdate, eos, areaNorm);

  // add contribution for viscous terms
  uncoupledScalar specRad(0.0, 0.0);
  // squareMatrix jac(5);
  if (inp.IsViscous()) {
    const auto offState = positive ? left : right;

    specRad.AddToFlowVariable(
        offState.ViscCellSpectralRadius(fAreaL, fAreaR, eos, suth, vol, turb));

    // jac = ApproxTSLJacobian(offState, eos, suth, areaNorm, 1.0, turb);

    if (inp.IsTurbulent()) {
      specRad.AddToTurbVariable(turb->ViscSpecRad(offState, fAreaL, fAreaR, eos,
                                                  suth, vol));
    }
  }

  // const auto invFluxJac = RusanovFluxJacobian(left, right, eos, areaNorm,
  //                                             positive);


  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  return positive ?
    fAreaL.Mag() * ((newFlux - oldFlux).ConvertToGenArray()) +
      0.5 * specRad.ArrayMult(update) :
    fAreaL.Mag() * ((newFlux - oldFlux).ConvertToGenArray()) -
      0.5 * specRad.ArrayMult(update);
    // fAreaL.Mag() * invFluxJac.VecMult(update) + 0.5 * specRadArr * update :
    // fAreaL.Mag() * invFluxJac.VecMult(update) - 0.5 * specRadArr * update;
}

// change of variable matrix going frim primative to conservative variables
// from Dwight
squareMatrix DelPrimativeDelConservative(const primVars &state,
                                         const idealGas &eos) {
  // state -- primative variables
  // eos -- equation of state

  const auto gammaMinusOne = eos.Gamma() - 1.0;
  const auto invRho = 1.0 / state.Rho();

  squareMatrix dPdC(5);

  // assign first column
  dPdC(0, 0) = 1.0;
  dPdC(1, 0) = -invRho * state.U();
  dPdC(2, 0) = -invRho * state.V();
  dPdC(3, 0) = -invRho * state.W();
  dPdC(4, 0) = 0.5 * gammaMinusOne * state.Velocity().DotProd(state.Velocity());

  // assign second column
  dPdC(1, 1) = invRho;
  dPdC(4, 1) = -gammaMinusOne * state.U();

  // assign third column
  dPdC(2, 2) = invRho;
  dPdC(4, 2) = -gammaMinusOne * state.V();

  // assign fourth column
  dPdC(3, 3) = invRho;
  dPdC(4, 3) = -gammaMinusOne * state.W();

  // assign fifth column
  dPdC(4, 4) = gammaMinusOne;

  return dPdC;
}

// approximate thin shear layer jacobian following implementation in Dwight.
// does not use any gradients
squareMatrix ApproxTSLJacobian(const primVars &state, const idealGas &eos,
                               const sutherland &suth,
                               const vector3d<double> &area, const double &dist,
                               const unique_ptr<turbModel> &turb) {
  // state -- primative variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // area -- face area unit vector
  // dist -- distance from cell center to cell center
  // turb --  turbulence model

  squareMatrix jacobian(5);

  const auto mu = suth.EffectiveViscosity(state.Temperature(eos));
  const auto mut = turb->EddyViscNoLim(state) * suth.NondimScaling();
  const auto velNorm = state.Velocity().DotProd(area);

  // assign first column
  jacobian(4, 0) = -eos.Conductivity(mu + mut) * state.Temperature(eos) /
    ((mu + mut) * state.Rho());

  // assign second column
  jacobian(1, 1) = (1.0 / 3.0) * area.X() * area.X() + 1.0;
  jacobian(1, 2) = (1.0 / 3.0) * area.X() * area.Y();
  jacobian(1, 3) = (1.0 / 3.0) * area.X() * area.Z();
  jacobian(1, 4) = (1.0 / 3.0) * area.X() * velNorm + state.U();

  // assign third column
  jacobian(2, 1) = (1.0 / 3.0) * area.Y() * area.X();
  jacobian(2, 2) = (1.0 / 3.0) * area.Y() * area.Y() + 1.0;
  jacobian(2, 3) = (1.0 / 3.0) * area.Y() * area.Z();
  jacobian(2, 4) = (1.0 / 3.0) * area.Y() * velNorm + state.V();

  // assign fourth column
  jacobian(3, 1) = (1.0 / 3.0) * area.Z() * area.X();
  jacobian(3, 2) = (1.0 / 3.0) * area.Z() * area.Y();
  jacobian(3, 3) = (1.0 / 3.0) * area.Z() * area.Z() + 1.0;
  jacobian(3, 4) = (1.0 / 3.0) * area.Z() * velNorm + state.W();

  // assign fifth column
  jacobian(4, 4) = eos.Conductivity(mu + mut) / ((mu + mut) * state.Rho());

  jacobian *= (mu + mut) / dist;

  const auto delPrimDelCons = DelPrimativeDelConservative(state, eos);

  return jacobian.MatMult(delPrimDelCons);
}

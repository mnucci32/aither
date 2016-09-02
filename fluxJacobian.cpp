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

#include <iostream>        // cout
#include <cmath>  // sqrt
#include <memory>
#include <string>
#include <algorithm>  // max
#include "fluxJacobian.hpp"
#include "turbulence.hpp"    // turbModel
#include "input.hpp"         // input
#include "primVars.hpp"      // primVars
#include "genArray.hpp"      // genArray
#include "inviscidFlux.hpp"  // ConvectiveFluxUpdate
#include "utility.hpp"       // TauNormal

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;

// constructor
// if constructed with two doubles, create scalar squareMatrix
fluxJacobian::fluxJacobian(const double &flow, const double &turb) {
  flowJacobian_ = squareMatrix(1);
  flowJacobian_ += flow;

  turbJacobian_ = squareMatrix(1);
  turbJacobian_ += turb;
}

// if constructed with two intss, create scalar squareMatrix with given size
fluxJacobian::fluxJacobian(const int &flowSize, const int &turbSize) {
  flowJacobian_ = squareMatrix(flowSize);
  turbJacobian_ = squareMatrix(turbSize);
}


// member functions
// member function to multiply the flux jacobians with a genArray
genArray fluxJacobian::ArrayMult(genArray arr) const {
  if (this->IsScalar()) {
    arr[0] *= flowJacobian_(0, 0);
    arr[1] *= flowJacobian_(0, 0);
    arr[2] *= flowJacobian_(0, 0);
    arr[3] *= flowJacobian_(0, 0);
    arr[4] *= flowJacobian_(0, 0);

    arr[5] *= turbJacobian_(0, 0);
    arr[6] *= turbJacobian_(0, 0);
  } else {
    arr = flowJacobian_.ArrayMult(arr);
    arr = turbJacobian_.ArrayMult(arr, NUMFLOWVARS);
  }
  return arr;
}

bool fluxJacobian::IsScalar() const {
  return (flowJacobian_.Size() > 1) ? false : true;
}

// function to take the inverse of a flux jacobian
void fluxJacobian::Inverse(const bool &isTurbulent) {
  flowJacobian_.Inverse();

  if (isTurbulent) {
    turbJacobian_.Inverse();
  }
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
void fluxJacobian::RusanovFluxJacobian(const primVars &state,
                                       const idealGas &eos,
                                       const unitVec3dMag<double> &area,
                                       const bool &positive,
                                       const input &inp,
                                       const unique_ptr<turbModel> &turb) {
  // state -- primative variables at face
  // eos -- ideal gas equation of state
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model

  // dot product of velocities with area vector (with magnitude)
  const auto specRad = 0.5 * area.Mag() *
      fabs(state.Velocity().DotProd(area.UnitVector())) + state.SoS(eos);

  // form dissipation matrix based on spectral radius
  fluxJacobian dissipation(inp.NumFlowEquations(), inp.NumTurbEquations());
  dissipation.flowJacobian_.Identity();
  dissipation.flowJacobian_ *= specRad;

  // begin jacobian calculation
  this->InvFluxJacobian(state, eos, area, inp, turb);

  // compute turbulent dissipation if necessary
  if (inp.IsTurbulent()) {
    dissipation.turbJacobian_ = 0.5 * turb->InviscidDissJacobian(state, area);
  }

  positive ? (*this) += dissipation : (*this) -= dissipation;
}

// function to calculate inviscid flux jacobian
void fluxJacobian::InvFluxJacobian(const primVars &state,
                                   const idealGas &eqnState,
                                   const unitVec3dMag<double> &area,
                                   const input &inp,
                                   const unique_ptr<turbModel> &turb) {
  // state -- primative variables at face
  // eqnState -- ideal gas equation of state
  // area -- face area vector
  // inp -- input variables
  // turb -- turbulence model

  const auto velNorm = state.Velocity().DotProd(area.UnitVector());
  const auto gammaMinusOne = eqnState.Gamma() - 1.0;
  const auto phi = 0.5 * gammaMinusOne * state.Velocity().MagSq();
  const auto a1 = eqnState.Gamma() * state.Energy(eqnState) - phi;
  const auto a3 = eqnState.Gamma() - 2.0;

  // begin jacobian calculation
  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  // calculate flux derivatives wrt left state
  // column zero
  flowJacobian_(0, 0) = 0.0;
  flowJacobian_(1, 0) = phi * area.UnitVector().X() - state.U() * velNorm;
  flowJacobian_(2, 0) = phi * area.UnitVector().Y() - state.V() * velNorm;
  flowJacobian_(3, 0) = phi * area.UnitVector().Z() - state.W() * velNorm;
  flowJacobian_(4, 0) = velNorm * (phi - a1);

  // column one
  flowJacobian_(0, 1) = area.UnitVector().X();
  flowJacobian_(1, 1) = velNorm - a3 * area.UnitVector().X() * state.U();
  flowJacobian_(2, 1) = state.V() * area.UnitVector().X() -
      gammaMinusOne * state.U() * area.UnitVector().Y();
  flowJacobian_(3, 1) = state.W() * area.UnitVector().X() -
      gammaMinusOne * state.U() * area.UnitVector().Z();
  flowJacobian_(4, 1) = a1 * area.UnitVector().X() - gammaMinusOne * state.U()
      * velNorm;

  // column two
  flowJacobian_(0, 2) = area.UnitVector().Y();
  flowJacobian_(1, 2) = state.U() * area.UnitVector().Y() -
      gammaMinusOne * state.V() * area.UnitVector().X();
  flowJacobian_(2, 2) = velNorm - a3 * area.UnitVector().Y() * state.V();
  flowJacobian_(3, 2) = state.W() * area.UnitVector().Y() -
      gammaMinusOne * state.V() * area.UnitVector().Z();
  flowJacobian_(4, 2) = a1 * area.UnitVector().Y() - gammaMinusOne * state.V()
      * velNorm;

  // column three
  flowJacobian_(0, 3) = area.UnitVector().Z();
  flowJacobian_(1, 3) = state.U() * area.UnitVector().Z() -
      gammaMinusOne * state.W() * area.UnitVector().X();
  flowJacobian_(2, 3) = state.V() * area.UnitVector().Z() -
      gammaMinusOne * state.W() * area.UnitVector().Y();
  flowJacobian_(3, 3) = velNorm - a3 * area.UnitVector().Z() * state.W();
  flowJacobian_(4, 3) = a1 * area.UnitVector().Z() - gammaMinusOne * state.W()
      * velNorm;

  // column four
  flowJacobian_(0, 4) = 0.0;
  flowJacobian_(1, 4) = gammaMinusOne * area.UnitVector().X();
  flowJacobian_(2, 4) = gammaMinusOne * area.UnitVector().Y();
  flowJacobian_(3, 4) = gammaMinusOne * area.UnitVector().Z();
  flowJacobian_(4, 4) = eqnState.Gamma() * velNorm;

  // multiply by 0.5 b/c 2 faces per i, j, k direction
  flowJacobian_ *= 0.5 * area.Mag();
  
  // turbulent jacobian here
  if (inp.IsTurbulent()) {
    turbJacobian_ = 0.5 * turb->InviscidConvJacobian(state, area);
  }
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
void fluxJacobian::ApproxRoeFluxJacobian(const primVars &left,
                                         const primVars &right,
                                         const idealGas &eos,
                                         const unitVec3dMag<double> &area,
                                         const bool &positive,
                                         const input &inp,
                                         const unique_ptr<turbModel> &turb) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eos -- ideal gas equation of state
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model
  
  // compute Roe averaged state
  const auto roeAvg = RoeAveragedState(left, right, eos);

  // compute Roe matrix
  fluxJacobian roeMatrix;
  roeMatrix.InvFluxJacobian(roeAvg, eos, area, inp, turb);

  // compute convective flux jacobian
  positive ? this->InvFluxJacobian(left, eos, area, inp, turb) :
      this->InvFluxJacobian(right, eos, area, inp, turb);

  positive ? (*this) += roeMatrix : (*this) -= roeMatrix;
}

// change of variable matrix going from primative to conservative variables
// from Dwight
void fluxJacobian::DelPrimativeDelConservative(const primVars &state,
                                               const idealGas &eos,
                                               const input &inp) {
  // state -- primative variables
  // eos -- equation of state
  // inp -- input variables

  const auto gammaMinusOne = eos.Gamma() - 1.0;
  const auto invRho = 1.0 / state.Rho();

  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  // assign column 0
  flowJacobian_(0, 0) = 1.0;
  flowJacobian_(1, 0) = -invRho * state.U();
  flowJacobian_(2, 0) = -invRho * state.V();
  flowJacobian_(3, 0) = -invRho * state.W();
  flowJacobian_(4, 0) = 0.5 * gammaMinusOne *
      state.Velocity().DotProd(state.Velocity());

  // assign column 1
  flowJacobian_(1, 1) = invRho;
  flowJacobian_(4, 1) = -gammaMinusOne * state.U();

  // assign column 2
  flowJacobian_(2, 2) = invRho;
  flowJacobian_(4, 2) = -gammaMinusOne * state.V();

  // assign column 3
  flowJacobian_(3, 3) = invRho;
  flowJacobian_(4, 3) = -gammaMinusOne * state.W();

  // assign column 4
  flowJacobian_(4, 4) = gammaMinusOne;

  // turbulent jacobian here
  if (inp.IsTurbulent()) {
    turbJacobian_(0, 0) = invRho;
    turbJacobian_(1, 1) = invRho;
  }
}

// approximate thin shear layer jacobian following implementation in Dwight.
// does not use any gradients
void fluxJacobian::ApproxTSLJacobian(const primVars &state,
                                     const double &lamVisc,
                                     const double &turbVisc, const double &f1,
                                     const idealGas &eos,
                                     const sutherland &suth,
                                     const unitVec3dMag<double> &area,
                                     const double &dist,
                                     const unique_ptr<turbModel> &turb,
                                     const input &inp, const bool &left,
                                     const tensor<double> &vGrad) {
  // state -- primative variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // area -- face area vector
  // dist -- distance from cell center to cell center
  // turb --  turbulence model
  // inp -- input variables
  // left -- flag that is negative if using left state
  // vGrad -- velocity gradient

  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  const auto mu = suth.NondimScaling() * lamVisc;
  const auto mut = suth.NondimScaling() * turbVisc;
  const auto velNorm = state.Velocity().DotProd(area.UnitVector());

  const auto tauNorm = TauNormal(vGrad, area.UnitVector(), mu, mut, suth);

  // DEBUG
  // const auto fac = left ? -1.0 : 1.0;
  auto fac = left ? -1.0 : 1.0;

  constexpr auto third = 1.0 / 3.0;

  // assign column 0
  flowJacobian_(4, 0) = -(eos.Conductivity(mu) +
                          eos.TurbConductivity(mut, turb->TurbPrandtlNumber()))
      * state.Temperature(eos) / ((mu + mut) * state.Rho());

  // assign column 1
  flowJacobian_(1, 1) = third * area.UnitVector().X() * area.UnitVector().X()
      + 1.0;
  flowJacobian_(2, 1) = third * area.UnitVector().X() * area.UnitVector().Y();
  flowJacobian_(3, 1) = third * area.UnitVector().X() * area.UnitVector().Z();
  flowJacobian_(4, 1) = fac * 0.5 * dist / (mu + mut) * tauNorm.X() +
      third * area.UnitVector().X() * velNorm + state.U();

  // assign column 2
  flowJacobian_(1, 2) = third * area.UnitVector().Y() * area.UnitVector().X();
  flowJacobian_(2, 2) = third * area.UnitVector().Y() * area.UnitVector().Y()
      + 1.0;
  flowJacobian_(3, 2) = third * area.UnitVector().Y() * area.UnitVector().Z();
  flowJacobian_(4, 2) = fac * 0.5 * dist / (mu + mut) * tauNorm.Y() +
      third * area.UnitVector().Y() * velNorm + state.V();

  // assign column 3
  flowJacobian_(1, 3) = third * area.UnitVector().Z() * area.UnitVector().X();
  flowJacobian_(2, 3) = third * area.UnitVector().Z() * area.UnitVector().Y();
  flowJacobian_(3, 3) = third * area.UnitVector().Z() * area.UnitVector().Z()
      + 1.0;
  flowJacobian_(4, 3) = fac * 0.5 * dist / (mu + mut) * tauNorm.Z() +
      third * area.UnitVector().Z() * velNorm + state.W();

  // assign column 4
  flowJacobian_(4, 4) = (eos.Conductivity(mu) +
                         eos.TurbConductivity(mut, turb->TurbPrandtlNumber()))
      / ((mu + mut) * state.Rho());

  // DEBUG
  flowJacobian_ *= area.Mag() * (mu + mut) / dist;
  // squareMatrix diag(5);
  // diag.Identity();
  // diag *= 1.0e3;
  // flowJacobian_ *= diag;

  auto specRad = std::max(4.0 / (3.0 * state.Rho()), eos.Gamma() / state.Rho());
  specRad *= (mu / eos.Prandtl() + mut / turb->TurbPrandtlNumber());
  specRad *= area.Mag() / dist;
  // flowJacobian_(0, 0) = 0.0;
  // fac = left ? 1.0 : 1.0;
  
  fluxJacobian prim2Cons;
  prim2Cons.DelPrimativeDelConservative(state, eos, inp);
  flowJacobian_ = fac * flowJacobian_.MatMult(prim2Cons.flowJacobian_);
  
  // flowJacobian_(0, 0) = fac * specRad;
  // flowJacobian_ *= 0.5;
  // flowJacobian_.Identity();
  // flowJacobian_ *= fac * specRad;

  // cout << "specRad: " << specRad << endl;
  // cout << "TSL Jacobian: " << flowJacobian_ << endl << endl;
  
  // calculate turbulent jacobian if necessary
  if (inp.IsTurbulent()) {
    turbJacobian_ = turb->ViscousJacobian(state, area, mu, suth, dist, mut, f1);
    // Don't need to multiply by prim2Cons b/c jacobian is already wrt
    // conservative variables
    turbJacobian_ *= fac;
  }
}


// non-member functions
// ----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}

genArray RusanovScalarOffDiagonal(const primVars &state,
                                  const genArray &update,
                                  const unitVec3dMag<double> &fArea,
                                  const double &mu,
                                  const double &mut, const double &f1,
                                  const double &dist,
                                  const idealGas &eos, const sutherland &suth,
                                  const unique_ptr<turbModel> &turb,
                                  const bool &isViscous, const bool &positive) {
  // state -- primative variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate = state.UpdateWithConsVars(eos, update, turb);

  // calculate updated convective flux
  // DEBUG -- removed const
  auto fluxChange = 0.5 * fArea.Mag() *
      ConvectiveFluxUpdate(state, stateUpdate, eos, fArea.UnitVector());
  // DEBUG
  fluxChange[5] = 0.0;
  fluxChange[6] = 0.0;

  // can't use stored cell spectral radius b/c it has contributions from i, j, k
  const uncoupledScalar specRad(state.FaceSpectralRadius(fArea, eos,
                                                         suth, dist, mu, mut,
                                                         turb, isViscous),
                                turb->FaceSpectralRadius(state, fArea, mu, suth,
                                                         dist, mut, f1,
                                                         positive));
  // DEBUG
  return positive ?
    fluxChange + specRad.ArrayMult(update) :
    fluxChange - specRad.ArrayMult(update);
}

genArray RusanovBlockOffDiagonal(const primVars &state,
                                 const genArray &update,
                                 const unitVec3dMag<double> &fArea,
                                 const double &mu,
                                 const double &mut, const double &f1,
                                 const double &dist,
                                 const idealGas &eos, const sutherland &suth,
                                 const unique_ptr<turbModel> &turb,
                                 const input &inp, const bool &positive,
                                 const tensor<double> &vGrad) {
  // state -- primative variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation
  // vGrad -- velocity gradient

  fluxJacobian jacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // calculate inviscid jacobian
  jacobian.RusanovFluxJacobian(state, eos, fArea, positive, inp, turb);

  // add viscous contribution
  fluxJacobian viscJac(inp.NumFlowEquations(), inp.NumTurbEquations());
  if (inp.IsViscous()) {
    viscJac.ApproxTSLJacobian(state, mu, mut, f1, eos, suth, fArea, dist, turb,
                              inp, positive, vGrad);
  }
  jacobian -= viscJac;

  return jacobian.ArrayMult(update);
}


genArray OffDiagonal(const primVars &offDiag, const primVars &diag,
                     const genArray &update,
                     const unitVec3dMag<double> &fArea, const double &mu,
                     const double &mut, const double &f1,
                     const double &dist, const tensor<double> &vGrad,
                     const idealGas &eos, const sutherland &suth,
                     const unique_ptr<turbModel> &turb,
                     const input &inp, const bool &positive) {
  // offDiag -- primative variables at off diagonal
  // diag -- primative variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // vGrad -- velocity gradient
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // input -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  genArray offDiagonal(0.0);

  if (inp.InvFluxJac() == "rusanov") {
    if (inp.IsBlockMatrix()) {
      offDiagonal = RusanovBlockOffDiagonal(offDiag, update, fArea, mu, mut, f1,
                                            dist, eos, suth, turb, inp,
                                            positive, vGrad);
    } else {
      offDiagonal = RusanovScalarOffDiagonal(offDiag, update, fArea, mu, mut,
                                             f1, dist, eos, suth, turb,
                                             inp.IsViscous(), positive);
    }
  } else if (inp.InvFluxJac() == "approximateRoe") {
    // always use block off diagonal for Roe method
    offDiagonal = RoeOffDiagonal(offDiag, diag, update, fArea, mu, mut, f1,
                                 dist, eos, suth, turb, inp.IsViscous(),
                                 inp.IsTurbulent(), positive);
  } else {
    cerr << "ERROR: Error in OffDiagonal(), inviscid flux jacobian method of "
         << inp.InvFluxJac() << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  return offDiagonal;
}


genArray RoeOffDiagonal(const primVars &offDiag, const primVars &diag,
                        const genArray &update,
                        const unitVec3dMag<double> &fArea,
                        const double &mu, const double &mut,
                        const double &dist, const double &f1,
                        const idealGas &eos,
                        const sutherland &suth,
                        const unique_ptr<turbModel> &turb,
                        const bool &isViscous, const bool &isTurbulent,
                        const bool &positive) {
  // offDiag -- primative variables at off diagonal
  // diag -- primative variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity at off diagonal
  // mut -- turbulent viscosity at off diagonal
  // f1 -- first blending coefficient at off diagonal
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // isTurbulent -- flag to determine if simulation is turbulent
  // positive -- flag to determine whether to add or subtract dissipation

  // DEBUG -- if positive, mu, mut, f1 should be from "left" side
  // REDO -- redo this whole function to not use the flux change and instead
  // calculate the matrix jacobian and multiply with the update

  
  const auto areaNorm = fArea.UnitVector();

  // calculate Roe flux with old variables
  const auto oldFlux = RoeFlux(offDiag, diag, eos, areaNorm);

  // calculate updated Roe flux on off diagonal
  const auto stateUpdate =  offDiag.UpdateWithConsVars(eos, update, turb);

  const auto newFlux = positive ?
    RoeFlux(stateUpdate, diag, eos, areaNorm) :
    RoeFlux(offDiag, stateUpdate, eos, areaNorm);

  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  const auto fluxChange = fArea.Mag() * (newFlux - oldFlux).ConvertToGenArray();
  
  // add contribution for viscous terms
  uncoupledScalar specRad(0.0, 0.0);
  if (isViscous) {
    specRad.AddToFlowVariable(
        offDiag.ViscFaceSpectralRadius(fArea, eos, suth, dist, mu,
                                        mut, turb));

    if (isTurbulent) {
      specRad.AddToTurbVariable(turb->ViscFaceSpecRad(offDiag, fArea, mu, suth,
                                                      dist, mut, f1));
    }
  }

  return positive ? fluxChange + specRad.ArrayMult(update) :
      fluxChange - specRad.ArrayMult(update);
}

void fluxJacobian::MultiplyOnDiagonal(const double &val,
                                      const bool &isTurbulent) {
  // val -- value to multiply along diagonal
  // isTurbulent -- flag identifiying if simulation is turbulent

  for (auto ii = 0; ii < flowJacobian_.Size(); ii++) {
    flowJacobian_(ii, ii) *= val;
  }

  if (isTurbulent) {
    for (auto ii = 0; ii < turbJacobian_.Size(); ii++) {
      turbJacobian_(ii, ii) *= val;
    }
  }
}

void fluxJacobian::AddOnDiagonal(const double &val, const bool &isTurbulent) {
  // val -- value to multiply along diagonal
  // isTurbulent -- flag identifiying if simulation is turbulent

  for (auto ii = 0; ii < flowJacobian_.Size(); ii++) {
    flowJacobian_(ii, ii) += val;
  }

  if (isTurbulent) {
    for (auto ii = 0; ii < turbJacobian_.Size(); ii++) {
      turbJacobian_(ii, ii) += val;
    }
  }
}

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
#include "transport.hpp"     // transport model
#include "eos.hpp"           // equation of state

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
    arr = turbJacobian_.ArrayMult(arr, flowJacobian_.Size());
  }
  return arr;
}

bool fluxJacobian::IsScalar() const {
  return (flowJacobian_.Size() > 1) ? false : true;
}

// function to take the inverse of a flux jacobian
void fluxJacobian::Inverse(const bool &isRANS) {
  flowJacobian_.Inverse();

  if (isRANS) {
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
                                       const unique_ptr<eos> &eqnState,
                                       const unitVec3dMag<double> &area,
                                       const bool &positive,
                                       const input &inp,
                                       const unique_ptr<turbModel> &turb) {
  // state -- primative variables at face
  // eqnState -- equation of state
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model

  // face inviscid spectral radius
  const auto specRad = state.InvFaceSpectralRadius(area, eqnState);

  // form dissipation matrix based on spectral radius
  fluxJacobian dissipation(inp.NumFlowEquations(), inp.NumTurbEquations());
  dissipation.flowJacobian_.Identity();
  dissipation.flowJacobian_ *= specRad;

  // begin jacobian calculation
  this->InvFluxJacobian(state, eqnState, area, inp, turb);

  // compute turbulent dissipation if necessary
  if (inp.IsRANS()) {
    // multiply by 0.5 b/c averaging with convection matrix
    dissipation.turbJacobian_ = 0.5 * turb->InviscidDissJacobian(state, area);
  }

  positive ? (*this) += dissipation : (*this) -= dissipation;
}

// function to calculate inviscid flux jacobian
void fluxJacobian::InvFluxJacobian(const primVars &state,
                                   const unique_ptr<eos> &eqnState,
                                   const unitVec3dMag<double> &area,
                                   const input &inp,
                                   const unique_ptr<turbModel> &turb) {
  // state -- primative variables at face
  // eqnState -- ideal gas equation of state
  // area -- face area vector
  // inp -- input variables
  // turb -- turbulence model

  const auto velNorm = state.Velocity().DotProd(area.UnitVector());
  const auto gammaMinusOne = eqnState->Gamma() - 1.0;
  const auto phi = 0.5 * gammaMinusOne * state.Velocity().MagSq();
  const auto a1 = eqnState->Gamma() * state.Energy(eqnState) - phi;
  const auto a3 = eqnState->Gamma() - 2.0;

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
  flowJacobian_(4, 4) = eqnState->Gamma() * velNorm;

  // multiply by 0.5 b/c averaging with dissipation matrix
  flowJacobian_ *= 0.5 * area.Mag();

  // turbulent jacobian here
  if (inp.IsRANS()) {
    // multiply by 0.5 b/c averaging with dissipation matrix
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
                                         const unique_ptr<eos> &eqnState,
                                         const unitVec3dMag<double> &area,
                                         const bool &positive,
                                         const input &inp,
                                         const unique_ptr<turbModel> &turb) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eqnState -- ideal gas equation of state
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model
  
  // compute Roe averaged state
  const auto roeAvg = RoeAveragedState(left, right);

  // compute Roe matrix
  fluxJacobian roeMatrix;
  roeMatrix.InvFluxJacobian(roeAvg, eqnState, area, inp, turb);

  // compute convective flux jacobian
  positive ? this->InvFluxJacobian(left, eqnState, area, inp, turb) :
      this->InvFluxJacobian(right, eqnState, area, inp, turb);

  positive ? (*this) += roeMatrix : (*this) -= roeMatrix;
}

// change of variable matrix going from primative to conservative variables
// from Dwight
void fluxJacobian::DelPrimativeDelConservative(const primVars &state,
                                               const unique_ptr<eos> &eqnState,
                                               const input &inp) {
  // state -- primative variables
  // eqnState -- equation of state
  // inp -- input variables

  const auto gammaMinusOne = eqnState->Gamma() - 1.0;
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
  if (inp.IsRANS()) {
    turbJacobian_(0, 0) = invRho;
    turbJacobian_(1, 1) = invRho;
  }
}

// approximate thin shear layer jacobian following implementation in Dwight.
// does not use any gradients
void fluxJacobian::ApproxTSLJacobian(const primVars &state,
                                     const double &lamVisc,
                                     const double &turbVisc, const double &f1,
                                     const unique_ptr<eos> &eqnState,
                                     const unique_ptr<transport> &trans,
                                     const unitVec3dMag<double> &area,
                                     const double &dist,
                                     const unique_ptr<turbModel> &turb,
                                     const input &inp, const bool &left,
                                     const tensor<double> &vGrad) {
  // state -- primative variables
  // eos -- equation of state
  // trans -- viscous transport model
  // area -- face area vector
  // dist -- distance from cell center to cell center
  // turb --  turbulence model
  // inp -- input variables
  // left -- flag that is negative if using left state
  // vGrad -- velocity gradient

  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  const auto mu = trans->NondimScaling() * lamVisc;
  const auto mut = trans->NondimScaling() * turbVisc;
  const auto velNorm = state.Velocity().DotProd(area.UnitVector());

  const auto tauNorm = TauNormal(vGrad, area.UnitVector(), mu, mut, trans);

  auto fac = left ? -1.0 : 1.0;

  constexpr auto third = 1.0 / 3.0;

  // assign column 0
  flowJacobian_(4, 0) =
      -(eqnState->Conductivity(mu) +
        eqnState->TurbConductivity(mut, turb->TurbPrandtlNumber())) *
      state.Temperature(eqnState) / ((mu + mut) * state.Rho());

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
  flowJacobian_(4, 4) =
      (eqnState->Conductivity(mu) +
       eqnState->TurbConductivity(mut, turb->TurbPrandtlNumber())) /
      ((mu + mut) * state.Rho());

  flowJacobian_ *= area.Mag() * (mu + mut) / dist;

  fluxJacobian prim2Cons;
  prim2Cons.DelPrimativeDelConservative(state, eqnState, inp);
  flowJacobian_ = flowJacobian_.MatMult(prim2Cons.flowJacobian_);

  // calculate turbulent jacobian if necessary
  if (inp.IsRANS()) {
    turbJacobian_ = fac * turb->ViscousJacobian(state, area, lamVisc, trans,
                                                dist, turbVisc, f1);
    // Don't need to multiply by prim2Cons b/c jacobian is already wrt
    // conservative variables
  }
}


// non-member functions
// ----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}

genArray RusanovScalarOffDiagonal(const primVars &state, const genArray &update,
                                  const unitVec3dMag<double> &fArea,
                                  const double &mu, const double &mut,
                                  const double &f1, const double &dist,
                                  const unique_ptr<eos> &eqnState,
                                  const unique_ptr<transport> &trans,
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
  // trans -- viscous transport model
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate = state.UpdateWithConsVars(eqnState, update, turb);

  // calculate updated convective flux
  auto fluxChange = 0.5 * fArea.Mag() *
      ConvectiveFluxUpdate(state, stateUpdate, eqnState, fArea.UnitVector());
  // zero out turbulence quantities b/c spectral radius is like full jacobian
  fluxChange[5] = 0.0;
  fluxChange[6] = 0.0;

  // can't use stored cell spectral radius b/c it has contributions from i, j, k
  const uncoupledScalar specRad(
      state.FaceSpectralRadius(fArea, eqnState, trans, dist, mu, mut, turb,
                               isViscous),
      turb->FaceSpectralRadius(state, fArea, mu, trans, dist, mut, f1,
                               positive));

  return positive ?
    fluxChange + specRad.ArrayMult(update) :
    fluxChange - specRad.ArrayMult(update);
}

genArray RusanovBlockOffDiagonal(
    const primVars &state, const genArray &update,
    const unitVec3dMag<double> &fArea, const double &mu, const double &mut,
    const double &f1, const double &dist, const unique_ptr<eos> &eqnState,
    const unique_ptr<transport> &trans, const unique_ptr<turbModel> &turb,
    const input &inp, const bool &positive, const tensor<double> &vGrad) {
  // state -- primative variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // eos -- equation of state
  // trans -- viscous transport model
  // turb -- turbulence model
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation
  // vGrad -- velocity gradient

  fluxJacobian jacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // calculate inviscid jacobian
  jacobian.RusanovFluxJacobian(state, eqnState, fArea, positive, inp, turb);

  // add viscous contribution
  if (inp.IsViscous()) {
    fluxJacobian viscJac(inp.NumFlowEquations(), inp.NumTurbEquations());
    viscJac.ApproxTSLJacobian(state, mu, mut, f1, eqnState, trans, fArea, dist,
                              turb, inp, positive, vGrad);
    positive ? jacobian -= viscJac : jacobian += viscJac;
  }
  return jacobian.ArrayMult(update);
}

genArray OffDiagonal(const primVars &offDiag, const primVars &diag,
                     const genArray &update, const unitVec3dMag<double> &fArea,
                     const double &mu, const double &mut, const double &f1,
                     const double &dist, const tensor<double> &vGrad,
                     const unique_ptr<eos> &eqnState,
                     const unique_ptr<transport> &trans,
                     const unique_ptr<turbModel> &turb, const input &inp,
                     const bool &positive) {
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
  // trans -- viscous transport model
  // turb -- turbulence model
  // input -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  genArray offDiagonal(0.0);

  if (inp.InvFluxJac() == "rusanov") {
    if (inp.IsBlockMatrix()) {
      offDiagonal = RusanovBlockOffDiagonal(offDiag, update, fArea, mu, mut, f1,
                                            dist, eqnState, trans, turb, inp,
                                            positive, vGrad);
    } else {
      offDiagonal = RusanovScalarOffDiagonal(offDiag, update, fArea, mu, mut,
                                             f1, dist, eqnState, trans, turb,
                                             inp.IsViscous(), positive);
    }
  } else if (inp.InvFluxJac() == "approximateRoe") {
    // always use flux change off diagonal with roe method
    offDiagonal = RoeOffDiagonal(offDiag, diag, update, fArea, mu, mut,
                                 f1, dist, eqnState, trans, turb,
                                 inp.IsViscous(), inp.IsRANS(),
                                 positive);
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
                        const unique_ptr<eos> &eqnState,
                        const unique_ptr<transport> &trans,
                        const unique_ptr<turbModel> &turb,
                        const bool &isViscous, const bool &isRANS,
                        const bool &positive) {
  // offDiag -- primative variables at off diagonal
  // diag -- primative variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity at off diagonal
  // mut -- turbulent viscosity at off diagonal
  // f1 -- first blending coefficient at off diagonal
  // eqnState -- equation of state
  // trans -- viscous transport model
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // isRANS -- flag to determine if simulation is turbulent
  // positive -- flag to determine whether to add or subtract dissipation

  // DEBUG -- redo this whole function to not use the flux change and instead
  // calculate the matrix jacobian and multiply with the update
  
  const auto areaNorm = fArea.UnitVector();

  // calculate Roe flux with old variables
  const auto oldFlux = RoeFlux(offDiag, diag, eqnState, areaNorm);

  // calculate updated Roe flux on off diagonal
  const auto stateUpdate =  offDiag.UpdateWithConsVars(eqnState, update, turb);

  const auto newFlux = positive ?
    RoeFlux(stateUpdate, diag, eqnState, areaNorm) :
    RoeFlux(diag, stateUpdate, eqnState, areaNorm);

  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  const auto fluxChange = fArea.Mag() * (newFlux - oldFlux).ConvertToGenArray();
  
  // add contribution for viscous terms
  uncoupledScalar specRad(0.0, 0.0);
  if (isViscous) {
    specRad.AddToFlowVariable(
        offDiag.ViscFaceSpectralRadius(fArea, eqnState, trans, dist, mu,
                                        mut, turb));

    if (isRANS) {
      specRad.AddToTurbVariable(turb->ViscFaceSpecRad(offDiag, fArea, mu, trans,
                                                      dist, mut, f1));
    }
  }

  return positive ? fluxChange + specRad.ArrayMult(update) :
      fluxChange - specRad.ArrayMult(update);
}

void fluxJacobian::MultiplyOnDiagonal(const double &val,
                                      const bool &isRANS) {
  // val -- value to multiply along diagonal
  // isRANS -- flag identifiying if simulation is turbulent

  for (auto ii = 0; ii < flowJacobian_.Size(); ii++) {
    flowJacobian_(ii, ii) *= val;
  }

  if (isRANS) {
    for (auto ii = 0; ii < turbJacobian_.Size(); ii++) {
      turbJacobian_(ii, ii) *= val;
    }
  }
}

void fluxJacobian::AddOnDiagonal(const double &val, const bool &isRANS) {
  // val -- value to multiply along diagonal
  // isRANS -- flag identifiying if simulation is turbulent

  for (auto ii = 0; ii < flowJacobian_.Size(); ii++) {
    flowJacobian_(ii, ii) += val;
  }

  if (isRANS) {
    for (auto ii = 0; ii < turbJacobian_.Size(); ii++) {
      turbJacobian_(ii, ii) += val;
    }
  }
}

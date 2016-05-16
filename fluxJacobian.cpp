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
#include "inviscidFlux.hpp"  // ConvectiveFluxUpdate

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std:: unique_ptr;

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
void fluxJacobian::RusanovFluxJacobian(const primVars &left,
                                       const primVars &right,
                                       const idealGas &eos,
                                       const vector3d<double> &areaNorm,
                                       const bool &positive,
                                       const input &inp) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eos -- ideal gas equation of state
  // areaNorm -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables

  // dot product of velocities with unit area vector
  const auto specRad = std::max(fabs(left.Velocity().DotProd(areaNorm)) +
                                left.SoS(eos),
                                fabs(right.Velocity().DotProd(areaNorm)) +
                                right.SoS(eos));

  // form dissipation matrix based on spectral radius
  fluxJacobian dissipation(inp.NumFlowEquations(), inp.NumTurbEquations());
  dissipation.flowJacobian_.Identity();
  dissipation.flowJacobian_ *= specRad;


  // begin jacobian calculation
  positive ? this->InvFluxJacobian(left, eos, areaNorm, inp) :
      this->InvFluxJacobian(right, eos, areaNorm, inp);


  // compute turbulent jacobian if necessary
  if (inp.IsTurbulent()) {
    dissipation.turbJacobian_.Identity();

    const auto turbSpecRad = std::max(fabs(left.Velocity().DotProd(areaNorm)),
                                      fabs(right.Velocity().DotProd(areaNorm)));

    dissipation.turbJacobian_ *= turbSpecRad;
  }

  positive ? (*this) += dissipation : (*this) -= dissipation;
  (*this) *= 0.5;
}

// function to calculate inviscid flux jacobian
void fluxJacobian::InvFluxJacobian(const primVars &state,
                                   const idealGas &eqnState,
                                   const vector3d<double> &areaNorm,
                                   const input &inp) {
  // state -- primative variables from left side
  // eqnState -- ideal gas equation of state
  // areaNorm -- face area vector
  // inp -- input variables

  const auto velNorm = state.Velocity().DotProd(areaNorm);
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
  flowJacobian_(1, 0) = phi * areaNorm.X() - state.U() * velNorm;
  flowJacobian_(2, 0) = phi * areaNorm.Y() - state.V() * velNorm;
  flowJacobian_(3, 0) = phi * areaNorm.Z() - state.W() * velNorm;
  flowJacobian_(4, 0) = velNorm * (phi - a1);

  // column one
  flowJacobian_(0, 1) = areaNorm.X();
  flowJacobian_(1, 1) = velNorm - a3 * areaNorm.X() * state.U();
  flowJacobian_(2, 1) = state.V() * areaNorm.X() -
      gammaMinusOne * state.U() * areaNorm.Y();
  flowJacobian_(3, 1) = state.W() * areaNorm.X() -
      gammaMinusOne * state.U() * areaNorm.Z();
  flowJacobian_(4, 1) = a1 * areaNorm.X() - gammaMinusOne * state.U() * velNorm;

  // column two
  flowJacobian_(0, 2) = areaNorm.Y();
  flowJacobian_(1, 2) = state.U() * areaNorm.Y() -
      gammaMinusOne * state.V() * areaNorm.X();
  flowJacobian_(2, 2) = velNorm - a3 * areaNorm.Y() * state.V();
  flowJacobian_(3, 2) = state.W() * areaNorm.Y() -
      gammaMinusOne * state.V() * areaNorm.Z();
  flowJacobian_(4, 2) = a1 * areaNorm.Y() - gammaMinusOne * state.V() * velNorm;

  // column three
  flowJacobian_(0, 3) = areaNorm.Z();
  flowJacobian_(1, 3) = state.U() * areaNorm.Z() -
      gammaMinusOne * state.W() * areaNorm.X();
  flowJacobian_(2, 3) = state.V() * areaNorm.Z() -
      gammaMinusOne * state.W() * areaNorm.Y();
  flowJacobian_(3, 3) = velNorm - a3 * areaNorm.Z() * state.W();
  flowJacobian_(4, 3) = a1 * areaNorm.Z() - gammaMinusOne * state.W() * velNorm;

  // column four
  flowJacobian_(0, 4) = 0.0;
  flowJacobian_(1, 4) = gammaMinusOne * areaNorm.X();
  flowJacobian_(2, 4) = gammaMinusOne * areaNorm.Y();
  flowJacobian_(3, 4) = gammaMinusOne * areaNorm.Z();
  flowJacobian_(4, 4) = eqnState.Gamma() * velNorm;

  // turbulent jacobian here
  if (inp.IsTurbulent()) {
    turbJacobian_(0, 0) = velNorm;
    turbJacobian_(1, 1) = velNorm;
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
                                         const vector3d<double> &areaNorm,
                                         const bool &positive,
                                         const input &inp) {
  // left -- primative variables from left side
  // right -- primative variables from right side
  // eos -- ideal gas equation of state
  // areaNorm -- face unit area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables

  // compute Roe averaged state
  const auto roeAvg = RoeAveragedState(left, right, eos);

  // compute Roe matrix
  fluxJacobian roeMatrix;
  roeMatrix.InvFluxJacobian(roeAvg, eos, areaNorm, inp);

  // compute convective flux jacobian
  positive ? this->InvFluxJacobian(left, eos, areaNorm, inp) :
      this->InvFluxJacobian(right, eos, areaNorm, inp);

  positive ? (*this) += roeMatrix : (*this) -= roeMatrix;
  (*this) *= 0.5;
}

// change of variable matrix going frim primative to conservative variables
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

  // assign first column
  flowJacobian_(0, 0) = 1.0;
  flowJacobian_(1, 0) = -invRho * state.U();
  flowJacobian_(2, 0) = -invRho * state.V();
  flowJacobian_(3, 0) = -invRho * state.W();
  flowJacobian_(4, 0) = 0.5 * gammaMinusOne *
      state.Velocity().DotProd(state.Velocity());

  // assign second column
  flowJacobian_(1, 1) = invRho;
  flowJacobian_(4, 1) = -gammaMinusOne * state.U();

  // assign third column
  flowJacobian_(2, 2) = invRho;
  flowJacobian_(4, 2) = -gammaMinusOne * state.V();

  // assign fourth column
  flowJacobian_(3, 3) = invRho;
  flowJacobian_(4, 3) = -gammaMinusOne * state.W();

  // assign fifth column
  flowJacobian_(4, 4) = gammaMinusOne;

  // turbulent jacobian here
  if (inp.IsTurbulent()) {
    turbJacobian_(0, 0) = invRho;
    turbJacobian_(1, 1) = invRho;
  }
}

// approximate thin shear layer jacobian following implementation in Dwight.
// does not use any gradients
void fluxJacobian::ApproxTSLJacobian(const primVars &state, const idealGas &eos,
                                     const sutherland &suth,
                                     const vector3d<double> &area,
                                     const double &dist,
                                     const unique_ptr<turbModel> &turb,
                                     const input &inp) {
  // state -- primative variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // area -- face area unit vector
  // dist -- distance from cell center to cell center
  // turb --  turbulence model
  // inp -- input variables

  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  const auto mu = suth.EffectiveViscosity(state.Temperature(eos));
  const auto mut = turb->EddyViscNoLim(state) * suth.NondimScaling();
  const auto velNorm = state.Velocity().DotProd(area);

  // assign first column
  flowJacobian_(4, 0) = -eos.Conductivity(mu + mut) * state.Temperature(eos) /
    ((mu + mut) * state.Rho());

  // assign second column
  flowJacobian_(1, 1) = (1.0 / 3.0) * area.X() * area.X() + 1.0;
  flowJacobian_(1, 2) = (1.0 / 3.0) * area.X() * area.Y();
  flowJacobian_(1, 3) = (1.0 / 3.0) * area.X() * area.Z();
  flowJacobian_(1, 4) = (1.0 / 3.0) * area.X() * velNorm + state.U();

  // assign third column
  flowJacobian_(2, 1) = (1.0 / 3.0) * area.Y() * area.X();
  flowJacobian_(2, 2) = (1.0 / 3.0) * area.Y() * area.Y() + 1.0;
  flowJacobian_(2, 3) = (1.0 / 3.0) * area.Y() * area.Z();
  flowJacobian_(2, 4) = (1.0 / 3.0) * area.Y() * velNorm + state.V();

  // assign fourth column
  flowJacobian_(3, 1) = (1.0 / 3.0) * area.Z() * area.X();
  flowJacobian_(3, 2) = (1.0 / 3.0) * area.Z() * area.Y();
  flowJacobian_(3, 3) = (1.0 / 3.0) * area.Z() * area.Z() + 1.0;
  flowJacobian_(3, 4) = (1.0 / 3.0) * area.Z() * velNorm + state.W();

  // assign fifth column
  flowJacobian_(4, 4) = eos.Conductivity(mu + mut) / ((mu + mut) * state.Rho());

  flowJacobian_ *= (mu + mut) / dist;

  fluxJacobian prim2Cons;
  prim2Cons.DelPrimativeDelConservative(state, eos, inp);

  flowJacobian_ = flowJacobian_.MatMult(prim2Cons.flowJacobian_);

  // calculate turbulent jacobian if necessary
  if (inp.IsTurbulent()) {
    turbJacobian_(0, 0) = mu + turb->SigmaK() * mut;
    turbJacobian_(1, 1) = mu + turb->SigmaW() * mut;
    turbJacobian_ *= 1.0 / dist;

    turbJacobian_ = turbJacobian_.MatMult(prim2Cons.turbJacobian_);
  }
}


// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}


genArray RusanovOffDiagonal(const primVars &state, const genArray &update,
                            const unitVec3dMag<double> &fAreaL,
                            const unitVec3dMag<double> &fAreaR,
                            const double &vol, const double &mu,
                            const idealGas &eos, const sutherland &suth,
                            const unique_ptr<turbModel> &turb,
                            const bool &isViscous, const bool &positive) {
  // state -- primative variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fAreaL -- face area vector on off diagonal boundary
  // fAreaR -- face area vector opposite off diagonal boundary
  // vol -- cell volume
  // mu -- laminar viscosity
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate = state.UpdateWithConsVars(eos, update, turb);

  // calculate updated convective flux
  const auto fluxChange = ConvectiveFluxUpdate(state, stateUpdate, eos,
                                               fAreaL.UnitVector());

  // can't use stored cell spectral radius b/c it has contribuitons from i, j, k
  const uncoupledScalar specRad(state.CellSpectralRadius(fAreaL, fAreaR, eos,
                                                         suth, vol, mu, turb,
                                                         isViscous),
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
                        const double &vol, const double &muL,
                        const double &muR, const idealGas &eos,
                        const sutherland &suth,
                        const unique_ptr<turbModel> &turb,
                        const bool &isViscous, const bool &isTurbulent,
                        const bool &positive) {
  // left -- primative variables at left side
  // right -- primative variables at right side
  // update -- conserved variable update at off diagonal
  // fAreaL -- face area vector on off diagonal boundary
  // fAreaR -- face area vector opposite off diagonal boundary
  // vol -- cell volume
  // muL -- laminar viscosity at left side
  // muR -- laminar viscosity at right side
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // isTurbulent -- flag to determine if simulation is turbulent
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
  if (isViscous) {
    const auto & offState = positive ? left : right;
    const auto & offMu = positive ? muL : muR;

    specRad.AddToFlowVariable(
        offState.ViscCellSpectralRadius(fAreaL, fAreaR, eos, suth, vol, offMu,
                                        turb));

    if (isTurbulent) {
      specRad.AddToTurbVariable(turb->ViscSpecRad(offState, fAreaL, fAreaR, eos,
                                                  suth, vol));
    }
  }

  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  return positive ?
    fAreaL.Mag() * ((newFlux - oldFlux).ConvertToGenArray()) +
      0.5 * specRad.ArrayMult(update) :
    fAreaL.Mag() * ((newFlux - oldFlux).ConvertToGenArray()) -
      0.5 * specRad.ArrayMult(update);
}


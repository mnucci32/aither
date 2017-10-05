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
#include "turbulence.hpp"     // turbModel
#include "input.hpp"          // input
#include "primitive.hpp"      // primitive
#include "inviscidFlux.hpp"   // ConvectiveFluxUpdate
#include "utility.hpp"        // TauNormal
#include "transport.hpp"      // transport model
#include "eos.hpp"            // equation of state
#include "thermodynamic.hpp"  // thermodynamic model
#include "conserved.hpp"      // conserved
#include "spectralRadius.hpp"

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



// non-member functions
// ----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}

varArray RusanovScalarOffDiagonal(
    const primitiveView &state, const varArrayView &update,
    const unitVec3dMag<double> &fArea, const double &mu, const double &mut,
    const double &f1, const double &dist, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, const bool &isViscous,
    const bool &positive) {
  // state -- primitive variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // eos -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate =
      state.UpdateWithConsVars(eqnState, thermo, update, turb);

  // calculate updated convective flux
  auto fluxChange =
      0.5 * fArea.Mag() * ConvectiveFluxUpdate(state, stateUpdate, eqnState,
                                               thermo, fArea.UnitVector());
  // zero out turbulence quantities b/c spectral radius is like full jacobian
  for (auto ii = 0; ii < fluxChange.NumTurbulence(); ++ii) {
    fluxChange[fluxChange.TurbulenceIndex() + ii] = 0.0;
  }

  // can't use stored cell spectral radius b/c it has contributions from i, j, k
  const uncoupledScalar specRad(
      FaceSpectralRadius(state, fArea, thermo, eqnState, trans, dist, mu, mut,
                         turb, isViscous),
      turb->FaceSpectralRadius(state, fArea, mu, trans, dist, mut, f1,
                               positive));

  return positive ?
    fluxChange + specRad.ArrayMult(update) :
    fluxChange - specRad.ArrayMult(update);
}

varArray RusanovBlockOffDiagonal(
    const primitiveView &state, const varArrayView &update,
    const unitVec3dMag<double> &fArea, const double &mu, const double &mut,
    const double &f1, const double &dist, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, const input &inp, const bool &positive,
    const tensor<double> &vGrad) {
  // state -- primitive variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // turb -- turbulence model
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation
  // vGrad -- velocity gradient

  fluxJacobian jacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // calculate inviscid jacobian
  jacobian.RusanovFluxJacobian(state, eqnState, thermo, fArea, positive, inp,
                               turb);

  // add viscous contribution
  if (inp.IsViscous()) {
    fluxJacobian viscJac(inp.NumFlowEquations(), inp.NumTurbEquations());
    viscJac.ApproxTSLJacobian(state, mu, mut, f1, eqnState, trans, thermo,
                              fArea, dist, turb, inp, positive, vGrad);
    positive ? jacobian -= viscJac : jacobian += viscJac;
  }
  return jacobian.ArrayMult(update);
}

varArray OffDiagonal(
    const primitiveView &offDiag, const primitiveView &diag,
    const varArrayView &update, const unitVec3dMag<double> &fArea,
    const double &mu, const double &mut, const double &f1, const double &dist,
    const tensor<double> &vGrad, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, const input &inp, const bool &positive) {
  // offDiag -- primitive variables at off diagonal
  // diag -- primitive variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // vGrad -- velocity gradient
  // eos -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // turb -- turbulence model
  // input -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  varArray offDiagonal(update.Size(), update.NumSpecies());

  if (inp.InvFluxJac() == "rusanov") {
    if (inp.IsBlockMatrix()) {
      offDiagonal = RusanovBlockOffDiagonal(offDiag, update, fArea, mu, mut, f1,
                                            dist, eqnState, thermo, trans, turb,
                                            inp, positive, vGrad);
    } else {
      offDiagonal = RusanovScalarOffDiagonal(offDiag, update, fArea, mu, mut,
                                             f1, dist, eqnState, thermo, trans,
                                             turb, inp.IsViscous(), positive);
    }
  } else if (inp.InvFluxJac() == "approximateRoe") {
    // always use flux change off diagonal with roe method
    offDiagonal = RoeOffDiagonal(offDiag, diag, update, fArea, mu, mut,
                                 f1, dist, eqnState, thermo, trans, turb,
                                 inp.IsViscous(), inp.IsRANS(),
                                 positive);
  } else {
    cerr << "ERROR: Error in OffDiagonal(), inviscid flux jacobian method of "
         << inp.InvFluxJac() << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  return offDiagonal;
}

varArray RoeOffDiagonal(
    const primitiveView &offDiag, const primitiveView &diag,
    const varArrayView &update, const unitVec3dMag<double> &fArea,
    const double &mu, const double &mut, const double &dist, const double &f1,
    const unique_ptr<eos> &eqnState, const unique_ptr<thermodynamic> &thermo,
    const unique_ptr<transport> &trans, const unique_ptr<turbModel> &turb,
    const bool &isViscous, const bool &isRANS, const bool &positive) {
  // offDiag -- primitive variables at off diagonal
  // diag -- primitive variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity at off diagonal
  // mut -- turbulent viscosity at off diagonal
  // f1 -- first blending coefficient at off diagonal
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // turb -- turbulence model
  // isViscous -- flag to determine if simulation is viscous
  // isRANS -- flag to determine if simulation is turbulent
  // positive -- flag to determine whether to add or subtract dissipation

  // DEBUG -- redo this whole function to not use the flux change and instead
  // calculate the matrix jacobian and multiply with the update
  
  const auto areaNorm = fArea.UnitVector();

  // calculate Roe flux with old variables
  const auto oldFlux = RoeFlux(offDiag, diag, eqnState, thermo, areaNorm);

  // calculate updated Roe flux on off diagonal
  const auto stateUpdate =
      offDiag.UpdateWithConsVars(eqnState, thermo, update, turb);

  const auto newFlux = positive ?
    RoeFlux(stateUpdate, diag, eqnState, thermo, areaNorm) :
    RoeFlux(diag, stateUpdate, eqnState, thermo, areaNorm);

  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  const auto fluxChange = fArea.Mag() * (newFlux - oldFlux);
  
  // add contribution for viscous terms
  uncoupledScalar specRad(0.0, 0.0);
  if (isViscous) {
    specRad.AddToFlowVariable(ViscFaceSpectralRadius(
        offDiag, fArea, thermo, eqnState, trans, dist, mu, mut, turb));

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

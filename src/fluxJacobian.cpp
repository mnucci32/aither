/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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
#include "matrix.hpp"
#include "physicsModels.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

// member functions
void fluxJacobian::AddToFlowJacobian(const squareMatrix &jac) {
  MSG_ASSERT(flowSize_ == jac.Size(), "matrix sizes don't match");
  std::transform(this->begin(), this->beginTurb(), jac.begin(), this->begin(),
                 std::plus<double>());
}
void fluxJacobian::SubtractFromFlowJacobian(const squareMatrix &jac) { 
  MSG_ASSERT(flowSize_ == jac.Size(), "matrix sizes don't match");
  std::transform(this->begin(), this->beginTurb(), jac.begin(), this->begin(),
                 std::minus<double>());
}

void fluxJacobian::AddToTurbJacobian(const squareMatrix &jac) {
  MSG_ASSERT(turbSize_ == jac.Size(), "matrix sizes don't match");
  std::transform(this->beginTurb(), this->end(), jac.begin(), this->beginTurb(),
                 std::plus<double>());
}
void fluxJacobian::SubtractFromTurbJacobian(const squareMatrix &jac) {
  MSG_ASSERT(turbSize_ == jac.Size(), "matrix sizes don't match");
  std::transform(this->beginTurb(), this->end(), jac.begin(), this->beginTurb(),
                 std::minus<double>());
}

void fluxJacobian::MultFlowJacobian(const double &fac) {
  std::for_each(this->begin(), this->beginTurb(),
                [&fac](auto &val) { val *= fac; });
}
void fluxJacobian::MultTurbJacobian(const double &fac) {
  std::for_each(this->beginTurb(), this->end(),
                [&fac](auto &val) { val *= fac; });
}

fluxJacobian fluxJacobian::FlowMatMult(const fluxJacobian &jac2) const {
  MSG_ASSERT(this->FlowSize() == jac2.FlowSize(),
             "Mismatch in flux jacobian size");
  fluxJacobian result(this->FlowSize(), this->TurbSize());
  MatrixMultiply(this->begin(), jac2.begin(), result.begin(), this->FlowSize());
  return result;
}

fluxJacobian fluxJacobian::TurbMatMult(const fluxJacobian &jac2) const {
  MSG_ASSERT(this->TurbSize() == jac2.TurbSize(),
             "Mismatch in flux jacobian size");
  fluxJacobian result(this->FlowSize(), this->TurbSize());
  MatrixMultiply(this->beginTurb(), jac2.beginTurb(), result.beginTurb(),
                 this->TurbSize());
  return result;
}



// non-member functions
// ----------------------------------------------------------------------------
ostream &operator<<(ostream &os, const fluxJacobian &jacobian) {
  // print flow jacobian
  for (auto rr = 0; rr < jacobian.FlowSize(); ++rr) {
    for (auto cc = 0; cc < jacobian.FlowSize(); ++cc) {
      os << jacobian.FlowJacobian(rr, cc);
      if (cc != (jacobian.FlowSize() - 1)) {
        os << ", ";
      } else {
        os << endl;
      }
    }
  }

  // print turbulence jacobian
  for (auto rr = 0; rr < jacobian.TurbSize(); ++rr) {
    for (auto cc = 0; cc < jacobian.TurbSize(); ++cc) {
      os << jacobian.TurbJacobian(rr, cc);
      if (cc != (jacobian.TurbSize() - 1)) {
        os << ", ";
      } else {
        os << endl;
      }
    }
  }

  return os;
}

varArray RusanovScalarOffDiagonal(const primitiveView &state,
                                  const varArrayView &update,
                                  const unitVec3dMag<double> &fArea,
                                  const double &mu, const double &mut,
                                  const double &f1, const double &dist,
                                  const physics &phys, const bool &isViscous,
                                  const bool &positive) {
  // state -- primitive variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // phys -- physics models
  // isViscous -- flag to determine if simulation is viscous
  // positive -- flag to determine whether to add or subtract dissipation

  // calculate updated state
  const auto stateUpdate =
      state.UpdateWithConsVars(phys, update);

  // calculate updated convective flux
  auto fluxChange =
      0.5 * fArea.Mag() *
      ConvectiveFluxUpdate(state, stateUpdate, phys, fArea.UnitVector());
  // zero out turbulence quantities b/c spectral radius is like full jacobian
  for (auto ii = 0; ii < fluxChange.NumTurbulence(); ++ii) {
    fluxChange[fluxChange.TurbulenceIndex() + ii] = 0.0;
  }

  // can't use stored cell spectral radius b/c it has contributions from i, j, k
  const uncoupledScalar specRad(
      FaceSpectralRadius(state, fArea, phys, dist, mu, mut, isViscous),
      phys.Turbulence()->FaceSpectralRadius(state, fArea, mu, phys.Transport(),
                                            dist, mut, f1, positive));

  return positive ?
    fluxChange + specRad.ArrayMult(update) :
    fluxChange - specRad.ArrayMult(update);
}

varArray RusanovBlockOffDiagonal(
    const primitiveView &state, const varArrayView &update,
    const unitVec3dMag<double> &fArea, const double &mu, const double &mut,
    const double &f1, const double &dist, const physics &phys, const input &inp,
    const bool &positive, const tensor<double> &vGrad) {
  // state -- primitive variables at off diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // phys -- physics models
  // inp -- input variables
  // positive -- flag to determine whether to add or subtract dissipation
  // vGrad -- velocity gradient

  fluxJacobian jacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // calculate inviscid jacobian
  jacobian.RusanovFluxJacobian(state, phys, fArea, positive, inp);

  // add viscous contribution
  if (inp.IsViscous()) {
    fluxJacobian viscJac(inp.NumFlowEquations(), inp.NumTurbEquations());
    viscJac.ApproxTSLJacobian(state, mu, mut, f1, phys, fArea, dist, inp,
                              positive, vGrad);
    positive ? jacobian -= viscJac : jacobian += viscJac;
  }
  return jacobian.ArrayMult(update);
}

varArray OffDiagonal(const primitiveView &offDiag, const primitiveView &diag,
                     const varArrayView &update,
                     const unitVec3dMag<double> &fArea, const double &mu,
                     const double &mut, const double &f1, const double &dist,
                     const tensor<double> &vGrad, const physics &phys,
                     const input &inp, const bool &positive) {
  // offDiag -- primitive variables at off diagonal
  // diag -- primitive variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // dist -- distance from cell center to cell center across face on diagonal
  // vGrad -- velocity gradient
  // phys -- physics models
  // input -- input variables
  // positive -- flag to determine whether to add or subtract dissipation

  varArray offDiagonal;

  if (inp.InvFluxJac() == "rusanov") {
    if (inp.IsBlockMatrix()) {
      offDiagonal = RusanovBlockOffDiagonal(offDiag, update, fArea, mu, mut, f1,
                                            dist, phys, inp, positive, vGrad);
    } else {
      offDiagonal =
          RusanovScalarOffDiagonal(offDiag, update, fArea, mu, mut, f1, dist,
                                   phys, inp.IsViscous(), positive);
    }
  } else if (inp.InvFluxJac() == "approximateRoe") {
    // always use flux change off diagonal with roe method
    offDiagonal =
        RoeOffDiagonal(offDiag, diag, update, fArea, mu, mut, f1, dist, phys,
                       inp.IsViscous(), inp.IsRANS(), positive);
  } else {
    cerr << "ERROR: Error in OffDiagonal(), inviscid flux jacobian method of "
         << inp.InvFluxJac() << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  return offDiagonal;
}

varArray RoeOffDiagonal(const primitiveView &offDiag, const primitiveView &diag,
                        const varArrayView &update,
                        const unitVec3dMag<double> &fArea, const double &mu,
                        const double &mut, const double &dist, const double &f1,
                        const physics &phys, const bool &isViscous,
                        const bool &isRANS, const bool &positive) {
  // offDiag -- primitive variables at off diagonal
  // diag -- primitive variables at diagonal
  // update -- conserved variable update at off diagonal
  // fArea -- face area vector on off diagonal boundary
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity at off diagonal
  // mut -- turbulent viscosity at off diagonal
  // f1 -- first blending coefficient at off diagonal
  // phys -- physics models
  // isViscous -- flag to determine if simulation is viscous
  // isRANS -- flag to determine if simulation is turbulent
  // positive -- flag to determine whether to add or subtract dissipation

  // DEBUG -- redo this whole function to not use the flux change and instead
  // calculate the matrix jacobian and multiply with the update
  
  const auto areaNorm = fArea.UnitVector();

  // calculate Roe flux with old variables
  const auto oldFlux = RoeFlux(offDiag, diag, phys, areaNorm);

  // calculate updated Roe flux on off diagonal
  const auto stateUpdate = offDiag.UpdateWithConsVars(phys, update);

  const auto newFlux = positive ? RoeFlux(stateUpdate, diag, phys, areaNorm)
                                : RoeFlux(diag, stateUpdate, phys, areaNorm);

  // don't need 0.5 factor on roe flux because RoeFlux function already does it
  const auto fluxChange = fArea.Mag() * (newFlux - oldFlux);
  
  // add contribution for viscous terms
  uncoupledScalar specRad(0.0, 0.0);
  if (isViscous) {
    specRad.AddToFlowVariable(
        ViscFaceSpectralRadius(offDiag, fArea, phys, dist, mu, mut));

    if (isRANS) {
      specRad.AddToTurbVariable(phys.Turbulence()->ViscFaceSpecRad(
          offDiag, fArea, mu, phys.Transport(), dist, mut, f1));
    }
  }

  return positive ? fluxChange + specRad.ArrayMult(update) :
      fluxChange - specRad.ArrayMult(update);
}


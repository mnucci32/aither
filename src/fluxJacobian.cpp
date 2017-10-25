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
#include "matrix.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;


// member functions
// member function to invert matrix using Gauss-Jordan elimination
void fluxJacobian::FlowInverse() {
  squareMatrix I(flowSize_);
  I.Identity();

  for (auto cPivot = 0, r = 0; r < flowSize_; ++r, ++cPivot) {
    // find pivot row
    auto rPivot = this->FlowFindMaxInCol(r, cPivot, flowSize_ - 1);

    // swap rows
    this->FlowSwapRows(r, rPivot);
    I.SwapRows(r, rPivot);

    if (r != 0) {  // if not first row, need to get rid entries ahead of pivot
      for (auto ii = 0; ii < cPivot; ++ii) {
        auto factor = this->FlowJacobian(r, ii) / this->FlowJacobian(ii, ii);
        this->FlowLinCombRow(ii, factor, r);
        I.LinCombRow(ii, factor, r);
      }
    }

    // normalize row by pivot
    if (this->FlowJacobian(r, cPivot) == 0.0) {
      cerr << "ERROR: Singular flow matrix in Gauss-Jordan elimination! Matrix "
              "(mid inversion) is" << endl << *this << endl;
      exit(EXIT_FAILURE);
    }
    // only normalize entries from pivot and to the right
    auto normFactor = 1.0 / this->FlowJacobian(r, cPivot);
    this->FlowRowMultiply(r, cPivot, normFactor);
    I.RowMultiply(r, 0, normFactor);  // multiply all entries
  }

  // matrix is now upper triangular, work way back up to identity matrix
  // start with second to last row
  for (auto cPivot = flowSize_ - 2, r = flowSize_ - 2; r >= 0; --r, --cPivot) {
    for (auto ii = flowSize_ - 1; ii > cPivot; --ii) {
      auto factor = this->FlowJacobian(r, ii);
      this->FlowLinCombRow(ii, factor, r);
      I.LinCombRow(ii, factor, r);
    }
  }

  // set this matrix equal to its inverse
  std::copy(I.begin(), I.end(), this->begin());
}

void fluxJacobian::TurbInverse() {
  MSG_ASSERT(turbSize_ <= 2, "expecting turbulence size < 2");
  if (turbSize_ == 1) {  // scalar turbulence
    this->TurbJacobian(0, 0) = 1.0 / this->TurbJacobian(0, 0);
  } else if (turbSize_ == 2) {  // block jacobian
    squareMatrix inv(turbSize_);
    auto det = this->TurbJacobian(0, 0) * this->TurbJacobian(1, 1) -
               this->TurbJacobian(0, 1) * this->TurbJacobian(1, 0);
    inv(0, 0) = this->TurbJacobian(1, 1);
    inv(0, 1) = -this->TurbJacobian(0, 1);
    inv(0, 1) = -this->TurbJacobian(1, 0);
    inv(1, 1) = this->TurbJacobian(0, 0);
    inv *= 1.0 / det;
    std::copy(inv.begin(), inv.end(), this->beginTurb());
  }
}


void fluxJacobian::FlowMultiplyOnDiagonal(const double &val) {
  // val -- value to multiply along diagonal
  for (auto ii = 0; ii < flowSize_; ++ii) {
    this->FlowJacobian(ii, ii) *= val;
  }
}
void fluxJacobian::TurbMultiplyOnDiagonal(const double &val) {
  // val -- value to multiply along diagonal
  for (auto ii = 0; ii < turbSize_; ++ii) {
    this->TurbJacobian(ii, ii) *= val;
  }
}

void fluxJacobian::FlowAddOnDiagonal(const double &val) {
  // val -- value to multiply along diagonal
  for (auto ii = 0; ii < flowSize_; ++ii) {
    this->FlowJacobian(ii, ii) += val;
  }
}
void fluxJacobian::TurbAddOnDiagonal(const double &val) {
  // val -- value to multiply along diagonal
  for (auto ii = 0; ii < turbSize_; ++ii) {
    this->TurbJacobian(ii, ii) += val;
  }
}

// member function to swap rows of matrix
void fluxJacobian::FlowSwapRows(const int &r1, const int &r2) {
  MSG_ASSERT(r1 < flowSize_ && r2 < flowSize_, "index outside of range");
  if (r1 != r2) {
    std::swap_ranges(this->begin() + this->GetFlowLoc(r1, 0),
                     this->begin() + this->GetFlowLoc(r1, flowSize_),
                     this->begin() + this->GetFlowLoc(r2, 0));
  }
}
void fluxJacobian::TurbSwapRows(const int &r1, const int &r2) {
  MSG_ASSERT(r1 < turbSize_ && r2 < turbSize_, "index outside of range");
  if (r1 != r2) {
    std::swap_ranges(this->begin() + this->GetTurbLoc(r1, 0),
                     this->begin() + this->GetTurbLoc(r1, flowSize_),
                     this->begin() + this->GetTurbLoc(r2, 0));
  }
}

// member function to multiply a row by a given factor
void fluxJacobian::FlowRowMultiply(const int &r, const int &c,
                                   const double &factor) {
  MSG_ASSERT(r < flowSize_ && c < flowSize_, "index outside of range");
  for_each(this->begin() + this->GetFlowLoc(r, c),
           this->begin() + this->GetFlowLoc(r, flowSize_),
           [&factor](auto &val) { val *= factor; });
}
void fluxJacobian::TurbRowMultiply(const int &r, const int &c,
                                   const double &factor) {
  MSG_ASSERT(r < turbSize_ && c < turbSize_, "index outside of range");
  for_each(this->begin() + this->GetTurbLoc(r, c),
           this->begin() + this->GetTurbLoc(r, turbSize_),
           [&factor](auto &val) { val *= factor; });
}

// member function to add a linear combination of one row to another
void fluxJacobian::FlowLinCombRow(const int &r1, const double &factor,
                                  const int &r2) {
  MSG_ASSERT(r1 < flowSize_ && r2 < flowSize_, "index outside of range");
  std::transform(
      this->begin() + this->GetFlowLoc(r2, 0),
      this->begin() + this->GetFlowLoc(r2, flowSize_),
      this->begin() + this->GetFlowLoc(r1, 0),
      this->begin() + this->GetFlowLoc(r2, 0),
      [&factor](const auto &v1, const auto &v2) { return v1 - factor * v2; });
}
void fluxJacobian::TurbLinCombRow(const int &r1, const double &factor,
                                  const int &r2) {
  MSG_ASSERT(r1 < turbSize_ && r2 < turbSize_, "index outside of range");
  std::transform(
      this->begin() + this->GetTurbLoc(r2, 0),
      this->begin() + this->GetTurbLoc(r2, turbSize_),
      this->begin() + this->GetTurbLoc(r1, 0),
      this->begin() + this->GetTurbLoc(r2, 0),
      [&factor](const auto &v1, const auto &v2) { return v1 - factor * v2; });
}

// member function to find maximum absolute value in a given column and range
// within that column and return the corresponding row indice
int fluxJacobian::FlowFindMaxInCol(const int &c, const int &start,
                                   const int &end) const {
  MSG_ASSERT(start < turbSize_ && end < turbSize_ && c < turbSize_,
             "index outside of range");
  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii <= end; ++ii) {
    if (fabs(this->FlowJacobian(ii, c)) > maxVal) {
      maxVal = fabs(this->FlowJacobian(ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}
int fluxJacobian::TurbFindMaxInCol(const int &c, const int &start,
                                   const int &end) const {
  MSG_ASSERT(start < turbSize_ && end < turbSize_ && c < turbSize_,
             "index outside of range");
  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii <= end; ++ii) {
    if (fabs(this->TurbJacobian(ii, c)) > maxVal) {
      maxVal = fabs(this->TurbJacobian(ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}

// member function to set matrix to Identity
void fluxJacobian::FlowIdentity() {
  for (auto rr = 0; rr < this->FlowSize(); ++rr) {
    for (auto cc = 0; cc < this->FlowSize(); ++cc) {
      this->FlowJacobian(rr, cc) = (rr == cc) ? 1.0 : 0.0;
    }
  }
}
void fluxJacobian::TurbIdentity() {
  for (auto rr = 0; rr < this->TurbSize(); ++rr) {
    for (auto cc = 0; cc < this->TurbSize(); ++cc) {
      this->TurbJacobian(rr, cc) = (rr == cc) ? 1.0 : 0.0;
    }
  }
}

// member function to find maximum absolute value on diagonal
// this can be used to find the spectral radius of a diagoanl matrix
double fluxJacobian::FlowMaxAbsValOnDiagonal() const {
  auto maxVal = 0.0;
  for (auto ii = 0; ii < flowSize_; ++ii) {
    maxVal = std::max(fabs(this->FlowJacobian(ii, ii)), maxVal);
  }
  return maxVal;
}
double fluxJacobian::TurbMaxAbsValOnDiagonal() const {
  auto maxVal = 0.0;
  for (auto ii = 0; ii < turbSize_; ++ii) {
    maxVal = std::max(fabs(this->TurbJacobian(ii, ii)), maxVal);
  }
  return maxVal;
}

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

  varArray offDiagonal;

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


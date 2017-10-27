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

#ifndef FLUXJACOBIANHEADERDEF  // only if the macro is not defined
                               // execute these lines of code
#define FLUXJACOBIANHEADERDEF  // define the macro

#include <vector>          // vector
#include <memory>          // unique_ptr
#include <type_traits>
#include <algorithm>
#include <functional>
#include "vector3d.hpp"
#include "tensor.hpp"
#include "uncoupledScalar.hpp"
#include "matrix.hpp"
#include "varArray.hpp"
#include "arrayView.hpp"
#include "utility.hpp"        // TauNormal
#include "inviscidFlux.hpp"   // ConvectiveFluxUpdate
#include "input.hpp"          // input

using std::vector;
using std::ostream;
using std::unique_ptr;

// forward class declarations
class primitive;
class conserved;
class eos;
class transport;
class thermodynamic;
class turbModel;

template <typename T1, typename T2,
          typename = std::enable_if_t<std::is_base_of<varArray, T1>::value ||
                                      std::is_same<varArrayView, T1>::value ||
                                      std::is_same<primitiveView, T1>::value ||
                                      std::is_same<conservedView, T1>::value ||
                                      std::is_same<residualView, T1>::value>,
          typename = std::enable_if_t<std::is_base_of<varArray, T2>::value>>
void ArrayMultiplication(const vector<double>::const_iterator &mat,
                         const int &flowSize, const int &turbSize,
                         const bool &isScalar, const T1 &orig, T2 &arr) {
  MSG_ASSERT(isScalar || (flowSize + turbSize) == orig.Size(),
             "matrix/vector size mismatch");
  MSG_ASSERT(isScalar || (flowSize + turbSize) == arr.Size(),
             "matrix/vector size mismatch");

  auto GetFlowVal = [&flowSize](const auto &mat, const int &r,
                                const int &c) -> decltype(auto) {
    return *(mat + r * flowSize + c);
  };
  auto GetTurbVal = [&flowSize, &turbSize](const auto &mat, const int &r,
                                           const int &c) -> decltype(auto) {
    return *(mat + flowSize * flowSize + r * turbSize + c);
  };

  if (isScalar) {
    for (auto ii = 0; ii < arr.TurbulenceIndex(); ++ii) {
      arr[ii] = orig[ii] * GetFlowVal(mat, 0, 0);
    }
    for (auto ii = arr.TurbulenceIndex(); ii < arr.Size(); ++ii) {
      arr[ii] = orig[ii] * GetTurbVal(mat, 0, 0);
    }
  } else {
    for (auto rr = 0; rr < flowSize; ++rr) {
      for (auto cc = 0; cc < flowSize; ++cc) {
        arr[rr] += GetFlowVal(mat, rr, cc) * orig[cc];
      }
    }

    for (auto rr = 0; rr < turbSize; ++rr) {
      for (auto cc = 0; cc < turbSize; ++cc) {
        arr[flowSize + rr] += GetTurbVal(mat, rr, cc) * orig[flowSize + cc];
      }
    }
  }
}

// This class holds the flux jacobians for the flow and turbulence equations.
// In the LU-SGS method the jacobians are scalars.

class fluxJacobian {
  vector<double> data_;
  int flowSize_;
  int turbSize_;

  // private member functions
  int GetFlowLoc(const int &r, const int &c) const {
    return r * flowSize_ + c;
  }
  int GetTurbLoc(const int &r, const int &c) const {
    return flowSize_ * flowSize_ + r * turbSize_ + c;
  }
  double & FlowJacobian(const int &r, const int &c) {
    return data_[this->GetFlowLoc(r, c)];
  }
  double & TurbJacobian(const int &r, const int &c) {
    return data_[this->GetTurbLoc(r, c)];
  }

 public:
  // constructors
  fluxJacobian(const int &flowSize, const int &turbSize)
      : data_(flowSize * flowSize + turbSize * turbSize, 0.0),
        flowSize_(flowSize),
        turbSize_(turbSize) {}
  fluxJacobian(const double &flow, const double &turb) : fluxJacobian(1, 1) {
    this->FlowJacobian(0, 0) = flow;
    this->TurbJacobian(0, 0) = turb;
  }
  fluxJacobian(const squareMatrix &flow, const squareMatrix &turb)
      : fluxJacobian(flow.Size(), turb.Size()) {
    std::copy(flow.begin(), flow.end(), data_.begin());
    std::copy(turb.begin(), turb.end(), data_.begin() + flowSize_ * flowSize_);
  }
  fluxJacobian() : fluxJacobian(0.0, 0.0) {}
  fluxJacobian(const uncoupledScalar &specRad, const bool &hasTurb)
      : fluxJacobian(1, hasTurb ? 1 : 0) {
    this->FlowJacobian(0, 0) = specRad.FlowVariable();
    if (hasTurb) {
      this->TurbJacobian(0, 0) = specRad.TurbVariable();
    }
  }

  // move constructor and assignment operator
  fluxJacobian(fluxJacobian&&) noexcept = default;
  fluxJacobian& operator=(fluxJacobian&&) noexcept = default;

  // copy constructor and assignment operator
  fluxJacobian(const fluxJacobian&) = default;
  fluxJacobian& operator=(const fluxJacobian&) = default;

  // member functions
  int Size() const {return data_.size();}
  int FlowSize() const { return flowSize_; }
  int TurbSize() const { return turbSize_; }

  const double & FlowJacobian(const int &r, const int &c) const {
    return data_[this->GetFlowLoc(r, c)];
  }
  const double & TurbJacobian(const int &r, const int &c) const {
    return data_[this->GetTurbLoc(r, c)];
  }

  // provide begin and end so std::begin and std::end can be used
  // use lower case to conform with std::begin, std::end
  auto begin() noexcept {return data_.begin();}
  const auto begin() const noexcept {return data_.begin();}
  auto end() noexcept {return data_.end();}
  const auto end() const noexcept {return data_.end();}
  auto beginTurb() noexcept { return data_.begin() + flowSize_ * flowSize_; }
  const auto beginTurb() const noexcept {
    return data_.begin() + flowSize_ * flowSize_;
  }

  // const squareMatrix & FlowJacobian() const {return flowJacobian_;}
  // const squareMatrix & TurbulenceJacobian() const {return turbJacobian_;}

  void AddToFlowJacobian(const squareMatrix &jac);
  void AddToTurbJacobian(const squareMatrix &jac);
  void SubtractFromFlowJacobian(const squareMatrix &jac);
  void SubtractFromTurbJacobian(const squareMatrix &jac);

  void MultFlowJacobian(const double &fac);
  void MultTurbJacobian(const double &fac);

  void MultiplyOnDiagonal(const double &fac) {
    this->FlowMultiplyOnDiagonal(fac);
    this->TurbMultiplyOnDiagonal(fac);
  }
  void FlowMultiplyOnDiagonal(const double &fac) {
    MultiplyFacOnDiagonal(this->begin(), flowSize_, fac);
  }
  void TurbMultiplyOnDiagonal(const double &fac) {
    MultiplyFacOnDiagonal(this->beginTurb(), turbSize_, fac);
  }

  void AddOnDiagonal(const double &fac) { 
    this->FlowAddOnDiagonal(fac);
    this->TurbAddOnDiagonal(fac);
  }
  void FlowAddOnDiagonal(const double &fac) {
    AddFacOnDiagonal(this->begin(), flowSize_, fac);
  }
  void TurbAddOnDiagonal(const double &fac) {
    AddFacOnDiagonal(this->beginTurb(), turbSize_, fac);
  }

  bool HasTurbulence() const { return turbSize_ > 0; }

  fluxJacobian FlowMatMult(const fluxJacobian &) const;
  fluxJacobian TurbMatMult(const fluxJacobian &) const;

  void FlowSwapRows(const int &r1, const int &r2) {
    SwapMatRows(this->begin(), flowSize_, r1, r2);
  }
  void TurbSwapRows(const int &r1, const int &r2) {
    SwapMatRows(this->beginTurb(), turbSize_, r1, r2);
  }

  void FlowRowMultiply(const int &r, const int &c, const double &fac) {
    RowMultiplyFactor(this->begin(), flowSize_, r, c, fac);
  }
  void TurbRowMultiply(const int &r, const int &c, const double &fac) {
    RowMultiplyFactor(this->beginTurb(), turbSize_, r, c, fac);
  }

  void FlowLinCombRow(const int &r1, const double &fac, const int &r2) {
    LinearCombRow(this->begin(), flowSize_, r1, fac, r2);
  }
  void TurbLinCombRow(const int &r1, const double &fac, const int &r2) {
    LinearCombRow(this->beginTurb(), turbSize_, r1, fac, r2);
  }

  int FlowFindMaxInCol(const int &c, const int &start, const int &end) const {
    return FindMaxInColumn(this->begin(), flowSize_, c, start, end);
  }
  int TurbFindMaxInCol(const int &c, const int &start, const int &end) const {
    return FindMaxInColumn(this->beginTurb(), turbSize_, c, start, end);
  }

  void FlowIdentity() { IdentityMatrix(this->begin(), flowSize_); }
  void TurbIdentity() { IdentityMatrix(this->beginTurb(), turbSize_); }

  double FlowMaxAbsValOnDiagonal() const {
    return MaximumAbsValOnDiagonal(this->begin(), flowSize_);
  }
  double TurbMaxAbsValOnDiagonal() const {
    return MaximumAbsValOnDiagonal(this->beginTurb(), turbSize_);
  }

  template <typename T>
  void RusanovFluxJacobian(const T &, const unique_ptr<eos> &,
                           const unique_ptr<thermodynamic> &,
                           const unitVec3dMag<double> &, const bool &,
                           const input &, const unique_ptr<turbModel> &);
  template <typename T>
  void InvFluxJacobian(const T &, const unique_ptr<eos> &,
                       const unique_ptr<thermodynamic> &,
                       const unitVec3dMag<double> &, const input &,
                       const unique_ptr<turbModel> &);
  template <typename T1, typename T2>
  void ApproxRoeFluxJacobian(const T1 &, const T2 &, const unique_ptr<eos> &,
                             const unique_ptr<thermodynamic> &,
                             const unitVec3dMag<double> &, const bool &,
                             const input &, const unique_ptr<turbModel> &);
  template <typename T>
  void DelprimitiveDelConservative(const T &, const unique_ptr<thermodynamic> &,
                                   const unique_ptr<eos> &, const input &);

  template <typename T>
  void ApproxTSLJacobian(const T &, const double &, const double &,
                         const double &, const unique_ptr<eos> &,
                         const unique_ptr<transport> &,
                         const unique_ptr<thermodynamic> &,
                         const unitVec3dMag<double> &, const double &,
                         const unique_ptr<turbModel> &, const input &,
                         const bool &, const tensor<double> &);

  void Zero() {
    std::fill(this->begin(), this->end(), 0.0);
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
  T ArrayMult(const T &orig) const {
    T arr(orig.Size(), orig.NumSpecies());
    ArrayMultiplication(this->begin(), flowSize_, turbSize_, this->IsScalar(),
                        orig, arr);
    return arr;
  }
  template <typename T,
            typename = std::enable_if_t<std::is_same<varArrayView, T>::value ||
                                        std::is_same<primitiveView, T>::value ||
                                        std::is_same<conservedView, T>::value ||
                                        std::is_same<residualView, T>::value>>
  auto ArrayMult(const T &arrView) const {
    auto arr = arrView.GetViewType();
    ArrayMultiplication(this->begin(), flowSize_, turbSize_, this->IsScalar(),
                        arrView, arr);
    return arr;
  }
  bool IsScalar() const {return flowSize_ == 1;}
  void Inverse() {
    this->FlowInverse();
    this->TurbInverse();
  }
  void FlowInverse() { MatrixInverse(this->begin(), flowSize_); }
  void TurbInverse() { MatrixInverse(this->beginTurb(), turbSize_); }

  inline fluxJacobian & operator+=(const fluxJacobian &);
  inline fluxJacobian & operator-=(const fluxJacobian &);
  inline fluxJacobian & operator*=(const fluxJacobian &);
  inline fluxJacobian & operator/=(const fluxJacobian &);

  inline fluxJacobian & operator+=(const double &);
  inline fluxJacobian & operator-=(const double &);
  inline fluxJacobian & operator*=(const double &);
  inline fluxJacobian & operator/=(const double &);

  inline fluxJacobian operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline fluxJacobian operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline fluxJacobian operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline fluxJacobian operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // destructor
  ~fluxJacobian() noexcept {}
};

// function definitions

// operator overload for addition
fluxJacobian & fluxJacobian::operator+=(const fluxJacobian &other) {
  MSG_ASSERT(this->Size() == other.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), other.begin(), this->begin(),
                 std::plus<double>());
  return *this;
}

// operator overload for subtraction with a scalar
fluxJacobian & fluxJacobian::operator-=(const fluxJacobian &other) {
  MSG_ASSERT(this->Size() == other.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), other.begin(), this->begin(),
                 std::minus<double>());
  return *this;
}

// operator overload for elementwise multiplication
fluxJacobian & fluxJacobian::operator*=(const fluxJacobian &other) {
  MSG_ASSERT(this->Size() == other.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), other.begin(), this->begin(),
                 std::multiplies<double>());
  return *this;
}

// operator overload for elementwise division
fluxJacobian & fluxJacobian::operator/=(const fluxJacobian &other) {
  MSG_ASSERT(this->Size() == other.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), other.begin(), this->begin(),
                 std::divides<double>());
  return *this;
}

inline const fluxJacobian operator+(fluxJacobian lhs, const fluxJacobian &rhs) {
  return lhs += rhs;
}

inline const fluxJacobian operator-(fluxJacobian lhs, const fluxJacobian &rhs) {
  return lhs -= rhs;
}

inline const fluxJacobian operator*(fluxJacobian lhs, const fluxJacobian &rhs) {
  return lhs *= rhs;
}

inline const fluxJacobian operator/(fluxJacobian lhs, const fluxJacobian &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
fluxJacobian & fluxJacobian::operator+=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val += scalar; });
  return *this;
}

// operator overload for subtraction with a scalar
fluxJacobian & fluxJacobian::operator-=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val -= scalar; });
  return *this;
}

// operator overload for elementwise multiplication
fluxJacobian & fluxJacobian::operator*=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val *= scalar; });
  return *this;
}

// operator overload for elementwise division
fluxJacobian & fluxJacobian::operator/=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val /= scalar; });
  return *this;
}

inline const fluxJacobian operator+(const double &lhs, fluxJacobian rhs) {
  return rhs += lhs;
}

inline const fluxJacobian operator-(const double &lhs, fluxJacobian rhs) {
  for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs - val; });
  return rhs;
}

inline const fluxJacobian operator*(const double &lhs, fluxJacobian rhs) {
  return rhs *= lhs;
}

inline const fluxJacobian operator/(const double &lhs, fluxJacobian rhs) {
  for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs / val; });
  return rhs;
}

ostream &operator<<(ostream &os, const fluxJacobian &jacobian);

// ---------------------------------------------------------------------------
// member functions

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
template <typename T>
void fluxJacobian::RusanovFluxJacobian(const T &state,
                                       const unique_ptr<eos> &eqnState,
                                       const unique_ptr<thermodynamic> &thermo,
                                       const unitVec3dMag<double> &area,
                                       const bool &positive, const input &inp,
                                       const unique_ptr<turbModel> &turb) {
  // state -- primitive variables at face
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  // face inviscid spectral radius
  const auto specRad = InvFaceSpectralRadius(state, area, thermo, eqnState);

  // form dissipation matrix based on spectral radius
  fluxJacobian dissipation(inp.NumFlowEquations(), inp.NumTurbEquations());
  dissipation.FlowIdentity();
  dissipation.FlowMultiplyOnDiagonal(specRad);

  // begin jacobian calculation
  this->InvFluxJacobian(state, eqnState, thermo, area, inp, turb);

  // compute turbulent dissipation if necessary
  if (inp.IsRANS()) {
    // multiply by 0.5 b/c averaging with convection matrix
    const auto tJac = 0.5 * turb->InviscidDissJacobian(state, area);
    std::copy(tJac.begin(), tJac.end(), dissipation.beginTurb());
  }

  positive ? (*this) += dissipation : (*this) -= dissipation;
}

// function to calculate inviscid flux jacobian
template <typename T>
void fluxJacobian::InvFluxJacobian(const T &state,
                                   const unique_ptr<eos> &eqnState,
                                   const unique_ptr<thermodynamic> &thermo,
                                   const unitVec3dMag<double> &area,
                                   const input &inp,
                                   const unique_ptr<turbModel> &turb) {
  // state -- primitive variables at face
  // eqnState -- ideal gas equation of state
  // thermo -- thermodynamic model
  // area -- face area vector
  // inp -- input variables
  // turb -- turbulence model
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  const auto t = state.Temperature(eqnState);
  const auto velNorm = state.Velocity().DotProd(area.UnitVector());
  const auto gammaMinusOne = thermo->Gamma(t) - 1.0;
  const auto phi = 0.5 * gammaMinusOne * state.Velocity().MagSq();
  const auto a1 = thermo->Gamma(t) * state.Energy(eqnState, thermo) - phi;
  const auto a3 = thermo->Gamma(t) - 2.0;

  // begin jacobian calculation
  *this = fluxJacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // calculate flux derivatives wrt left state
  // column zero
  this->FlowJacobian(0, 0) = 0.0;
  this->FlowJacobian(1, 0) = phi * area.UnitVector().X() - state.U() * velNorm;
  this->FlowJacobian(2, 0) = phi * area.UnitVector().Y() - state.V() * velNorm;
  this->FlowJacobian(3, 0) = phi * area.UnitVector().Z() - state.W() * velNorm;
  this->FlowJacobian(4, 0) = velNorm * (phi - a1);

  // column one
  this->FlowJacobian(0, 1) = area.UnitVector().X();
  this->FlowJacobian(1, 1) = velNorm - a3 * area.UnitVector().X() * state.U();
  this->FlowJacobian(2, 1) = state.V() * area.UnitVector().X() -
      gammaMinusOne * state.U() * area.UnitVector().Y();
  this->FlowJacobian(3, 1) = state.W() * area.UnitVector().X() -
      gammaMinusOne * state.U() * area.UnitVector().Z();
  this->FlowJacobian(4, 1) =
      a1 * area.UnitVector().X() - gammaMinusOne * state.U() * velNorm;

  // column two
  this->FlowJacobian(0, 2) = area.UnitVector().Y();
  this->FlowJacobian(1, 2) = state.U() * area.UnitVector().Y() -
      gammaMinusOne * state.V() * area.UnitVector().X();
  this->FlowJacobian(2, 2) = velNorm - a3 * area.UnitVector().Y() * state.V();
  this->FlowJacobian(3, 2) = state.W() * area.UnitVector().Y() -
      gammaMinusOne * state.V() * area.UnitVector().Z();
  this->FlowJacobian(4, 2) =
      a1 * area.UnitVector().Y() - gammaMinusOne * state.V() * velNorm;

  // column three
  this->FlowJacobian(0, 3) = area.UnitVector().Z();
  this->FlowJacobian(1, 3) = state.U() * area.UnitVector().Z() -
      gammaMinusOne * state.W() * area.UnitVector().X();
  this->FlowJacobian(2, 3) = state.V() * area.UnitVector().Z() -
      gammaMinusOne * state.W() * area.UnitVector().Y();
  this->FlowJacobian(3, 3) = velNorm - a3 * area.UnitVector().Z() * state.W();
  this->FlowJacobian(4, 3) =
      a1 * area.UnitVector().Z() - gammaMinusOne * state.W() * velNorm;

  // column four
  this->FlowJacobian(0, 4) = 0.0;
  this->FlowJacobian(1, 4) = gammaMinusOne * area.UnitVector().X();
  this->FlowJacobian(2, 4) = gammaMinusOne * area.UnitVector().Y();
  this->FlowJacobian(3, 4) = gammaMinusOne * area.UnitVector().Z();
  this->FlowJacobian(4, 4) = thermo->Gamma(t) * velNorm;

  // multiply by 0.5 b/c averaging with dissipation matrix
  this->MultFlowJacobian(0.5 * area.Mag());

  // turbulent jacobian here
  if (inp.IsRANS()) {
    // multiply by 0.5 b/c averaging with dissipation matrix
    const auto tJac = 0.5 * turb->InviscidConvJacobian(state, area);
    std::copy(tJac.begin(), tJac.end(), this->beginTurb());
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
template <typename T1, typename T2>
void fluxJacobian::ApproxRoeFluxJacobian(
    const T1 &left, const T2 &right,
    const unique_ptr<eos> &eqnState, const unique_ptr<thermodynamic> &thermo,
    const unitVec3dMag<double> &area, const bool &positive, const input &inp,
    const unique_ptr<turbModel> &turb) {
  // left -- primitive variables from left side
  // right -- primitive variables from right side
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // area -- face area vector
  // positive -- flag to determine whether to add or subtract dissipation
  // inp -- input variables
  // turb -- turbulence model
  static_assert(std::is_same<primitive, T1>::value ||
                    std::is_same<primitiveView, T1>::value,
                "T1 requires primitive or primativeView type");
  static_assert(std::is_same<primitive, T2>::value ||
                    std::is_same<primitiveView, T2>::value,
                "T2 requires primitive or primativeView type");

  // compute Roe averaged state
  const auto roeAvg = RoeAveragedState(left, right);

  // compute Roe matrix
  fluxJacobian roeMatrix;
  roeMatrix.InvFluxJacobian(roeAvg, eqnState, thermo, area, inp, turb);

  // compute convective flux jacobian
  positive ? this->InvFluxJacobian(left, eqnState, thermo, area, inp, turb) :
      this->InvFluxJacobian(right, eqnState, thermo, area, inp, turb);

  positive ? (*this) += roeMatrix : (*this) -= roeMatrix;
}

// change of variable matrix going from primitive to conservative variables
// from Dwight
template <typename T>
void fluxJacobian::DelprimitiveDelConservative(
    const T &state, const unique_ptr<thermodynamic> &thermo,
    const unique_ptr<eos> &eqnState, const input &inp) {
  // state -- primitive variables
  // thermo -- thermodynamic model
  // inp -- input variables
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  const auto t = state.Temperature(eqnState);
  const auto gammaMinusOne = thermo->Gamma(t) - 1.0;
  const auto invRho = 1.0 / state.Rho();

  *this = fluxJacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  // assign column 0
  this->FlowJacobian(0, 0) = 1.0;
  this->FlowJacobian(1, 0) = -invRho * state.U();
  this->FlowJacobian(2, 0) = -invRho * state.V();
  this->FlowJacobian(3, 0) = -invRho * state.W();
  this->FlowJacobian(4, 0) = 0.5 * gammaMinusOne *
      state.Velocity().DotProd(state.Velocity());

  // assign column 1
  this->FlowJacobian(1, 1) = invRho;
  this->FlowJacobian(4, 1) = -gammaMinusOne * state.U();

  // assign column 2
  this->FlowJacobian(2, 2) = invRho;
  this->FlowJacobian(4, 2) = -gammaMinusOne * state.V();

  // assign column 3
  this->FlowJacobian(3, 3) = invRho;
  this->FlowJacobian(4, 3) = -gammaMinusOne * state.W();

  // assign column 4
  this->FlowJacobian(4, 4) = gammaMinusOne;

  // turbulent jacobian here
  if (inp.IsRANS()) {
    this->TurbJacobian(0, 0) = invRho;
    this->TurbJacobian(1, 1) = invRho;
  }
}

// approximate thin shear layer jacobian following implementation in Dwight.
// does not use any gradients
template <typename T>
void fluxJacobian::ApproxTSLJacobian(
    const T &state, const double &lamVisc, const double &turbVisc,
    const double &f1, const unique_ptr<eos> &eqnState,
    const unique_ptr<transport> &trans, const unique_ptr<thermodynamic> &thermo,
    const unitVec3dMag<double> &area, const double &dist,
    const unique_ptr<turbModel> &turb, const input &inp, const bool &left,
    const tensor<double> &vGrad) {
  // state -- primitive variables
  // eos -- equation of state
  // trans -- viscous transport model
  // area -- face area vector
  // dist -- distance from cell center to cell center
  // turb --  turbulence model
  // inp -- input variables
  // left -- flag that is negative if using left state
  // vGrad -- velocity gradient
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");
  *this = fluxJacobian(inp.NumFlowEquations(), inp.NumTurbEquations());

  const auto t = state.Temperature(eqnState);
  const auto mu = trans->NondimScaling() * lamVisc;
  const auto mut = trans->NondimScaling() * turbVisc;
  const auto velNorm = state.Velocity().DotProd(area.UnitVector());

  const auto tauNorm = TauNormal(vGrad, area.UnitVector(), mu, mut, trans);

  auto fac = left ? -1.0 : 1.0;

  constexpr auto third = 1.0 / 3.0;

  // assign column 0
  this->FlowJacobian(4, 0) =
      -(trans->Conductivity(mu, t, thermo) +
        trans->TurbConductivity(mut, turb->TurbPrandtlNumber(), t, thermo)) *
      state.Temperature(eqnState) / ((mu + mut) * state.Rho());

  // assign column 1
  this->FlowJacobian(1, 1) =
      third * area.UnitVector().X() * area.UnitVector().X() + 1.0;
  this->FlowJacobian(2, 1) =
      third * area.UnitVector().X() * area.UnitVector().Y();
  this->FlowJacobian(3, 1) =
      third * area.UnitVector().X() * area.UnitVector().Z();
  this->FlowJacobian(4, 1) = fac * 0.5 * dist / (mu + mut) * tauNorm.X() +
      third * area.UnitVector().X() * velNorm + state.U();

  // assign column 2
  this->FlowJacobian(1, 2) =
      third * area.UnitVector().Y() * area.UnitVector().X();
  this->FlowJacobian(2, 2) =
      third * area.UnitVector().Y() * area.UnitVector().Y() + 1.0;
  this->FlowJacobian(3, 2) =
      third * area.UnitVector().Y() * area.UnitVector().Z();
  this->FlowJacobian(4, 2) = fac * 0.5 * dist / (mu + mut) * tauNorm.Y() +
                             third * area.UnitVector().Y() * velNorm +
                             state.V();

  // assign column 3
  this->FlowJacobian(1, 3) =
      third * area.UnitVector().Z() * area.UnitVector().X();
  this->FlowJacobian(2, 3) =
      third * area.UnitVector().Z() * area.UnitVector().Y();
  this->FlowJacobian(3, 3) =
      third * area.UnitVector().Z() * area.UnitVector().Z() + 1.0;
  this->FlowJacobian(4, 3) = fac * 0.5 * dist / (mu + mut) * tauNorm.Z() +
                             third * area.UnitVector().Z() * velNorm +
                             state.W();

  // assign column 4
  this->FlowJacobian(4, 4) =
      (trans->Conductivity(mu, t, thermo) +
       trans->TurbConductivity(mut, turb->TurbPrandtlNumber(), t, thermo)) /
      ((mu + mut) * state.Rho());

  this->MultFlowJacobian(area.Mag() * (mu + mut) / dist);

  fluxJacobian prim2Cons;
  prim2Cons.DelprimitiveDelConservative(state, thermo, eqnState, inp);
  const auto product = this->FlowMatMult(prim2Cons);
  std::copy(product.begin(), product.beginTurb(), this->begin());

  // calculate turbulent jacobian if necessary
  if (inp.IsRANS()) {
    const auto turbProd =
        fac *
        turb->ViscousJacobian(state, area, lamVisc, trans, dist, turbVisc, f1);
    std::copy(turbProd.begin(), turbProd.end(), this->beginTurb());
    // Don't need to multiply by prim2Cons b/c jacobian is already wrt
    // conservative variables
  }
}

// ---------------------------------------------------------------------------
// non member functions
varArray RusanovScalarOffDiagonal(const primitiveView &, const varArrayView &,
                                  const unitVec3dMag<double> &, const double &,
                                  const double &, const double &,
                                  const double &, const unique_ptr<eos> &,
                                  const unique_ptr<thermodynamic> &,
                                  const unique_ptr<transport> &,
                                  const unique_ptr<turbModel> &, const bool &,
                                  const bool &);
varArray RusanovBlockOffDiagonal(const primitiveView &, const varArrayView &,
                                 const unitVec3dMag<double> &, const double &,
                                 const double &, const double &, const double &,
                                 const unique_ptr<eos> &,
                                 const unique_ptr<thermodynamic> &,
                                 const unique_ptr<transport> &,
                                 const unique_ptr<turbModel> &, const input &,
                                 const bool &, const tensor<double> &);

varArray RoeOffDiagonal(const primitiveView &, const primitiveView &,
                        const varArrayView &, const unitVec3dMag<double> &,
                        const double &, const double &, const double &,
                        const double &, const unique_ptr<eos> &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &, const bool &,
                        const bool &, const bool &);

varArray OffDiagonal(const primitiveView &, const primitiveView &,
                     const varArrayView &, const unitVec3dMag<double> &,
                     const double &, const double &, const double &,
                     const double &, const tensor<double> &,
                     const unique_ptr<eos> &, const unique_ptr<thermodynamic> &,
                     const unique_ptr<transport> &,
                     const unique_ptr<turbModel> &, const input &,
                     const bool &);

#endif

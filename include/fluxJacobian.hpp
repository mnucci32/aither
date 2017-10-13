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

// This class holds the flux jacobians for the flow and turbulence equations.
// In the LU-SGS method the jacobians are scalars.

class fluxJacobian {
  squareMatrix flowJacobian_;
  squareMatrix turbJacobian_;

 public:
  // constructors
  fluxJacobian(const double &flow, const double &turb);
  fluxJacobian(const int &flowSize, const int &turbSize);
  fluxJacobian(const squareMatrix &flow, const squareMatrix &turb)
      : flowJacobian_(flow), turbJacobian_(turb) {}
  fluxJacobian() : fluxJacobian(0.0, 0.0) {}
  explicit fluxJacobian(const uncoupledScalar &specRad) :
      fluxJacobian(specRad.FlowVariable(), specRad.TurbVariable()) {}

  // move constructor and assignment operator
  fluxJacobian(fluxJacobian&&) noexcept = default;
  fluxJacobian& operator=(fluxJacobian&&) noexcept = default;

  // copy constructor and assignment operator
  fluxJacobian(const fluxJacobian&) = default;
  fluxJacobian& operator=(const fluxJacobian&) = default;

  // member functions
  squareMatrix FlowJacobian() const {return flowJacobian_;}
  squareMatrix TurbulenceJacobian() const {return turbJacobian_;}

  void AddToFlowJacobian(const squareMatrix &jac) {flowJacobian_ += jac;}
  void AddToTurbJacobian(const squareMatrix &jac) {turbJacobian_ += jac;}
  void SubtractFromFlowJacobian(const squareMatrix &jac) {flowJacobian_ -= jac;}
  void SubtractFromTurbJacobian(const squareMatrix &jac) {turbJacobian_ -= jac;}

  void MultiplyOnDiagonal(const double &, const bool &);
  void AddOnDiagonal(const double &, const bool &);
  
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
    flowJacobian_.Zero();
    turbJacobian_.Zero();
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
  T ArrayMult(const T &orig) const {
    auto arr = orig;
    if (this->IsScalar()) {
      for (auto ii = 0; ii < arr.TurbulenceIndex(); ++ii) {
        arr[ii] *= flowJacobian_(0, 0);
      }
      for (auto ii = arr.TurbulenceIndex(); ii < arr.Size(); ++ii) {
        arr[ii] *= turbJacobian_(0, 0);
      }
    } else {
      arr = flowJacobian_.ArrayMult(arr);
      arr = turbJacobian_.ArrayMult(arr, flowJacobian_.Size());
    }
    return arr;
  }
  template <typename T,
            typename = std::enable_if_t<std::is_same<varArrayView, T>::value ||
                                        std::is_same<primitiveView, T>::value ||
                                        std::is_same<conservedView, T>::value ||
                                        std::is_same<residualView, T>::value>>
  auto ArrayMult(const T &arrView) const {
    auto arr = arrView.CopyData();
    return this->ArrayMult(arr);
  }
  bool IsScalar() const;
  void Inverse(const bool &);

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

  friend inline const fluxJacobian operator-(const double &lhs,
                                             fluxJacobian rhs);
  friend inline const fluxJacobian operator/(const double &lhs,
                                             fluxJacobian rhs);

  // destructor
  ~fluxJacobian() noexcept {}
};

// function definitions

// operator overload for addition
fluxJacobian & fluxJacobian::operator+=(const fluxJacobian &other) {
  flowJacobian_ += other.flowJacobian_;
  turbJacobian_ += other.turbJacobian_;
  return *this;
}

// operator overload for subtraction with a scalar
fluxJacobian & fluxJacobian::operator-=(const fluxJacobian &other) {
  flowJacobian_ -= other.flowJacobian_;
  turbJacobian_ -= other.turbJacobian_;
  return *this;
}

// operator overload for elementwise multiplication
fluxJacobian & fluxJacobian::operator*=(const fluxJacobian &other) {
  flowJacobian_ *= other.flowJacobian_;
  turbJacobian_ *= other.turbJacobian_;
  return *this;
}

// operator overload for elementwise division
fluxJacobian & fluxJacobian::operator/=(const fluxJacobian &other) {
  flowJacobian_ /= other.flowJacobian_;
  turbJacobian_ /= other.turbJacobian_;
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
  flowJacobian_ += scalar;
  turbJacobian_ += scalar;
  return *this;
}

// operator overload for subtraction with a scalar
fluxJacobian & fluxJacobian::operator-=(const double &scalar) {
  flowJacobian_ -= scalar;
  turbJacobian_ -= scalar;
  return *this;
}

// operator overload for elementwise multiplication
fluxJacobian & fluxJacobian::operator*=(const double &scalar) {
  flowJacobian_ *= scalar;
  turbJacobian_ *= scalar;
  return *this;
}

// operator overload for elementwise division
fluxJacobian & fluxJacobian::operator/=(const double &scalar) {
  flowJacobian_ /= scalar;
  turbJacobian_ /= scalar;
  return *this;
}

inline const fluxJacobian operator+(const double &lhs, fluxJacobian rhs) {
  return rhs += lhs;
}

inline const fluxJacobian operator-(const double &lhs, fluxJacobian rhs) {
  rhs.flowJacobian_ = lhs - rhs.flowJacobian_;
  rhs.turbJacobian_ = lhs - rhs.turbJacobian_;
  return rhs;
}

inline const fluxJacobian operator*(const double &lhs, fluxJacobian rhs) {
  return rhs *= lhs;
}

inline const fluxJacobian operator/(const double &lhs, fluxJacobian rhs) {
  rhs.flowJacobian_ = lhs / rhs.flowJacobian_;
  rhs.turbJacobian_ = lhs / rhs.turbJacobian_;
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
  dissipation.flowJacobian_.Identity();
  dissipation.flowJacobian_ *= specRad;

  // begin jacobian calculation
  this->InvFluxJacobian(state, eqnState, thermo, area, inp, turb);

  // compute turbulent dissipation if necessary
  if (inp.IsRANS()) {
    // multiply by 0.5 b/c averaging with convection matrix
    dissipation.turbJacobian_ = 0.5 * turb->InviscidDissJacobian(state, area);
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
  flowJacobian_(4, 4) = thermo->Gamma(t) * velNorm;

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

  flowJacobian_ = squareMatrix(inp.NumFlowEquations());
  turbJacobian_ = squareMatrix(inp.NumTurbEquations());

  const auto t = state.Temperature(eqnState);
  const auto mu = trans->NondimScaling() * lamVisc;
  const auto mut = trans->NondimScaling() * turbVisc;
  const auto velNorm = state.Velocity().DotProd(area.UnitVector());

  const auto tauNorm = TauNormal(vGrad, area.UnitVector(), mu, mut, trans);

  auto fac = left ? -1.0 : 1.0;

  constexpr auto third = 1.0 / 3.0;

  // assign column 0
  flowJacobian_(4, 0) =
      -(trans->Conductivity(mu, t, thermo) +
        trans->TurbConductivity(mut, turb->TurbPrandtlNumber(), t, thermo)) *
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
      (trans->Conductivity(mu, t, thermo) +
       trans->TurbConductivity(mut, turb->TurbPrandtlNumber(), t, thermo)) /
      ((mu + mut) * state.Rho());

  flowJacobian_ *= area.Mag() * (mu + mut) / dist;

  fluxJacobian prim2Cons;
  prim2Cons.DelprimitiveDelConservative(state, thermo, eqnState, inp);
  flowJacobian_ = flowJacobian_.MatMult(prim2Cons.flowJacobian_);

  // calculate turbulent jacobian if necessary
  if (inp.IsRANS()) {
    turbJacobian_ = fac * turb->ViscousJacobian(state, area, lamVisc, trans,
                                                dist, turbVisc, f1);
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

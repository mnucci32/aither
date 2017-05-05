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
#include "vector3d.hpp"
#include "tensor.hpp"
#include "uncoupledScalar.hpp"
#include "matrix.hpp"

using std::vector;
using std::ostream;
using std::unique_ptr;

// forward class declarations
class primVars;
class eos;
class transport;
class thermodynamic;
class turbModel;
class input;
class genArray;

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

  void RusanovFluxJacobian(const primVars &, const unique_ptr<eos> &,
                           const unique_ptr<thermodynamic> &,
                           const unitVec3dMag<double> &, const bool &,
                           const input &, const unique_ptr<turbModel> &);
  void InvFluxJacobian(const primVars &, const unique_ptr<eos> &,
                       const unique_ptr<thermodynamic> &,
                       const unitVec3dMag<double> &, const input &,
                       const unique_ptr<turbModel> &);
  void ApproxRoeFluxJacobian(const primVars &, const primVars &,
                             const unique_ptr<eos> &,
                             const unique_ptr<thermodynamic> &,
                             const unitVec3dMag<double> &, const bool &,
                             const input &, const unique_ptr<turbModel> &);
  void DelPrimativeDelConservative(const primVars &,
                                   const unique_ptr<thermodynamic> &,
                                   const input &);

  void ApproxTSLJacobian(
      const primVars &, const double &, const double &, const double &,
      const unique_ptr<eos> &, const unique_ptr<transport> &,
      const unique_ptr<thermodynamic> &, const unitVec3dMag<double> &,
      const double &, const unique_ptr<turbModel> &, const input &,
      const bool &, const tensor<double> &);

  void Zero() {
    flowJacobian_.Zero();
    turbJacobian_.Zero();
  }

  genArray ArrayMult(genArray) const;
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

genArray RusanovScalarOffDiagonal(const primVars &, const genArray &,
                                  const unitVec3dMag<double> &,
                                  const double &, const double &,
                                  const double &, const double &,
                                  const unique_ptr<eos> &,
                                  const unique_ptr<thermodynamic> &,
                                  const unique_ptr<transport> &,
                                  const unique_ptr<turbModel> &,
                                  const bool &, const bool &);
genArray RusanovBlockOffDiagonal(const primVars &, const genArray &,
                                 const unitVec3dMag<double> &,
                                 const double &, const double &,
                                 const double &, const double &,
                                 const unique_ptr<eos> &,
                                 const unique_ptr<thermodynamic> &,
                                 const unique_ptr<transport> &,
                                 const unique_ptr<turbModel> &,
                                 const input &, const bool &,
                                 const tensor<double> &);

genArray RoeOffDiagonal(const primVars &, const primVars &,
                        const genArray &,
                        const unitVec3dMag<double> &, const double &,
                        const double &, const double &,
                        const double &, const unique_ptr<eos> &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &, const bool &,
                        const bool &, const bool &);

genArray OffDiagonal(const primVars &, const primVars &, const genArray &,
                     const unitVec3dMag<double> &, const double &,
                     const double &, const double &, const double &,
                     const tensor<double> &, const unique_ptr<eos> &,
                     const unique_ptr<thermodynamic> &,
                     const unique_ptr<transport> &,
                     const unique_ptr<turbModel> &, const input &,
                     const bool &);

#endif

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

#ifndef FLUXJACOBIANHEADERDEF  // only if the macro is not defined
                               // execute these lines of code
#define FLUXJACOBIANHEADERDEF  // define the macro

#include <iostream>        // cout
#include <vector>          // vector
#include <string>          // string
#include <memory>          // unique_ptr
#include "vector3d.hpp"
#include "uncoupledScalar.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std:: unique_ptr;

// forward class declarations
class primVars;
class idealGas;
class sutherland;
class turbModel;
class input;
class genArray;
class squareMatrix;

// This class holds the flux jacobians for the flow and turbulence equations.
// In the LU-SGS method the jacobians are scalars.

class fluxJacobian {
  double flowJacobian_;
  double turbJacobian_;

 public:
  // constructors
  fluxJacobian(const double &flow, const double &turb) :
      flowJacobian_(flow), turbJacobian_(turb) {}
  fluxJacobian() : fluxJacobian(0.0, 0.0) {}
    fluxJacobian(const uncoupledScalar &specRad) :
      fluxJacobian(specRad.FlowVariable(), specRad.TurbVariable()) {}  
  fluxJacobian(const primVars &, const unitVec3dMag<double> &,
               const unitVec3dMag<double> &, const idealGas &,
               const sutherland &, const double &,
               const unique_ptr<turbModel> &, const input &, const bool &);

  // move constructor and assignment operator
  fluxJacobian(fluxJacobian&&) noexcept = default;
  fluxJacobian& operator=(fluxJacobian&&) noexcept = default;

  // copy constructor and assignment operator
  fluxJacobian(const fluxJacobian&) = default;
  fluxJacobian& operator=(const fluxJacobian&) = default;

  // member functions
  double FlowJacobian() const { return flowJacobian_; }
  double TurbulenceJacobian() const { return turbJacobian_; }

  void AddInviscidJacobian(const primVars &, const unitVec3dMag<double> &,
                           const unitVec3dMag<double> &, const idealGas &,
                           const unique_ptr<turbModel> &, const bool &);
  void AddViscousJacobian(const primVars &, const unitVec3dMag<double> &,
                          const unitVec3dMag<double> &, const idealGas &,
                          const sutherland &, const double &,
                          const unique_ptr<turbModel> &, const bool &);
  void AddTurbSourceJacobian(const primVars &, const sutherland &,
                             const double &, const unique_ptr<turbModel> &);

  void Zero() {
    flowJacobian_ = 0.0;
    turbJacobian_ = 0.0;
  }

  genArray ArrayMult(genArray) const;

  fluxJacobian Inverse(const bool &) const;

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

ostream &operator<<(ostream &os, const fluxJacobian &);

squareMatrix RusanovFluxJacobian(const primVars &, const primVars &,
                                 const idealGas &, const vector3d<double> &,
				 const bool &);
squareMatrix InvFluxJacobian(const primVars &, const idealGas &,
                             const vector3d<double> &);
squareMatrix ApproxRoeFluxJacobian(const primVars &, const primVars &,
				   const idealGas &, const vector3d<double> &,
				   const bool &);
genArray RusanovOffDiagonal(const primVars &, const genArray &,
			    const unitVec3dMag<double> &, const unitVec3dMag<double> &,
			    const double &, const idealGas &, const sutherland &,
			    const unique_ptr<turbModel> &, const input &, const bool &);
genArray RoeOffDiagonal(const primVars &, const primVars &, const genArray &,
			const unitVec3dMag<double> &, const unitVec3dMag<double> &,
			const double &, const idealGas &, const sutherland &,
			const unique_ptr<turbModel> &, const input &, const bool &);


#endif

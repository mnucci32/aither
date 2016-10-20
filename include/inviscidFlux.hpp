/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef INVFLUXHEADERDEF  // only if the macro INVFLUXHEADERDEF is not defined
                          // execute these lines of code
#define INVFLUXHEADERDEF  // define the macro

#include <vector>        // vector
#include <string>        // string
#include <iostream>      // cout
#include <memory>        // unique_ptr
#include "vector3d.hpp"  // vector3d
#include "macros.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declaration
class idealGas;
class primVars;
class genArray;
class squareMatrix;
class turbModel;

class inviscidFlux {
  double data_[NUMVARS];  // rho dot velocity vector
  // rho dot velocity vector * u-velocity + pressure * i-dir-vector
  // rho dot velocity vector * v-velocity + pressure * j-dir-vector
  // rho dot velocity vector * w-velocity + pressure * k-dir-vector
  // rho dot velocity vector * enthalpy

  // private member functions
  void ConstructFromPrim(const primVars&, const idealGas&,
                         const vector3d<double>&);

 public:
  // constructors
  inviscidFlux() : data_{0.0} {}
  inviscidFlux(const primVars&, const idealGas&, const vector3d<double>&);
  inviscidFlux(const genArray&, const idealGas&, const unique_ptr<turbModel>&,
               const vector3d<double>&);

  // move constructor and assignment operator
  inviscidFlux(inviscidFlux&&) noexcept = default;
  inviscidFlux& operator=(inviscidFlux&&) noexcept = default;

  // copy constructor and assignment operator
  inviscidFlux(const inviscidFlux&) = default;
  inviscidFlux& operator=(const inviscidFlux&) = default;

  // member functions
  double RhoVel() const { return data_[0]; }
  double RhoVelU() const { return data_[1]; }
  double RhoVelV() const { return data_[2]; }
  double RhoVelW() const { return data_[3]; }
  double RhoVelH() const { return data_[4]; }
  double RhoVelK() const { return data_[5]; }
  double RhoVelO() const { return data_[6]; }

  void RoeFlux(const inviscidFlux&, const genArray&);

  inline inviscidFlux & operator+=(const inviscidFlux &);
  inline inviscidFlux & operator-=(const inviscidFlux &);
  inline inviscidFlux & operator*=(const inviscidFlux &);
  inline inviscidFlux & operator/=(const inviscidFlux &);

  inline inviscidFlux & operator+=(const double &);
  inline inviscidFlux & operator-=(const double &);
  inline inviscidFlux & operator*=(const double &);
  inline inviscidFlux & operator/=(const double &);

  inline inviscidFlux operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline inviscidFlux operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline inviscidFlux operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline inviscidFlux operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  friend inline const inviscidFlux operator-(const double &lhs,
                                             inviscidFlux rhs);
  friend inline const inviscidFlux operator/(const double &lhs,
                                             inviscidFlux rhs);

  genArray ConvertToGenArray() const;

  // destructor
  ~inviscidFlux() noexcept {}
};

// function definitions
// function to calculate Roe flux with entropy fix
inviscidFlux RoeFlux(const primVars&, const primVars&, const idealGas&,
                     const vector3d<double>&);
inviscidFlux RusanovFlux(const primVars&, const primVars&, const idealGas&,
                         const vector3d<double>&, const bool&);

// function to calculate Roe flux with entropy fix for implicit methods
void ApproxRoeFluxJacobian(const primVars&, const primVars&, const idealGas&,
                           const vector3d<double>&, double&, squareMatrix&,
                           squareMatrix&);

genArray ConvectiveFluxUpdate(const primVars&, const primVars&, const idealGas&,
                              const vector3d<double>&);

// operator overload for addition
inviscidFlux & inviscidFlux::operator+=(const inviscidFlux &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] += arr.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
inviscidFlux & inviscidFlux::operator-=(const inviscidFlux &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] -= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
inviscidFlux & inviscidFlux::operator*=(const inviscidFlux &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] *= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
inviscidFlux & inviscidFlux::operator/=(const inviscidFlux &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] /= arr.data_[rr];
  }
  return *this;
}

inline const inviscidFlux operator+(inviscidFlux lhs, const inviscidFlux &rhs) {
  return lhs += rhs;
}

inline const inviscidFlux operator-(inviscidFlux lhs, const inviscidFlux &rhs) {
  return lhs -= rhs;
}

inline const inviscidFlux operator*(inviscidFlux lhs, const inviscidFlux &rhs) {
  return lhs *= rhs;
}

inline const inviscidFlux operator/(inviscidFlux lhs, const inviscidFlux &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
inviscidFlux & inviscidFlux::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
inviscidFlux & inviscidFlux::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
inviscidFlux & inviscidFlux::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
inviscidFlux & inviscidFlux::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const inviscidFlux operator+(const double &lhs, inviscidFlux rhs) {
  return rhs += lhs;
}

inline const inviscidFlux operator-(const double &lhs, inviscidFlux rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs - rhs.data_[rr];
  }
  return rhs;
}

inline const inviscidFlux operator*(const double &lhs, inviscidFlux rhs) {
  return rhs *= lhs;
}

inline const inviscidFlux operator/(const double &lhs, inviscidFlux rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs / rhs.data_[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const inviscidFlux &);

#endif

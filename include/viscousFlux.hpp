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

#ifndef VISCFLUXHEADERDEF  // only if the macro VISCFLUXHEADERDEF is not defined
                           // execute these lines of code
#define VISCFLUXHEADERDEF  // define the macro

#include <iostream>        // cout
#include <memory>          // unique_ptr
#include "vector3d.hpp"    // vector3d
#include "tensor.hpp"      // tensor
#include "macros.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declarations
class primVars;
class eos;
class transport;
class thermodynamic;
class turbModel;
class squareMatrix;
struct wallVars;

class viscousFlux {
  double data_[NUMVARS - 1];  // viscous flux for x-momentum equation
  // viscous flux for y-momentum equation
  // viscous flux for z-momentum equation
  // viscous flux for energy equation

 public:
  // constructors
  viscousFlux() : data_{0.0} {}

  // move constructor and assignment operator
  viscousFlux(viscousFlux&&) noexcept = default;
  viscousFlux& operator=(viscousFlux&&) noexcept = default;

  // copy constructor and assignment operator
  viscousFlux(const viscousFlux&) = default;
  viscousFlux& operator=(const viscousFlux&) = default;

  // member functions
  double MomX() const { return data_[0]; }
  double MomY() const { return data_[1]; }
  double MomZ() const { return data_[2]; }
  double Engy() const { return data_[3]; }
  double MomK() const { return data_[4]; }
  double MomO() const { return data_[5]; }

  void CalcFlux(const tensor<double> &, const unique_ptr<transport> &,
                const unique_ptr<thermodynamic> &,
                const unique_ptr<eos> &eqnState, const vector3d<double> &,
                const vector3d<double> &, const vector3d<double> &,
                const vector3d<double> &, const unique_ptr<turbModel> &,
                const primVars &, const double &, const double &,
                const double &);
  wallVars CalcWallFlux(const tensor<double> &, const unique_ptr<transport> &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<eos> &, const vector3d<double> &,
                        const vector3d<double> &, const vector3d<double> &,
                        const vector3d<double> &, const unique_ptr<turbModel> &,
                        const primVars &, const double &, const double &,
                        const double &);
  void CalcWallLawFlux(const vector3d<double> &, const double &, const double &,
                       const double &, const vector3d<double> &,
                       const vector3d<double> &, const vector3d<double> &,
                       const vector3d<double> &, const unique_ptr<turbModel> &);

  inline viscousFlux &operator+=(const viscousFlux &);
  inline viscousFlux & operator-=(const viscousFlux &);
  inline viscousFlux & operator*=(const viscousFlux &);
  inline viscousFlux & operator/=(const viscousFlux &);

  inline viscousFlux & operator+=(const double &);
  inline viscousFlux & operator-=(const double &);
  inline viscousFlux & operator*=(const double &);
  inline viscousFlux & operator/=(const double &);

  inline viscousFlux operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline viscousFlux operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline viscousFlux operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline viscousFlux operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  friend inline const viscousFlux operator-(const double &lhs, viscousFlux rhs);
  friend inline const viscousFlux operator/(const double &lhs, viscousFlux rhs);

  // destructor
  ~viscousFlux() noexcept {}
};

// function definitions
// operator overload for addition
viscousFlux & viscousFlux::operator+=(const viscousFlux &arr) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    data_[rr] += arr.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
viscousFlux & viscousFlux::operator-=(const viscousFlux &arr) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    data_[rr] -= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
viscousFlux & viscousFlux::operator*=(const viscousFlux &arr) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    data_[rr] *= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
viscousFlux & viscousFlux::operator/=(const viscousFlux &arr) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    data_[rr] /= arr.data_[rr];
  }
  return *this;
}

inline const viscousFlux operator+(viscousFlux lhs, const viscousFlux &rhs) {
  return lhs += rhs;
}

inline const viscousFlux operator-(viscousFlux lhs, const viscousFlux &rhs) {
  return lhs -= rhs;
}

inline const viscousFlux operator*(viscousFlux lhs, const viscousFlux &rhs) {
  return lhs *= rhs;
}

inline const viscousFlux operator/(viscousFlux lhs, const viscousFlux &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
viscousFlux & viscousFlux::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
viscousFlux & viscousFlux::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
viscousFlux & viscousFlux::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
viscousFlux & viscousFlux::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const viscousFlux operator+(const double &lhs, viscousFlux rhs) {
  return rhs += lhs;
}

inline const viscousFlux operator-(const double &lhs, viscousFlux rhs) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    rhs.data_[rr] = lhs - rhs.data_[rr];
  }
  return rhs;
}

inline const viscousFlux operator*(const double &lhs, viscousFlux rhs) {
  return rhs *= lhs;
}

inline const viscousFlux operator/(const double &lhs, viscousFlux rhs) {
  for (auto rr = 0; rr < NUMVARS - 1; rr++) {
    rhs.data_[rr] = lhs / rhs.data_[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const viscousFlux &);

#endif

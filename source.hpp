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

#ifndef SOURCEHEADERDEF  // only if the macro SOURCEHEADERDEF is not defined
                         // execute these lines of code
#define SOURCEHEADERDEF  // define the macro

/* This header contains the source class.

   The source class stores the source terms for the Euler and Navier-Stokes
   equations. */

#include <vector>  // vector
#include <string>  // string
#include <iostream>
#include <memory>
#include "macros.hpp"
#include "vector3d.hpp"
#include "tensor.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declaration
class primVars;
class turbModel;
class sutherland;

class source {
  double data_[NUMVARS];  // source variables at cell center

 public:
  // constructors
  source() : data_{0.0} {}

  // move constructor and assignment operator
  source(source&&) noexcept = default;
  source& operator=(source&&) noexcept = default;

  // copy constructor and assignment operator
  source(const source&) = default;
  source& operator=(const source&) = default;

  // member functions
  double SrcMass() const { return data_[0]; }
  double SrcMomX() const { return data_[1]; }
  double SrcMomY() const { return data_[2]; }
  double SrcMomZ() const { return data_[3]; }
  double SrcEngy() const { return data_[4]; }
  double SrcTke() const { return data_[5]; }
  double SrcOmg() const { return data_[6]; }

  double CalcTurbSrc(const unique_ptr<turbModel> &, const primVars &,
                     const tensor<double> &, const vector3d<double> &,
                     const vector3d<double> &, const vector3d<double> &,
                     const sutherland &, const double &, const double &,
                     const double &);

  inline source & operator+=(const source &);
  inline source & operator-=(const source &);
  inline source & operator*=(const source &);
  inline source & operator/=(const source &);

  inline source & operator+=(const double &);
  inline source & operator-=(const double &);
  inline source & operator*=(const double &);
  inline source & operator/=(const double &);

  inline source operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline source operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline source operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline source operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  friend inline const source operator-(const double &lhs, source rhs);
  friend inline const source operator/(const double &lhs, source rhs);

  // destructor
  ~source() noexcept {}
};

// function definitions -------------------------------------

// operator overload for addition
source & source::operator+=(const source &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] += arr.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
source & source::operator-=(const source &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] -= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
source & source::operator*=(const source &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] *= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
source & source::operator/=(const source &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] /= arr.data_[rr];
  }
  return *this;
}

inline const source operator+(source lhs, const source &rhs) {
  return lhs += rhs;
}

inline const source operator-(source lhs, const source &rhs) {
  return lhs -= rhs;
}

inline const source operator*(source lhs, const source &rhs) {
  return lhs *= rhs;
}

inline const source operator/(source lhs, const source &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
source & source::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
source & source::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
source & source::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
source & source::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const source operator+(const double &lhs, source rhs) {
  return rhs += lhs;
}

inline const source operator-(const double &lhs, source rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs - rhs.data_[rr];
  }
  return rhs;
}

inline const source operator*(const double &lhs, source rhs) {
  return rhs *= lhs;
}

inline const source operator/(const double &lhs, source rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs / rhs.data_[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const source &);

#endif

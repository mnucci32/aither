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

#ifndef ARRAYHEADERDEF  // only if the macro ARRAYHEADERDEF is not defined
                         // execute these lines of code
#define ARRAYHEADERDEF  // define the macro

#include <iostream>
#include "macros.hpp"

using std::ostream;

// forward class declarations
class uncoupledScalar;

/*class to store an array of a fixed size equal to the number of variables
being solved for. This is useful because a vector of these will be
contiguous in memory. */
class genArray {
  double data_[NUMVARS];

 public:
  // constructor
  genArray(const double &a, const double &b, const double &c, const double &d,
           const double &e, const double &f, const double &g)
      : data_{a, b, c, d, e, f, g} {}
  genArray(const double &a, const double &b, const double &c, const double &d,
           const double &e)
      : genArray(a, b, c, d, e, 0.0, 0.0) {}
  genArray() : genArray(0.0, 0.0, 0.0, 0.0, 0.0) {}
  explicit genArray(const double &a) : genArray(a, a, a, a, a, a, a) {}
  explicit genArray(const uncoupledScalar &a);

  // member functions
  void Zero();
  double Sum();
  void SquareRoot();

  // move constructor and assignment operator
  genArray(genArray&&) noexcept = default;
  genArray& operator=(genArray&&) noexcept = default;

  // copy constructor and assignment operator
  genArray(const genArray&) = default;
  genArray& operator=(const genArray&) = default;

  // operator overloads
  const double & operator[](const int &r) const { return data_[r]; }
  double & operator[](const int &r) { return data_[r]; }

  inline genArray & operator+=(const genArray &);
  inline genArray & operator-=(const genArray &);
  inline genArray & operator*=(const genArray &);
  inline genArray & operator/=(const genArray &);

  inline genArray & operator+=(const double &);
  inline genArray & operator-=(const double &);
  inline genArray & operator*=(const double &);
  inline genArray & operator/=(const double &);

  inline genArray operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline genArray operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline genArray operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline genArray operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  void GlobalReduceMPI(const int &, const int &);

  // destructor
  ~genArray() noexcept {}
};

// function declarations --------------------------------------
// operator overload for addition
genArray & genArray::operator+=(const genArray &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] += arr[rr];
  }
  return *this;
}

// operator overload for subtraction
genArray & genArray::operator-=(const genArray &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] -= arr[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
genArray & genArray::operator*=(const genArray &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] *= arr[rr];
  }
  return *this;
}

// operator overload for elementwise division
genArray & genArray::operator/=(const genArray &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] /= arr[rr];
  }
  return *this;
}

inline const genArray operator+(genArray lhs, const genArray &rhs) {
  return lhs += rhs;
}

inline const genArray operator-(genArray lhs, const genArray &rhs) {
  return lhs -= rhs;
}

inline const genArray operator*(genArray lhs, const genArray &rhs) {
  return lhs *= rhs;
}

inline const genArray operator/(genArray lhs, const genArray &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
genArray & genArray::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
genArray & genArray::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
genArray & genArray::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
genArray & genArray::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const genArray operator+(const double &lhs, genArray rhs) {
  return rhs += lhs;
}

inline const genArray operator-(const double &lhs, genArray rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs[rr] = lhs - rhs[rr];
  }
  return rhs;
}

inline const genArray operator*(const double &lhs, genArray rhs) {
  return rhs *= lhs;
}

inline const genArray operator/(const double &lhs, genArray rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs[rr] = lhs / rhs[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const genArray &);

#endif

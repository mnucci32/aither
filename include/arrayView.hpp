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

#ifndef ARRAYVIEWHEADERDEF
#define ARRAYVIEWHEADERDEF

#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include "macros.hpp"
#include "varArray.hpp"

using std::ostream;
using std::vector;
using std::endl;

/* class to store a view of an array. This is useful to slice out data from a
 * std::vector.
 */

// forward class declarations
class primitive;
class conserved;
class residual;

template <typename T1, typename T2>
class arrayView {
  static_assert(std::is_base_of<varArray, T1>::value,
                "arrayView<T1, T2> requires T1 to be varArray type!");
  static_assert(std::is_arithmetic<T2>::value,
                "arrayView<T1, T2> requires T2 to be an arithmetic type!");

  typename vector<const T2>::iterator begin_;
  typename vector<const T2>::iterator end_;
  int momentumIndex_;
  int energyIndex_;
  int turbulenceIndex_;

 public:
  // constructor
  arrayView(const typename vector<const T2>::iterator &b,
            const typename vector<const T2>::iterator &e, const int &numSpecies)
      : begin_(b),
        end_(e),
        momentumIndex_(numSpecies),
        energyIndex_(momentumIndex_ + 3),
        turbulenceIndex_(energyIndex_ + 1) {}

  // member functions
  T2 Sum() const { return std::accumulate(begin_, end_, T2(0)); }
  auto Size() const { return std::distance(begin_, end_); }
  int NumSpecies() const { return momentumIndex_; }
  T1 CopyData() const { return T1{begin_, end_, this->NumSpecies()}; }

  // move constructor and assignment operator
  arrayView(arrayView&&) noexcept = default;
  arrayView& operator=(arrayView&&) noexcept = default;

  // copy constructor and assignment operator
  arrayView(const arrayView &) = default;
  arrayView& operator=(const arrayView &assign) {
    std::copy(assign.begin_, assign.end_, begin_);
  }

  // operator overloads
  const T2 & operator[](const int &r) const { return *(begin_ + r); }

  inline T1 operator+(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs += s;
  }
  inline T1 operator-(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs -= s;
  }
  inline T1 operator*(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs *= s;
  }
  inline T1 operator/(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs /= s;
  }

  // destructor
  ~arrayView() noexcept {}
};

// function declarations --------------------------------------
template <typename T1, typename T2>
inline const T1 operator+(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll += rhs;
}

template <typename T1, typename T2>
inline const T1 operator-(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll -= rhs;
}

template <typename T1, typename T2>
inline const T1 operator*(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll *= rhs;
}

template <typename T1, typename T2>
inline const T1 operator/(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll /= rhs;
}

// operator overloads for type T -------------------------------------
template <typename T1, typename T2>
inline const T1 operator+(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return result += lhs;
}

template <typename T1, typename T2>
inline const T1 operator-(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return lhs - result;
}

template <typename T1, typename T2>
inline const T1 operator*(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return result *= lhs;
}

template <typename T1, typename T2>
inline const T1 operator/(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return lhs / rhs;
}

// operation overload for << - allows use of cout, cerr, etc.
template <typename T1, typename T2>
ostream &operator<<(ostream &os, const arrayView<T1, T2> &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// using typedefs
using primitiveView = arrayView<primitive, double>;
using conservedView = arrayView<conserved, double>;
using residualView = arrayView<residual, double>;

#endif

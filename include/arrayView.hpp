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
#include <cmath>
#include "macros.hpp"

using std::ostream;
using std::vector;
using std::endl;

/* class to store a view of an array. This is useful to slice out data from a
 * std::vector.
 */

template <typename T>
class arrayView {
  static_assert(std::is_arithmetic<T>::value,
                "arrayView<T> requires an arithmetic type!");

  typename vector<T>::iterator begin_;
  typename vector<T>::iterator end_;
  int momentumIndex_;
  int energyIndex_;
  int turbulenceIndex_;

 public:
  // constructor
  arrayView(const typename vector<T>::iterator &b,
            const typename vector<T>::iterator &e, const int &numSpecies)
      : begin_(b),
        end_(e),
        momentumIndex_(numSpecies),
        energyIndex_(momentumIndex_ + 3),
        turbulenceIndex_(energyIndex_ + 1) {}

  // member functions
  void Zero() { std::fill(begin_, end_, T(0)); }
  T Sum() const { return std::accumulate(begin_, end_, T(0)); }
  void SquareRoot() {
    std::for_each(begin_, end_, [](T &val) { val = sqrt(val); });
  }
  auto Size() const { return std::distance(begin_, end_); }

  // move constructor and assignment operator
  arrayView(arrayView&&) noexcept = default;
  arrayView& operator=(arrayView&&) noexcept = default;

  // copy constructor and assignment operator
  arrayView(const arrayView &) = default;
  arrayView& operator=(const arrayView &assign) {
    std::copy(assign.begin_, assign.end_, begin_);
  }

  // operator overloads
  const T & operator[](const int &r) const { return *(begin_ + r); }
  T & operator[](const int &r) { return *(begin_ + r); }

  inline arrayView & operator+=(const arrayView &);
  inline arrayView & operator-=(const arrayView &);
  inline arrayView & operator*=(const arrayView &);
  inline arrayView & operator/=(const arrayView &);

  inline arrayView & operator+=(const T &);
  inline arrayView & operator-=(const T &);
  inline arrayView & operator*=(const T &);
  inline arrayView & operator/=(const T &);

  inline arrayView operator+(const T &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline arrayView operator-(const T &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline arrayView operator*(const T &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline arrayView operator/(const T &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // destructor
  ~arrayView() noexcept {}
};

// function declarations --------------------------------------
// operator overload for addition
template <typename T>
arrayView<T> & arrayView<T>::operator+=(const arrayView<T> &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "arrayViews must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    (*this)[rr] += arr[rr];
  }
  return *this;
}

// operator overload for subtraction
template <typename T>
arrayView<T> & arrayView<T>::operator-=(const arrayView<T> &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "arrayViews must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    (*this)[rr] -= arr[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
arrayView<T> & arrayView<T>::operator*=(const arrayView<T> &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "arrayViews must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    (*this)[rr] *= arr[rr];
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
arrayView<T> & arrayView<T>::operator/=(const arrayView<T> &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "arrayViews must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    (*this)[rr] /= arr[rr];
  }
  return *this;
}

template <typename T>
inline const arrayView<T> operator+(arrayView<T> lhs, const arrayView<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const arrayView<T> operator-(arrayView<T> lhs, const arrayView<T> &rhs) {
  return lhs -= rhs;
}

template <typename T>
inline const arrayView<T> operator*(arrayView<T> lhs, const arrayView<T> &rhs) {
  return lhs *= rhs;
}

template <typename T>
inline const arrayView<T> operator/(arrayView<T> lhs, const arrayView<T> &rhs) {
  return lhs /= rhs;
}

// operator overloads for type T -------------------------------------
// operator overload for addition
template <typename T>
arrayView<T> & arrayView<T>::operator+=(const T &scalar) {
  std::for_each(begin_, end_, [&](T &val) { val += scalar; });
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
arrayView<T> & arrayView<T>::operator-=(const T &scalar) {
  std::for_each(begin_, end_, [&](T &val) { val -= scalar; });
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
arrayView<T> & arrayView<T>::operator*=(const T &scalar) {
  std::for_each(begin_, end_, [&](T &val) { val *= scalar; });
  return *this;
}

// operator overload for elementwise division
template <typename T>
arrayView<T> & arrayView<T>::operator/=(const T &scalar) {
  std::for_each(begin_, end_, [&](T &val) { val /= scalar; });
  return *this;
}

template <typename T>
inline const arrayView<T> operator+(const T &lhs, arrayView<T> rhs) {
  return rhs += lhs;
}

template <typename T>
inline const arrayView<T> operator-(const T &lhs, arrayView<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs[rr] = lhs - rhs[rr];
  }
  return rhs;
}

template <typename T>
inline const arrayView<T> operator*(const T &lhs, arrayView<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const arrayView<T> operator/(const T &lhs, arrayView<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs[rr] = lhs / rhs[rr];
  }
  return rhs;
}

// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const arrayView<T> &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}


#endif

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

#ifndef MATRIXVIEWHEADERDEF
#define MATRIXVIEWHEADERDEF

#include <iostream>
#include <vector>
#include "macros.hpp"
#include "matrix.hpp"

using std::ostream;
using std::vector;

// forward class declaration
class varArray;

// class to store a square matrix
class squareMatrixView {
  int size_;
  typename vector<double>::iterator begin_;
  typename vector<double>::iterator end_;

  // private member functions
  int GetLoc(const int &r, const int &c) const {
    return r * size_ + c;
  }

 public:
  // constructor
  squareMatrixView(const typename vector<double>::iterator &b,
                   const typename vector<double>::iterator &e,
                   const int &size)
      : size_(size), begin_(b), end_(e) {}

  // move constructor and assignment operator
  squareMatrixView(squareMatrixView &&) noexcept = default;
  squareMatrixView& operator=(squareMatrixView &&) = default;

  // copy constructor and assignment operator
  squareMatrixView(const squareMatrixView &) = default;
  squareMatrixView& operator=(const squareMatrixView &) = default;

  // member functions
  int Size() const {return size_;}
  void SwapRows(const int &, const int &);
  void Inverse();
  int FindMaxInCol(const int &, const int &, const int &) const;
  void RowMultiply(const int &, const int &, const double &);
  void LinCombRow(const int &, const double &, const int &);
  void Zero();
  void Identity();
  template <typename T>
  squareMatrix MatMult(const T &) const;
  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
  T ArrayMult(const T &, const int = 0) const;
  double MaxAbsValOnDiagonal() const;

  // provide begin and end so std::begin and std::end can be used
  // use lower case to conform with std::begin, std::end
  auto begin() noexcept {return begin_;}
  const auto begin() const noexcept {return begin_;}
  auto end() noexcept {return end_;}
  const auto end() const noexcept {return end_;}

  // operator overloads
  double & operator()(const int &r, const int &c) {
    return *(begin_ + this->GetLoc(r, c));
  }
  const double & operator()(const int &r, const int &c) const {
    return *(begin_ + this->GetLoc(r, c));
  }

  inline squareMatrixView & operator+=(const squareMatrixView &);
  inline squareMatrixView & operator-=(const squareMatrixView &);
  inline squareMatrixView & operator*=(const squareMatrixView &);
  inline squareMatrixView & operator/=(const squareMatrixView &);

  inline squareMatrixView & operator+=(const double &);
  inline squareMatrixView & operator-=(const double &);
  inline squareMatrixView & operator*=(const double &);
  inline squareMatrixView & operator/=(const double &);

  inline squareMatrixView operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline squareMatrixView operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline squareMatrixView operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline squareMatrixView operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // destructor
  ~squareMatrixView() noexcept {}
};

// ---------------------------------------------------------------------------
// function declarations
// matrix-matrix multiplication
// using cache efficient implimentation
template <typename T>
squareMatrix squareMatrixView::MatMult(const T &s2) const {
  static_assert(std::is_same<T, squareMatrix>::value ||
                    std::is_same<T, squareMatrixView>::value,
                "T should be squareMatrix type");
  squareMatrix s1(s2.Size());
  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      for (auto ii = 0; ii < s2.Size(); ii++) {
        s1(rr, ii) += (*this)(rr, cc) * s2(cc, ii);
      }
    }
  }
  return s1;
}

// member function to do matrix/vector multplication with varArray type
template <typename T, typename TT>
T squareMatrixView::ArrayMult(const T &vec, const int pos) const {
  // vec -- vector to multiply with
  auto product = vec;

  // zero out portion of array that will be written over
  if (pos == 0) {
    for (auto ii = 0; ii < vec.TurbulenceIndex(); ii++) {
      product[ii] = 0.0;
    }
  } else {
    for (auto ii = pos; ii < vec.Size(); ii++) {
      product[ii] = 0.0;
    }
  }

  for (auto rr = 0; rr < size_; rr++) {
    for (auto cc = 0; cc < size_; cc++) {
      product[pos + rr] += (*this)(rr, cc) * vec[pos + cc];
    }
  }
  return product;
}

ostream &operator<<(ostream &os, const squareMatrixView &);

// operator overload for addition
squareMatrixView & squareMatrixView::operator+=(const squareMatrixView &mat) {
  std::transform(begin_, end_, mat.begin_, begin_,
                 [](auto &lhs, auto &rhs) { return lhs + rhs; });
  return *this;
}

// operator overload for subtraction
squareMatrixView & squareMatrixView::operator-=(const squareMatrixView &mat) {
  std::transform(begin_, end_, mat.begin_, begin_,
                 [](auto &lhs, auto &rhs) { return lhs - rhs; });
  return *this;
}

// operator overload for elementwise multiplication
squareMatrixView & squareMatrixView::operator*=(const squareMatrixView &mat) {
  std::transform(begin_, end_, mat.begin_, begin_,
                 [](auto &lhs, auto &rhs) { return lhs * rhs; });
  return *this;
}

// operator overload for elementwise multiplication
squareMatrixView & squareMatrixView::operator/=(const squareMatrixView &mat) {
  std::transform(begin_, end_, mat.begin_, begin_,
                 [](auto &lhs, auto &rhs) { return lhs / rhs; });
  return *this;
}

inline const squareMatrixView operator+(squareMatrixView lhs,
                                        const squareMatrixView &rhs) {
  return lhs += rhs;
}

inline const squareMatrixView operator-(squareMatrixView lhs,
                                        const squareMatrixView &rhs) {
  return lhs -= rhs;
}

inline const squareMatrixView operator*(squareMatrixView lhs,
                                        const squareMatrixView &rhs) {
  return lhs *= rhs;
}

inline const squareMatrixView operator/(squareMatrixView lhs,
                                        const squareMatrixView &rhs) {
  return lhs /= rhs;
}

// operator overloads for double --------------------------------------------
// operator overload for addition
squareMatrixView & squareMatrixView::operator+=(const double &scalar) {
  for_each(begin_, end_, [&scalar](auto &val) { val += scalar; });
  return *this;
}

// operator overload for subtraction
squareMatrixView & squareMatrixView::operator-=(const double &scalar) {
  for_each(begin_, end_, [&scalar](auto &val) { val -= scalar; });
  return *this;
}

// operator overload for multiplication
squareMatrixView & squareMatrixView::operator*=(const double &scalar) {
  for_each(begin_, end_, [&scalar](auto &val) { val *= scalar; });
  return *this;
}

// operator overload for division
squareMatrixView & squareMatrixView::operator/=(const double &scalar) {
  for_each(begin_, end_, [&scalar](auto &val) { val /= scalar; });
  return *this;
}

inline const squareMatrixView operator+(const double &lhs,
                                        squareMatrixView rhs) {
  return rhs += lhs;
}

inline const squareMatrixView operator-(const double &lhs,
                                        squareMatrixView rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    for (auto cc = 0; cc < rhs.Size(); cc++) {
      rhs(rr, cc) = lhs - rhs(rr, cc);
    }
  }
  return rhs;
}

inline const squareMatrixView operator*(const double &lhs,
                                        squareMatrixView rhs) {
  return rhs *= lhs;
}

inline const squareMatrixView operator/(const double &lhs,
                                        squareMatrixView rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    for (auto cc = 0; cc < rhs.Size(); cc++) {
      rhs(rr, cc) = lhs / rhs(rr, cc);
    }
  }
  return rhs;
}

#endif

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

#ifndef MATRIXHEADERDEF  // only if the macro MATRIXHEADERDEF is not defined
                         // execute these lines of code
#define MATRIXHEADERDEF  // define the macro

#include <iostream>
#include <vector>
#include "macros.hpp"

using std::ostream;
using std::vector;

// forward class declarations
class genArray;

// class to store a square matrix
class squareMatrix {
  int size_;
  vector<double> data_;

  // private member functions
  int GetLoc(const int &r, const int &c) const {
    return r * size_ + c;
  }

 public:
  // constructor
  explicit squareMatrix(const int &a) : size_(a), data_(a * a, 0.0) {}
  squareMatrix() : squareMatrix(0) {}

  // move constructor and assignment operator
  squareMatrix(squareMatrix &&) noexcept = default;
  squareMatrix& operator=(squareMatrix &&) = default;

  // copy constructor and assignment operator
  squareMatrix(const squareMatrix &) = default;
  squareMatrix& operator=(const squareMatrix &) = default;

  // member functions
  int Size() const {return size_;}
  void SwapRows(const int &, const int &);
  void Inverse();
  int FindMaxInCol(const int &, const int &, const int &) const;
  void RowMultiply(const int &, const int &, const double &);
  void LinCombRow(const int &, const double &, const int &);
  void Zero();
  void Identity();
  squareMatrix MatMult(const squareMatrix &) const;
  genArray ArrayMult(const genArray &, const int = 0) const;
  double MaxAbsValOnDiagonal() const;

  // operator overloads
  double & operator()(const int &r, const int &c) {
    return data_[this->GetLoc(r, c)];
  }
  const double & operator()(const int &r, const int &c) const {
    return data_[this->GetLoc(r, c)];
  }

  inline squareMatrix & operator+=(const squareMatrix &);
  inline squareMatrix & operator-=(const squareMatrix &);
  inline squareMatrix & operator*=(const squareMatrix &);
  inline squareMatrix & operator/=(const squareMatrix &);

  inline squareMatrix & operator+=(const double &);
  inline squareMatrix & operator-=(const double &);
  inline squareMatrix & operator*=(const double &);
  inline squareMatrix & operator/=(const double &);

  inline squareMatrix operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline squareMatrix operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline squareMatrix operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline squareMatrix operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // destructor
  ~squareMatrix() noexcept {}
};


// function declarations
ostream &operator<<(ostream &os, const squareMatrix &);

// operator overload for addition
squareMatrix & squareMatrix::operator+=(const squareMatrix &mat) {
  for (auto ii = 0U; ii < mat.data_.size(); ii++) {
    data_[ii] += mat.data_[ii];
  }
  return *this;
}

// operator overload for subtraction
squareMatrix & squareMatrix::operator-=(const squareMatrix &mat) {
  for (auto ii = 0U; ii < mat.data_.size(); ii++) {
    data_[ii] -= mat.data_[ii];
  }
  return *this;
}

// operator overload for elementwise multiplication
squareMatrix & squareMatrix::operator*=(const squareMatrix &mat) {
  for (auto ii = 0U; ii < mat.data_.size(); ii++) {
    data_[ii] *= mat.data_[ii];
  }
  return *this;
}

// operator overload for elementwise multiplication
squareMatrix & squareMatrix::operator/=(const squareMatrix &mat) {
  for (auto ii = 0U; ii < mat.data_.size(); ii++) {
    data_[ii] /= mat.data_[ii];
  }
  return *this;
}

inline const squareMatrix operator+(squareMatrix lhs, const squareMatrix &rhs) {
  return lhs += rhs;
}

inline const squareMatrix operator-(squareMatrix lhs, const squareMatrix &rhs) {
  return lhs -= rhs;
}

inline const squareMatrix operator*(squareMatrix lhs, const squareMatrix &rhs) {
  return lhs *= rhs;
}

inline const squareMatrix operator/(squareMatrix lhs, const squareMatrix &rhs) {
  return lhs /= rhs;
}

// operator overloads for double --------------------------------------------
// operator overload for addition
squareMatrix & squareMatrix::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction
squareMatrix & squareMatrix::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for multiplication
squareMatrix & squareMatrix::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for division
squareMatrix & squareMatrix::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const squareMatrix operator+(const double &lhs, squareMatrix rhs) {
  return rhs += lhs;
}

inline const squareMatrix operator-(const double &lhs, squareMatrix rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    for (auto cc = 0; cc < rhs.Size(); cc++) {
      rhs(rr, cc) = lhs - rhs(rr, cc);
    }
  }
  return rhs;
}

inline const squareMatrix operator*(const double &lhs, squareMatrix rhs) {
  return rhs *= lhs;
}

inline const squareMatrix operator/(const double &lhs, squareMatrix rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    for (auto cc = 0; cc < rhs.Size(); cc++) {
      rhs(rr, cc) = lhs / rhs(rr, cc);
    }
  }
  return rhs;
}


#endif

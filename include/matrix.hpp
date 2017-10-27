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

#ifndef MATRIXHEADERDEF  // only if the macro MATRIXHEADERDEF is not defined
                         // execute these lines of code
#define MATRIXHEADERDEF  // define the macro

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "macros.hpp"

using std::ostream;
using std::vector;

// forward class declaration
class varArray;

// ---------------------------------------------------------------------------
// matrix functions

// function for matrix multiplication
// using cache efficient implimentation
void MatrixMultiply(const vector<double>::const_iterator &matL,
                    const vector<double>::const_iterator &matR,
                    const vector<double>::iterator &result, const int &size);

double MaximumAbsValOnDiagonal(const vector<double>::const_iterator &mat,
                               const int &size);

void IdentityMatrix(const vector<double>::iterator &mat, const int &size);

int FindMaxInColumn(const vector<double>::const_iterator &mat, const int &size,
                    const int &c, const int &start, const int &end);

void RowMultiplyFactor(const vector<double>::iterator &mat, const int &size,
                       const int &r, const int &c, const double &factor);

void LinearCombRow(const vector<double>::iterator &mat, const int &size,
                   const int &r1, const double &factor, const int &r2);

void SwapMatRows(const vector<double>::iterator &mat, const int &size,
                 const int &r1, const int &r2);

void MatrixInverse(const vector<double>::iterator &mat, const int &size);

void MultiplyFacOnDiagonal(const vector<double>::iterator &mat, const int &size,
                           const double &val);

void AddFacOnDiagonal(const vector<double>::iterator &mat, const int &size,
                      const double &val);

// ---------------------------------------------------------------------------                      

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
  // provide begin and end so std::begin and std::end can be used
  // use lower case to conform with std::begin, std::end
  auto begin() noexcept {return data_.begin();}
  const auto begin() const noexcept {return data_.begin();}
  auto end() noexcept {return data_.end();}
  const auto end() const noexcept {return data_.end();}

  int Size() const {return size_;}
  void SwapRows(const int &r1, const int &r2) {
    SwapMatRows(this->begin(), size_, r1, r2);
  }
  void Inverse() { MatrixInverse(this->begin(), size_); }
  int FindMaxInCol(const int &c, const int &start, const int &end) const {
    return FindMaxInColumn(this->begin(), size_, c, start, end);
  }
  void RowMultiply(const int &r, const int &c, const double &fac) {
    RowMultiplyFactor(this->begin(), size_, r, c, fac);
  }
  void LinCombRow(const int &r1, const double &fac, const int &r2) {
    LinearCombRow(this->begin(), size_, r1, fac, r2);
  }
  void Zero() { std::fill(this->begin(), this->end(), 0.0); }
  void Identity() { IdentityMatrix(this->begin(), size_); }
  squareMatrix MatMult(const squareMatrix &) const;
  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
  T ArrayMult(const T &, const int = 0) const;
  double MaxAbsValOnDiagonal() const {
    return MaximumAbsValOnDiagonal(this->begin(), size_);
  }

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
// member function to do matrix/vector multplication with varArray type
template <typename T, typename TT>
T squareMatrix::ArrayMult(const T &vec, const int pos) const {
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

ostream &operator<<(ostream &os, const squareMatrix &);

// operator overload for addition
squareMatrix & squareMatrix::operator+=(const squareMatrix &mat) {
  MSG_ASSERT(this->Size() == mat.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), mat.begin(), this->begin(),
                 std::plus<double>());
  return *this;
}

// operator overload for subtraction
squareMatrix & squareMatrix::operator-=(const squareMatrix &mat) {
  MSG_ASSERT(this->Size() == mat.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), mat.begin(), this->begin(),
                 std::minus<double>());
  return *this;
}

// operator overload for elementwise multiplication
squareMatrix & squareMatrix::operator*=(const squareMatrix &mat) {
  MSG_ASSERT(this->Size() == mat.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), mat.begin(), this->begin(),
                 std::multiplies<double>());
  return *this;
}

// operator overload for elementwise multiplication
squareMatrix & squareMatrix::operator/=(const squareMatrix &mat) {
  MSG_ASSERT(this->Size() == mat.Size(), "matrix sizes must be equal");
  std::transform(this->begin(), this->end(), mat.begin(), this->begin(),
                 std::divides<double>());
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
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val += scalar; });
  return *this;
}

// operator overload for subtraction
squareMatrix & squareMatrix::operator-=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val -= scalar; });
  return *this;
}

// operator overload for multiplication
squareMatrix & squareMatrix::operator*=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val *= scalar; });
  return *this;
}

// operator overload for division
squareMatrix & squareMatrix::operator/=(const double &scalar) {
  for_each(this->begin(), this->end(), [&scalar](auto &val) { val /= scalar; });
  return *this;
}

inline const squareMatrix operator+(const double &lhs, squareMatrix rhs) {
  return rhs += lhs;
}

inline const squareMatrix operator-(const double &lhs, squareMatrix rhs) {
  for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs - val; });
  return rhs;
}

inline const squareMatrix operator*(const double &lhs, squareMatrix rhs) {
  return rhs *= lhs;
}

inline const squareMatrix operator/(const double &lhs, squareMatrix rhs) {
  for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs / val; });
  return rhs;
}


#endif

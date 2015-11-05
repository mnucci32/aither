/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <cstdlib>  // exit()
#include <iostream>  // cout
#include <cmath>
#include <vector>
#include <algorithm>  // swap
#include "matrix.hpp"
#include "plot3d.hpp"  // matrix location functions

using std::cout;
using std::endl;
using std::cerr;
using std::copy;
using std::swap_ranges;
using std::fabs;

// copy constructor
squareMatrix::squareMatrix(const squareMatrix &cp) {
  size_ = cp.Size();
  data_ = new double[cp.Size() * cp.Size()];
  copy(&cp.data_[0], &cp.data_[0] + cp.Size() * cp.Size(), &data_[0]);
}

// copy assignment operator
squareMatrix &squareMatrix::operator=(squareMatrix other) {
  swap(*this, other);
  return *this;
}

// friend function to allow for swap functionality
void swap(squareMatrix &first, squareMatrix &second) noexcept {
  std::swap(first.size_, second.size_);
  std::swap(first.data_, second.data_);
}

// member function to get the data_ from the matrix
double squareMatrix::Data(const int &r, const int &c) const {
  // test to see that row and column inputs are within bounds
  if ((r >= (size_)) || (c >= (size_))) {
    cerr << "ERROR: The requested data_, does not lie within the matrix "
            "bounds. Check row and column inputs." << endl;
    exit(1);
  }
  return data_[c + r * size_];
}

// member function to set the data_ in the matrix
void squareMatrix::SetData(const int &r, const int &c, const double &val) {
  // test to see that row and column inputs are within bounds
  if ((r >= size_) || (c >= size_)) {
    cerr << "ERROR: Cannot assign data_ to given location because it does not "
            "lie within the matrix bounds. Check row and column inputs."
         << endl;
    exit(1);
  }
  data_[c + r * size_] = val;
}

// member function to swap rows of matrix
void squareMatrix::SwapRows(const int &r1, const int &r2) {
  if (r1 != r2) {
    for (auto ii = 0; ii < size_; ii++) {
      auto ind1 = ii + r1 * size_;
      auto ind2 = ii + r2 * size_;
      std::swap(data_[ind1], data_[ind2]);
    }
  }
}

// member function to invert matrix using Gauss-Jordan elimination
void squareMatrix::Inverse() {
  squareMatrix I(size_);
  auto r = 0;

  I.Identity();

  auto cPivot = 0;
  for (cPivot = 0, r = 0; r < size_; r++, cPivot++) {
    // find pivot row
    auto rPivot = this->FindMaxInCol(r, cPivot, size_ - 1);

    // swap rows
    this->SwapRows(r, rPivot);
    I.SwapRows(r, rPivot);

    if (r != 0) {  // if not on first row, need to get rid entries ahead of
                   // pivot
      for (auto ii = 0; ii < cPivot; ii++) {
        auto factor = this->Data(r, ii) / this->Data(ii, ii);
        this->LinCombRow(ii, factor, r);
        I.LinCombRow(ii, factor, r);
      }
    }

    // normalize row by pivot
    if (this->Data(r, cPivot) == 0.0) {
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination! Matrix (mid "
              "inversion) is" << endl << *this << endl;
      exit(1);
    }
    auto normFactor = 1.0 / this->Data(r, cPivot);
    this->RowMultiply(r, cPivot, normFactor);  // only multiply entries from
                                                 // pivot and to the right
    I.RowMultiply(r, 0, normFactor);  // multiply all entries
  }

  // matrix is now upper triangular, work way back up to identity matrix
  cPivot = size_ - 2;  // start with second to last row
  for (cPivot = size_ - 2, r = size_ - 2; r >= 0; r--, cPivot--) {
    for (auto ii = size_ - 1; ii > cPivot; ii--) {
      auto factor = this->Data(r, ii);
      this->LinCombRow(ii, factor, r);
      I.LinCombRow(ii, factor, r);
    }
  }

  // set this matrix equal to its inverse
  (*this) = I;
}

// member function to add a linear combination of one row to another
void squareMatrix::LinCombRow(const int &r1, const double &factor,
                              const int &r2) {
  for (auto ii = 0; ii < size_; ii++) {
    this->SetData(r2, ii, this->Data(r2, ii) - this->Data(r1, ii) * factor);
  }
}

// member function to multiply a row by a given factor
void squareMatrix::RowMultiply(const int &r, const int &c,
                               const double &factor) {
  for (auto ii = c; ii < size_; ii++) {
    this->SetData(r, ii, this->Data(r, ii) * factor);
  }
}

// member function to find maximum absolute value in a given column and range
// within that column and return the corresponding row indice
int squareMatrix::FindMaxInCol(const int &c, const int &start,
                               const int &end) const {
  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii < end + 1; ii++) {
    if (fabs(this->Data(ii, c)) > maxVal) {
      maxVal = fabs(this->Data(ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}

// operator overload for addition
squareMatrix squareMatrix::operator+(const squareMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot add matrices, dimensions do not agree." << endl;
  }

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) + s2.Data(rr, cc));
    }
  }
  return s1;
}

// operator overload for addition with a scalar
squareMatrix squareMatrix::operator+(const double &scalar) const {
  auto s1 = *this;

  for (auto cc = 0; cc < s1.Size(); cc++) {
    for (auto rr = 0; rr < s1.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) + scalar);
    }
  }
  return s1;
}

// operator overload for addition with a scalar
squareMatrix operator+(const double &scalar, const squareMatrix &s2) {
  squareMatrix s1(s2.Size());

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, s2.Data(rr, cc) + scalar);
    }
  }
  return s1;
}

// operator overload for subtraction
squareMatrix squareMatrix::operator-(const squareMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot subtract matrices, dimensions do not agree." << endl;
  }

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) - s2.Data(rr, cc));
    }
  }
  return s1;
}

// operator overload for subtraction with a scalar
squareMatrix squareMatrix::operator-(const double &scalar) const {
  auto s1 = *this;

  for (auto cc = 0; cc < s1.Size(); cc++) {
    for (auto rr = 0; rr < s1.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) - scalar);
    }
  }
  return s1;
}

// operator overload for subtraction with a scalar
squareMatrix operator-(const double &scalar, const squareMatrix &s2) {
  squareMatrix s1(s2.Size());

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, scalar - s2.Data(rr, cc));
    }
  }
  return s1;
}

// operator overload for multiplication
squareMatrix squareMatrix::operator*(const squareMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot multiply matrices, dimensions do not agree." << endl;
  }

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      auto newVal = 0.0;
      for (auto ii = 0; ii < s2.Size(); ii++) {
        newVal += (this->Data(rr, ii) * s2.Data(ii, cc));
      }
      s1.SetData(rr, cc, newVal);
    }
  }
  return s1;
}

// operator overload for multiplication with a scalar
squareMatrix squareMatrix::operator*(const double &scalar) const {
  auto s1 = *this;

  for (auto cc = 0; cc < s1.Size(); cc++) {
    for (auto rr = 0; rr < s1.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) * scalar);
    }
  }
  return s1;
}

// operator overload for multiplication with a scalar
squareMatrix operator*(const double &scalar, const squareMatrix &s2) {
  squareMatrix s1(s2.Size());

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, s2.Data(rr, cc) * scalar);
    }
  }
  return s1;
}

// operator overload for division with a scalar
squareMatrix squareMatrix::operator/(const double &scalar) const {
  auto s1 = *this;

  for (auto cc = 0; cc < s1.Size(); cc++) {
    for (auto rr = 0; rr < s1.Size(); rr++) {
      s1.SetData(rr, cc, s1.Data(rr, cc) / scalar);
    }
  }
  return s1;
}

// operator overload for division with a scalar
squareMatrix operator/(const double &scalar, const squareMatrix &s2) {
  squareMatrix s1(s2.Size());

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      s1.SetData(rr, cc, scalar / s2.Data(rr, cc));
    }
  }
  return s1;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const squareMatrix &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    for (auto cc = 0; cc < m.Size(); cc++) {
      cout << m.Data(rr, cc);
      if (cc != (m.Size() - 1)) {
        cout << ", ";
      } else {
        cout << endl;
      }
    }
  }
  return os;
}

// member function to zero the matrix
void squareMatrix::Zero() {
  for (auto cc = 0; cc < size_; cc++) {
    for (auto rr = 0; rr < size_; rr++) {
      this->SetData(rr, cc, 0.0);
    }
  }
}

// member function to set matrix to Identity
void squareMatrix::Identity() {
  for (auto rr = 0; rr < this->Size(); rr++) {
    for (auto cc = 0; cc < this->Size(); cc++) {
      if (rr == cc) {
        this->SetData(rr, cc, 1.0);
      } else {
        this->SetData(rr, cc, 0.0);
      }
    }
  }
}

colMatrix squareMatrix::Multiply(const colMatrix &X) const {
  // Test to see that column matrix can be multiplied with square matrix
  if (this->Size() != X.Size()) {
    cerr << "ERROR: Column matrix cannot be multiplied with square matrix. "
            "Sizes do not agree!" << endl;
    exit(1);
  }

  colMatrix B(X.Size());
  B.Zero();

  for (auto rr = 0; rr < X.Size(); rr++) {
    auto tempData = 0.0;
    for (auto cc = 0; cc < X.Size(); cc++) {
      tempData += this->Data(rr, cc) * X.Data(cc);
    }
    B.SetData(rr, tempData);
  }

  return B;
}

// -------------------------------------------------------------
// function definitions for matrixDiagonal class

// copy constructor
matrixDiagonal::matrixDiagonal(const matrixDiagonal &cp) {
  size_ = cp.Size();
  data_ = new squareMatrix[cp.Size()];
  copy(&cp.data_[0], &cp.data_[0] + cp.Size(), &data_[0]);
}

// copy assignment operator
matrixDiagonal &matrixDiagonal::operator=(matrixDiagonal other) {
  swap(*this, other);
  return *this;
}

// friend function to allow for swap functionality
void swap(matrixDiagonal &first, matrixDiagonal &second) noexcept {
  std::swap(first.size_, second.size_);
  std::swap(first.data_, second.data_);
}

// member function to get the data_ from the matrix
squareMatrix matrixDiagonal::Data(const int &ind) const {
  // test to see that the index input is within bounds
  if (ind >= size_) {
    cerr << "ERROR: The requested data_, does not lie within the matrix "
            "bounds. Check index input." << endl;
    exit(1);
  }
  return data_[ind];
}

// member function to set the data_ in the matrix
void matrixDiagonal::SetData(const int &ind, const squareMatrix &val) {
  // test to see that the index input is within bounds
  if (ind >= size_) {
    cerr << "ERROR: The requested data_, does not lie within the matrix "
            "bounds. Check index input." << endl;
    exit(1);
  }
  data_[ind] = val;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const matrixDiagonal &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    cout << "In index " << rr << " matrix block is:" << endl;
    cout << m.Data(rr) << endl;
  }
  return os;
}

// member function to zero the matrix
void matrixDiagonal::Zero(const int &s) {
  squareMatrix mZero(s);
  mZero.Zero();

  for (auto cc = 0; cc < size_; cc++) {
    this->SetData(cc, mZero);
  }
}

// member function to delete the contents of the data_ structure and resize it
void matrixDiagonal::CleanResizeZero(const int &s, const int &m) {
  squareMatrix mZero(m);
  mZero.Zero();

  delete[] data_;
  data_ = new squareMatrix[s];
  size_ = s;

  for (auto cc = 0; cc < size_; cc++) {
    this->SetData(cc, mZero);
  }
}

// member function to invert each matrix stored
void matrixDiagonal::Inverse() {
  for (auto ii = 0; ii < this->Size(); ii++) {
    squareMatrix temp = this->Data(ii);
    temp.Inverse();
    this->SetData(ii, temp);
  }
}

// functions
// ------------------------------------------------------------------------

// functions for colMatrix class
// copy constructor
colMatrix::colMatrix(const colMatrix &cp) {
  size_ = cp.Size();
  data_ = new double[cp.Size()];
  copy(&cp.data_[0], &cp.data_[0] + cp.Size(), &data_[0]);
}

// copy assignment operator
colMatrix &colMatrix::operator=(colMatrix other) {
  swap(*this, other);
  return *this;
}

// friend function to allow for swap functionality
void swap(colMatrix &first, colMatrix &second) noexcept {
  std::swap(first.size_, second.size_);
  std::swap(first.data_, second.data_);
}

// member function to get the data_ from the matrix
double colMatrix::Data(const int &r) const {
  // test to see that row and column inputs are within bounds
  if (r >= size_) {
    cerr << "ERROR: The requested data_, does not lie within the column matrix "
            "bounds. Check row input." << endl;
    exit(1);
  }
  return data_[r];
}

// member function to set the data_ in the matrix
void colMatrix::SetData(const int &r, const double &val) {
  // test to see that row and column inputs are within bounds
  if (r >= size_) {
    cerr << "ERROR: Cannot assign data_ to given location because it does not "
            "lie within the column matrix bounds. Check row input." << endl;
    exit(1);
  }
  data_[r] = val;
}

// operator overload for addition
colMatrix colMatrix::operator+(const colMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot add column matrices, dimensions do not agree."
         << endl;
    cerr << "Dimensions are " << s1.size_ << " and " << s2.size_ << endl;
  }

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) + s2.Data(rr));
  }
  return s1;
}

// operator overload for addition
colMatrix colMatrix::operator+(const vector<double> &v1) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != static_cast<int> (v1.size())) {
    cerr << "ERROR: Cannot add column matrix of size_ " << s1.size_
         << " to vector class of size_ " << v1.size()
         << ", dimensions do not agree." << endl;
  }

  for (auto rr = 0; rr < s1.size_; rr++) {
    s1.SetData(rr, s1.Data(rr) + v1[rr]);
  }
  return s1;
}

// operator overload for addition with a scalar
colMatrix colMatrix::operator+(const double &scalar) const {
  auto s1 = *this;

  for (auto rr = 0; rr < s1.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) + scalar);
  }
  return s1;
}

// operator overload for addition with a scalar
colMatrix operator+(const double &scalar, const colMatrix &s2) {
  colMatrix s1(s2.Size());

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s2.Data(rr) + scalar);
  }
  return s1;
}

// operator overload for subtraction
colMatrix colMatrix::operator-(const colMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot subtract column matrices, dimensions do not agree."
         << endl;
  }

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) - s2.Data(rr));
  }
  return s1;
}

// operator overload for addition
colMatrix colMatrix::operator-(const vector<double> &v1) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != static_cast<int> (v1.size())) {
    cerr << "ERROR: Cannot subtract column matrix of size_ " << s1.size_
         << " to vector class of size_ " << v1.size()
         << ", dimensions do not agree." << endl;
  }

  for (auto rr = 0; rr < s1.size_; rr++) {
    s1.SetData(rr, s1.Data(rr) - v1[rr]);
  }
  return s1;
}

// operator overload for subtraction with a scalar
colMatrix colMatrix::operator-(const double &scalar) const {
  auto s1 = *this;

  for (auto rr = 0; rr < s1.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) - scalar);
  }
  return s1;
}

// operator overload for subtraction with a scalar
colMatrix operator-(const double &scalar, const colMatrix &s2) {
  colMatrix s1(s2.Size());

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, scalar - s2.Data(rr));
  }
  return s1;
}

// operator overload for elementwise multiplication
colMatrix colMatrix::operator*(const colMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot elementwise multiply column matrices, dimensions do "
            "not agree." << endl;
  }

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) * s2.Data(rr));
  }
  return s1;
}

// operator overload for multiplication with a scalar
colMatrix colMatrix::operator*(const double &scalar) const {
  auto s1 = *this;

  for (auto rr = 0; rr < s1.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) * scalar);
  }
  return s1;
}

// operator overload for multiplication with a scalar
colMatrix operator*(const double &scalar, const colMatrix &s2) {
  colMatrix s1(s2.Size());

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s2.Data(rr) * scalar);
  }
  return s1;
}

// operator overload for elementwise division
colMatrix colMatrix::operator/(const colMatrix &s2) const {
  auto s1 = *this;

  // check to see that matrix dimensions are the same
  if (s1.size_ != s2.size_) {
    cerr << "ERROR: Cannot elementwise divide column matrices, dimensions do "
            "not agree." << endl;
  }

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) / s2.Data(rr));
  }
  return s1;
}

// operator overload for division with a scalar
colMatrix colMatrix::operator/(const double &scalar) const {
  auto s1 = *this;

  for (auto rr = 0; rr < s1.Size(); rr++) {
    s1.SetData(rr, s1.Data(rr) / scalar);
  }
  return s1;
}

// operator overload for division with a scalar
colMatrix operator/(const double &scalar, const colMatrix &s2) {
  colMatrix s1(s2.Size());

  for (auto rr = 0; rr < s2.Size(); rr++) {
    s1.SetData(rr, scalar / s2.Data(rr));
  }
  return s1;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const colMatrix &m) {

  for (auto rr = 0; rr < m.Size(); rr++) {
    cout << m.Data(rr) << endl;
  }

  return os;
}

// member function to zero the matrix
void colMatrix::Zero() {
  for (auto rr = 0; rr < size_; rr++) {
    this->SetData(rr, 0.0);
  }
}

// member function to sum column matrix
double colMatrix::Sum() {
  auto sum = 0.0;
  for (auto ii = 0; ii < this->Size(); ii++) {
    sum += this->Data(ii);
  }
  return sum;
}

// member function to delete the contents of the data_ structure and resize it
void colMatrix::CleanResizeZero(const int &s) {

  delete[] data_;
  data_ = new double[s];
  size_ = s;

  for (auto cc = 0; cc < size_; cc++) {
    this->SetData(cc, 0.0);
  }
}

// ------------------------------------------------------------------
// functions for genArray class

// operator overload for addition
genArray genArray::operator+(const genArray &s2) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] += s2[rr];
  }
  return s1;
}

// operator overload for addition
genArray genArray::operator+(const vector<double> &v1) const {
  // check to see that vector is appropriate size_
  if (v1.size() != NUMVARS) {
    cerr << "ERROR: Cannot add vector and genArray because vector not "
            "appropriate size_! Vector is of size_ " << v1.size() << endl;
    exit(0);
  }

  genArray s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] += v1[rr];
  }
  return s1;
}

// operator overload for addition with a scalar
genArray genArray::operator+(const double &scalar) const {
  auto s1 = *this;

  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] += scalar;
  }
  return s1;
}

// operator overload for addition with a scalar
genArray operator+(const double &scalar, const genArray &s2) {
  genArray s1;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] = scalar + s2[rr];
  }
  return s1;
}

// operator overload for subtraction
genArray genArray::operator-(const genArray &s2) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] -= s2[rr];
  }
  return s1;
}

// operator overload for addition
genArray genArray::operator-(const vector<double> &v1) const {
  // check to see that vector is appropriate size_
  if (v1.size() != NUMVARS) {
    cerr << "ERROR: Cannot subtract vector and genArray because vector not "
            "appropriate size_! Vector is of size_ " << v1.size() << endl;
    exit(0);
  }

  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] -= v1[rr];
  }
  return s1;
}

// operator overload for subtraction with a scalar
genArray genArray::operator-(const double &scalar) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] -= scalar;
  }
  return s1;
}

// operator overload for subtraction with a scalar
genArray operator-(const double &scalar, const genArray &s2) {
  genArray s1;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] = scalar - s2[rr];
  }
  return s1;
}

// operator overload for elementwise multiplication
genArray genArray::operator*(const genArray &s2) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] *= s2[rr];
  }
  return s1;
}

// operator overload for multiplication with a scalar
genArray genArray::operator*(const double &scalar) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] *= scalar;
  }
  return s1;
}

// operator overload for multiplication with a scalar
genArray operator*(const double &scalar, const genArray &s2) {
  genArray s1;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] = s2[rr] * scalar;
  }
  return s1;
}

// operator overload for elementwise division
genArray genArray::operator/(const genArray &s2) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] /= s2[rr];
  }
  return s1;
}

// operator overload for division with a scalar
genArray genArray::operator/(const double &scalar) const {
  auto s1 = *this;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] /= scalar;
  }
  return s1;
}

// operator overload for division with a scalar
genArray operator/(const double &scalar, const genArray &s2) {
  genArray s1;
  for (auto rr = 0; rr < NUMVARS; rr++) {
    s1[rr] = scalar / s2[rr];
  }
  return s1;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const genArray &m) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    cout << m[rr] << endl;
  }
  return os;
}

// member function to zero the matrix
void genArray::Zero() {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    (*this)[rr] = 0.0;
  }
}

// member function to sum column matrix
double genArray::Sum() {
  auto sum = 0.0;
  for (auto ii = 0; ii < NUMVARS; ii++) {
    sum += (*this)[ii];
  }
  return sum;
}

// member function to sum the residuals from all processors
void genArray::GlobalReduceMPI(const int &rank, const int &numEqns) {
  // Get residuals from all processors
  if (rank == ROOTP) {
    MPI_Reduce(MPI_IN_PLACE, &data_[0], numEqns, MPI_DOUBLE, MPI_SUM,
               ROOTP, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&data_[0], &data_[0], numEqns, MPI_DOUBLE, MPI_SUM, ROOTP,
               MPI_COMM_WORLD);
  }
}

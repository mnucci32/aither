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

#include <cstdlib>  // exit()
#include <iostream>  // cout
#include <algorithm>  // swap
#include "matrix.hpp"
#include "genArray.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::copy;
using std::swap_ranges;
using std::fabs;

// member function to swap rows of matrix
void squareMatrix::SwapRows(const int &r1, const int &r2) {
  if (r1 != r2) {
    for (auto ii = 0; ii < size_; ii++) {
      std::swap(data_[this->GetLoc(ii, r1)], data_[this->GetLoc(ii, r2)]);
    }
  }
}

// member function to invert matrix using Gauss-Jordan elimination
void squareMatrix::Inverse() {
  squareMatrix I(size_);
  I.Identity();

  for (auto cPivot = 0, r = 0; r < size_; r++, cPivot++) {
    // find pivot row
    auto rPivot = this->FindMaxInCol(r, cPivot, size_ - 1);

    // swap rows
    this->SwapRows(r, rPivot);
    I.SwapRows(r, rPivot);

    if (r != 0) {  // if not on first row, need to get rid entries ahead of
                   // pivot
      for (auto ii = 0; ii < cPivot; ii++) {
        auto factor = (*this)(r, ii) / (*this)(ii, ii);
        this->LinCombRow(ii, factor, r);
        I.LinCombRow(ii, factor, r);
      }
    }

    // normalize row by pivot
    if ((*this)(r, cPivot) == 0.0) {
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination! Matrix (mid "
              "inversion) is" << endl << *this << endl;
      exit(1);
    }
    auto normFactor = 1.0 / (*this)(r, cPivot);
    this->RowMultiply(r, cPivot, normFactor);  // only multiply entries from
                                                 // pivot and to the right
    I.RowMultiply(r, 0, normFactor);  // multiply all entries
  }

  // matrix is now upper triangular, work way back up to identity matrix
  // start with second to last row
  for (auto cPivot = size_ - 2, r = size_ - 2; r >= 0; r--, cPivot--) {
    for (auto ii = size_ - 1; ii > cPivot; ii--) {
      auto factor = (*this)(r, ii);
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
    (*this)(r2, ii) = (*this)(r2, ii) - (*this)(r1, ii) * factor;
  }
}

// member function to multiply a row by a given factor
void squareMatrix::RowMultiply(const int &r, const int &c,
                               const double &factor) {
  for (auto ii = c; ii < size_; ii++) {
    (*this)(r, ii) = (*this)(r, ii) * factor;
  }
}

// member function to find maximum absolute value in a given column and range
// within that column and return the corresponding row indice
int squareMatrix::FindMaxInCol(const int &c, const int &start,
                               const int &end) const {
  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii < end + 1; ii++) {
    if (fabs((*this)(ii, c)) > maxVal) {
      maxVal = fabs((*this)(ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}


// operator overload for multiplication
squareMatrix squareMatrix::MatMult(const squareMatrix &s2) const {
  // check to see that matrix dimensions are the same
  if (size_ != s2.size_) {
    cerr << "ERROR: Cannot multiply matrices, dimensions do not agree." << endl;
  }

  auto s1 = *this;

  for (auto cc = 0; cc < s2.Size(); cc++) {
    for (auto rr = 0; rr < s2.Size(); rr++) {
      auto newVal = 0.0;
      for (auto ii = 0; ii < s2.Size(); ii++) {
        newVal += (*this)(rr, ii) * s2(ii, cc);
      }
      s1(rr, cc) = newVal;
    }
  }
  return s1;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const squareMatrix &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    for (auto cc = 0; cc < m.Size(); cc++) {
      cout << m(rr, cc);
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
  for (auto &val : data_) {
    val = 0.0;
  }
}

// member function to set matrix to Identity
void squareMatrix::Identity() {
  for (auto rr = 0; rr < this->Size(); rr++) {
    for (auto cc = 0; cc < this->Size(); cc++) {
      if (rr == cc) {
        (*this)(rr, cc) = 1.0;
      } else {
        (*this)(rr, cc) = 0.0;
      }
    }
  }
}

// member function to do matrix/vector multplication
genArray squareMatrix::ArrayMult(const genArray &vec, const int pos) const {
  // vec -- vector to multiply with

  auto product = vec;

  // zero out portion of genArray that will be written over
  if (pos == 0) {
    for (auto ii = 0; ii < NUMFLOWVARS; ii++) {
      product[ii] = 0.0;
    }
  } else {
    for (auto ii = pos; ii < NUMVARS; ii++) {
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

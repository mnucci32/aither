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

#include <cstdlib>  // exit()
#include <cmath>  // fabs
#include <iostream>  // cout
#include <algorithm>  // swap
#include "matrix.hpp"
#include "varArray.hpp"

using std::cout;
using std::endl;
using std::cerr;

// member function to swap rows of matrix
void squareMatrix::SwapRows(const int &r1, const int &r2) {
  MSG_ASSERT(r1 < size_ && r2 < size_, "index outside of range");
  if (r1 != r2) {
    std::swap_ranges(this->begin() + this->GetLoc(r1, 0),
                     this->begin() + this->GetLoc(r1, size_),
                     this->begin() + this->GetLoc(r2, 0));
  }
}

// member function to invert matrix using Gauss-Jordan elimination
void squareMatrix::Inverse() {
  squareMatrix I(size_);
  I.Identity();

  for (auto cPivot = 0, r = 0; r < size_; ++r, ++cPivot) {
    // find pivot row
    auto rPivot = this->FindMaxInCol(r, cPivot, size_ - 1);

    // swap rows
    this->SwapRows(r, rPivot);
    I.SwapRows(r, rPivot);

    if (r != 0) {  // if not first row, need to get rid entries ahead of pivot
      for (auto ii = 0; ii < cPivot; ++ii) {
        auto factor = (*this)(r, ii) / (*this)(ii, ii);
        this->LinCombRow(ii, factor, r);
        I.LinCombRow(ii, factor, r);
      }
    }

    // normalize row by pivot
    if ((*this)(r, cPivot) == 0.0) {
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination! Matrix (mid "
              "inversion) is" << endl << *this << endl;
      exit(EXIT_FAILURE);
    }
    // only normalize entries from pivot and to the right
    auto normFactor = 1.0 / (*this)(r, cPivot);
    this->RowMultiply(r, cPivot, normFactor);
    I.RowMultiply(r, 0, normFactor);  // multiply all entries
  }

  // matrix is now upper triangular, work way back up to identity matrix
  // start with second to last row
  for (auto cPivot = size_ - 2, r = size_ - 2; r >= 0; --r, --cPivot) {
    for (auto ii = size_ - 1; ii > cPivot; --ii) {
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
  MSG_ASSERT(r1 < size_ && r2 < size_, "index outside of range");
  std::transform(
      this->begin() + this->GetLoc(r2, 0),
      this->begin() + this->GetLoc(r2, size_),
      this->begin() + this->GetLoc(r1, 0), this->begin() + this->GetLoc(r2, 0),
      [&factor](const auto &v1, const auto &v2) { return v1 - factor * v2; });
}

// member function to multiply a row by a given factor
void squareMatrix::RowMultiply(const int &r, const int &c,
                               const double &factor) {
  MSG_ASSERT(r < size_ && c < size_, "index outside of range");
  for_each(this->begin() + this->GetLoc(r, c),
           this->begin() + this->GetLoc(r, size_),
           [&factor](auto &val) { val *= factor; });
}

// member function to find maximum absolute value in a given column and range
// within that column and return the corresponding row indice
int squareMatrix::FindMaxInCol(const int &c, const int &start,
                               const int &end) const {
  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii <= end; ii++) {
    if (fabs((*this)(ii, c)) > maxVal) {
      maxVal = fabs((*this)(ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}


// operator overload for multiplication
// using cache efficient implimentation
squareMatrix squareMatrix::MatMult(const squareMatrix &s2) const {
  MSG_ASSERT(this->Size() == s2.Size(), "matrix size mismatch");
  squareMatrix result(s2.Size());
  auto m1 = [&](const int &r, const int &c) -> decltype(auto) {
    return (*this)(r, c);
  };
  auto m2 = [&](const int &r, const int &c) -> decltype(auto) {
    return s2(r, c);
  };
  auto res = [&](const int &r, const int &c) -> decltype(auto) {
    return result(r, c);
  };
  MatrixMultiply(m1, m2, this->Size(), res);
  return result;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const squareMatrix &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    for (auto cc = 0; cc < m.Size(); cc++) {
      os << m(rr, cc);
      if (cc != (m.Size() - 1)) {
        os << ", ";
      } else {
        os << endl;
      }
    }
  }
  return os;
}

// member function to set matrix to Identity
void squareMatrix::Identity() {
  for (auto rr = 0; rr < this->Size(); ++rr) {
    for (auto cc = 0; cc < this->Size(); ++cc) {
      (*this)(rr, cc) = (rr == cc) ? 1.0 : 0.0;
    }
  }
}

// member function to find maximum absolute value on diagonal
// this can be used to find the spectral radius of a diagoanl matrix
double squareMatrix::MaxAbsValOnDiagonal() const {
  auto maxVal = 0.0;
  for (auto ii = 0; ii < size_; ++ii) {
    maxVal = std::max(fabs((*this)(ii, ii)), maxVal);
  }
  return maxVal;
}

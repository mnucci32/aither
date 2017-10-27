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

// operator overload for multiplication
// using cache efficient implimentation
squareMatrix squareMatrix::MatMult(const squareMatrix &s2) const {
  MSG_ASSERT(this->Size() == s2.Size(), "matrix size mismatch");
  squareMatrix result(s2.Size());
  MatrixMultiply(this->begin(), s2.begin(), result.begin(), this->Size());
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


// ---------------------------------------------------------------------------
// generic matrix functions

// function to invert matrix using Gauss-Jordan elimination
void MatrixInverse(const vector<double>::iterator &mat, const int &size) {
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };

  squareMatrix I(size);
  I.Identity();

  for (auto cPivot = 0, r = 0; r < size; ++r, ++cPivot) {
    // find pivot row
    auto rPivot = FindMaxInColumn(mat, size, r, cPivot, size - 1);

    // swap rows
    SwapMatRows(mat, size, r, rPivot);
    I.SwapRows(r, rPivot);

    if (r != 0) {  // if not first row, need to get rid entries ahead of pivot
      for (auto ii = 0; ii < cPivot; ++ii) {
        auto factor = GetVal(mat, r, ii) / GetVal(mat, ii, ii);
        LinearCombRow(mat, size, ii, factor, r);
        I.LinCombRow(ii, factor, r);
      }
    }

    // normalize row by pivot
    if (GetVal(mat, r, cPivot) == 0.0) {
      cerr << "ERROR: Singular matrix in Gauss-Jordan elimination!" << endl;
      exit(EXIT_FAILURE);
    }
    // only normalize entries from pivot and to the right
    auto normFactor = 1.0 / GetVal(mat, r, cPivot);
    RowMultiplyFactor(mat, size, r, cPivot, normFactor);
    I.RowMultiply(r, 0, normFactor);  // multiply all entries
  }

  // matrix is now upper triangular, work way back up to identity matrix
  // start with second to last row
  for (auto cPivot = size - 2, r = size - 2; r >= 0; --r, --cPivot) {
    for (auto ii = size - 1; ii > cPivot; --ii) {
      auto factor = GetVal(mat, r, ii);
      LinearCombRow(mat, size, ii, factor, r);
      I.LinCombRow(ii, factor, r);
    }
  }

  // set matrix equal to its inverse
  std::copy(I.begin(), I.end(), mat);
}

// member function to swap rows of matrix
void SwapMatRows(const vector<double>::iterator &mat, const int &size,
                 const int &r1, const int &r2) {
  MSG_ASSERT(r1 < size && r2 < size, "index outside of range");
  if (r1 != r2) {
    auto GetLoc = [&size](const int &r, const int &c) -> decltype(auto) {
      return r * size + c;
    };
    std::swap_ranges(mat + GetLoc(r1, 0),
                     mat + GetLoc(r1, size),
                     mat + GetLoc(r2, 0));
  }
}

// function to add a linear combination of one row to another
void LinearCombRow(const vector<double>::iterator &mat, const int &size,
                   const int &r1, const double &factor, const int &r2) {
  MSG_ASSERT(r1 < size && r2 < size, "index outside of range");
  auto GetLoc = [&size](const int &r, const int &c) -> decltype(auto) {
    return r * size + c;
  };
  std::transform(
      mat + GetLoc(r2, 0), mat + GetLoc(r2, size), mat + GetLoc(r1, 0),
      mat + GetLoc(r2, 0),
      [&factor](const auto &v1, const auto &v2) { return v1 - factor * v2; });
}

// function to multiply a row by a given factor
void RowMultiplyFactor(const vector<double>::iterator &mat, const int &size,
                       const int &r, const int &c, const double &factor) {
  MSG_ASSERT(r < size && c < size, "index outside of range");
  auto GetLoc = [&size](const int &r, const int &c) -> decltype(auto) {
    return r * size + c;
  };
  for_each(mat + GetLoc(r, c), mat + GetLoc(r, size),
           [&factor](auto &val) { val *= factor; });
}

// function to find maximum absolute value in a given column and range
// within that column and return the corresponding row indice
int FindMaxInColumn(const vector<double>::const_iterator &mat, const int &size,
                    const int &c, const int &start, const int &end) {
  MSG_ASSERT(start < size && end < size && c < size,
             "index outside of range");
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };

  auto maxVal = 0.0;
  auto maxRow = 0;
  for (auto ii = start; ii <= end; ++ii) {
    if (fabs(GetVal(mat, ii, c)) > maxVal) {
      maxVal = fabs(GetVal(mat, ii, c));
      maxRow = ii;
    }
  }
  return maxRow;
}

// member function to set matrix to Identity
void IdentityMatrix(const vector<double>::iterator &mat, const int &size) {
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };
  for (auto rr = 0; rr < size; ++rr) {
    for (auto cc = 0; cc < size; ++cc) {
      GetVal(mat, rr, cc) = (rr == cc) ? 1.0 : 0.0;
    }
  }
}

// function to find maximum absolute value on diagonal
// this can be used to find the spectral radius of a diagoanl matrix
double MaximumAbsValOnDiagonal(const vector<double>::const_iterator &mat,
                               const int &size) {
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };

  auto maxVal = 0.0;
  for (auto ii = 0; ii < size; ++ii) {
    maxVal = std::max(fabs(GetVal(mat, ii, ii)), maxVal);
  }
  return maxVal;
}

// function to multiply two matrices
void MatrixMultiply(const vector<double>::const_iterator &matL,
                    const vector<double>::const_iterator &matR,
                    const vector<double>::iterator &result, const int &size) {
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };

  for (auto cc = 0; cc < size; ++cc) {
    for (auto rr = 0; rr < size; ++rr) {
      for (auto ii = 0; ii < size; ++ii) {
        GetVal(result, rr, ii) += GetVal(matL, rr, cc) * GetVal(matR, cc, ii);
      }
    }
  }
}

void MultiplyFacOnDiagonal(const vector<double>::iterator &mat, const int &size,
                           const double &val) {
  // val -- value to multiply along diagonal
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };
  for (auto ii = 0; ii < size; ++ii) {
    GetVal(mat, ii, ii) *= val;
  }
}

void AddFacOnDiagonal(const vector<double>::iterator &mat, const int &size,
                      const double &val) {
  // val -- value to multiply along diagonal
  auto GetVal = [&size](const auto &mat, const int &r,
                        const int &c) -> decltype(auto) {
    return *(mat + r * size + c);
  };
  for (auto ii = 0; ii < size; ++ii) {
    GetVal(mat, ii, ii) += val;
  }
}


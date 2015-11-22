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

#ifndef MATRIXHEADERDEF  // only if the macro MATRIXHEADERDEF is not defined
                         // execute these lines of code
#define MATRIXHEADERDEF  // define the macro

#include <iostream>
#include "macros.hpp"

using std::ostream;

// class to store a square matrix
class squareMatrix {
  int size_;
  double *data_;

 public:
  // constructor
  explicit squareMatrix(const int &a) : size_(a) { data_ = new double[a * a]; }
  squareMatrix() : size_(0), data_(nullptr) {}

  // copy constructor
  squareMatrix(const squareMatrix &cp);

  // copy assignment operator
  squareMatrix& operator=(squareMatrix other);

  // move constructor
  explicit squareMatrix(squareMatrix &&other) noexcept : squareMatrix() {
    swap(*this, other);
    other.data_ = nullptr;
  }

  // member functions
  double Data(const int &, const int &) const;
  void SetData(const int &, const int &, const double &);
  int Size() const { return size_; }
  void SwapRows(const int &, const int &);
  void Inverse();
  int FindMaxInCol(const int &, const int &, const int &) const;
  void RowMultiply(const int &, const int &, const double &);
  void LinCombRow(const int &, const double &, const int &);
  void Zero();
  void Identity();

  // operator overloads
  squareMatrix operator+(const squareMatrix &) const;
  squareMatrix operator-(const squareMatrix &) const;
  squareMatrix operator*(const squareMatrix &) const;

  squareMatrix operator+(const double &) const;
  squareMatrix operator-(const double &) const;
  squareMatrix operator*(const double &) const;
  squareMatrix operator/(const double &) const;

  friend squareMatrix operator+(const double &, const squareMatrix &);
  friend squareMatrix operator-(const double &, const squareMatrix &);
  friend squareMatrix operator*(const double &, const squareMatrix &);
  friend squareMatrix operator/(const double &, const squareMatrix &);
  friend ostream &operator<<(ostream &os, const squareMatrix &);

  friend void swap(squareMatrix &first, squareMatrix &second) noexcept;

  // destructor
  ~squareMatrix() noexcept {
    delete[] data_;
    data_ = nullptr;
  }
};

// Class to store the implicit flux jacobians for the entire mesh. Only values
// on populated diagonals are stored.
class matrixDiagonal {
  int size_;
  squareMatrix *data_;

 public:
  // constructor
  explicit matrixDiagonal(const int &a) : size_(a) {
    data_ = new squareMatrix[a];}
  matrixDiagonal() : size_(0), data_(nullptr) {}

  // copy constructor
  matrixDiagonal(const matrixDiagonal &cp);

  // copy assignment operator
  matrixDiagonal &operator=(matrixDiagonal other);

  // move constructor
  explicit matrixDiagonal(matrixDiagonal &&other) noexcept : matrixDiagonal() {
    swap(*this, other);
    other.data_ = nullptr;
  }

  // member functions
  squareMatrix Data(const int &) const;
  void SetData(const int &, const squareMatrix &);
  int Size() const { return size_; }
  void Zero(const int &);
  void CleanResizeZero(const int &, const int &);
  void Inverse();

  // operator overloads
  friend ostream &operator<<(ostream &os, const matrixDiagonal &);

  friend void swap(matrixDiagonal &first, matrixDiagonal &second) noexcept;

  // destructor
  ~matrixDiagonal() noexcept {
    delete[] data_;
    data_ = nullptr;
  }
};

// function declarations

#endif

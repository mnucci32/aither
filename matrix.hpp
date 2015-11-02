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

#ifndef MATRIXHEADERDEF  // only if the macro MATRIXHEADERDEF is not defined
                         // execute these lines of code
#define MATRIXHEADERDEF  // define the macro

#include <math.h>  // sqrt
#include <vector>  // vector
#include <string>  // string
#include <iostream>
#include <algorithm>   // swap
#include "mpi.h"       // parallelism
#include "macros.hpp"

using std::vector;
using std::string;
using std::ostream;

// class to store a column matrix
class colMatrix {
  int size_;
  double *data_;

 public:
  // constructor
  explicit colMatrix(const int &a) : size_(a) { data_ = new double[a]; }
  colMatrix() : size_(0), data_(nullptr) {}

  // copy constructor
  colMatrix(const colMatrix &cp);

  // copy assignment operator
  colMatrix &operator=(colMatrix other);

  // move constructor
  explicit colMatrix(colMatrix &&other) noexcept : colMatrix() {
    swap(*this, other);
    other.data_ = nullptr;
  }

  // member functions
  double Data(const int &) const;
  void SetData(const int &, const double &);
  int Size() const { return size_; }
  void Zero();
  double Sum();
  void CleanResizeZero(const int &);

  // operator overloads
  colMatrix operator+(const colMatrix &) const;
  colMatrix operator-(const colMatrix &) const;
  colMatrix operator*(const colMatrix &) const;
  colMatrix operator/(const colMatrix &) const;

  colMatrix operator+(const vector<double> &) const;
  colMatrix operator-(const vector<double> &) const;

  colMatrix operator+(const double &) const;
  colMatrix operator-(const double &) const;
  colMatrix operator*(const double &) const;
  colMatrix operator/(const double &) const;

  friend colMatrix operator+(const double &, const colMatrix &);
  friend colMatrix operator-(const double &, const colMatrix &);
  friend colMatrix operator*(const double &, const colMatrix &);
  friend colMatrix operator/(const double &, const colMatrix &);
  friend ostream &operator<<(ostream &os, const colMatrix &);

  friend void swap(colMatrix &first, colMatrix &second);

  // destructor
  ~colMatrix() noexcept {
    delete[] data_;
    data_ = nullptr;
  }
};

/*class to store an array of a fixed size_ equal to the number of variables
being solved for. This is useful because a vector of these will be
contiguous in memory. */
class genArray {
  double data_[NUMVARS];

 public:
  // constructor
  genArray() : data_{0.0} {}
  explicit genArray(const double &);
  genArray(const double &a, const double &b, const double &c, const double &d,
           const double &e)
      : data_{a, b, c, d, e} {}
  genArray(const double &a, const double &b, const double &c, const double &d,
           const double &e, const double &f, const double &g)
      : data_{a, b, c, d, e, f, g} {}

  // member functions
  void Zero();
  double Sum();

  // move constructor and assignment operator
  genArray(genArray&&) noexcept = default;
  genArray& operator=(genArray&&) noexcept = default;

  // copy constructor and assignment operator
  genArray(const genArray&) = default;
  genArray& operator=(const genArray&) = default;

  // operator overloads
  genArray operator+(const genArray &) const;
  genArray operator-(const genArray &) const;
  genArray operator*(const genArray &) const;
  genArray operator/(const genArray &) const;

  genArray operator+(const vector<double> &) const;
  genArray operator-(const vector<double> &) const;

  genArray operator+(const double &) const;
  genArray operator-(const double &) const;
  genArray operator*(const double &) const;
  genArray operator/(const double &) const;

  const double & operator[](const int &r) const { return data_[r]; }
  double & operator[](const int &r) { return data_[r]; }

  friend genArray operator+(const double &, const genArray &);
  friend genArray operator-(const double &, const genArray &);
  friend genArray operator*(const double &, const genArray &);
  friend genArray operator/(const double &, const genArray &);
  friend ostream &operator<<(ostream &os, const genArray &);

  void GlobalReduceMPI(const int &, const int &);

  // destructor
  ~genArray() noexcept {}
};

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
  colMatrix Multiply(const colMatrix &) const;

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

  friend void swap(squareMatrix &first, squareMatrix &second);

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

  friend void swap(matrixDiagonal &first, matrixDiagonal &second);

  // destructor
  ~matrixDiagonal() noexcept {
    delete[] data_;
    data_ = nullptr;
  }
};

// function declarations

#endif

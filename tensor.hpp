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

#ifndef TENSORHEADERDEF  // only if the macro TENSORHEADERDEF is not defined
                         // execute these lines of code

#define TENSORHEADERDEF  // define the macro

// This file contains the header and implementation for the vector3d templated
// class. The implementation is included in
// this file because the class is templated. If the implementation were not
// included, it would also have to be included in
// any files that depend on the header. Leaving the implementation in
// streamlines the compiling process.

#include <math.h>    // sqrt()
#include <iostream>  // ostream
#include "vector3d.hpp"

using std::ostream;
using std::endl;

// Templated class for a 2D tensor holding 9 elements
template <typename T>
class tensor {
  T data_[9];

 public:
  // constructor
  tensor(const T &a, const T &b, const T &c, const T &d, const T &e,
         const T &f, const T &g, const T &h, const T &i)
      : data_{a, b, c, d, e, f, g, h, i} {}
  tensor() : data_{0, 0, 0, 0, 0, 0, 0, 0, 0} {}
  explicit tensor(const T &i) : data_{i, 0, 0, 0, i, 0, 0, 0, i} {}
  tensor(const vector3d<T> &v1, const vector3d<T> &v2, const vector3d<T> &v3)
      : data_{v1.X(), v1.Y(), v1.Z(), v2.X(), v2.Y(),
              v2.Z(), v3.X(), v3.Y(), v3.Z()} {}

  // move constructor and assignment operator
  tensor(tensor<T>&&) noexcept = default;
  tensor& operator=(tensor<T>&&) noexcept = default;

  // copy constructor and assignment operator
  tensor(const tensor<T>&) = default;
  tensor& operator=(const tensor<T>&) = default;

  // member functions
  // operator overloads
  tensor<T> operator+(const tensor &) const;
  tensor<T> operator-(const tensor &) const;
  tensor<T> operator*(const tensor &) const;
  tensor<T> operator+(const T &) const;
  tensor<T> operator-(const T &) const;
  tensor<T> operator*(const T &) const;
  tensor<T> operator/(const T &) const;
  template <typename TT>
  friend tensor<TT> operator*(const TT &, const tensor<TT> &);
  template <typename TT>
  friend tensor<TT> operator/(const TT &, const tensor<TT> &);
  template <typename TT>
  friend tensor<TT> operator+(const TT &, const tensor<TT> &);
  template <typename TT>
  friend tensor<TT> operator-(const TT &, const tensor<TT> &);
  template <typename TT>
  friend ostream &operator<<(ostream &os, const tensor<TT> &);
  // assignment of data_ members
  void SetXX(const T &val) { data_[0] = val; }
  void SetXY(const T &val) { data_[1] = val; }
  void SetXZ(const T &val) { data_[2] = val; }
  void SetYX(const T &val) { data_[3] = val; }
  void SetYY(const T &val) { data_[4] = val; }
  void SetYZ(const T &val) { data_[5] = val; }
  void SetZX(const T &val) { data_[6] = val; }
  void SetZY(const T &val) { data_[7] = val; }
  void SetZZ(const T &val) { data_[8] = val; }
  // access of data_ members
  T XX() const { return data_[0]; }
  T XY() const { return data_[1]; }
  T XZ() const { return data_[2]; }
  T YX() const { return data_[3]; }
  T YY() const { return data_[4]; }
  T YZ() const { return data_[5]; }
  T ZX() const { return data_[6]; }
  T ZY() const { return data_[7]; }
  T ZZ() const { return data_[8]; }
  // math functions
  T Trace() const { return data_[0] + data_[4] + data_[8]; }
  tensor<T> Transpose() const;
  vector3d<T> MatMult(const vector3d<T> &) const;
  T DoubleDot(const tensor<T> &) const;
  T DoubleDotTrans(const tensor<T> &) const;
  void Identity();
  void Zero();

  // destructor
  ~tensor() noexcept {}
};

// operator overload for addition - element wise addition
template <typename T>
tensor<T> tensor<T>::operator+(const tensor &v2) const {
  tensor<T> temp = (*this);
  temp.data_[0] += v2.data_[0];
  temp.data_[1] += v2.data_[1];
  temp.data_[2] += v2.data_[2];

  temp.data_[3] += v2.data_[3];
  temp.data_[4] += v2.data_[4];
  temp.data_[5] += v2.data_[5];

  temp.data_[6] += v2.data_[6];
  temp.data_[7] += v2.data_[7];
  temp.data_[8] += v2.data_[8];

  return temp;
}

// operator overload for subtraction - element wise subtraction
template <typename T>
tensor<T> tensor<T>::operator-(const tensor &v2) const {
  tensor<T> temp = *this;
  temp.data_[0] -= v2.data_[0];
  temp.data_[1] -= v2.data_[1];
  temp.data_[2] -= v2.data_[2];

  temp.data_[3] -= v2.data_[3];
  temp.data_[4] -= v2.data_[4];
  temp.data_[5] -= v2.data_[5];

  temp.data_[6] -= v2.data_[6];
  temp.data_[7] -= v2.data_[7];
  temp.data_[8] -= v2.data_[8];

  return temp;
}

// operator overload for multiplication - matrix multiplication
template <typename T>
tensor<T> tensor<T>::operator*(const tensor &v2) const {
  tensor<T> temp;

  temp.data_[0] = (*this).data_[0] * v2.data_[0] +
                  (*this).data_[1] * v2.data_[3] +
                  (*this).data_[2] * v2.data_[6];
  temp.data_[1] = (*this).data_[0] * v2.data_[1] +
                  (*this).data_[1] * v2.data_[4] +
                  (*this).data_[2] * v2.data_[7];
  temp.data_[2] = (*this).data_[0] * v2.data_[2] +
                  (*this).data_[1] * v2.data_[5] +
                  (*this).data_[2] * v2.data_[8];

  temp.data_[3] = (*this).data_[3] * v2.data_[0] +
                  (*this).data_[4] * v2.data_[3] +
                  (*this).data_[5] * v2.data_[6];
  temp.data_[4] = (*this).data_[3] * v2.data_[1] +
                  (*this).data_[4] * v2.data_[4] +
                  (*this).data_[5] * v2.data_[7];
  temp.data_[5] = (*this).data_[3] * v2.data_[2] +
                  (*this).data_[4] * v2.data_[5] +
                  (*this).data_[5] * v2.data_[8];

  temp.data_[6] = (*this).data_[6] * v2.data_[0] +
                  (*this).data_[7] * v2.data_[3] +
                  (*this).data_[8] * v2.data_[6];
  temp.data_[7] = (*this).data_[6] * v2.data_[1] +
                  (*this).data_[7] * v2.data_[4] +
                  (*this).data_[8] * v2.data_[7];
  temp.data_[8] = (*this).data_[6] * v2.data_[2] +
                  (*this).data_[7] * v2.data_[5] +
                  (*this).data_[8] * v2.data_[8];

  return temp;
}

// operator overload for addition with a scalar - element wise
// addition
template <typename T>
tensor<T> tensor<T>::operator+(const T &scalar) const {
  tensor<T> temp = *this;
  temp.data_[0] += scalar;
  temp.data_[1] += scalar;
  temp.data_[2] += scalar;

  temp.data_[3] += scalar;
  temp.data_[4] += scalar;
  temp.data_[5] += scalar;

  temp.data_[6] += scalar;
  temp.data_[7] += scalar;
  temp.data_[8] += scalar;

  return temp;
}

// operator overload for addition with a scalar - element wise
// addition
// this function is a friend function of the class so that double + tensor
// behaves as tensor + double
template <typename TT>
tensor<TT> operator+(const TT &scalar, const tensor<TT> &v1) {
  tensor<TT> temp;
  temp.data_[0] = v1.data_[0] + scalar;
  temp.data_[1] = v1.data_[1] + scalar;
  temp.data_[2] = v1.data_[2] + scalar;

  temp.data_[3] = v1.data_[3] + scalar;
  temp.data_[4] = v1.data_[4] + scalar;
  temp.data_[5] = v1.data_[5] + scalar;

  temp.data_[6] = v1.data_[6] + scalar;
  temp.data_[7] = v1.data_[7] + scalar;
  temp.data_[8] = v1.data_[8] + scalar;

  return temp;
}

// operator overload for subtraction with a scalar - element wise
// subtraction
template <typename T>
tensor<T> tensor<T>::operator-(const T &scalar) const {
  tensor<T> temp = *this;
  temp.data_[0] -= scalar;
  temp.data_[1] -= scalar;
  temp.data_[2] -= scalar;

  temp.data_[3] -= scalar;
  temp.data_[4] -= scalar;
  temp.data_[5] -= scalar;

  temp.data_[6] -= scalar;
  temp.data_[7] -= scalar;
  temp.data_[8] -= scalar;

  return temp;
}

// operator overload for subtraction with a scalar - element wise
// subtraction
// this function is a friend function of the class so that double - tensor
// behaves as tensor - double
template <typename TT>
tensor<TT> operator-(const TT &scalar, const tensor<TT> &v1) {
  tensor<TT> temp;
  temp.data_[0] = v1.data_[0] - scalar;
  temp.data_[1] = v1.data_[1] - scalar;
  temp.data_[2] = v1.data_[2] - scalar;

  temp.data_[3] = v1.data_[3] - scalar;
  temp.data_[4] = v1.data_[4] - scalar;
  temp.data_[5] = v1.data_[5] - scalar;

  temp.data_[6] = v1.data_[6] - scalar;
  temp.data_[7] = v1.data_[7] - scalar;
  temp.data_[8] = v1.data_[8] - scalar;

  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
template <typename T>
tensor<T> tensor<T>::operator*(const T &scalar) const {
  tensor<T> temp = *this;
  temp.data_[0] *= scalar;
  temp.data_[1] *= scalar;
  temp.data_[2] *= scalar;

  temp.data_[3] *= scalar;
  temp.data_[4] *= scalar;
  temp.data_[5] *= scalar;

  temp.data_[6] *= scalar;
  temp.data_[7] *= scalar;
  temp.data_[8] *= scalar;

  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
// this function is a friend function of the class so that double * tensor
// behaves as tensor * double
template <typename TT>
tensor<TT> operator*(const TT &scalar, const tensor<TT> &v1) {
  tensor<TT> temp;
  temp.data_[0] = v1.data_[0] * scalar;
  temp.data_[1] = v1.data_[1] * scalar;
  temp.data_[2] = v1.data_[2] * scalar;

  temp.data_[3] = v1.data_[3] * scalar;
  temp.data_[4] = v1.data_[4] * scalar;
  temp.data_[5] = v1.data_[5] * scalar;

  temp.data_[6] = v1.data_[6] * scalar;
  temp.data_[7] = v1.data_[7] * scalar;
  temp.data_[8] = v1.data_[8] * scalar;

  return temp;
}

// operator overload for division with a scalar - element wise division
template <typename T>
tensor<T> tensor<T>::operator/(const T &scalar) const {
  tensor<T> temp = *this;
  temp.data_[0] /= scalar;
  temp.data_[1] /= scalar;
  temp.data_[2] /= scalar;

  temp.data_[3] /= scalar;
  temp.data_[4] /= scalar;
  temp.data_[5] /= scalar;

  temp.data_[6] /= scalar;
  temp.data_[7] /= scalar;
  temp.data_[8] /= scalar;

  return temp;
}

// operator overload for division with a scalar - element wise division
// this function is a friend function of the class so that double / tensor works
template <typename TT>
tensor<TT> operator/(const TT &scalar, const tensor<TT> &v1) {
  tensor<TT> temp;
  temp.data_[0] = scalar / v1.data_[0];
  temp.data_[1] = scalar / v1.data_[1];
  temp.data_[2] = scalar / v1.data_[2];

  temp.data_[3] = scalar / v1.data_[3];
  temp.data_[4] = scalar / v1.data_[4];
  temp.data_[5] = scalar / v1.data_[5];

  temp.data_[6] = scalar / v1.data_[6];
  temp.data_[7] = scalar / v1.data_[7];
  temp.data_[8] = scalar / v1.data_[8];

  return temp;
}

// operator overload for << - allows use of cout, cerr, etc.
template <typename TT>
ostream &operator<<(ostream &os, const tensor<TT> &v1) {
  os << v1.data_[0] << ", " << v1.data_[1] << ", " << v1.data_[2] << endl;
  os << v1.data_[3] << ", " << v1.data_[4] << ", " << v1.data_[5] << endl;
  os << v1.data_[6] << ", " << v1.data_[7] << ", " << v1.data_[8] << endl;
  return os;
}

// Function to return the transpose of the given tensor
template <typename T>
tensor<T> tensor<T>::Transpose() const {
  tensor<T> temp;

  temp.data_[0] = (*this).data_[0];
  temp.data_[1] = (*this).data_[3];
  temp.data_[2] = (*this).data_[6];

  temp.data_[3] = (*this).data_[1];
  temp.data_[4] = (*this).data_[4];
  temp.data_[5] = (*this).data_[7];

  temp.data_[6] = (*this).data_[2];
  temp.data_[7] = (*this).data_[5];
  temp.data_[8] = (*this).data_[8];

  return temp;
}

// Function to return the matrix multiplication of a tensor and vector3d
template <typename T>
vector3d<T> tensor<T>::MatMult(const vector3d<T> &vec) const {
  vector3d<T> temp;

  temp.SetX((*this).data_[0] * vec.X() + (*this).data_[1] * vec.Y() +
            (*this).data_[2] * vec.Z());
  temp.SetY((*this).data_[3] * vec.X() + (*this).data_[4] * vec.Y() +
            (*this).data_[5] * vec.Z());
  temp.SetZ((*this).data_[6] * vec.X() + (*this).data_[7] * vec.Y() +
            (*this).data_[8] * vec.Z());

  return temp;
}

// Function to convert the matrix to the identity matrix
template <typename T>
void tensor<T>::Identity() {
  T one = 1;
  T zero = 0;

  data_[0] = one;
  data_[1] = zero;
  data_[2] = zero;

  data_[3] = zero;
  data_[4] = one;
  data_[5] = zero;

  data_[6] = zero;
  data_[7] = zero;
  data_[8] = one;
}

// Function to zero a tensor
template <typename T>
void tensor<T>::Zero() {
  T var = 0;
  data_[0] = var;
  data_[1] = var;
  data_[2] = var;
  data_[3] = var;
  data_[4] = var;
  data_[5] = var;
  data_[6] = var;
  data_[7] = var;
  data_[8] = var;
}

// Function to return the double dot product of two tensors
// Aij Bij
template <typename T>
T tensor<T>::DoubleDotTrans(const tensor<T> &temp) const {
  return data_[0] * temp.data_[0] + data_[1] * temp.data_[1] +
      data_[2] * temp.data_[2] + data_[3] * temp.data_[3] +
      data_[4] * temp.data_[4] + data_[5] * temp.data_[5] +
      data_[6] * temp.data_[6] + data_[7] * temp.data_[7] +
      data_[8] * temp.data_[8];
}

// Function to return the double dot product of two tensors
// Aij Bji
template <typename T>
T tensor<T>::DoubleDot(const tensor<T> &temp) const {
  return data_[0] * temp.data_[0] + data_[1] * temp.data_[3] +
      data_[2] * temp.data_[6] + data_[3] * temp.data_[1] +
      data_[4] * temp.data_[4] + data_[5] * temp.data_[7] +
      data_[6] * temp.data_[2] + data_[7] * temp.data_[5] +
      data_[8] * temp.data_[8];
}


#endif

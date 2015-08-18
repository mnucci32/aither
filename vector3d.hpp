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

#ifndef VECTOR3DHEADERDEF  // only if the macro VECTOR3DHEADERDEF is not defined
                           // execute these lines of code

#define VECTOR3DHEADERDEF  // define the macro

// This file contains the header and implementation for the vector3d templated
// class. The implementation is included in this file because the class is
// templated. If the implementation were not included, it would also have to be
// included in any files that depend on the header. Leaving the implementation
// in streamlines the compiling process.

#include <math.h>    // sqrt()
#include <iostream>  // ostream

using std::ostream;

// Templated class for a vector holding 3 entries
template <class T>
class vector3d {
  T data_[3];

 public:
  // constructor
  vector3d(const T &a, const T &b, const T &c) : data_{a, b, c} {}
  vector3d() : data_{0, 0, 0} {}
  // copy constructor
  vector3d(const vector3d<T> &a) : data_{a.data_[0], a.data_[1], a.data_[2]} {}

  // member functions
  // operator overloads
  vector3d<T> operator+(const vector3d &) const;
  vector3d<T> operator-(const vector3d &) const;
  vector3d<T> operator*(const T &) const;
  vector3d<T> operator/(const T &) const;
  T operator[](const int &a) const { return data_[a]; }
  T &operator[](const int &a) { return data_[a]; }
  bool operator==(const vector3d &) const;
  template <class TT>
  friend vector3d<TT> operator*(const TT &, const vector3d<TT> &);
  template <class TT>
  friend vector3d<TT> operator/(const TT &, const vector3d<TT> &);
  template <class TT>
  friend ostream &operator<<(ostream &os, const vector3d<TT> &);
  // assignment of data_ members
  void SetX(const T &val) { data_[0] = val; }
  void SetY(const T &val) { data_[1] = val; }
  void SetZ(const T &val) { data_[2] = val; }
  // access of data_ members
  T X() const { return data_[0]; }
  T Y() const { return data_[1]; }
  T Z() const { return data_[2]; }
  // math functions
  T DotProd(const vector3d &) const;
  vector3d<T> CrossProd(const vector3d &) const;
  inline T Mag() const;
  inline T MagSq() const;
  T SumElem() const;
  T Distance(const vector3d &) const;

  // destructor
  ~vector3d() {}
};

template <class T>
class unitVec3dMag {
  vector3d<T> unitVec_;
  T mag_;

 public:
  // constructor
  unitVec3dMag() : unitVec_(), mag_(0) {}
  explicit unitVec3dMag(const vector3d<T> &a) : unitVec_{a / a.Mag()},
    mag_(a.Mag()) {}
  // copy constructor
  unitVec3dMag(const unitVec3dMag<T> &a) : unitVec_(a.unitVec_), mag_(a.mag_) {}

  // member functions
  // operator overloads
  unitVec3dMag<T> operator+(const unitVec3dMag &) const;
  unitVec3dMag<T> operator-(const unitVec3dMag &) const;
  unitVec3dMag<T> operator*(const T &) const;
  unitVec3dMag<T> operator/(const T &) const;
  bool operator==(const unitVec3dMag &) const;
  template <class TT>
  friend unitVec3dMag<TT> operator*(const TT &, const unitVec3dMag<TT> &);
  template <class TT>
  friend unitVec3dMag<TT> operator/(const TT &, const unitVec3dMag<TT> &);
  template <class TT>
  friend ostream &operator<<(ostream &os, const unitVec3dMag<TT> &);
  // access of data_ members
  vector3d<T> UnitVector() const { return unitVec_; }
  T X() const { return unitVec_.X(); }
  T Y() const { return unitVec_.Y(); }
  T Z() const { return unitVec_.Z(); }
  T Mag() const { return mag_; }

  T MagSq() const { return mag_ * mag_; }
  vector3d<T> Vector() const { return unitVec_ * mag_; }

  // math functions
  T DotProd(const unitVec3dMag &) const;
  unitVec3dMag<T> CrossProd(const unitVec3dMag &) const;

  // destructor
  ~unitVec3dMag() {}
};


// operator overload for addition - element wise addition
template <class T>
vector3d<T> vector3d<T>::operator+(const vector3d &v2) const {
  vector3d<T> temp = *this;
  temp.data_[0] += v2.data_[0];
  temp.data_[1] += v2.data_[1];
  temp.data_[2] += v2.data_[2];
  return temp;
}

// operator overload for subtraction - element wise subtraction
template <class T>
vector3d<T> vector3d<T>::operator-(const vector3d &v2) const {
  vector3d<T> temp = *this;
  temp.data_[0] -= v2.data_[0];
  temp.data_[1] -= v2.data_[1];
  temp.data_[2] -= v2.data_[2];
  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
template <class T>
vector3d<T> vector3d<T>::operator*(const T &scalar) const {
  vector3d<T> temp = *this;
  temp.data_[0] *= scalar;
  temp.data_[1] *= scalar;
  temp.data_[2] *= scalar;
  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
// this function is a friend function of the class so that double * vector3d
// behaves as vector3d * double
template <class TT>
vector3d<TT> operator*(const TT &scalar, const vector3d<TT> &v1) {
  vector3d<TT> temp;
  temp.data_[0] = v1.data_[0] * scalar;
  temp.data_[1] = v1.data_[1] * scalar;
  temp.data_[2] = v1.data_[2] * scalar;
  return temp;
}

// operator overload for divisioin with a scalar - element wise division
template <class T>
vector3d<T> vector3d<T>::operator/(const T &scalar) const {
  vector3d<T> temp = *this;
  temp.data_[0] /= scalar;
  temp.data_[1] /= scalar;
  temp.data_[2] /= scalar;
  return temp;
}

// operator overload for division with a scalar - element wise division
// this function is a friend function of the class so that double / vector3d
// works
template <class TT>
vector3d<TT> operator/(const TT &scalar, const vector3d<TT> &v1) {
  vector3d<TT> temp;
  temp.data_[0] = scalar / v1.data_[0];
  temp.data_[1] = scalar / v1.data_[1];
  temp.data_[2] = scalar / v1.data_[2];
  return temp;
}

// operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream &operator<<(ostream &os, const vector3d<TT> &v1) {
  os << v1.data_[0] << ", " << v1.data_[1] << ", " << v1.data_[2];
  return os;
}

// Function to calculate the dot product of two vectors
template <class T>
T vector3d<T>::DotProd(const vector3d &v2) const {
  return data_[0] * v2.data_[0] + data_[1] * v2.data_[1] +
         data_[2] * v2.data_[2];
}

// operator overload for comparison
template <class T>
bool vector3d<T>::operator==(const vector3d &v2) const {
  bool test = false;
  if (data_[0] == v2.data_[0] && data_[1] == v2.data_[1] &&
      data_[2] == v2.data_[2]) {
    test = true;
  }
  return test;
}

// Function to calculate the cross product of two vectors
template <class T>
vector3d<T> vector3d<T>::CrossProd(const vector3d &v2) const {
  vector3d<T> crossProd;

  crossProd.data_[0] = data_[1] * v2.data_[2] - data_[2] * v2.data_[1];
  crossProd.data_[1] = -1.0 * (data_[0] * v2.data_[2] - data_[2] * v2.data_[0]);
  crossProd.data_[2] = data_[0] * v2.data_[1] - data_[1] * v2.data_[0];

  return crossProd;
}

// Function to calculate the magnitude of the vector
template <class T>
T vector3d<T>::Mag() const {
  return sqrt(data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2]);
}

// Function to calculate the square of the magnitude of the vector
template <class T>
T vector3d<T>::MagSq() const {
  return data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2];
}

// Function to sum the elements in the vector
template <class T>
T vector3d<T>::SumElem() const {
  return data_[0] + data_[1] + data_[2];
}

// Function to calculate the distance between two vector3ds
template <class T>
T vector3d<T>::Distance(const vector3d &v2) const {
  return sqrt(pow(data_[0] - v2.data_[0], 2) + pow(data_[1] - v2.data_[1], 2) +
              pow(data_[2] - v2.data_[2], 2));
}

// -----------------------------------------------------------------
// Functions for unitVec3dMag class
// operator overload for addition - element wise addition
template <class T>
unitVec3dMag<T> unitVec3dMag<T>::operator+(const unitVec3dMag &v2) const {
  vector3d<T> newVec = (*this).unitVec_ * (*this).mag_ +
      v2.unitVec_ * v2.mag_;
  unitVec3dMag<T> temp(newVec);
  return temp;
}

// operator overload for subtraction - element wise subtraction
template <class T>
unitVec3dMag<T> unitVec3dMag<T>::operator-(const unitVec3dMag &v2) const {
  vector3d<T> newVec = (*this).unitVec_ * (*this).mag_ -
      v2.unitVec_ * v2.mag_;
  unitVec3dMag<T> temp(newVec);
  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
template <class T>
unitVec3dMag<T> unitVec3dMag<T>::operator*(const T &scalar) const {
  unitVec3dMag<T> temp = *this;
  temp.mag_ *= scalar;
  return temp;
}

// operator overload for multiplication with a scalar - element wise
// multiplication
// this function is a friend function of the class so that double * unitVec3dMag
// behaves as unitVec3dMag * double
template <class TT>
unitVec3dMag<TT> operator*(const TT &scalar, const unitVec3dMag<TT> &v1) {
  unitVec3dMag<TT> temp = v1;
  temp.mag_ *= scalar;
  return temp;
}

// operator overload for divisioin with a scalar - element wise division
template <class T>
unitVec3dMag<T> unitVec3dMag<T>::operator/(const T &scalar) const {
  unitVec3dMag<T> temp = *this;
  temp.mag_ /= scalar;
  return temp;
}

// operator overload for division with a scalar - element wise division
// this function is a friend function of the class so that double / unitVec3dMag
// works
template <class TT>
unitVec3dMag<TT> operator/(const TT &scalar, const unitVec3dMag<TT> &v1) {
  unitVec3dMag<TT> temp = v1;
  temp.mag_ = scalar / temp.mag_;
  return temp;
}

// operator overload for << - allows use of cout, cerr, etc.
template <class TT>
ostream &operator<<(ostream &os, const unitVec3dMag<TT> &v1) {
  os << v1.unitVec_.X() << ", " << v1.unitVec_.Y() << ", " << v1.unitVec_.Z()
     << ", " << v1.mag_;
  return os;
}

// Function to calculate the dot product of two vectors
template <class T>
T unitVec3dMag<T>::DotProd(const unitVec3dMag &v2) const {
  return mag_ * v2.mag_ * (unitVec_.dotProd(v2.unitVec_));
}

// operator overload for comparison
template <class T>
bool unitVec3dMag<T>::operator==(const unitVec3dMag &v2) const {
  bool test = false;
  if (unitVec_ == v2.unitVec_ && mag_ == v2.mag_) {
    test = true;
  }
  return test;
}

// Function to calculate the cross product of two vectors
template <class T>
unitVec3dMag<T> unitVec3dMag<T>::CrossProd(const unitVec3dMag &v2) const {
  vector3d<T> newVec = (unitVec_ * mag_).CrossProd(v2.unitVec * v2.mag_);
  unitVec3dMag<T> crossProd(newVec);
  return crossProd;
}

#endif

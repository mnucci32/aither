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

#ifndef VECTOR3DHEADERDEF  // only if the macro VECTOR3DHEADERDEF is not defined
                           // execute these lines of code

#define VECTOR3DHEADERDEF  // define the macro

// This file contains the header and implementation for the vector3d templated
// class. The implementation is included in this file because the class is
// templated. If the implementation were not included, it would also have to be
// included in any files that depend on the header. Leaving the implementation
// in streamlines the compiling process.

#include <cmath>        // sqrt()
#include <iostream>     // ostream
#include <type_traits>  // is_arithmetic

#define LENGTH 3

using std::ostream;

// Templated class for a vector holding 3 entries
template <typename T>
class vector3d {
  static_assert(std::is_arithmetic<T>::value,
                "vector3d<T> requires an arithmetic type!");

  T data_[LENGTH];

 public:
  // constructor
  vector3d(const T &a, const T &b, const T &c) : data_{a, b, c} {}
  vector3d() : data_{0, 0, 0} {}

  // move constructor and assignment operator
  vector3d(vector3d<T>&&) noexcept = default;
  vector3d& operator=(vector3d<T>&&) noexcept = default;

  // copy constructor and assignment operator
  vector3d(const vector3d<T>&) = default;
  vector3d& operator=(const vector3d<T>&) = default;

  // member functions
  // operator overloads
  const T& operator[](const int &a) const { return data_[a]; }
  T& operator[](const int &a) { return data_[a]; }

  bool operator==(const vector3d<T>&) const;

  inline vector3d<T> & operator+=(const vector3d<T> &);
  inline vector3d<T> & operator-=(const vector3d<T> &);
  inline vector3d<T> & operator*=(const vector3d<T> &);
  inline vector3d<T> & operator/=(const vector3d<T> &);

  inline vector3d<T> & operator+=(const T &);
  inline vector3d<T> & operator-=(const T &);
  inline vector3d<T> & operator*=(const T &);
  inline vector3d<T> & operator/=(const T &);

  inline vector3d<T> operator+(const T &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline vector3d<T> operator-(const T &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline vector3d<T> operator*(const T &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline vector3d<T> operator/(const T &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // assignment of data members
  void SetX(const T &val) { data_[0] = val; }
  void SetY(const T &val) { data_[1] = val; }
  void SetZ(const T &val) { data_[2] = val; }

  // access of data members
  T X() const { return data_[0]; }
  T Y() const { return data_[1]; }
  T Z() const { return data_[2]; }

  // math functions
  T DotProd(const vector3d<T>&) const;
  vector3d<T> CrossProd(const vector3d<T>&) const;
  inline T Mag() const;
  inline T MagSq() const;
  T SumElem() const;
  T Distance(const vector3d<T>&) const;
  T DistSq(const vector3d<T>&) const;
  vector3d<T> Normalize() const;
  void Zero();

  // destructor
  ~vector3d() noexcept {}
};

template <typename T>
class unitVec3dMag {
  static_assert(std::is_arithmetic<T>::value,
                "unitVect3dMag<T> requires an arithmetic type!");

  vector3d<T> unitVec_;
  T mag_;

 public:
  // constructor
  unitVec3dMag() : unitVec_(), mag_(0) {}
  explicit unitVec3dMag(const vector3d<T> &a) : unitVec_(a.Normalize()),
    mag_(a.Mag()) {}

  // move constructor and assignment operator
  unitVec3dMag(unitVec3dMag<T>&&) noexcept = default;
  unitVec3dMag& operator=(unitVec3dMag<T>&&) noexcept = default;

  // copy constructor and assignment operator
  unitVec3dMag(const unitVec3dMag<T>&) = default;
  unitVec3dMag& operator=(const unitVec3dMag<T>&) = default;

  // member functions
  // operator overloads
  bool operator==(const unitVec3dMag<T>&) const;

  inline unitVec3dMag<T> & operator+=(const unitVec3dMag<T> &);
  inline unitVec3dMag<T> & operator-=(const unitVec3dMag<T> &);

  inline unitVec3dMag<T> & operator*=(const T &);
  inline unitVec3dMag<T> & operator/=(const T &);

  inline unitVec3dMag<T> operator*(const T &s) const {
    auto lhs = *this;
    lhs.mag_ *= s;
    return *this;
  }
  inline unitVec3dMag<T> operator/(const T &s) const {
    auto lhs = *this;
    lhs.mag_ /= s;
    return *this;
  }

  // access of data_ members
  vector3d<T> UnitVector() const { return unitVec_; }
  T X() const { return unitVec_.X(); }
  T Y() const { return unitVec_.Y(); }
  T Z() const { return unitVec_.Z(); }
  T Mag() const { return mag_; }

  T MagSq() const { return mag_ * mag_; }
  vector3d<T> Vector() const { return unitVec_ * mag_; }

  // math functions
  T DotProd(const unitVec3dMag<T>&) const;
  unitVec3dMag<T> CrossProd(const unitVec3dMag<T>&) const;

  // destructor
  ~unitVec3dMag() noexcept {}
};

// operator overload for addition
template <typename T>
vector3d<T> & vector3d<T>::operator+=(const vector3d<T> &vec) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    data_[rr] += vec[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
vector3d<T> & vector3d<T>::operator-=(const vector3d<T> &vec) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    data_[rr] -= vec[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
vector3d<T> & vector3d<T>::operator*=(const vector3d<T> &vec) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    data_[rr] *= vec[rr];
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
vector3d<T> & vector3d<T>::operator/=(const vector3d<T> &vec) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    data_[rr] /= vec[rr];
  }
  return *this;
}

template <typename T>
inline const vector3d<T> operator+(vector3d<T> lhs, const vector3d<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const vector3d<T> operator-(vector3d<T> lhs, const vector3d<T> &rhs) {
  return lhs -= rhs;
}

template <typename T>
inline const vector3d<T> operator*(vector3d<T> lhs, const vector3d<T> &rhs) {
  return lhs *= rhs;
}

template <typename T>
inline const vector3d<T> operator/(vector3d<T> lhs, const vector3d<T> &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
template <typename T>
vector3d<T> & vector3d<T>::operator+=(const T &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
vector3d<T> & vector3d<T>::operator-=(const T &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
vector3d<T> & vector3d<T>::operator*=(const T &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
vector3d<T> & vector3d<T>::operator/=(const T &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

template <typename T>
inline const vector3d<T> operator+(const T &lhs, vector3d<T> rhs) {
  return rhs += lhs;
}

template <typename T>
inline const vector3d<T> operator-(const T &lhs, vector3d<T> rhs) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    rhs[rr] = lhs - rhs[rr];
  }
  return rhs;
}

template <typename T>
inline const vector3d<T> operator*(const T &lhs, vector3d<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const vector3d<T> operator/(const T &lhs, vector3d<T> rhs) {
  for (auto rr = 0; rr < LENGTH; rr++) {
    rhs[rr] = lhs / rhs[rr];
  }
  return rhs;
}

template <typename T>
ostream &operator<<(ostream &os, const vector3d<T> &v) {
  os << v[0] << ", " << v[1] << ", " << v[2];
  return os;
}

// Function to calculate the dot product of two vectors
template <typename T>
T vector3d<T>::DotProd(const vector3d<T>&v2) const {
  return data_[0] * v2.data_[0] + data_[1] * v2.data_[1] +
         data_[2] * v2.data_[2];
}

// operator overload for comparison
template <typename T>
bool vector3d<T>::operator==(const vector3d<T>&v2) const {
  auto test = false;
  if (data_[0] == v2.data_[0] && data_[1] == v2.data_[1] &&
      data_[2] == v2.data_[2]) {
    test = true;
  }
  return test;
}

// Function to calculate the cross product of two vectors
template <typename T>
vector3d<T> vector3d<T>::CrossProd(const vector3d<T>&v2) const {
  vector3d<T> crossProd;

  crossProd.data_[0] = data_[1] * v2.data_[2] - data_[2] * v2.data_[1];
  crossProd.data_[1] = -1.0 * (data_[0] * v2.data_[2] - data_[2] * v2.data_[0]);
  crossProd.data_[2] = data_[0] * v2.data_[1] - data_[1] * v2.data_[0];

  return crossProd;
}

// Function to calculate the magnitude of the vector
template <typename T>
T vector3d<T>::Mag() const {
  return sqrt(this->MagSq());
}

// Function to calculate the square of the magnitude of the vector
template <typename T>
T vector3d<T>::MagSq() const {
  return this->DotProd(*this);
}

// Function to sum the elements in the vector
template <typename T>
T vector3d<T>::SumElem() const {
  return data_[0] + data_[1] + data_[2];
}

// Function to calculate the distance between two vector3ds
template <typename T>
T vector3d<T>::Distance(const vector3d<T>&v2) const {
  return sqrt(this->DistSq(v2));
}

// Function to calculate the distance squared between two vector3ds
template <typename T>
T vector3d<T>::DistSq(const vector3d<T>&v2) const {
  return pow(data_[0] - v2.data_[0], 2) + pow(data_[1] - v2.data_[1], 2) +
      pow(data_[2] - v2.data_[2], 2);
}

// Function to normalize a vector3d into a unit vector
template <typename T>
vector3d<T> vector3d<T>::Normalize() const {
  return (*this) / this->Mag();
}

// Function to zero out a vector
template <typename T>
void vector3d<T>::Zero() {
  for (auto &val : data_) {
    val = 0;
  }
}


// -----------------------------------------------------------------
// Functions for unitVec3dMag class
// Function to calculate the dot product of two vectors
template <typename T>
T unitVec3dMag<T>::DotProd(const unitVec3dMag<T>&v2) const {
  return mag_ * v2.mag_ * (unitVec_.DotProd(v2.unitVec_));
}

// operator overload for comparison
template <typename T>
bool unitVec3dMag<T>::operator==(const unitVec3dMag<T>&v2) const {
  auto test = false;
  if (unitVec_ == v2.unitVec_ && mag_ == v2.mag_) {
    test = true;
  }
  return test;
}

// Function to calculate the cross product of two vectors
template <typename T>
unitVec3dMag<T> unitVec3dMag<T>::CrossProd(const unitVec3dMag<T>&v2) const {
  auto newVec = this->Vector().CrossProd(v2.Vector());
  unitVec3dMag<T> crossProd(newVec);
  return crossProd;
}

// operator overload for addition
template <typename T>
unitVec3dMag<T> & unitVec3dMag<T>::operator+=(const unitVec3dMag<T> &vec) {
  auto newVec = this->Vector() + vec.Vector();
  *this = unitVec3dMag<T>(newVec);
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
unitVec3dMag<T> & unitVec3dMag<T>::operator-=(const unitVec3dMag<T> &vec) {
  auto newVec = this->Vector() - vec.Vector();
  *this = unitVec3dMag<T>(newVec);
  return *this;
}

template <typename T>
inline const unitVec3dMag<T> operator+(unitVec3dMag<T> lhs,
                                       const unitVec3dMag<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const unitVec3dMag<T> operator-(unitVec3dMag<T> lhs,
                                       const unitVec3dMag<T> &rhs) {
  return lhs -= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for elementwise multiplication
template <typename T>
unitVec3dMag<T> & unitVec3dMag<T>::operator*=(const T &scalar) {
  if (scalar < 0.0) {
    unitVec_ *= -1.0;
  }
  mag_ *= fabs(scalar);
  return *this;
}

// operator overload for elementwise division
template <typename T>
unitVec3dMag<T> & unitVec3dMag<T>::operator/=(const T &scalar) {
  if (scalar < 0.0) {
    unitVec_ *= -1.0;
  }
  mag_ /= fabs(scalar);
  return *this;
}

template <typename T>
inline const unitVec3dMag<T> operator*(const T &lhs, unitVec3dMag<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const unitVec3dMag<T> operator/(const T &lhs, unitVec3dMag<T> rhs) {
  if (lhs < 0.0) {
    rhs.unitVec_ *= -1.0;
  }
  rhs.mag_ = fabs(lhs) / rhs.mag_;
  return rhs;
}

template <typename T>
ostream &operator<<(ostream &os, const unitVec3dMag<T> &v) {
  os << v.Vector() << ", " << v.Mag();
  return os;
}


#endif

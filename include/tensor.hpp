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

#ifndef TENSORHEADERDEF  // only if the macro TENSORHEADERDEF is not defined
                         // execute these lines of code

#define TENSORHEADERDEF  // define the macro

// This file contains the header and implementation for the vector3d templated
// class. The implementation is included in
// this file because the class is templated. If the implementation were not
// included, it would also have to be included in
// any files that depend on the header. Leaving the implementation in
// streamlines the compiling process.

#include <cmath>        // sqrt()
#include <iostream>     // ostream
#include <type_traits>  // is_arithmetic
#include <numeric>      // accumulate
#include "vector3d.hpp"

#define TENSORSIZE 9

using std::ostream;
using std::endl;

// Templated class for a 2D tensor holding 9 elements
template <typename T>
class tensor {
  static_assert(std::is_arithmetic<T>::value,
                "tensor<T> requires an arithmetic type!");

  T data_[TENSORSIZE];

  // private member functions
  int GetLoc(const int &rr, const int &cc) const {
    return cc + rr * this->Size();
  }

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
  T& operator()(const int &rr, const int &cc) {
    return data_[this->GetLoc(rr, cc)];
  }
  const T& operator()(const int &rr, const int &cc) const {
    return data_[this->GetLoc(rr, cc)];
  }


  inline tensor<T> & operator+=(const tensor<T> &);
  inline tensor<T> & operator-=(const tensor<T> &);
  inline tensor<T> & operator*=(const tensor<T> &);
  inline tensor<T> & operator/=(const tensor<T> &);

  inline tensor<T> & operator+=(const T &);
  inline tensor<T> & operator-=(const T &);
  inline tensor<T> & operator*=(const T &);
  inline tensor<T> & operator/=(const T &);

  inline tensor<T> operator+(const T &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline tensor<T> operator-(const T &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline tensor<T> operator*(const T &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline tensor<T> operator/(const T &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  template <typename TT>
  friend inline const tensor<TT> operator-(const TT &lhs, tensor<TT> rhs);
  template <typename TT>
  friend inline const tensor<TT> operator/(const TT &lhs, tensor<TT> rhs);

  // assignment of data members
  void SetXX(const T &val) { data_[0] = val; }
  void SetXY(const T &val) { data_[1] = val; }
  void SetXZ(const T &val) { data_[2] = val; }
  void SetYX(const T &val) { data_[3] = val; }
  void SetYY(const T &val) { data_[4] = val; }
  void SetYZ(const T &val) { data_[5] = val; }
  void SetZX(const T &val) { data_[6] = val; }
  void SetZY(const T &val) { data_[7] = val; }
  void SetZZ(const T &val) { data_[8] = val; }

  // access of data members
  T XX() const { return data_[0]; }
  T XY() const { return data_[1]; }
  T XZ() const { return data_[2]; }
  T YX() const { return data_[3]; }
  T YY() const { return data_[4]; }
  T YZ() const { return data_[5]; }
  T ZX() const { return data_[6]; }
  T ZY() const { return data_[7]; }
  T ZZ() const { return data_[8]; }

  vector3d<T> X() const { return {data_[0], data_[1], data_[2]}; }
  vector3d<T> Y() const { return {data_[3], data_[4], data_[5]}; }
  vector3d<T> Z() const { return {data_[6], data_[7], data_[8]}; }

  // math functions
  T Trace() const { return data_[0] + data_[4] + data_[8]; }
  tensor<T> Transpose() const;
  vector3d<T> MatMult(const vector3d<T> &) const;
  tensor<T> MatMult(const tensor<T> &) const;
  T DoubleDot(const tensor<T> &) const;
  T DoubleDotTrans(const tensor<T> &) const;
  void Identity();
  void Zero();
  int Size() const {return 3;}
  tensor<T> RemoveComponent(const vector3d<T> &) const;
  vector3d<T> LinearCombination(const vector3d<T> &) const;
  T Sum() const {
    return std::accumulate(std::begin(data_), std::end(data_), T(0));
  }

  // destructor
  ~tensor() noexcept {}
};

// operator overload for addition
template <typename T>
tensor<T> & tensor<T>::operator+=(const tensor<T> &ten) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    data_[rr] += ten.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
tensor<T> & tensor<T>::operator-=(const tensor<T> &ten) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    data_[rr] -= ten.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
tensor<T> & tensor<T>::operator*=(const tensor<T> &ten) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    data_[rr] *= ten.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
tensor<T> & tensor<T>::operator/=(const tensor<T> &ten) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    data_[rr] /= ten.data_[rr];
  }
  return *this;
}

template <typename T>
inline const tensor<T> operator+(tensor<T> lhs, const tensor<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const tensor<T> operator-(tensor<T> lhs, const tensor<T> &rhs) {
  return lhs -= rhs;
}

template <typename T>
inline const tensor<T> operator*(tensor<T> lhs, const tensor<T> &rhs) {
  return lhs *= rhs;
}

template <typename T>
inline const tensor<T> operator/(tensor<T> lhs, const tensor<T> &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
template <typename T>
tensor<T> & tensor<T>::operator+=(const T &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
tensor<T> & tensor<T>::operator-=(const T &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
tensor<T> & tensor<T>::operator*=(const T &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
tensor<T> & tensor<T>::operator/=(const T &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

template <typename T>
inline const tensor<T> operator+(const T &lhs, tensor<T> rhs) {
  return rhs += lhs;
}

template <typename T>
inline const tensor<T> operator-(const T &lhs, tensor<T> rhs) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    rhs.data_[rr] = lhs - rhs.data_[rr];
  }
  return rhs;
}

template <typename T>
inline const tensor<T> operator*(const T &lhs, tensor<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const tensor<T> operator/(const T &lhs, tensor<T> rhs) {
  for (auto rr = 0; rr < TENSORSIZE; rr++) {
    rhs.data_[rr] = lhs / rhs.data_[rr];
  }
  return rhs;
}

// function for matrix multiplication
// using more cache efficient method
template <typename T>
tensor<T> tensor<T>::MatMult(const tensor &v2) const {
  tensor<T> temp;
  for (auto cc = 0; cc < v2.Size(); cc++) {
    for (auto rr = 0; rr < v2.Size(); rr++) {
      for (auto ii = 0; ii < v2.Size(); ii++) {
        temp(rr, ii) += (*this)(rr, cc) * v2(cc, ii);
      }
    }
  }
  return temp;
}

// operator overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const tensor<T> &v1) {
  os << v1.XX() << ", " << v1.XY() << ", " << v1.XZ() << endl;
  os << v1.YX() << ", " << v1.YY() << ", " << v1.YZ() << endl;
  os << v1.ZX() << ", " << v1.ZY() << ", " << v1.ZZ() << endl;
  return os;
}

// Function to return the transpose of the given tensor
template <typename T>
tensor<T> tensor<T>::Transpose() const {
  tensor<T> temp;

  temp.data_[0] = data_[0];
  temp.data_[1] = data_[3];
  temp.data_[2] = data_[6];

  temp.data_[3] = data_[1];
  temp.data_[4] = data_[4];
  temp.data_[5] = data_[7];

  temp.data_[6] = data_[2];
  temp.data_[7] = data_[5];
  temp.data_[8] = data_[8];

  return temp;
}

// Function to return the matrix multiplication of a tensor and vector3d
template <typename T>
vector3d<T> tensor<T>::MatMult(const vector3d<T> &vec) const {
  vector3d<T> temp;

  temp.SetX(data_[0] * vec.X() + data_[1] * vec.Y() + data_[2] * vec.Z());
  temp.SetY(data_[3] * vec.X() + data_[4] * vec.Y() + data_[5] * vec.Z());
  temp.SetZ(data_[6] * vec.X() + data_[7] * vec.Y() + data_[8] * vec.Z());

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

// function to remove the components that are aligned with a given direction
template <typename T>
tensor<T> tensor<T>::RemoveComponent(const vector3d<T> &dir) const {
  auto x = this->X();
  x -= x.DotProd(dir) * dir;
  auto y = this->Y();
  y -= y.DotProd(dir) * dir;
  auto z = this->Z();
  z -= z.DotProd(dir) * dir;
  return {x, y, z};
}

// function to scale the rows of the tensor by a scale factor and sum them
// together
template <typename T>
vector3d<T> tensor<T>::LinearCombination(const vector3d<T> &vec) const {
  auto comb = this->X() * vec.X();
  comb += this->Y() * vec.Y();
  comb += this->Z() * vec.Z();
  return comb;
}

#endif

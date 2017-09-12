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

#ifndef VARARRAYHEADERDEF
#define VARARRAYHEADERDEF

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include "macros.hpp"

using std::ostream;
using std::vector;

// Class to hold an array of variables. Length is equal to number of equations
// being solved for.
class varArray {
  vector<double> data_;
  int momentumIndex_;
  int energyIndex_;
  int turbulenceIndex_;

 public:
  // constructor
  varArray(const int &numEqns, const int &numSpecies) : data_(numEqns, 0.0) {
    MSG_ASSERT(numEqns > numSpecies && numEqns >= 5,
               "number of equations should be greater than number of species");
    momentumIndex_ = numSpecies;
    energyIndex_ = momentumIndex_ + 3;
    turbulenceIndex_ = energyIndex_ + 1;
  }

  // member functions
  int Size() const { return data_.size(); }
  int NumSpecies() const { return momentumIndex_; }
  bool IsMultiSpecies() const { return this->NumSpecies() > 1; }
  bool HasTurbulenceData() const { return this->Size() != turbulenceIndex_; }
  int MomentumIndex() const { return momentumIndex_; }
  int EnergyIndex() const { return energyIndex_; }
  int TurbulenceIndex() const { return turbulenceIndex_; }

  const double &SpeciesN(const int &ii) const {
    MSG_ASSERT(ii < momentumIndex_, "requesting species variable out of range");
    return (*this)[ii];
  }
  const double &MomentumX() const { return (*this)[momentumIndex_]; }
  const double &MomentumY() const { return (*this)[momentumIndex_ + 1]; }
  const double &MomentumZ() const { return (*this)[momentumIndex_ + 2]; }
  const double &Energy() const { return (*this)[this->EnergyIndex()]; }
  const double &TurbulenceN(const int &ii) const {
    MSG_ASSERT(tubulenceIndex_ + ii >= this->Size(),
               "requesting turbulence variable out of range");
    return (*this)[turbulenceIndex_ + ii];
  }

  void Zero() { std::fill(std::begin(data_), std::end(data_), 0.0); }
  double Sum() {
    return std::accumulate(std::begin(data_), std::end(data_), 0.0);
  }
  void SquareRoot() {
    std::for_each(std::begin(data_), std::end(data_),
                  [](double &val) { val = sqrt(val); });
  }

  // move constructor and assignment operator
  varArray(varArray&&) noexcept = default;
  varArray& operator=(varArray&&) noexcept = default;

  // copy constructor and assignment operator
  varArray(const varArray&) = default;
  varArray& operator=(const varArray&) = default;

  // operator overloads
  const double & operator[](const int &r) const { return data_[r]; }
  double & operator[](const int &r) { return data_[r]; }

  template <typename T>
  inline varArray & operator+=(const T &);
  template <typename T>
  inline varArray & operator-=(const T &);
  template <typename T>
  inline varArray & operator*=(const T &);
  template <typename T>
  inline varArray & operator/=(const T &);

  inline varArray & operator+=(const double &);
  inline varArray & operator-=(const double &);
  inline varArray & operator*=(const double &);
  inline varArray & operator/=(const double &);

  inline varArray operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline varArray operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline varArray operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline varArray operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // destructor
  virtual ~varArray() noexcept {}
};

// function declarations --------------------------------------
// operator overload for addition
template <typename T>
varArray & varArray::operator+=(const T &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "array types must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] += arr[rr];
  }
  return *this;
}

// operator overload for subtraction
template <typename T>
varArray & varArray::operator-=(const T &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "array types must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] -= arr[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
varArray & varArray::operator*=(const T &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "array types must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] *= arr[rr];
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
varArray & varArray::operator/=(const T &arr) {
  MSG_ASSERT(this->Size() == arr.Size(), "array types must be same size");
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] /= arr[rr];
  }
  return *this;
}

inline const varArray operator+(varArray lhs, const varArray &rhs) {
  return lhs += rhs;
}

inline const varArray operator-(varArray lhs, const varArray &rhs) {
  return lhs -= rhs;
}

inline const varArray operator*(varArray lhs, const varArray &rhs) {
  return lhs *= rhs;
}

inline const varArray operator/(varArray lhs, const varArray &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
varArray & varArray::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
varArray & varArray::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
varArray & varArray::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
varArray & varArray::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const varArray operator+(const double &lhs, varArray rhs) {
  return rhs += lhs;
}

inline const varArray operator-(const double &lhs, varArray rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs[rr] = lhs - rhs[rr];
  }
  return rhs;
}

inline const varArray operator*(const double &lhs, varArray rhs) {
  return rhs *= lhs;
}

inline const varArray operator/(const double &lhs, varArray rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs[rr] = lhs / rhs[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const varArray &);

// --------------------------------------------------------------------------
// Derived wrapper classes for different variable types
// --------------------------------------------------------------------------
class primative : public varArray {
 public:
  // constructor
  primative(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}

  // member functions
  const double & RhoN(const int &ii) const { return this->SpeciesN(ii); }
  const double & U() const { return this->MomentumX(); }
  const double & V() const { return this->MomentumY(); }
  const double & W() const { return this->MomentumZ(); }
  const double & P() const { return this->Energy(); }
  const double & Tke() const { return this->TurbulenceN(0); }
  const double & Omega() const { return this->TurbulenceN(1); }

  // move constructor and assignment operator
  primative(primative&&) noexcept = default;
  primative& operator=(primative&&) noexcept = default;

  // copy constructor and assignment operator
  primative(const primative&) = default;
  primative& operator=(const primative&) = default;

  // destructor
  ~primative() noexcept {}
};

ostream &operator<<(ostream &os, const primative &);

// --------------------------------------------------------------------------
class conserved : public varArray {
 public:
  // constructor
  conserved(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}

  // member functions
  const double & RhoN(const int &ii) const { return this->SpeciesN(ii); }
  const double & RhoU() const { return this->MomentumX(); }
  const double & RhoV() const { return this->MomentumY(); }
  const double & RhoW() const { return this->MomentumZ(); }
  const double & RhoE() const { return this->Energy(); }
  const double & RhoTke() const { return this->TurbulenceN(0); }
  const double & RhoOmega() const { return this->TurbulenceN(1); }

  // move constructor and assignment operator
  conserved(conserved&&) noexcept = default;
  conserved& operator=(conserved&&) noexcept = default;

  // copy constructor and assignment operator
  conserved(const conserved&) = default;
  conserved& operator=(const conserved&) = default;

  // destructor
  ~conserved() noexcept {}
};

ostream &operator<<(ostream &os, const conserved &);

// --------------------------------------------------------------------------
class residual : public varArray {
 public:
  // constructor
  residual(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}

  // member functions
  const double & MassN(const int &ii) const { return this->SpeciesN(ii); }
  const double & MomentumX() const { return this->MomentumX(); }
  const double & MomentumY() const { return this->MomentumY(); }
  const double & MomentumZ() const { return this->MomentumZ(); }
  const double & Energy() const { return this->Energy(); }
  const double & TurbulenceN(const int &ii) const { 
    return this->TurbulenceN(ii); 
  }

  // move constructor and assignment operator
  residual(residual&&) noexcept = default;
  residual& operator=(residual&&) noexcept = default;

  // copy constructor and assignment operator
  residual(const residual&) = default;
  residual& operator=(const residual&) = default;

  // destructor
  virtual ~residual() noexcept {}
};

ostream &operator<<(ostream &os, const residual &);

// --------------------------------------------------------------------------
class source : public residual {
 public:
  // constructor
  source(const int &numEqns, const int &numSpecies)
      : residual(numEqns, numSpecies) {}

  // move constructor and assignment operator
  source(source&&) noexcept = default;
  source& operator=(source&&) noexcept = default;

  // copy constructor and assignment operator
  source(const source&) = default;
  source& operator=(const source&) = default;

  // destructor
  ~source() noexcept {}
};

ostream &operator<<(ostream &os, const source &);

// --------------------------------------------------------------------------
class flux : public residual {
 public:
  // constructor
  flux(const int &numEqns, const int &numSpecies)
      : residual(numEqns, numSpecies) {}

  // move constructor and assignment operator
  flux(flux&&) noexcept = default;
  flux& operator=(flux&&) noexcept = default;

  // copy constructor and assignment operator
  flux(const flux&) = default;
  flux& operator=(const flux&) = default;

  // destructor
  ~flux() noexcept {}
};

ostream &operator<<(ostream &os, const flux &);


#endif

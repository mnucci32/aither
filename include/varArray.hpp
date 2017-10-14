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
#include <iterator>
#include <cmath>
#include <type_traits>
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
  varArray(const int &numEqns, const int &numSpecies, const double &val)
      : data_(numEqns, val),
        momentumIndex_(numSpecies),
        energyIndex_(momentumIndex_ + 3),
        turbulenceIndex_(energyIndex_ + 1) {
    MSG_ASSERT(numEqns > numSpecies && numEqns >= 5,
               "number of equations should be greater than number of species");
  }
  varArray(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies, 0.0) {}
  varArray() {}
  varArray(const vector<double>::const_iterator &b,
           const vector<double>::const_iterator &e, const int &numSpecies)
      : data_(b, e),
        momentumIndex_(numSpecies),
        energyIndex_(momentumIndex_ + 3),
        turbulenceIndex_(energyIndex_ + 1) {
    MSG_ASSERT(data_.size() > numSpecies && data_.size() >= 5,
               "number of equations should be greater than number of species");
  }

  // move constructor and assignment operator
  varArray(varArray&&) noexcept = default;
  varArray& operator=(varArray&&) noexcept = default;

  // copy constructor and assignment operator
  varArray(const varArray&) = default;
  varArray& operator=(const varArray&) = default;

  // member functions
  int Size() const { return data_.size(); }
  int NumSpecies() const { return momentumIndex_; }
  int NumTurbulence() const { return this->Size() - turbulenceIndex_; }
  bool IsMultiSpecies() const { return this->NumSpecies() > 1; }
  bool HasTurbulenceData() const { return this->Size() != turbulenceIndex_; }
  int MomentumXIndex() const { return momentumIndex_; }
  int MomentumYIndex() const { return momentumIndex_ + 1; }
  int MomentumZIndex() const { return momentumIndex_ + 2; }
  int EnergyIndex() const { return energyIndex_; }
  int TurbulenceIndex() const { return turbulenceIndex_; }
  double SpeciesSum() const {
    return std::accumulate(std::begin(data_),
                           std::begin(data_) + this->NumSpecies(), 0.0);
  }

  const double &SpeciesN(const int &ii) const {
    MSG_ASSERT(ii < momentumIndex_, "requesting species variable out of range");
    return (*this)[ii];
  }
  const double &MomentumX() const { return (*this)[momentumIndex_]; }
  const double &MomentumY() const { return (*this)[momentumIndex_ + 1]; }
  const double &MomentumZ() const { return (*this)[momentumIndex_ + 2]; }
  const double &Energy() const { return (*this)[energyIndex_]; }
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
  bool IsZero() const {
    return std::all_of(std::begin(data_), std::end(data_),
                       [](const double &val) { return val == 0.0; });
  }
  varArray Squared() const;

  // provide begin and end so std::begin and std::end can be used
  // use lower case to conform with std::begin, std::end
  auto begin() noexcept {return data_.begin();}
  const auto begin() const noexcept {return data_.begin();}
  auto end() noexcept {return data_.end();}
  const auto end() const noexcept {return data_.end();}

  // operator overloads
  const double & operator[](const int &r) const { return data_[r]; }
  double & operator[](const int &r) { return data_[r]; }

  // destructor
  virtual ~varArray() noexcept {}
};

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// function declarations ----------------------------------------------------
// operator overload for addition
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator+=(T &lhs, const T &rhs) {
  MSG_ASSERT(lhs.Size() == rhs.Size(), "array types must be same size");
  for (auto rr = 0; rr < lhs.Size(); ++rr) {
    lhs[rr] += rhs[rr];
  }
  return lhs;
}

// operator overload for different type addition - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
& operator+=(T1 &lhs, const T2 &rhs) {
  return lhs += dynamic_cast<const T1 &>(rhs);
}

// operator overload for different type addition - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
& operator+=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<T2 &>(lhs) += rhs;
}

// operator overload for different type addition - T2/T1 not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
& operator+=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<varArray &>(lhs) += dynamic_cast<const varArray &>(rhs);
}

// --------------------------------------------------------------------------
// operator overload for subtraction
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator-=(T &lhs, const T &rhs) {
  MSG_ASSERT(lhs.Size() == rhs.Size(), "array types must be same size");
  for (auto rr = 0; rr < lhs.Size(); ++rr) {
    lhs[rr] -= rhs[rr];
  }
  return lhs;
}

// operator overload for different type subtraction - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
& operator-=(T1 &lhs, const T2 &rhs) {
  return lhs -= dynamic_cast<const T1 &>(rhs);
}

// operator overload for different type subtraction - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
& operator-=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<T2 &>(lhs) -= rhs;
}

// operator overload for different type subtraction - T2/T1 not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
& operator-=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<varArray &>(lhs) -= dynamic_cast<const varArray &>(rhs);
}

// --------------------------------------------------------------------------
// operator overload for elementwise multiplication
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator*=(T &lhs, const T &rhs) {
  MSG_ASSERT(lhs.Size() == rhs.Size(), "array types must be same size");
  for (auto rr = 0; rr < lhs.Size(); ++rr) {
    lhs[rr] *= rhs[rr];
  }
  return lhs;
}

// operator overload for different type multiplication - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
& operator*=(T1 &lhs, const T2 &rhs) {
  return lhs *= dynamic_cast<const T1 &>(rhs);
}

// operator overload for different type multiplication - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
& operator*=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<T2 &>(lhs) *= rhs;
}

// operator overload for different type addition - T2/T1 not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
& operator*=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<varArray &>(lhs) *= dynamic_cast<const varArray &>(rhs);
}

// --------------------------------------------------------------------------
// operator overload for elementwise division
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator/=(T &lhs, const T &rhs) {
  MSG_ASSERT(lhs.Size() == rhs.Size(), "array types must be same size");
  for (auto rr = 0; rr < lhs.Size(); ++rr) {
    lhs[rr] /= rhs[rr];
  }
  return lhs;
}

// operator overload for different type division - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
& operator/=(T1 &lhs, const T2 &rhs) {
  return lhs /= dynamic_cast<const T1 &>(rhs);
}

// operator overload for different type division - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
& operator/=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<T2 &>(lhs) /= rhs;
}

// operator overload for different type division - T2/T1 not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
& operator/=(T1 &lhs, const T2 &rhs) {
  return dynamic_cast<varArray &>(lhs) /= dynamic_cast<const varArray &>(rhs);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// operator overload for same type addition
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator+(T lhs, const T &rhs) {
  return lhs += rhs;
}

// operator overload for different type addition - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
operator+(T1 lhs, const T2 &rhs) {
  return lhs += rhs;
}

// operator overload for different type addition - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
operator+(T1 lhs, const T2 &rhs) {
  return lhs += rhs;
}

// operator overload for different type addition - T1/T2 are not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
operator+(T1 lhs, const T2 &rhs) {
  return lhs += rhs;
}

// operator overload for same type subtraction
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator-(T lhs, const T &rhs) {
  return lhs -= rhs;
}

// operator overload for different type subtraction - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
operator-(T1 lhs, const T2 &rhs) {
  return lhs -= rhs;
}

// operator overload for different type subtraction - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
operator-(T1 lhs, const T2 &rhs) {
  return lhs -= rhs;
}

// operator overload for different type subtraction - T1/T2 are not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
operator-(T1 lhs, const T2 &rhs) {
  return lhs -= rhs;
}


// operator overload for same type multiplication
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator*(T lhs, const T &rhs) {
  return lhs *= rhs;
}

// operator overload for different type multiplication - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
operator*(T1 lhs, const T2 &rhs) {
  return lhs *= rhs;
}

// operator overload for different type multiplication - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
operator*(T1 lhs, const T2 &rhs) {
  return lhs *= rhs;
}

// operator overload for different type multiplication - T1/T2 are not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
operator*(T1 lhs, const T2 &rhs) {
  return lhs *= rhs;
}

// operator overload for same type division
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator/(T lhs, const T &rhs) {
  return lhs /= rhs;
}

// operator overload for different type division - T1 is base of T2
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T1, T2>::value, T1>
operator/(T1 lhs, const T2 &rhs) {
  return lhs /= rhs;
}

// operator overload for different type division - T2 is base of T1
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T1>::value, T2>
operator/(T1 lhs, const T2 &rhs) {
  return lhs /= rhs;
}

// operator overload for different type division - T1/T2 are not base/derived
template <typename T1, typename T2>
inline const typename std::enable_if_t<std::is_base_of<varArray, T1>::value && 
                                       std::is_base_of<varArray, T2>::value &&
                                       !(std::is_base_of<T1, T2>::value ||
                                         std::is_base_of<T2, T1>::value), varArray>
operator/(T1 lhs, const T2 &rhs) {
  return lhs /= rhs;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// operator overloads for double --------------------------------------------
// operator overload for addition
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator+=(T &lhs, const double &scalar) {
  for (auto ii = 0; ii < lhs.Size(); ++ii) {
    lhs[ii] += scalar;
  }
  return lhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline T operator+(T lhs, const double &s) {
  return lhs += s;
}

// operator overload for subtraction with a scalar
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator-=(T &lhs, const double &scalar) {
    for (auto ii = 0; ii < lhs.Size(); ++ii) { 
    lhs[ii] -= scalar;
  }
  return lhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline T operator-(T lhs, const double &s) {
  return lhs -= s;
}

// operator overload for elementwise multiplication
template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator*=(T & lhs, const double &scalar) {
  for (auto ii = 0; ii < lhs.Size(); ++ii) { 
    lhs[ii] *= scalar;
  }
  return lhs;
}

template <typename T, 
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline T operator*(T lhs, const double &s) {
  return lhs *= s;
}

// operator overload for elementwise division
template <typename T, 
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
T & operator/=(T &lhs, const double &scalar) {
  for (auto ii = 0; ii < lhs.Size(); ++ii) {
    lhs[ii] /= scalar;
  }
  return lhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline T operator/(T lhs, const double &s) {
  return lhs /= s;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator+(const double &lhs, T rhs) {
  return rhs += lhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator-(const double &lhs, T rhs) {
  for (auto rr = 0; rr < rhs.Size(); ++rr) {
    rhs[rr] = lhs - rhs[rr];
  }
  return rhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator*(const double &lhs, T rhs) {
  return rhs *= lhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
inline const T operator/(const double &lhs, T rhs) {
  for (auto rr = 0; rr < rhs.Size(); ++rr) {
    rhs[rr] = lhs / rhs[rr];
  }
  return rhs;
}

template <typename T,
          typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
ostream &operator<<(ostream &os, const T &m) {
  for (auto rr = 0; rr < m.Size(); ++rr) {
    os << m[rr] << std::endl;
  }
  return os;
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// Derived wrapper classes for different variable types
// --------------------------------------------------------------------------
class residual : public varArray {
 public:
  // constructor
  residual(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}
  residual(const vector<double>::const_iterator &b, 
           const vector<double>::const_iterator &e, const int &numSpecies)
      : varArray(b, e, numSpecies) {}

  // move constructor and assignment operator
  residual(residual&&) noexcept = default;
  residual& operator=(residual&&) noexcept = default;

  // copy constructor and assignment operator
  residual(const residual&) = default;
  residual& operator=(const residual&) = default;

  // member functions
  const double & MassN(const int &ii) const { return this->SpeciesN(ii); }
  void GlobalReduceMPI(const int &);

  // destructor
  virtual ~residual() noexcept {}
};


#endif

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

#ifndef ARRAYVIEWHEADERDEF
#define ARRAYVIEWHEADERDEF

#include <iostream>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include <memory>
#include "macros.hpp"
#include "varArray.hpp"
#include "vector3d.hpp"
#include "thermodynamic.hpp"
#include "eos.hpp"
#include "conserved.hpp"

using std::ostream;
using std::vector;
using std::endl;
using std::unique_ptr;

/* class to store a view of an array. This is useful to slice out data from a
 * std::vector.
 */

// forward class declarations
class primitive;
class residual;
class turbModel;

template <typename T1, typename T2>
class arrayView {
  static_assert(std::is_base_of<varArray, T1>::value,
                "arrayView<T1, T2> requires T1 to be varArray type!");
  static_assert(std::is_arithmetic<T2>::value,
                "arrayView<T1, T2> requires T2 to be an arithmetic type!");

  typename vector<T2>::const_iterator begin_;
  typename vector<T2>::const_iterator end_;
  int momentumIndex_;
  int energyIndex_;
  int turbulenceIndex_;

 public:
  // constructor
  arrayView(const typename vector<T2>::const_iterator &b,
            const typename vector<T2>::const_iterator &e, const int &numSpecies)
      : begin_(b),
        end_(e),
        momentumIndex_(numSpecies),
        energyIndex_(momentumIndex_ + 3),
        turbulenceIndex_(energyIndex_ + 1) {}

  // move constructor and assignment operator
  arrayView(arrayView&&) noexcept = default;
  arrayView& operator=(arrayView&&) noexcept = default;

  // copy constructor and assignment operator
  arrayView(const arrayView &) = default;
  arrayView &operator=(const arrayView &) = default;

  // member functions
  T2 Sum() const { return std::accumulate(begin_, end_, T2(0)); }
  T1 CopyData() const { return T1{begin_, end_, this->NumSpecies()}; }
  T1 GetViewType() const { return T1(this->Size(), this->NumSpecies(), 0.0); }
  T1 Squared() const { return (*this) * (*this); }
  auto Size() const { return std::distance(begin_, end_); }
  auto begin() const { return begin_; }
  auto end() const { return end_; }
  int NumSpecies() const { return momentumIndex_; }
  int NumTurbulence() const { return this->Size() - turbulenceIndex_; }
  bool IsMultiSpecies() const { return this->NumSpecies() > 1; }
  bool HasTurbulenceData() const { return this->Size() != turbulenceIndex_; }
  int MomentumXIndex() const { return momentumIndex_; }
  int MomentumYIndex() const { return momentumIndex_ + 1; }
  int MomentumZIndex() const { return momentumIndex_ + 2; }
  int EnergyIndex() const { return energyIndex_; }
  int TurbulenceIndex() const { return turbulenceIndex_; }
  T2 SpeciesSum() const {
    return std::accumulate(begin_, begin_ + this->NumSpecies(), T2(0));
  }
  const T2 &SpeciesN(const int &ii) const {
    MSG_ASSERT(ii < momentumIndex_, "requesting species variable out of range");
    return (*this)[ii];
  }
  const T2 &MomentumX() const { return (*this)[momentumIndex_]; }
  const T2 &MomentumY() const { return (*this)[momentumIndex_ + 1]; }
  const T2 &MomentumZ() const { return (*this)[momentumIndex_ + 2]; }
  const T2 &Energy() const { return (*this)[energyIndex_]; }
  const T2 &TurbulenceN(const int &ii) const {
    MSG_ASSERT(tubulenceIndex_ + ii >= this->Size(),
               "requesting turbulence variable out of range");
    return (*this)[turbulenceIndex_ + ii];
  }

  arrayView<T1, T2> GetView() const { return *this; }

  // --------------------------------------------------------------------------
  // getters for primitives ---------------------------------------------------
  const T2 &RhoN(const int &ii) const {
    static_assert(std::is_same<primitive, T1>::value ||
                      std::is_same<conserved, T1>::value,
                  "getter only valid for primitive/conserved type!");
    return this->SpeciesN(ii); 
  }
  T2 Rho() const {
    static_assert(std::is_same<primitive, T1>::value ||
                      std::is_same<conserved, T1>::value,
                  "getter only valid for primitive/conserved type!");
    return this->SpeciesSum();
  }
  T2 MassFractionN(const int &ii) const {
    static_assert(std::is_same<primitive, T1>::value ||
                      std::is_same<conserved, T1>::value,
                  "getter only valid for primitive/conserved type!");
    return this->RhoN(ii) / this->Rho(); 
  }
  vector<T2> RhoVec() const {
    static_assert(std::is_same<primitive, T1>::value ||
                      std::is_same<conserved, T1>::value,
                  "getter only valid for primitive/conserved type!");
    return {this->begin(), this->begin() + this->NumSpecies()};
  }
  vector<T2> MassFractions() const {
    static_assert(std::is_same<primitive, T1>::value ||
                      std::is_same<conserved, T1>::value,
                  "getter only valid for primitive/conserved type!");
    vector<T2> mf(this->NumSpecies());
    const auto totalRho = this->Rho();
    for (auto ii = 0U; ii < mf.size(); ++ii) {
      mf[ii] = this->RhoN(ii) / totalRho;
    }
    return mf;
  }
  const T2 &U() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->MomentumX(); 
  }
  const T2 &V() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->MomentumY(); 
  }
  const T2 &W() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->MomentumZ(); 
  }
  const T2 &P() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->Energy(); 
  }
  const T2 &Tke() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->TurbulenceN(0); 
  }
  const T2 &Omega() const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->TurbulenceN(1); 
  }
  const T2 &TurbN(const int &ii) const { 
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return this->TurbulenceN(ii); 
  }
  vector3d<T2> Velocity() const {
    static_assert(std::is_same<primitive, T1>::value,
                "getter only valid for primitive type!");
    return {this->U(), this->V(), this->W()};
  }
  T2 SoS(const unique_ptr<thermodynamic> &thermo,
         const unique_ptr<eos> &eqnState) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "getter only valid for primitive type!");
    return sqrt(
        thermo->Gamma(this->Temperature(eqnState), this->MassFractions()) *
        this->P() / this->Rho());
  }
  T2 Temperature(const unique_ptr<eos> &eqnState) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "getter only valid for primitive type!");
    return eqnState->Temperature(this->P(), this->RhoVec());
  }
  T2 Energy(const unique_ptr<eos> &eqnState,
            const unique_ptr<thermodynamic> &thermo) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "getter only valid for primitive type!");
    const auto t = this->Temperature(eqnState);
    return eqnState->Energy(
        eqnState->SpecEnergy(thermo, t, this->MassFractions()),
        this->Velocity().Mag());
  }
  T2 Enthalpy(const unique_ptr<eos> &eqnState,
              const unique_ptr<thermodynamic> &thermo) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "getter only valid for primitive type!");
    const auto t = this->Temperature(eqnState);
    return eqnState->Enthalpy(thermo, t, this->Velocity().Mag(),
                              this->MassFractions());
  }
  conserved ConsVars(const unique_ptr<eos> &eqnState,
                     const unique_ptr<thermodynamic> &thermo) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "function only valid for primitive type!");
    return PrimToCons((*this), eqnState, thermo);
  }
  T1 UpdateWithConsVars(const unique_ptr<eos> &eqnState,
                        const unique_ptr<thermodynamic> &thermo,
                        const arrayView<varArray, double> &du,
                        const unique_ptr<turbModel> &turb) const {
    static_assert(std::is_same<primitive, T1>::value,
                  "function only valid for primitive type!");
    return UpdatePrimWithCons((*this), eqnState, thermo, du, turb);
  }

  // --------------------------------------------------------------------------
  // operator overloads
  const T2 & operator[](const int &r) const { return *(begin_ + r); }

  inline T1 operator+(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs += s;
  }
  inline T1 operator-(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs -= s;
  }
  inline T1 operator*(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs *= s;
  }
  inline T1 operator/(const T2 &s) const {
    auto lhs = this->CopyData();
    return lhs /= s;
  }

  // destructor
  ~arrayView() noexcept {}
};

// function declarations --------------------------------------
template <typename T1, typename T2>
inline const T1 operator+(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll + rhs;
}

template <typename T1, typename T2>
inline const T1 operator-(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll - rhs;
}

template <typename T1, typename T2>
inline const T1 operator*(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll * rhs;
}

template <typename T1, typename T2>
inline const T1 operator/(const arrayView<T1, T2> &lhs,
                          const arrayView<T1, T2> &rhs) {
  auto ll = lhs.CopyData();
  return ll / rhs;
}

// operator overloads for type T -------------------------------------
template <typename T1, typename T2>
inline const T1 operator+(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return result += lhs;
}

template <typename T1, typename T2>
inline const T1 operator-(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return lhs - result;
}

template <typename T1, typename T2>
inline const T1 operator*(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return result *= lhs;
}

template <typename T1, typename T2>
inline const T1 operator/(const T2 &lhs, const arrayView<T1, T2> &rhs) {
  auto result = rhs.CopyData();
  return lhs / rhs;
}

// operator overloads for varArray type T -------------------------------------
template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator+(const T3 &lhs, const arrayView<T1, T2> &rhs) {
  T3 result(rhs.begin(), rhs.end(), rhs.NumSpecies());
  return lhs + result;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator-(const T3 &lhs, const arrayView<T1, T2> &rhs) {
  T3 result(rhs.begin(), rhs.end(), rhs.NumSpecies());
  return lhs - result;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator*(const T3 &lhs, const arrayView<T1, T2> &rhs) {
  T3 result(rhs.begin(), rhs.end(), rhs.NumSpecies());
  return lhs * result;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator/(const T3 &lhs, const arrayView<T1, T2> &rhs) {
  T3 result(rhs.begin(), rhs.end(), rhs.NumSpecies());
  return lhs / result;
}

// ---------------------------------------------------------------------------
template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator+(const arrayView<T1, T2> &lhs, const T3 &rhs) {
  T3 result(lhs.begin(), lhs.end(), lhs.NumSpecies());
  return result += rhs;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator-(const arrayView<T1, T2> &lhs, const T3 &rhs) {
  T3 result(lhs.begin(), lhs.end(), lhs.NumSpecies());
  return result -= rhs;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator*(const arrayView<T1, T2> &lhs, const T3 &rhs) {
  T3 result(lhs.begin(), lhs.end(), lhs.NumSpecies());
  return result *= rhs;
}

template <typename T1, typename T2, typename T3>
inline const typename std::enable_if_t<
    std::is_base_of<varArray, T3>::value, T3>
operator/(const arrayView<T1, T2> &lhs, const T3 &rhs) {
  T3 result(lhs.begin(), lhs.end(), lhs.NumSpecies());
  return result /= rhs;
}


// ---------------------------------------------------------------------------
// operation overload for << - allows use of cout, cerr, etc.
template <typename T1, typename T2>
ostream &operator<<(ostream &os, const arrayView<T1, T2> &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// using typedefs
using varArrayView = arrayView<varArray, double>;
using primitiveView = arrayView<primitive, double>;
using conservedView = arrayView<conserved, double>;
using residualView = arrayView<residual, double>;

#endif

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

#ifndef MATMULTIARRAY3DHEADERDEF
#define MATMULTIARRAY3DHEADERDEF

/* This file contains the header and implementation for a multidimensional (3D)
   array class. The class is to act as a container to store flux jacobian
   matrices and provide easy access to elements with i, j, k indexing.
 */

#include <iostream>  // ostream
#include <vector>    // vector
#include <algorithm>
#include <functional>
#include "multiArray3d.hpp"
#include "fluxJacobian.hpp"
#include "varArray.hpp"
#include "arrayView.hpp"

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::string;

class matMultiArray3d : public multiArray3d<double> {
  int flowSize_;
  int turbSize_;

  // private member functions
  auto BeginFlow(const int &ii, const int &jj, const int &kk) noexcept {
    return this->begin() + this->GetLoc1D(ii, jj, kk);
  }
  const auto BeginFlow(const int &ii, const int &jj, const int &kk) const
      noexcept {
    return this->begin() + this->GetLoc1D(ii, jj, kk);
  }
  auto BeginTurb(const int &ii, const int &jj, const int &kk) noexcept {
    return this->begin() + this->GetLoc1D(ii, jj, kk) +
           flowSize_ * flowSize_;
  }
  const auto BeginTurb(const int &ii, const int &jj, const int &kk) const
      noexcept {
    return this->begin() + this->GetLoc1D(ii, jj, kk) +
           flowSize_ * flowSize_;
  }

 public:
  // constructor
  matMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const fluxJacobian &init)
      : multiArray3d<double>(ii, jj, kk, ng, init.Size()),
        flowSize_(init.FlowSize()), turbSize_(init.TurbSize()) {
    for (auto kk = this->StartK(); kk < this->EndK(); ++kk) {
      for (auto jj = this->StartJ(); jj < this->EndJ(); ++jj) {
        for (auto ii = this->StartI(); ii < this->EndI(); ++ii) {
          std::copy(init.begin(), init.end(), this->BeginFlow(ii, jj, kk));
        }
      }
    }
  }
  matMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const int &fs, const int &ts)
      : multiArray3d<double>(ii, jj, kk, ng, bs),
        flowSize_(fs),
        turbSize_(ts) {}
  matMultiArray3d() : multiArray3d<double>(), flowSize_(0), turbSize_(0) {}

  // move constructor and assignment operator
  matMultiArray3d(matMultiArray3d &&) noexcept = default;
  matMultiArray3d &operator=(matMultiArray3d &&) noexcept = default;

  // copy constructor and assignment operator
  matMultiArray3d(const matMultiArray3d &) = default;
  matMultiArray3d &operator=(const matMultiArray3d &) = default;

  // member functions
  int FlowSize() const { return flowSize_; }
  int TurbSize() const { return turbSize_; }
  bool IsScalar() const { return flowSize_ == 1; }
  void Zero() { std::fill(this->begin(), this->end(), 0.0); }

  void MultiplyOnDiagonal(const int &ii, const int &jj, const int &kk,
                          const double &fac) {
    MultiplyFacOnDiagonal(this->BeginFlow(ii, jj, kk), flowSize_, fac);
    MultiplyFacOnDiagonal(this->BeginTurb(ii, jj, kk), turbSize_, fac);
  }
  void AddOnDiagonal(const int &ii, const int &jj, const int &kk,
                     const double &fac) {
    AddFacOnDiagonal(this->BeginFlow(ii, jj, kk), flowSize_, fac);
    AddFacOnDiagonal(this->BeginTurb(ii, jj, kk), turbSize_, fac);
  }
  void Inverse(const int &ii, const int &jj, const int &kk) {
    MatrixInverse(this->BeginFlow(ii, jj, kk), flowSize_);
    MatrixInverse(this->BeginTurb(ii, jj, kk), turbSize_);
  }

  void Add(const int &ii, const int &jj, const int &kk,
           const fluxJacobian &jac) {
    MSG_ASSERT(jac.Size() == this->BlockSize(), "block size must match");
    std::transform(this->BeginFlow(ii, jj, kk),
                   this->BeginFlow(ii, jj, kk) + this->BlockSize(), jac.begin(),
                   this->BeginFlow(ii, jj, kk), std::plus<double>());
  }
  void Subtract(const int &ii, const int &jj, const int &kk,
           const fluxJacobian &jac) {
    MSG_ASSERT(jac.Size() == this->BlockSize(), "block size must match");
    std::transform(this->BeginFlow(ii, jj, kk),
                   this->BeginFlow(ii, jj, kk) + this->BlockSize(), jac.begin(),
                   this->BeginFlow(ii, jj, kk), std::minus<double>());
  }
  void SubtractFromTurb(const int &ii, const int &jj, const int &kk,
                        const squareMatrix &jac) {
    MSG_ASSERT(jac.Size() == this->TurbSize(), "block size must match");
    std::transform(
        this->BeginTurb(ii, jj, kk),
        this->BeginTurb(ii, jj, kk) + this->TurbSize() * this->TurbSize(),
        jac.begin(), this->BeginTurb(ii, jj, kk), std::minus<double>());
  }

  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value>>
  T ArrayMult(const int &ii, const int &jj, const int &kk,
              const T &orig) const {
    T arr(orig.Size(), orig.NumSpecies());
    ArrayMultiplication(this->BeginFlow(ii, jj, kk), flowSize_, turbSize_,
                        this->IsScalar(), orig, arr);
    return arr;
  }
  template <typename T,
            typename = std::enable_if_t<std::is_same<varArrayView, T>::value ||
                                        std::is_same<primitiveView, T>::value ||
                                        std::is_same<conservedView, T>::value ||
                                        std::is_same<residualView, T>::value>>
  auto ArrayMult(const int &ii, const int &jj, const int &kk,
                 const T &arrView) const {
    auto arr = arrView.GetViewType();
    ArrayMultiplication(this->BeginFlow(ii, jj, kk), flowSize_, turbSize_,
                        this->IsScalar(), arrView, arr);
    return arr;
  }

  void ClearResize(const int &ii, const int &jj, const int &kk,
                   const int &ng, const int &bs, const int &fs, const int &ts) {
    *this = matMultiArray3d(ii, jj, kk, ng, bs, fs, ts);
  }
  void ClearResize(const int &ii, const int &jj, const int &kk, const int &ng,
                   const fluxJacobian &val) {
    *this = matMultiArray3d(ii, jj, kk, ng, val);
  }

  // destructor
  ~matMultiArray3d() noexcept {}
};

// ---------------------------------------------------------------------------
// member function definitions
ostream &operator<<(ostream &os, const matMultiArray3d &arr);

#endif

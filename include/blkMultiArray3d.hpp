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

#ifndef BLKMULTIARRAY3DHEADERDEF
#define BLKMULTIARRAY3DHEADERDEF

/* This file contains the header and implementation for a multidimensional (3D)
   array class. The class is to act as a container to store all data types,
   and provide easy access to elements with i, j, k indexing.
 */

#include <iostream>  // ostream
#include <vector>    // vector
#include <string>    // string
#include <memory>    // unique_ptr
#include "mpi.h"
#include "vector3d.hpp"
#include "boundaryConditions.hpp"  // connection
#include "range.hpp"  // range
#include "multiArray3d.hpp"
#include "varArray.hpp"
#include "arrayView.hpp"

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;

template <typename T>
class blkMultiArray3d : public multiArray3d<double> {
  static_assert(std::is_base_of<varArray, T>::value,
                "blkMultiArray3d<T> requires a varArray based type!");

  int numSpecies_;

  // private member functions
  int GetBlkLoc1D(const int &ii, const int &jj, const int &kk) const {
    return this->BlockSize() * this->GetLoc1D(ii, jj, kk);
  }

 public:
  // constructor
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const int &ns, const double &init)
      : multiArray3d<double>(ii, jj, kk, ng, bs, init), numSpecies_(ns) {}
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const int &ns)
      : multiArray3d<double>(ii, jj, kk, ng, bs), numSpecies_(ns) {}
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs)
      : blkMultiArray3d(ii, jj, kk, ng, bs, 0) {}
  blkMultiArray3d() : multiArray3d<double>(), numSpecies_(0) {}

  // move constructor and assignment operator
  blkMultiArray3d(blkMultiArray3d&&) noexcept = default;
  blkMultiArray3d& operator=(blkMultiArray3d&&) noexcept = default;

  // copy constructor and assignment operator
  blkMultiArray3d(const blkMultiArray3d&) = default;
  blkMultiArray3d& operator=(const blkMultiArray3d&) = default;

  // member functions
  // operator overloads
  arrayView<T, double> operator()(const int &ii, const int &jj, const int &kk) const {
    auto start = this->GetBlkLoc1D(ii, jj, kk);
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }
    
  double &operator()(const int &ii, const int &jj, const int &kk,
                const int &bb) {
    MSG_ASSERT(bb <= this->BlockSize(), "accessing data out of block index");
    auto ind = this->GetBlkLoc1D(ii, jj, kk) + bb;
    return *ind;
  }
  const double &operator()(const int &ii, const int &jj, const int &kk,
                      const int &bb) const {
    MSG_ASSERT(bb <= this->BlockSize(), "accessing data out of block index");
    auto ind = this->GetBlkLoc1D(ii, jj, kk) + bb;
    return *ind;
  }

  arrayView<T, double> operator()(const int &ind) const {
    auto start = ind * this->BlockSize();
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }
  //const T &operator()(const int &ind) const {
  //  auto start = ind * this->BlockSize();
  //  return {this->begin() + start, this->begin() + start + this->BlockSize(),
  //          numSpecies_};
  //}

  // destructor
  ~blkMultiArray3d() noexcept {}
};

// ---------------------------------------------------------------------------
// member function definitions

// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const blkMultiArray3d<T> &arr) {
  os << "Size: " << arr.NumI() << ", " << arr.NumJ() << ", " << arr.NumK()
     << endl;
  os << "Number of ghost layers: " << arr.GhostLayers() << endl;

  for (auto kk = arr.StartK(); kk < arr.EndK(); kk++) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); jj++) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ii++) {
        os << ii << ", " << jj << ", " << kk << ": " << arr(ii, jj, kk) << endl;
      }
    }
  }
  return os;
}

#endif

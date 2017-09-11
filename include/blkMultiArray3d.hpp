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
#include "arrayView.hpp"
#include "multiArray3d.hpp"

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;

template <typename T>
class blkMultiArray3d : public multiArray3d<T> {

  // private member functions
  int GetBlkLoc1D(const int &ii, const int &jj, const int &kk) const {
    return this->BlockSize() * this->GetLoc1D(ii, jj, kk);
  }

 public:
  // constructor
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const T &init)
      : multiArray3d<T>(ii, jj, kk, ng, bs, init) {}
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs)
      : multiArray3d<T>(ii, jj, kk, ng, bs) {}
  blkMultiArray3d() : multiArray3d<T>() {}

  // move constructor and assignment operator
  blkMultiArray3d(blkMultiArray3d&&) noexcept = default;
  blkMultiArray3d& operator=(blkMultiArray3d&&) noexcept = default;

  // copy constructor and assignment operator
  blkMultiArray3d(const blkMultiArray3d&) = default;
  blkMultiArray3d& operator=(const blkMultiArray3d&) = default;

  // member functions
  arrayView<T> GetElem(const int &ii, const int &jj, const int &kk) const;

  // operator overloads
  arrayView<T> &operator()(const int &ii, const int &jj, const int &kk) {
    auto start = this->GetBlkLoc1D(ii, jj, kk);
    return {this->begin() + start, this->begin() + start + this->BlockSize()};
  }
  const arrayView<T> &operator()(const int &ii, const int &jj,
                                 const int &kk) const {
    auto start = this->GetBlkLoc1D(ii, jj, kk);
    return {this->begin() + start, this->begin() + start + this->BlockSize()};
  }

  arrayView<T> &operator()(const int &ind) {
    auto start = ind * this->BlockSize();
    return {this->begin() + start, this->begin() + start + this->BlockSize()};
  }
  const arrayView<T> &operator()(const int &ind) const { 
    auto start = ind * this->BlockSize();
    return {this->begin() + start, this->begin() + start + this->BlockSize()};
  }

  arrayView<T> &operator()(const string &dir, const int &d1, const int &d2,
                           const int &d3) {
    if (dir == "i") {  // direction 1 is i
      return (*this)(d1, d2, d3);
    } else if (dir == "j") {  // direction 1 is j
      return (*this)(d3, d1, d2);
    } else if (dir == "k") {  // direction 1 is k
      rreturn (*this)(d2, d3, d1);
    } else {
      cerr << "ERROR: Direction " << dir << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  const arrayView<T> &operator()(const string &dir, const int &d1,
                                 const int &d2, const int &d3) const {
    if (dir == "i") {  // direction 1 is i
      return (*this)(d1, d2, d3);
    } else if (dir == "j") {  // direction 1 is j
      return (*this)(d3, d1, d2);
    } else if (dir == "k") {  // direction 1 is k
      return (*this)(d2, d3, d1);
    } else {
      cerr << "ERROR: Direction " << dir << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // destructor
  ~blkMultiArray3d() noexcept {}
};

// ---------------------------------------------------------------------------
// member function definitions

template <typename T>
arrayView<T> blkMultiArray3d<T>::GetElem(const int &ii, const int &jj,
                                         const int &kk) const {
  if (ii < this->EndI() && jj < this->EndJ() && kk < this->Endk() &&
      ii >= this->StartI() && jj >= this->StartJ() && kk >= this->StartK()) {
    return (*this)(ii, jj, kk);
  } else {
    cerr << "ERROR: Tried to access location outside of bounds of "
         << "multiArray3d" << endl;
    cerr << "Tried to access " << ii << ", " << jj << ", " << kk << endl;
    cerr << "Maximum locations are " << this->EndI() - 1 << ", "
         << this->EndJ() - 1 << ", " << this->EndK() - 1 << endl;
    exit(EXIT_FAILURE);
  }
}

// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const blkMultiArray3d<T> &arr) {
  os << "Size: " << arr.NumI() << ", " << arr.NumJ() << ", "
     << arr.NumK() << endl;
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

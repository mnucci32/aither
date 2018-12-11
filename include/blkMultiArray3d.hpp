/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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
#include <type_traits>
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
                "blkMultiArray3d<T> requires a varArray based type");

  int numSpecies_;

 public:
  // constructor
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const T &init)
      : multiArray3d<double>(ii, jj, kk, ng, init.Size()),
        numSpecies_(init.NumSpecies()) {
    for (auto kk = this->StartK(); kk < this->EndK(); kk++) {
      for (auto jj = this->StartJ(); jj < this->EndJ(); jj++) {
        for (auto ii = this->StartI(); ii < this->EndI(); ii++) {
          for (auto bb = 0; bb < this->BlockSize(); ++bb) {
            (*this)(ii, jj, kk, bb) = init[bb];
          }
        }
      }
    }
  }
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const int &ns, const double &init)
      : multiArray3d<double>(ii, jj, kk, ng, bs, init), numSpecies_(ns) {}
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const int &bs, const int &ns)
      : multiArray3d<double>(ii, jj, kk, ng, bs), numSpecies_(ns) {}
  blkMultiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
                  const std::pair<int, int> &info)
      : multiArray3d<double>(ii, jj, kk, ng, info),
        numSpecies_(info.second) {}
  blkMultiArray3d() : multiArray3d<double>(), numSpecies_(0) {}

  // move constructor and assignment operator
  blkMultiArray3d(blkMultiArray3d &&) noexcept = default;
  blkMultiArray3d &operator=(blkMultiArray3d &&) noexcept = default;

  // copy constructor and assignment operator
  blkMultiArray3d(const blkMultiArray3d &) = default;
  blkMultiArray3d &operator=(const blkMultiArray3d &) = default;

  // member functions
  std::pair<int, int> BlockInfo() const override {
    return std::make_pair(this->BlockSize(), numSpecies_);
  }

  template <typename T2>
  void InsertBlock(const int &ii, const int &jj, const int &kk, const T2 &arr) {
    static_assert(std::is_same<T, T2>::value ||
                      std::is_same<arrayView<T, double>, T2>::value,
                  "blkMultiArray3d<T> requires a varArray based type!");
    MSG_ASSERT(arr.Size() == this->BlockSize(), "block size must match");
    auto start = this->GetLoc1D(ii, jj, kk);
    std::copy(arr.begin(), arr.end(), this->begin() + start);
  }

  template <typename T2>
  void InsertBlock(const string &dir, const int &d1, const int &d2,
                   const int &d3, const T2 &arr) {
    static_assert(std::is_same<T, T2>::value ||
                      std::is_same<arrayView<T, double>, T2>::value,
                  "blkMultiArray3d<T> requires a varArray based type!");
    MSG_ASSERT(arr.Size() == this->BlockSize(), "block size must match");

    auto start = 0;
    if (dir == "i") {  // direction 1 is i
      start = this->GetLoc1D(d1, d2, d3);
    } else if (dir == "j") {  // direction 1 is j
      start = this->GetLoc1D(d3, d1, d2);
    } else if (dir == "k") {  // direction 1 is k
      start = this->GetLoc1D(d2, d3, d1);
    } else {
      cerr << "ERROR: Direction " << dir << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
    std::copy(arr.begin(), arr.end(), this->begin() + start);
  }

  blkMultiArray3d<T> Slice(const range &, const range &, const range &) const;
  blkMultiArray3d<T> Slice(const string &, const range &,
                           const bool = false) const;
  blkMultiArray3d<T> Slice(const string &, int, int, const bool = false,
                           const string = "cell", const bool = false,
                           const bool = false) const;
  blkMultiArray3d<T> Slice(const string &, int, range, range,
                           const string = "cell", const int = 0) const;

  template <typename TT>
  void Insert(const range &, const range &, const range &, const TT &);
  template <typename TT>
  void Insert(const string &, const range &, const TT &, const bool = false);
  template <typename TT>
  void Insert(const string &, int, int, const TT &, const bool = false,
              const string = "cell", const bool = false, const bool = false);
  template <typename TT>
  void Insert(const string &, int, range, range, const TT &,
              const string = "cell", const int = 0);

  template <typename TT>
  void SwapSlice(const connection &conn, TT &array);
  void SwapSliceMPI(const connection &conn, const int &rank,
                    const MPI_Datatype &MPI_arrData, const int tag = 1);
  template <typename TT>
  void PutSlice(const TT &, const connection &, const int &);

  T GetCopy(const int &ii, const int &jj, const int &kk) const {
    auto start = this->GetLoc1D(ii, jj, kk);
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }

  // operator overloads
  arrayView<T, double> operator()(const int &ii, const int &jj,
                                  const int &kk) const {
    auto start = this->GetLoc1D(ii, jj, kk);
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }
  arrayView<T, double> operator()(const int &ind) const {
    auto start = ind * this->BlockSize();
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }
  arrayView<T, double> operator()(const string &dir, const int &d1,
                                  const int &d2, const int &d3) const {
    auto start = 0;
    if (dir == "i") {  // direction 1 is i
      start = this->GetLoc1D(d1, d2, d3);
    } else if (dir == "j") {  // direction 1 is j
      start = this->GetLoc1D(d3, d1, d2);
    } else if (dir == "k") {  // direction 1 is k
      start = this->GetLoc1D(d2, d3, d1);
    } else {
      cerr << "ERROR: Direction " << dir << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
    return {this->begin() + start, this->begin() + start + this->BlockSize(),
            numSpecies_};
  }
  double &operator()(const int &ii, const int &jj, const int &kk,
                     const int &bb) {
    MSG_ASSERT(bb <= this->BlockSize(), "accessing data out of block index");
    return (*this)[this->GetLoc1D(ii, jj, kk) + bb];
  }
  const double &operator()(const int &ii, const int &jj, const int &kk,
                           const int &bb) const {
    MSG_ASSERT(bb <= this->BlockSize(), "accessing data out of block index");
    return (*this)[this->GetLoc1D(ii, jj, kk) + bb];
  }

  void ClearResize(const int &ii, const int &jj, const int &kk,
                   const int &ng, const int &bs, const int &ns) {
    *this = blkMultiArray3d<T>(ii, jj, kk, ng, bs, ns);
  }
  void ClearResize(const int &ii, const int &jj, const int &kk, const int &ng,
                   const int &bs, const int &ns, const T &val) {
    *this = blkMultiArray3d<T>(ii, jj, kk, ng, bs, ns, val);
  }

  blkMultiArray3d<T> GrowI() const;
  blkMultiArray3d<T> GrowJ() const;
  blkMultiArray3d<T> GrowK() const;

  // arithmetic with same type
  template <typename T2>
  inline typename std::enable_if_t<std::is_base_of<varArray, T>::value &&
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T>::value,
                                   blkMultiArray3d<T>>
      &operator+=(const blkMultiArray3d<T2> &arr) {
    for (auto rr = 0; rr < this->Size(); ++rr) {
      (*this)[rr] += arr[rr];
    }
    return *this;
  }

  template <typename T2>
  inline typename std::enable_if_t<std::is_base_of<varArray, T>::value &&
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T>::value,
                                   blkMultiArray3d<T>>
      &operator-=(const blkMultiArray3d<T2> &arr) {
    for (auto rr = 0; rr < this->Size(); rr++) {
      (*this)[rr] -= arr[rr];
    }
    return *this;
  }

  template <typename T2>
  inline typename std::enable_if_t<std::is_base_of<varArray, T>::value &&
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T>::value,
                                   blkMultiArray3d<T>>
      &operator*=(const blkMultiArray3d<T2> &arr) {
    for (auto rr = 0; rr < this->Size(); rr++) {
      (*this)[rr] *= arr[rr];
    }
    return *this;
  }

  template <typename T2>
  inline typename std::enable_if_t<std::is_base_of<varArray, T>::value &&
                                       std::is_base_of<varArray, T2>::value &&
                                       std::is_base_of<T2, T>::value,
                                   blkMultiArray3d<T>>
      &operator/=(const blkMultiArray3d<T2> &arr) {
    for (auto rr = 0; rr < this->Size(); rr++) {
      (*this)[rr] /= arr[rr];
    }
    return *this;
  }

  // arithmetic with scalar
  inline blkMultiArray3d<T> & operator+=(const double &);
  inline blkMultiArray3d<T> & operator-=(const double &);
  inline blkMultiArray3d<T> & operator*=(const double &);
  inline blkMultiArray3d<T> & operator/=(const double &);


  // destructor
  ~blkMultiArray3d() noexcept {}
};

// ---------------------------------------------------------------------------
// Same type operator overloads
template <typename T>
inline const blkMultiArray3d<T> operator+(blkMultiArray3d<T> lhs,
                                          const blkMultiArray3d<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator-(blkMultiArray3d<T> lhs,
                                          const blkMultiArray3d<T> &rhs) {
  return lhs -= rhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator*(blkMultiArray3d<T> lhs,
                                          const blkMultiArray3d<T> &rhs) {
  return lhs *= rhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator/(blkMultiArray3d<T> lhs,
                                          const blkMultiArray3d<T> &rhs) {
  return lhs /= rhs;
}

// -------------------------------------------------------------------------
// operator overload for addition
template <typename T>
blkMultiArray3d<T> & blkMultiArray3d<T>::operator+=(const double &scalar) {
  for (auto rr = 0; rr < this->Size(); ++rr) {
    (*this)[rr] += scalar;
  }
  return *this;
}

// operator overload for subtraction
template <typename T>
blkMultiArray3d<T> & blkMultiArray3d<T>::operator-=(const double &scalar) {
  for (auto rr = 0; rr < this->Size(); ++rr) {
    (*this)[rr] -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
blkMultiArray3d<T> & blkMultiArray3d<T>::operator*=(const double &scalar) {
  for (auto rr = 0; rr < this->Size(); ++rr) {
    (*this)[rr] *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
blkMultiArray3d<T> & blkMultiArray3d<T>::operator/=(const double &scalar) {
  for (auto rr = 0; rr < this->Size(); ++rr) {
    (*this)[rr] /= scalar;
  }
  return *this;
}

template <typename T>
inline const blkMultiArray3d<T> operator+(const double &lhs,
                                          blkMultiArray3d<T> rhs) {
  return rhs += lhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator-(const double &lhs,
                                          blkMultiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); ++rr) {
    rhs(rr) = lhs - rhs(rr);
  }
  return rhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator*(const double &lhs,
                                          blkMultiArray3d<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const blkMultiArray3d<T> operator/(const double &lhs,
                                          blkMultiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); ++rr) {
    rhs(rr) = lhs / rhs(rr);
  }
  return rhs;
}

// ----------------------------------------------------------------


// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const blkMultiArray3d<T> &arr) {
  os << "Size: " << arr.NumI() << ", " << arr.NumJ() << ", " << arr.NumK()
     << endl;
  os << "Number of ghost layers: " << arr.GhostLayers() << endl;

  for (auto kk = arr.StartK(); kk < arr.EndK(); ++kk) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); ++jj) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ++ii) {
        os << ii << ", " << jj << ", " << kk << ": " << arr(ii, jj, kk) << endl;
      }
    }
  }
  return os;
}

// member function definitions
// member function to return a slice of the array
// this is the main slice function that all other overloaded slice functions call
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::Slice(const range &ir, const range &jr,
                                             const range &kr) const {
  // ir -- i-index range to take slice [inclusive, exclusive)
  // jr -- j-index range to take slice [inclusive, exclusive)
  // kr -- k-index range to take slice [inclusive, exclusive)
  return SliceArray((*this), ir, jr, kr);
}

// member function to return a slice of the array
// Overload to slice only in one direction. Given a 3D array, this slice returns
// a plane with normal direction dir, or a smaller 3D array where the direction
// dir is sliced over dirRange. It also has the ability to include or ignore
// ghost cells in its planar slices
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::Slice(const string &dir,
                                             const range &dirRange,
                                             const bool physOnly) const {
  // dir -- direction of slice
  // dirRange -- range of slice in direction given
  // phsOnly -- flag to only include physical cells in the two directions that
  //            are not specified as dir
  return SliceArray((*this), dir, dirRange, physOnly);
}

// member function to return a slice of the array
// overload to slice line out of array
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::Slice(const string &dir, int d2Ind,
                                             int d3Ind, const bool physOnly,
                                             const string id, const bool upper2,
                                             const bool upper3) const {
  // dir -- direction of line slice (direction 1)
  // d2Ind -- index of direction 2
  // d3Ind -- index of direction 3
  // physOnly -- flag to only include physical cells in line slice
  // id -- type of multiArray3d being sliced: cell, i, j, or k
  //       d2Ind and d3Ind are supplied as cell indices, but may need to be
  //       altered if the array is storing i, j, or k face data
  // upper2 -- flag to determine if direction 2 is at upper index
  // upper3 -- flag to determine if direction 3 is at upper index
  return SliceArray((*this), dir, d2Ind, d3Ind, physOnly, id, upper2, upper3);
}

// overload to slice plane out of array
// Identical to previous slice overload, but more general in that in can slice
// over a subset of direction 2 & 3. This is useful to slice out a plane that
// borders a boundary condition patch.
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::Slice(const string &dir, int dirInd,
                                             range dir1, range dir2,
                                             const string id,
                                             const int type) const {
  // dir -- normal direction of planar slice
  // dirInd -- index in normal direction
  // dir1 -- range of direction 1 (direction 3 is normal to slice)
  // dir2 -- range of direction 2 (direction 3 is normal to slice)
  // id -- id of array being sliced (i, j, k for faces, cell for cells)
  // type -- surface type of dir
  return SliceArray((*this), dir, dirInd, dir1, dir2, id, type);
}

// member function to insert an array into this one
// this is the main insert funciton that all other overloaded insert functions
// call
template <typename T>
template <typename TT>
void blkMultiArray3d<T>::Insert(const range &ir, const range &jr,
                                const range &kr, const TT &arr) {
  // ir -- i-index range to take slice [inclusive, exclusive)
  // jr -- j-index range to take slice [inclusive, exclusive)
  // kr -- k-index range to take slice [inclusive, exclusive)
  // arr -- array to insert into this one
  InsertArray((*this), ir, jr, kr, arr);
}

// Overload to insert only in one direction. Given a 3D array, this inserts a
// plane with normal direction dir, or a smaller 3D array where the direction
// dir is inserted over dirRange. It also has the ability to include or ignore
// ghost cells in its planar inserts
template <typename T>
template <typename TT>
void blkMultiArray3d<T>::Insert(const string &dir, const range &dirRange,
                                const TT &arr, const bool physOnly) {
  // dir -- direction of slice to insert
  // dirRange -- range to insert slice into in direction given
  // arr -- array to insert
  // phsOnly -- flag to only include physical cells in the two directions that
  //            are not specified as dir
  InsertArray((*this), dir, dirRange, arr, physOnly);
}

// overload to insert line into array
template <typename T>
template <typename TT>
void blkMultiArray3d<T>::Insert(const string &dir, int d2Ind, int d3Ind,
                                const TT &arr, const bool physOnly,
                                const string id, const bool upper2,
                                const bool upper3) {
  // dir -- direction of line slice to insert (direction 1)
  // d2Ind -- index of direction 2 to insert into
  // d3Ind -- index of direction 3 to insert into
  // physOnly -- flag to only include physical cells in line insert
  // id -- type of multiArray3d being sliced: cell, i, j, or k
  //       d2Ind and d3Ind are supplied as cell indices, but may need to be
  //       altered if the array is storing i, j, or k face data
  // upper2 -- flag to determine if direction 2 is at upper index
  // upper3 -- flag to determine if direction 3 is at upper index
  InsertArray((*this), dir, d2Ind, d3Ind, arr, physOnly, id, upper2, upper3);
}

// overload to insert plane into array
// Identical to previous insert overload, but more general in that in can insert
// over a subset of direction 2 & 3. This is useful to insert into a plane that
// borders a boundary condition patch.
template <typename T>
template <typename TT>
void blkMultiArray3d<T>::Insert(const string &dir, int dirInd, range dir1,
                                range dir2, const TT &arr, const string id,
                                const int type) {
  // dir -- normal direction of planar slice
  // dirInd -- index in normal direction
  // dir1 -- range of direction 1 (direction 3 is normal to slice)
  // dir2 -- range of direction 2 (direction 3 is normal to slice)
  // arr -- array to insert
  // id -- id of array being inserted into (i, j, k for faces, cell for cells)
  // type -- surface type of dir
  InsertArray((*this), dir, dirInd, dir1, dir2, arr, id, type);
}

template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::GrowI() const {
  return GrowInI(*this);
}
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::GrowJ() const {
  return GrowInJ(*this);
}
template <typename T>
blkMultiArray3d<T> blkMultiArray3d<T>::GrowK() const {
  return GrowInK(*this);
}

/* Function to swap slice using MPI. This is similar to the SwapSlice
   function, but is called when the neighboring procBlocks are on different
   processors.
*/
template <typename T>
void blkMultiArray3d<T>::SwapSliceMPI(const connection &conn, const int &rank,
                                      const MPI_Datatype &MPI_arrData,
                                      const int tag) {
  // conn -- connection boundary information
  // rank -- processor rank
  // MPI_arrData -- MPI datatype for passing data in *this
  // tag -- id for MPI swap (default 1)
  SwapSliceParallel((*this), conn, rank, MPI_arrData, tag);
}

/* Function to swap ghost cells between two blocks at an connection
boundary. Slices are removed from the physical cells (extending into ghost cells
at the edges) of one block and inserted into the ghost cells of its partner
block. The reverse is also true. The slices are taken in the coordinate system
orientation of their parent block.

   Interior Cells    Ghost Cells               Ghost Cells   Interior Cells
   ________ ______|________ _________       _______________|_______ _________
Ui-3/2   Ui-1/2   |    Uj+1/2    Uj+3/2  Ui-3/2    Ui-1/2  |    Uj+1/2    Uj+3/2
  |        |      |        |         |     |        |      |       |         |
  | Ui-1   |  Ui  |  Uj    |  Uj+1   |     |  Ui-1  |   Ui |  Uj   |  Uj+1   |
  |        |      |        |         |     |        |      |       |         |
  |________|______|________|_________|     |________|______|_______|_________|
                  |                                        |

The above diagram shows the resulting values after the ghost cell swap. The
logic ensures that the ghost cells at the connection boundary exactly match
their partner block as if there were no separation in the grid.
*/
template <typename T>
template <typename TT>
void blkMultiArray3d<T>::SwapSlice(const connection &conn, TT &array) {
  // conn -- connection boundary information
  // array -- second array involved in connection boundary
  SwapSliceLocal((*this), conn, array);
}

template <typename T>
template <typename TT>
void blkMultiArray3d<T>::PutSlice(const TT &array, const connection &inter,
                                  const int &d3) {
  // array -- array to insert into *this
  // inter -- connection data structure defining patches and orientation
  // d3 -- distance of direction normal to patch to insert
  InsertSlice((*this), array, inter, d3);
}


#endif

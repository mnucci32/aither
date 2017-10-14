/*  This file is part of aither.
 
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

#ifndef MULTIARRAY3DHEADERDEF  // only if the macro MULTIARRAY3DHEADERDEF
                               // is not defined execute these lines of code

#define MULTIARRAY3DHEADERDEF  // define the macro

/* This file contains the header and implementation for a multidimensional (3D)
   array class. The class is to act as a container to store all data types,
   and provide easy access to elements with i, j, k indexing.
 */

#include <iostream>  // ostream
#include <vector>    // vector
#include <string>    // string
#include <memory>    // unique_ptr
#include <utility>   // pair
#include <type_traits>
#include "mpi.h"
#include "vector3d.hpp"
#include "boundaryConditions.hpp"  // connection
#include "range.hpp"  // range
#include "macros.hpp"  // MSG_ASSERT

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;

template <typename T>
class multiArray3d {
  vector<T> data_;
  int numI_;
  int numJ_;
  int numK_;
  int numGhosts_;
  int blkSize_;

 public:
  // constructor
  multiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
               const int &bs, const T &init) :
      data_(bs * (ii + 2 * ng) * (jj + 2 * ng) * (kk + 2 * ng), init),
      numI_(ii + 2 * ng), numJ_(jj + 2 * ng), numK_(kk + 2 * ng),
      numGhosts_(ng), blkSize_(bs) {}
  multiArray3d(const int &ii, const int &jj, const int &kk, const int &ng, 
               const int &bs=1) :
      data_(bs * (ii + 2 * ng) * (jj + 2 * ng) * (kk + 2 * ng)),
      numI_(ii + 2 * ng), numJ_(jj + 2 * ng), numK_(kk + 2 * ng),
      numGhosts_(ng), blkSize_(bs) {}
  multiArray3d(const int &ii, const int &jj, const int &kk, const int &ng,
               const std::pair<int, int> &info)
      : multiArray3d(ii, jj, kk, ng, info.first) {}
  multiArray3d() : multiArray3d(0, 0, 0, 0, 1) {}

  // move constructor and assignment operator
  multiArray3d(multiArray3d&&) noexcept = default;
  multiArray3d& operator=(multiArray3d&&) noexcept = default;

  // copy constructor and assignment operator
  multiArray3d(const multiArray3d&) = default;
  multiArray3d& operator=(const multiArray3d&) = default;

  // member functions
  int Size() const {return data_.size();}
  int BlockSize() const { return blkSize_; }
  int NumBlocks() const {return this->Size() / this->BlockSize();}
  virtual std::pair<int, int> BlockInfo() const {
    return std::make_pair(this->BlockSize(), 0);
  }
  bool IsEmpty() const { return data_.empty(); }
  int NumI() const {return numI_;}
  int NumJ() const {return numJ_;}
  int NumK() const {return numK_;}
  int NumINoGhosts() const {return numI_ - 2 * numGhosts_;}
  int NumJNoGhosts() const {return numJ_ - 2 * numGhosts_;}
  int NumKNoGhosts() const {return numK_ - 2 * numGhosts_;}
  int GhostLayers() const {return numGhosts_;}
  int StartI() const {return -numGhosts_;}
  int StartJ() const {return -numGhosts_;}
  int StartK() const {return -numGhosts_;}
  int EndI() const {return numI_ - numGhosts_;}
  int EndJ() const {return numJ_ - numGhosts_;}
  int EndK() const {return numK_ - numGhosts_;}
  int Start(const string &) const;
  int End(const string &) const;
  int GetLoc1D(const int &ii, const int &jj, const int &kk) const {
    return (ii + numGhosts_) + (jj + numGhosts_) * numI_ +
        (kk + numGhosts_) * numI_ * numJ_;
  }

  int PhysStartI() const {return 0;}
  int PhysStartJ() const {return 0;}
  int PhysStartK() const {return 0;}
  int PhysEndI() const {return this->NumINoGhosts();}
  int PhysEndJ() const {return this->NumJNoGhosts();}
  int PhysEndK() const {return this->NumKNoGhosts();}
  int PhysStart(const string &) const;
  int PhysEnd(const string &) const;

  int PhysicalSize() const {
    return this->NumINoGhosts() * this->NumJNoGhosts() * this->NumKNoGhosts();
  }
  range RangeI() const {return {this->StartI(), this->EndI()};}
  range RangeJ() const {return {this->StartJ(), this->EndJ()};}
  range RangeK() const {return {this->StartK(), this->EndK()};}
  range PhysRangeI() const {return {this->PhysStartI(), this->PhysEndI()};}
  range PhysRangeJ() const {return {this->PhysStartJ(), this->PhysEndJ()};}
  range PhysRangeK() const {return {this->PhysStartK(), this->PhysEndK()};}

  // provide begin and end so std::begin and std::end can be used
  // use lower case to conform with std::begin, std::end
  auto begin() noexcept {return data_.begin();}
  const auto begin() const noexcept {return data_.begin();}
  auto end() noexcept {return data_.end();}
  const auto end() const noexcept {return data_.end();}

  multiArray3d<T> Slice(const range &, const range &, const range &) const;
  multiArray3d<T> Slice(const string &, const range &,
                        const bool = false) const;
  multiArray3d<T> Slice(const string &, int, int, const bool = false,
                        const string = "cell", const bool = false,
                        const bool = false) const;
  multiArray3d<T> Slice(const string &, int, range, range,
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
  void Fill(const TT &);
  template <typename TT>
  void PutSlice(const TT &, const connection &, const int &);
  void SwapSliceMPI(const connection &, const int &, const MPI_Datatype &,
                    const int = 1);
  template <typename TT>
  void SwapSlice(const connection &, TT &);

  void Zero(const T &z) { std::fill(this->begin(), this->end(), z); }
  void Zero() { this->Zero(T()); }

  multiArray3d<T> GrowI() const;
  multiArray3d<T> GrowJ() const;
  multiArray3d<T> GrowK() const;

  void PackSwapUnpackMPI(const connection &, const MPI_Datatype &, const int &,
                         const int = 1);

  T GetElem(const int &ii, const int &jj, const int &kk) const;

  void InsertBlock(const int &ii, const int &jj, const int &kk, const T &val) {
    (*this)(ii, jj, kk) = val;
  }
  void InsertBlock(const string &dir, const int &d1, const int &d2,
                   const int &d3, const T &val) {
    (*this)(dir, d1, d2, d3) = val;
  }

  // operator overloads
  T &operator()(const int &ii, const int &jj, const int &kk) {
    MSG_ASSERT(this->BlockSize() == 1,
               "shouldn't be used with blkMultiArray3d");
    return data_[this->GetLoc1D(ii, jj, kk)];
  }
  const T &operator()(const int &ii, const int &jj, const int &kk) const {
    MSG_ASSERT(this->BlockSize() == 1,
               "shouldn't be used with blkMultiArray3d");
    return data_[this->GetLoc1D(ii, jj, kk)];
  }
  T &operator()(const int &ind) {
    MSG_ASSERT(this->BlockSize() == 1,
               "shouldn't be used with blkMultiArray3d");
    return data_[ind];
  }
  const T &operator()(const int &ind) const {
    MSG_ASSERT(this->BlockSize() == 1,
               "shouldn't be used with blkMultiArray3d");
    return data_[ind];
  }
  T &operator()(const string &dir, const int &d1, const int &d2,
                const int &d3) {
    MSG_ASSERT(this->BlockSize() == 1, "shouldn't be used with blkMultiArray3d");
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
  const T& operator()(const string &dir, const int &d1, const int &d2,
                      const int &d3) const {
    MSG_ASSERT(this->BlockSize() == 1, "shouldn't be used with blkMultiArray3d");
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

  T &operator[](const int &ind) {
    return data_[ind];
  }
  const T &operator[](const int &ind) const {
    return data_[ind];
  }


  // arithmetic with same type
  inline multiArray3d<T> & operator+=(const multiArray3d<T> &);
  inline multiArray3d<T> & operator-=(const multiArray3d<T> &);
  inline multiArray3d<T> & operator*=(const multiArray3d<T> &);
  inline multiArray3d<T> & operator/=(const multiArray3d<T> &);

  // arithmetic with type that *this is holding
  inline multiArray3d<T> & operator+=(const T &);
  inline multiArray3d<T> & operator-=(const T &);
  inline multiArray3d<T> & operator*=(const T &);
  inline multiArray3d<T> & operator/=(const T &);

  inline multiArray3d<T> operator+(const T &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline multiArray3d<T> operator-(const T &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline multiArray3d<T> operator*(const T &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline multiArray3d<T> operator/(const T &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  // arithmetic with a type that is not *this or the type *this is holding
  // used for example to multiply with a double if *this is a
  // multiArray3d<vector3d<double>>

  template <typename TT>
  multiArray3d<T> & operator+=(const TT &scalar) {
    for (auto &val : data_) {
      val += scalar;
    }
    return *this;
  }

  template <typename TT>
  multiArray3d<T> & operator-=(const TT &scalar) {
    for (auto &val : data_) {
      val -= scalar;
    }
    return *this;
  }

  template <typename TT>
  multiArray3d<T> & operator*=(const TT &scalar) {
    for (auto &val : data_) {
      val *= scalar;
    }
    return *this;
  }

  template <typename TT>
  multiArray3d<T> & operator/=(const TT &scalar) {
    for (auto &val : data_) {
      val /= scalar;
    }
    return *this;
  }

  template <typename TT>
  multiArray3d<T> operator+(
      std::enable_if_t<!std::is_same<T, TT>::value, TT> const &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  template <typename TT>
  multiArray3d<T> operator-(
      std::enable_if_t<!std::is_same<T, TT>::value, TT> const &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  template <typename TT>
  multiArray3d<T> operator*(
      std::enable_if_t<!std::is_same<T, TT>::value, TT> const &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  template <typename TT>
  multiArray3d<T> operator/(
      std::enable_if_t<!std::is_same<T, TT>::value, TT> const &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  void ClearResize(const int &ii, const int &jj, const int &kk, const int &ng) {
    *this = multiArray3d<T>(ii, jj, kk, ng, blkSize_);
  }
  void ClearResize(const int &ii, const int &jj, const int &kk, const int &ng,
                   const T &val) {
    *this = multiArray3d<T>(ii, jj, kk, ng, blkSize_, val);
  }

  void SameSizeResize(const int &ii, const int &jj, const int &kk);
  void SameSizeResizeGhosts(const int &ii, const int &jj, const int &kk,
                            const int &ng);

  // destructor
  virtual ~multiArray3d() noexcept {}
};

// ---------------------------------------------------------------------------
// non member functions
// main slice function that all other overloaded slice functions call
template <typename T>
auto SliceArray(const T &parent, const range &ir, const range &jr,
                const range &kr) {
  // ir -- i-index range to take slice [inclusive, exclusive)
  // jr -- j-index range to take slice [inclusive, exclusive)
  // kr -- k-index range to take slice [inclusive, exclusive)
  // T should be multiArray3d or blkMultiArray3d type

  // check that slice bounds are within parent array and that end bounds are
  // greater than or equal to start bounds
  if (!ir.IsInside(parent.RangeI()) || !jr.IsInside(parent.RangeJ()) ||
      !kr.IsInside(parent.RangeK()) || !ir.IsValid() || !jr.IsValid() ||
      !kr.IsValid()) {
    cerr << "ERROR: Error in SliceArray. Cannot take slice with "
         << "boundaries " << ir << ", " << jr << ", " << kr << endl
         << "from array with ranges " << parent.RangeI() << ", "
         << parent.RangeJ() << ", " << parent.RangeK() << endl;
    exit(EXIT_FAILURE);
  }

  // slices always have 0 ghost cells
  T arr(ir.Size(), jr.Size(), kr.Size(), 0, parent.BlockInfo());

  // s is for index of sliced array, p is for index of parent array
  for (int ks = arr.StartK(), kp = kr.Start(); ks < arr.EndK(); ks++, kp++) {
    for (int js = arr.StartJ(), jp = jr.Start(); js < arr.EndJ(); js++, jp++) {
      for (int is = arr.StartI(), ip = ir.Start(); is < arr.EndI(); is++, ip++) {
        arr.InsertBlock(is, js, ks, parent(ip, jp, kp));
      }
    }
  }
  return arr;
}

// Overload to slice only in one direction. Given a 3D array, this slice returns
// a plane with normal direction dir, or a smaller 3D array where the direction
// dir is sliced over dirRange. It also has the ability to include or ignore
// ghost cells in its planar slices
template <typename T>
auto SliceArray(const T &parent, const string &dir, const range &dirRange,
                const bool physOnly) {
  // dir -- direction of slice
  // dirRange -- range of slice in direction given
  // phsOnly -- flag to only include physical cells in the two directions that
  //            are not specified as dir

  if (dir == "i") {
    if (physOnly) {
      return parent.Slice(dirRange, parent.PhysRangeJ(), parent.PhysRangeK());
    } else {
      return parent.Slice(dirRange, parent.RangeJ(), parent.RangeK());
    }
  } else if (dir == "j") {
    if (physOnly) {
      return parent.Slice(parent.PhysRangeI(), dirRange, parent.PhysRangeK());
    } else {
      return parent.Slice(parent.RangeI(), dirRange, parent.RangeK());
    }
  } else if (dir == "k") {
    if (physOnly) {
      return parent.Slice(parent.PhysRangeI(), parent.PhysRangeJ(), dirRange);
    } else {
      return parent.Slice(parent.RangeI(), parent.RangeJ(), dirRange);
    }
  } else {
    cerr << "ERROR: Error in multiArray3d::Slice, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// function to return a slice of the array
// overload to slice line out of array
template <typename T>
auto SliceArray(const T &parent, const string &dir, int d2Ind, int d3Ind,
                const bool physOnly, const string id, const bool upper2,
                const bool upper3) {
  // dir -- direction of line slice (direction 1)
  // d2Ind -- index of direction 2
  // d3Ind -- index of direction 3
  // physOnly -- flag to only include physical cells in line slice
  // id -- type of multiArray3d being sliced: cell, i, j, or k
  //       d2Ind and d3Ind are supplied as cell indices, but may need to be
  //       altered if the array is storing i, j, or k face data
  // upper2 -- flag to determine if direction 2 is at upper index
  // upper3 -- flag to determine if direction 3 is at upper index

  if (dir == "i") {  // d2 = j, d3 = k
    if (upper2 && id == "j") {
      d2Ind++;
    } else if (upper3 && id == "k") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Slice(parent.PhysRangeI(), d2Ind, d3Ind);
    } else {
      return parent.Slice(parent.RangeI(), d2Ind, d3Ind);
    }
  } else if (dir == "j") {  // d2 = k, d3 = i
    if (upper2 && id == "k") {
      d2Ind++;
    } else if (upper3 && id == "i") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Slice(d3Ind, parent.PhysRangeJ(), d2Ind);
    } else {
      return parent.Slice(d3Ind, parent.RangeJ(), d2Ind);
    }
  } else if (dir == "k") {  // d2 = i, d3 = j
    if (upper2 && id == "i") {
      d2Ind++;
    } else if (upper3 && id == "j") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Slice(d2Ind, d3Ind, parent.PhysRangeK());
    } else {
      return parent.Slice(d2Ind, d3Ind, parent.RangeK());
    }
  } else {
    cerr << "ERROR: Error in multiArray3d::Slice, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// overload to slice plane out of array
// Identical to previous slice overload, but more general in that in can slice
// over a subset of direction 2 & 3. This is useful to slice out a plane that
// borders a boundary condition patch.
template <typename T>
auto SliceArray(const T &parent, const string &dir, int dirInd, range dir1,
                range dir2, const string id, const int type) {
  // dir -- normal direction of planar slice
  // dirInd -- index in normal direction
  // dir1 -- range of direction 1 (direction 3 is normal to slice)
  // dir2 -- range of direction 2 (direction 3 is normal to slice)
  // id -- id of array being sliced (i, j, k for faces, cell for cells)
  // type -- surface type of dir

  if (dir == "i") {  // d1 = j, d2 = k
    if (type == 2 && id == "i") {  // upper i-surface & i face data
      dirInd++;
    }
    if (id == "j") {
      dir1.GrowEnd();
    } else if (id == "k") {
      dir2.GrowEnd();
    }
    return parent.Slice(dirInd, dir1, dir2);
  } else if (dir == "j") {  // d1 = k, d2 = i
    if (type == 4 && id == "j") {  // upper j-surface & j face data
      dirInd++;
    }
    if (id == "k") {
      dir1.GrowEnd();
    } else if (id == "i") {
      dir2.GrowEnd();
    }
    return parent.Slice(dir2, dirInd, dir1);
  } else if (dir == "k") {  // d1 = i, d2 = j
    if (type == 6 && id == "k") {  // upper k-surface & k face data
      dirInd++;
    }
    if (id == "i") {
      dir1.GrowEnd();
    } else if (id == "j") {
      dir2.GrowEnd();
    }
    return parent.Slice(dir1, dir2, dirInd);
  } else {
    cerr << "ERROR: Error in multiArray3d::Slice, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// main insert funciton that all other overloaded insert functions call
template <typename T>
void InsertArray(T &parent, const range &ir, const range &jr, const range &kr,
                 const T &arr) {
  // ir -- i-index range to take slice [inclusive, exclusive)
  // jr -- j-index range to take slice [inclusive, exclusive)
  // kr -- k-index range to take slice [inclusive, exclusive)
  // arr -- array to insert into this one

#ifndef NDEBUG
  // check that given bounds fit in this, sizes match, and that bounds are valid
  if (!ir.IsInside(parent.RangeI()) || !jr.IsInside(parent.RangeJ()) ||
      !kr.IsInside(parent.RangeK()) || ir.Size() != arr.RangeI().Size() ||
      jr.Size() != arr.RangeJ().Size() || kr.Size() != arr.RangeK().Size() ||
      !ir.IsValid() || !jr.IsValid() || !kr.IsValid()) {
    cerr << "ERROR: Error in InsertArray. Given array does not fit in "
         << "given bounds" << endl
         << "Given bounds: " << ir << ", " << jr << ", " << kr << endl
         << "Bounds of array being inserted: " << arr.RangeI() << ", "
         << arr.RangeJ() << ", " << arr.RangeK() << endl;
    cerr << "Bounds of array accepting data: " << parent.RangeI() << ", "
         << parent.RangeJ() << ", " << parent.RangeK() << endl;
    exit(EXIT_FAILURE);
  }
#endif

  // s is for index of sliced array, p is for index of parent array
  for (int ks = arr.StartK(), kp = kr.Start(); ks < arr.EndK(); ks++, kp++) {
    for (int js = arr.StartJ(), jp = jr.Start(); js < arr.EndJ(); js++, jp++) {
      for (int is = arr.StartI(), ip = ir.Start(); is < arr.EndI(); is++, ip++) {
        parent.InsertBlock(ip, jp, kp, arr(is, js, ks));
      }
    }
  }
}

// Overload to insert only in one direction. Given a 3D array, this inserts a
// plane with normal direction dir, or a smaller 3D array where the direction
// dir is inserted over dirRange. It also has the ability to include or ignore
// ghost cells in its planar inserts
template <typename T>
void InsertArray(T &parent, const string &dir, const range &dirRange,
                 const T &arr, const bool physOnly) {
  // dir -- direction of slice to insert
  // dirRange -- range to insert slice into in direction given
  // arr -- array to insert
  // phsOnly -- flag to only include physical cells in the two directions that
  //            are not specified as dir

  if (dir == "i") {
    if (physOnly) {
      return parent.Insert(dirRange, parent.PhysRangeJ(), parent.PhysRangeK(),
                           arr);
    } else {
      return parent.Insert(dirRange, parent.RangeJ(), parent.RangeK(), arr);
    }
  } else if (dir == "j") {
    if (physOnly) {
      return parent.Insert(parent.PhysRangeI(), dirRange, parent.PhysRangeK(),
                           arr);
    } else {
      return parent.Insert(parent.RangeI(), dirRange, parent.RangeK(), arr);
    }
  } else if (dir == "k") {
    if (physOnly) {
      return parent.Insert(parent.PhysRangeI(), parent.PhysRangeJ(), dirRange,
                           arr);
    } else {
      return parent.Insert(parent.RangeI(), parent.RangeJ(), dirRange, arr);
    }
  } else {
    cerr << "ERROR: Error in InsertArray, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// overload to insert line into array
template <typename T>
void InsertArray(T &parent, const string &dir, int d2Ind, int d3Ind,
                 const T &arr, const bool physOnly, const string id,
                 const bool upper2, const bool upper3) {
  // dir -- direction of line slice to insert (direction 1)
  // d2Ind -- index of direction 2 to insert into
  // d3Ind -- index of direction 3 to insert into
  // physOnly -- flag to only include physical cells in line insert
  // id -- type of multiArray3d being sliced: cell, i, j, or k
  //       d2Ind and d3Ind are supplied as cell indices, but may need to be
  //       altered if the array is storing i, j, or k face data
  // upper2 -- flag to determine if direction 2 is at upper index
  // upper3 -- flag to determine if direction 3 is at upper index

  if (dir == "i") {  // d2 = j, d3 = k
    if (upper2 && id == "j") {
      d2Ind++;
    } else if (upper3 && id == "k") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Insert(parent.PhysRangeI(), d2Ind, d3Ind, arr);
    } else {
      return parent.Insert(parent.RangeI(), d2Ind, d3Ind, arr);
    }
  } else if (dir == "j") {  // d2 = k, d3 = i
    if (upper2 && id == "k") {
      d2Ind++;
    } else if (upper3 && id == "i") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Insert(d3Ind, parent.PhysRangeJ(), d2Ind, arr);
    } else {
      return parent.Insert(d3Ind, parent.RangeJ(), d2Ind, arr);
    }
  } else if (dir == "k") {  // d2 = i, d3 = j
    if (upper2 && id == "i") {
      d2Ind++;
    } else if (upper3 && id == "j") {
      d3Ind++;
    }

    if (physOnly) {
      return parent.Insert(d2Ind, d3Ind, parent.PhysRangeK(), arr);
    } else {
      return parent.Insert(d2Ind, d3Ind, parent.RangeK(), arr);
    }
  } else {
    cerr << "ERROR: Error in InsertArray, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// overload to insert plane into array
// Identical to previous insert overload, but more general in that in can insert
// over a subset of direction 2 & 3. This is useful to insert into a plane that
// borders a boundary condition patch.
template <typename T>
void InsertArray(T &parent, const string &dir, int dirInd, range dir1,
                 range dir2, const T &arr, const string id, const int type) {
  // dir -- normal direction of planar slice
  // dirInd -- index in normal direction
  // dir1 -- range of direction 1 (direction 3 is normal to slice)
  // dir2 -- range of direction 2 (direction 3 is normal to slice)
  // arr -- array to insert
  // id -- id of array being inserted into (i, j, k for faces, cell for cells)
  // type -- surface type of dir

  if (dir == "i") {  // d1 = j, d2 = k
    if (type == 2 && id == "i") {  // upper i-surface & i normal
      dirInd++;
    }
    if (id == "j") {
      dir1.GrowEnd();
    } else if (id == "k") {
      dir2.GrowEnd();
    }
    return parent.Insert(dirInd, dir1, dir2, arr);
  } else if (dir == "j") {  // d1 = k, d2 = i
    if (type == 4 && id == "j") {  // upper j-surface & j normal
      dirInd++;
    }
    if (id == "k") {
      dir1.GrowEnd();
    } else if (id == "i") {
      dir2.GrowEnd();
    }
    return parent.Insert(dir2, dirInd, dir1, arr);
  } else if (dir == "k") {  // d1 = i, d2 = j
    if (type == 6 && id == "k") {  // upper k-surface & k normal
      dirInd++;
    }
    if (id == "i") {
      dir1.GrowEnd();
    } else if (id == "j") {
      dir2.GrowEnd();
    }
    return parent.Insert(dir1, dir2, dirInd, arr);
  } else {
    cerr << "ERROR: Error in InsertArray, direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
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
void SwapSliceLocal(T &array1, const connection &conn, T &array2) {
  // array1 -- first array involved in connection boundary
  // conn -- connection boundary information
  // array2 -- second array involved in connection boundary

  // Get indices for slice coming from first block to swap
  auto is1 = 0, ie1 = 0;
  auto js1 = 0, je1 = 0;
  auto ks1 = 0, ke1 = 0;

  conn.FirstSliceIndices(is1, ie1, js1, je1, ks1, ke1, array1.GhostLayers());

  // Get indices for slice coming from second block to swap
  auto is2 = 0, ie2 = 0;
  auto js2 = 0, je2 = 0;
  auto ks2 = 0, ke2 = 0;

  conn.SecondSliceIndices(is2, ie2, js2, je2, ks2, ke2, array2.GhostLayers());

  // get slices to swap
  auto slice1 = array1.Slice({is1, ie1}, {js1, je1}, {ks1, ke1});
  auto slice2 = array2.Slice({is2, ie2}, {js2, je2}, {ks2, ke2});

  // change connections to work with slice and ghosts
  connection conn1 = conn;
  connection conn2 = conn;
  conn1.AdjustForSlice(false, array1.GhostLayers());
  conn2.AdjustForSlice(true, array2.GhostLayers());

  // put slices in proper blocks
  array1.PutSlice(slice2, conn2, array2.GhostLayers());
  array2.PutSlice(slice1, conn1, array1.GhostLayers());
}

/* Function to swap slice using MPI. This is similar to the SwapSlice
   function, but is called when the neighboring procBlocks are on different
   processors.
*/
template <typename T>
void SwapSliceParallel(T &array, const connection &conn, const int &rank,
                       const MPI_Datatype &MPI_arrData, const int tag) {
  // array -- array on local processor to swap
  // conn -- connection boundary information
  // rank -- processor rank
  // MPI_arrData -- MPI datatype for passing data in *this
  // tag -- id for MPI swap (default 1)

  // Get indices for slice coming from block to swap
  auto is = 0, ie = 0;
  auto js = 0, je = 0;
  auto ks = 0, ke = 0;

  if (rank == conn.RankFirst()) {  // local block first in connection
    conn.FirstSliceIndices(is, ie, js, je, ks, ke, array.GhostLayers());
  } else if (rank == conn.RankSecond()) {  // local block second in connection
    conn.SecondSliceIndices(is, ie, js, je, ks, ke, array.GhostLayers());
  } else {
    cerr << "ERROR: Error in SwapSliceParallel(). Processor rank does "
            "not match either of connection ranks!" << endl;
    exit(EXIT_FAILURE);
  }

  // get local state slice to swap
  auto slice = array.Slice({is, ie}, {js, je}, {ks, ke});

  // swap state slices with partner block
  slice.PackSwapUnpackMPI(conn, MPI_arrData, rank, tag);

  // change connections to work with slice and ghosts
  auto connAdj = conn;

  // change connections to work with slice and ghosts
  // block to insert into is first in connection
  if (rank == conn.RankFirst()) {
    connAdj.AdjustForSlice(true, array.GhostLayers());
  } else {  // block to insert into is second in connection, so pass swapped
            // version
    connAdj.AdjustForSlice(false, array.GhostLayers());
  }

  // insert state slice into procBlock
  array.PutSlice(slice, connAdj, array.GhostLayers());
}

template <typename T>
void InsertSlice(T &array1, const T &array2, const connection &inter,
                 const int &d3) {
  // array1 -- array to accept data
  // array2 -- array to insert into array1
  // inter -- connection data structure defining patches and orientation
  // d3 -- distance of direction normal to patch to insert

#ifndef NDEBUG
  // check that number of cells to insert matches
  auto blkCell = inter.Dir1LenFirst() * inter.Dir2LenFirst() * d3;
  if (blkCell != array2.NumBlocks()) {
    cerr << "ERROR: Error in InsertSlice(). Number of cells "
            "being inserted does not match designated space to insert."
         << endl;
    cerr << "Direction 1, 2, 3 of array to insert into: "
         << inter.Dir1LenFirst() << ", " << inter.Dir2LenFirst() << ", " << d3
         << endl;
    cerr << "Direction I, J, K of array3d to insert: " << array2.NumI() << ", "
         << array2.NumJ() << ", " << array2.NumK() << endl;
    exit(EXIT_FAILURE);
  }
#endif

  // adjust insertion indices if patch borders another connection on the same
  // surface of the block
  const auto adjS1 =
      (inter.Dir1StartInterBorderFirst()) ? array1.GhostLayers() : 0;
  const auto adjE1 =
      (inter.Dir1EndInterBorderFirst()) ? array1.GhostLayers() : 0;
  const auto adjS2 =
      (inter.Dir2StartInterBorderFirst()) ? array1.GhostLayers() : 0;
  const auto adjE2 =
      (inter.Dir2EndInterBorderFirst()) ? array1.GhostLayers() : 0;

  // loop over cells to insert
  for (auto l3 = 0; l3 < d3; l3++) {
    for (auto l2 = adjS2; l2 < inter.Dir2LenFirst() - adjE2; l2++) {
      for (auto l1 = adjS1; l1 < inter.Dir1LenFirst() - adjE1; l1++) {
        // get acceptor and inserter indices
        auto indA =
            GetSwapLoc(l1, l2, l3, array1.GhostLayers(), inter, d3, true);
        auto indI =
            GetSwapLoc(l1, l2, l3, array2.GhostLayers(), inter, d3, false);

        // assign cell data
        array1.InsertBlock(indA[0], indA[1], indA[2],
                           array2(indI[0], indI[1], indI[2]));
      }
    }
  }
}

template <typename T>
auto GrowInI(const T &orig) {
  T arr(orig.NumINoGhosts() + 1, orig.NumJNoGhosts(), orig.NumKNoGhosts(),
        orig.GhostLayers(), orig.BlockInfo());
  for (auto kk = arr.StartK(); kk < arr.EndK(); kk++) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); jj++) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ii++) {
        auto val = (ii == arr.EndI() - 1) ? orig(ii - 1, jj, kk)
                                          : orig(ii, jj, kk);
        arr.InsertBlock(ii, jj, kk, val);
      }
    }
  }
  return arr;
}

template <typename T>
auto GrowInJ(const T &orig) {
  T arr(orig.NumINoGhosts(), orig.NumJNoGhosts() + 1, orig.NumKNoGhosts(),
        orig.GhostLayers(), orig.BlockInfo());
  for (auto kk = arr.StartK(); kk < arr.EndK(); kk++) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); jj++) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ii++) {
        auto val = (jj == arr.EndJ() - 1) ? orig(ii, jj - 1, kk)
                                          : orig(ii, jj, kk);
        arr.InsertBlock(ii, jj, kk, val);
      }
    }
  }
  return arr;
}

template <typename T>
auto GrowInK(const T &orig) {
  T arr(orig.NumINoGhosts(), orig.NumJNoGhosts(), orig.NumKNoGhosts() + 1,
        orig.GhostLayers(), orig.BlockInfo());
  for (auto kk = arr.StartK(); kk < arr.EndK(); kk++) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); jj++) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ii++) {
        auto val = (kk == arr.EndK() - 1) ? orig(ii, jj, kk - 1)
                                          : orig(ii, jj, kk);
        arr.InsertBlock(ii, jj, kk, val);
      }
    }
  }
  return arr;
}


// ---------------------------------------------------------------------------
// member function definitions

template <typename T>
int multiArray3d<T>::Start(const string &dir) const {
  if (dir == "i") {
    return this->StartI();
  } else if (dir == "j") {
    return this->StartJ();
  } else if (dir == "k") {
    return this->StartK();
  } else {
    cerr << "ERROR: Error in multiArray3d::Start. Direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

template <typename T>
int multiArray3d<T>::End(const string &dir) const {
  if (dir == "i") {
    return this->EndI();
  } else if (dir == "j") {
    return this->EndJ();
  } else if (dir == "k") {
    return this->EndK();
  } else {
    cerr << "ERROR: Error in multiArray3d::End. Direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

template <typename T>
int multiArray3d<T>::PhysStart(const string &dir) const {
  if (dir == "i") {
    return this->PhysStartI();
  } else if (dir == "j") {
    return this->PhysStartJ();
  } else if (dir == "k") {
    return this->PhysStartK();
  } else {
    cerr << "ERROR: Error in multiArray3d::PhysStart. Direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

template <typename T>
int multiArray3d<T>::PhysEnd(const string &dir) const {
  if (dir == "i") {
    return this->PhysEndI();
  } else if (dir == "j") {
    return this->PhysEndJ();
  } else if (dir == "k") {
    return this->PhysEndK();
  } else {
    cerr << "ERROR: Error in multiArray3d::PhysEnd. Direction " << dir
         << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}


template <typename T>
T multiArray3d<T>::GetElem(const int &ii, const int &jj, const int &kk) const {
  if (ii < this->EndI() && jj < this->EndJ() && kk < this->EndK() &&
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

// Same type operator overloads
// operator overload for addition
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator+=(const multiArray3d<T> &arr) {
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] += arr.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator-=(const multiArray3d<T> &arr) {
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] -= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator*=(const multiArray3d<T> &arr) {
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] *= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator/=(const multiArray3d<T> &arr) {
  for (auto rr = 0; rr < this->Size(); rr++) {
    data_[rr] /= arr.data_[rr];
  }
  return *this;
}

template <typename T>
inline const multiArray3d<T> operator+(multiArray3d<T> lhs,
                                       const multiArray3d<T> &rhs) {
  return lhs += rhs;
}

template <typename T>
inline const multiArray3d<T> operator-(multiArray3d<T> lhs,
                                       const multiArray3d<T> &rhs) {
  return lhs -= rhs;
}

template <typename T>
inline const multiArray3d<T> operator*(multiArray3d<T> lhs,
                                       const multiArray3d<T> &rhs) {
  return lhs *= rhs;
}

template <typename T>
inline const multiArray3d<T> operator/(multiArray3d<T> lhs,
                                       const multiArray3d<T> &rhs) {
  return lhs /= rhs;
}

// ----------------------------------------------------------------
// operator overloads for type *this is holding

// operator overload for addition
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator+=(const T &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator-=(const T &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator*=(const T &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
template <typename T>
multiArray3d<T> & multiArray3d<T>::operator/=(const T &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

template <typename T>
inline const multiArray3d<T> operator+(const T &lhs, multiArray3d<T> rhs) {
  return rhs += lhs;
}

template <typename T>
inline const multiArray3d<T> operator-(const T &lhs, multiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs(rr) = lhs - rhs(rr);
  }
  return rhs;
}

template <typename T>
inline const multiArray3d<T> operator*(const T &lhs, multiArray3d<T> rhs) {
  return rhs *= lhs;
}

template <typename T>
inline const multiArray3d<T> operator/(const T &lhs, multiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs.data_(rr) = lhs / rhs(rr);
  }
  return rhs;
}

// ----------------------------------------------------------------
// operator overloads type that is not *this, or the type *this is holding

template <typename T, typename TT>
inline const multiArray3d<T> operator+(const TT &lhs, multiArray3d<T> rhs) {
  return rhs += lhs;
}

template <typename T, typename TT>
inline const multiArray3d<T> operator-(const TT &lhs, multiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs(rr) = lhs - rhs(rr);
  }
  return rhs;
}

template <typename T, typename TT>
inline const multiArray3d<T> operator*(const TT &lhs, multiArray3d<T> rhs) {
  return rhs *= lhs;
}

template <typename T, typename TT>
inline const multiArray3d<T> operator/(const TT &lhs, multiArray3d<T> rhs) {
  for (auto rr = 0; rr < rhs.Size(); rr++) {
    rhs(rr) = lhs / rhs(rr);
  }
  return rhs;
}

// member function to return a slice of the array
// this is the main slice function that all other overloaded slice functions call
template <typename T>
multiArray3d<T> multiArray3d<T>::Slice(const range &ir, const range &jr,
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
multiArray3d<T> multiArray3d<T>::Slice(const string &dir, const range &dirRange,
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
multiArray3d<T> multiArray3d<T>::Slice(const string &dir, int d2Ind, int d3Ind,
                                       const bool physOnly, const string id,
                                       const bool upper2,
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
multiArray3d<T> multiArray3d<T>::Slice(const string &dir, int dirInd,
                                       range dir1, range dir2, const string id,
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
void multiArray3d<T>::Insert(const range &ir, const range &jr, const range &kr,
                             const TT &arr) {
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
void multiArray3d<T>::Insert(const string &dir, const range &dirRange,
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
void multiArray3d<T>::Insert(const string &dir, int d2Ind, int d3Ind,
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
void multiArray3d<T>::Insert(const string &dir, int dirInd, range dir1,
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

// member function to fill this array with data from a provided one
// the size of the provided array must be identical to this array
// the data can now be accessed with the i, j, k indices from *this
// if this is not wanted SameSizeResize or SameSizeResizeGhosts can be used
template <typename T>
template <typename TT>
void multiArray3d<T>::Fill(const TT &arr) {
  // arr -- array to insert into this one

#ifndef NDEBUG
  // check that given array is same size
  if (this->Size() != arr.Size()) {
    cerr << "ERROR: Error in multiArray3d::Fill. Size of given array "
         << "does not match size of array to fill!" << endl;
    cerr << "Size of given array is " << arr.NumI() << ", " << arr.NumJ()
         << ", " << arr.NumK() << endl;
    cerr << "With " << arr.GhostLayers() << " ghost cell layers" << endl;
    cerr << "Resulting in a total size of " << arr.Size() << endl;
    cerr << "Size of array to fill is " << numI_ << ", " << numJ_ << ", "
         << numK_ << endl;
    cerr << "With " << numGhosts_ << " ghost cell layers" << endl;
    cerr << "Resulting in a total size of " << this->Size() << endl;
    exit(EXIT_FAILURE);
  }
#endif

  data_ = arr.data_;
}

template <typename T>
void multiArray3d<T>::SameSizeResize(const int &ii, const int&jj,
                                     const int &kk) {
  if (this->PhysicalSize() != (ii * jj * kk)) {
    cerr << "ERROR: Error in multiArray3d<T>::SameSizeResize. Attempting to "
         << "resize array of " << this->PhysicalSize() << " cells to " <<
        ii * jj * kk << " cells." << endl;
    exit(EXIT_FAILURE);
  }
  numI_ = ii + 2 * numGhosts_;
  numJ_ = jj + 2 * numGhosts_;
  numK_ = kk + 2 * numGhosts_;
}

template <typename T>
void multiArray3d<T>::SameSizeResizeGhosts(const int &ii, const int&jj,
                                           const int &kk, const int &ng) {
  if (this->Size() != ((ii + 2 * ng) * (jj + 2 * ng) * (kk + 2 * ng))) {
    cerr << "ERROR: Error in multiArray3d<T>::SameSizeResizeGhosts. Attempting "
         << "to resize array of " << this->Size() << " cells to " <<
        (ii + 2 * ng) * (jj + 2 * ng) * (kk + 2 * ng) << " cells." << endl;
    exit(EXIT_FAILURE);
  }
  numI_ = ii;
  numJ_ = jj;
  numK_ = kk;
  numGhosts_ = ng;
}

// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const multiArray3d<T> &arr) {
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

template <typename T>
multiArray3d<T> multiArray3d<T>::GrowI() const {
  return GrowInI(*this);
}
template <typename T>
multiArray3d<T> multiArray3d<T>::GrowJ() const {
  return GrowInJ(*this);
}
template <typename T>
multiArray3d<T> multiArray3d<T>::GrowK() const {
  return GrowInK(*this);
}

template <typename T>
template <typename TT>
void multiArray3d<T>::PutSlice(const TT &array, const connection &inter,
                               const int &d3) {
  // array -- array to insert into *this
  // inter -- connection data structure defining patches and orientation
  // d3 -- distance of direction normal to patch to insert
  InsertSlice((*this), array, inter, d3);
}

/*Member function to pack an array into a buffer, swap it with its
  connection partner, and then unpack it into an array.*/
template <typename T>
void multiArray3d<T>::PackSwapUnpackMPI(const connection &inter,
                                        const MPI_Datatype &MPI_arrData,
                                        const int &rank, const int tag) {
  // inter -- connection boundary for the swap
  // MPI_arrData -- MPI datatype to pass data type in array
  // rank -- processor rank
  // tag -- id to send data with (default 1)

  // swap with mpi_send_recv_replace
  // pack data into buffer, but first get size
  auto bufSize = 0;
  auto tempSize = 0;
  // add size for states
  MPI_Pack_size(this->Size(), MPI_arrData, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  // add size for 5 ints for multiArray3d dims and num ghosts
  MPI_Pack_size(5, MPI_INT, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;

  // allocate buffer to pack data into
  // use unique_ptr to manage memory; use underlying pointer for MPI calls
  auto unqBuffer = unique_ptr<char>(new char[bufSize]);
  auto *buffer = unqBuffer.get();  

  // pack data into buffer
  auto numI = this->NumI();
  auto numJ = this->NumJ();
  auto numK = this->NumK();
  auto numGhosts = this->GhostLayers();
  auto blkSize = this->BlockSize();
  auto position = 0;
  MPI_Pack(&numI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numGhosts, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&blkSize, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(data_)), this->Size(), MPI_arrData, buffer, bufSize,
           &position, MPI_COMM_WORLD);

  MPI_Status status;
  if (rank == inter.RankFirst()) {  // send/recv with second entry in connection
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankSecond(), tag,
                         inter.RankSecond(), tag, MPI_COMM_WORLD, &status);
  } else {  // send/recv with first entry in connection
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankFirst(), tag,
                         inter.RankFirst(), tag, MPI_COMM_WORLD, &status);
  }

  // put slice back into multiArray3d
  position = 0;
  MPI_Unpack(buffer, bufSize, &position, &numI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numGhosts, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &blkSize, 1, MPI_INT, MPI_COMM_WORLD);

  // resize slice
  this->SameSizeResize(numI, numJ, numK);

  MPI_Unpack(buffer, bufSize, &position, &(*std::begin(data_)),
             this->Size(), MPI_arrData, MPI_COMM_WORLD);
}

/* Function to swap slice using MPI. This is similar to the SwapSlice
   function, but is called when the neighboring procBlocks are on different
   processors.
*/
template <typename T>
void multiArray3d<T>::SwapSliceMPI(const connection &conn, const int &rank,
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
void multiArray3d<T>::SwapSlice(const connection &conn, TT &array) {
  // conn -- connection boundary information
  // array -- second array involved in connection boundary
  SwapSliceLocal((*this), conn, array);
}


#endif

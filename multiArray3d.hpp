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

#ifndef MULTIARRAY3DHEADERDEF  // only if the macro MULTIARRAY3DHEADERDEF
                               // is not defined execute these lines of code

#define MULTIARRAY3DHEADERDEF  // define the macro

/* This file contains the header and implementation for a multidimensional (3D)
   array class. The class is to act as a container to store all data types,
   and provide easy access to elements with i, j, k indexing.
 */

#include <iostream>  // ostream
#include <vector>    // vector
#include "mpi.h"
#include "vector3d.hpp"
#include "boundaryConditions.hpp"  // interblock

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;
using std::vector;

template <typename T>
class multiArray3d {
  vector<T> data_;
  int numI_;
  int numJ_;
  int numK_;

  // private member functions
  int GetLoc1D(const int &ii, const int &jj, const int &kk) const {
    return ii + jj * numI_ + kk * numI_ * numJ_;
  }

 public:
  // constructor
  multiArray3d(const int &ii, const int &jj, const int &kk, const T &init) :
      data_(vector<T>(ii * jj * kk, init)),
      numI_(ii), numJ_(jj), numK_(kk) {}
  multiArray3d(const int &ii, const int &jj, const int &kk) :
      data_(vector<T>(ii * jj * kk)),
      numI_(ii), numJ_(jj), numK_(kk) {}
  multiArray3d() : multiArray3d(1, 1, 1) {}

  // move constructor and assignment operator
  multiArray3d(multiArray3d&&) noexcept = default;
  multiArray3d& operator=(multiArray3d&&) noexcept = default;

  // copy constructor and assignment operator
  multiArray3d(const multiArray3d&) = default;
  multiArray3d& operator=(const multiArray3d&) = default;

  // member functions
  int Size() const { return data_.size(); }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }

  multiArray3d<T> Slice(const int&, const int&, const int&, const int&,
                        const int&, const int&) const;
  void Insert(const int&, const int&, const int&, const int&, const int&,
              const int&, const multiArray3d<T>&);
  void PutSlice(const multiArray3d<T> &, const interblock &, const int &,
                const int &);
  void SwapSliceMPI(const interblock&, const int&, const MPI_Datatype&,
                    const int&);
  void SwapSlice(const interblock&, multiArray3d<T>&, const int&,
                      const int&);

  void Zero(const T&);

  void GrowI();
  void GrowJ();
  void GrowK();

  void PackSwapUnpackMPI(const interblock &, const MPI_Datatype &, const int &);

  T GetElem(const int &ii, const int &jj, const int &kk) const;

  // operator overloads
  T& operator()(const int &ii, const int &jj, const int &kk) {
    return data_[this->GetLoc1D(ii, jj, kk)];
  }
  const T& operator()(const int &ii, const int &jj, const int &kk) const {
    return data_[this->GetLoc1D(ii, jj, kk)];
  }
  T& operator()(const int &ind) {
    return data_[ind];
  }
  const T& operator()(const int &ind) const {
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
  inline multiArray3d<T> & operator+=(const TT &);
  template <typename TT>
  inline multiArray3d<T> & operator-=(const TT &);
  template <typename TT>
  inline multiArray3d<T> & operator*=(const TT &);
  template <typename TT>
  inline multiArray3d<T> & operator/=(const TT &);

  template <typename TT>
  inline multiArray3d<T> operator+(const TT &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  template <typename TT>
  inline multiArray3d<T> operator-(const TT &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  template <typename TT>
  inline multiArray3d<T> operator*(const TT &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  template <typename TT>
  inline multiArray3d<T> operator/(const TT &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  void ClearResize(const int &ii, const int &jj, const int &kk) {
    *this = multiArray3d<T>(ii, jj, kk);
  }
  void ClearResize(const int &ii, const int &jj, const int &kk, const T &val) {
    *this = multiArray3d<T>(ii, jj, kk, val);
  }

  void SameSizeResize(const int &ii, const int&jj, const int &kk);

  // destructor
  ~multiArray3d() noexcept {}
};

template <typename T>
T multiArray3d<T>::GetElem(const int &ii, const int &jj, const int &kk) const {
  if (ii < numI_ && jj < numJ_ && kk < numK_ && ii >= 0 && jj >= 0 && kk >= 0) {
    return data_[this->GetLoc1D(ii, jj, kk)];
  } else {
    cerr << "ERROR: Tried to access location outside of bounds of "
         << "multiArray3d" << endl;
    cerr << "Tried to access " << ii << ", " << jj << ", " << kk << endl;
    cerr << "Maximum locations are " << numI_ - 1 << ", " << numJ_ - 1 << ", "
         << numK_ - 1 << endl;
    exit(0);
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

// operator overload for addition
template <typename T> template <typename TT>
multiArray3d<T> & multiArray3d<T>::operator+=(const TT &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
template <typename T> template <typename TT>
multiArray3d<T> & multiArray3d<T>::operator-=(const TT &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
template <typename T> template <typename TT>
multiArray3d<T> & multiArray3d<T>::operator*=(const TT &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
template <typename T> template <typename TT>
multiArray3d<T> & multiArray3d<T>::operator/=(const TT &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

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
template <typename T>
multiArray3d<T> multiArray3d<T>::Slice(const int &is, const int &ie,
                                       const int &js, const int &je,
                                       const int &ks, const int &ke) const {
  // is -- starting i-index to take slice (inclusive)
  // ie -- ending i-index to take slice (inclusive)
  // js -- starting j-index to take slice (inclusive)
  // je -- ending j-index to take slice (inclusive)
  // ks -- starting k-index to take slice (inclusive)
  // ke -- ending k-index to take slice (inclusive)

  // check that slice bounds are within parent array and that end bounds are
  // greater than or equal to start bounds
  if (!(ie <= numI_ && ie >= 0 && is <= numI_ && is >= 0 &&
        je <= numJ_ && je >= 0 && js <= numJ_ && js >= 0 &&
        ke <= numK_ && ke >= 0 && ks <= numK_ && ks >= 0 &&
        ie >= is && je >= js && ke >= ks)) {
    cerr << "ERROR: Error in multiArray3d::Slice. Cannot take slice with "
         << "boundaries " << is << ", " << ie << ", " << js << ", " << je
         << ", " << ks << ", " << ke << " from array with dimensions " <<
        numI_ << ", " << numJ_ << ", " << numK_ << endl;
    exit(0);
  }

  multiArray3d<T> array(ie - is + 1, je - js + 1, ke - ks + 1);

// in struct s is for index of sliced array, p is for index of parent array
  for (struct {int s; int p;} kk = {0, ks}; kk.s < array.numK_;
       kk.s++, kk.p++) {
    for (struct {int s; int p;} jj = {0, js}; jj.s < array.numJ_;
         jj.s++, jj.p++) {
      for (struct {int s; int p;} ii = {0, is}; ii.s < array.numI_;
           ii.s++, ii.p++) {
        array(ii.s, jj.s, kk.s) = (*this)(ii.p, jj.p, kk.p);
      }
    }
  }
  return array;
}

// member function to insert an array into this one
template <typename T>
void multiArray3d<T>::Insert(const int &is, const int &ie, const int &js,
                             const int &je, const int &ks, const int &ke,
                             const multiArray3d<T> &array) {
  // is -- starting i-index to insert slice (inclusive)
  // ie -- ending i-index to insert slice (inclusive)
  // js -- starting j-index to insert slice (inclusive)
  // je -- ending j-index to insert slice (inclusive)
  // ks -- starting k-index to insert slice (inclusive)
  // ke -- ending k-index to insert slice (inclusive)
  // array -- array to insert into this one

  // check that given array fits in given bounds and that end bounds are
  // greater than or equal to start bounds
  if (ie - is + 1 != array.numI_ && je - js + 1 != array.numJ_ &&
      ke - ks +1 != array.numK_ && ie < is && je < js && ke < ks) {
    cerr << "ERROR: Error in multiArray3d::Insert. Bounds of given array " <<
        "does not match size of location given to insert array into!" << endl;
    cerr << "Size of given array is " << array.numI_ << ", " << array.numJ_ <<
        ", " << array.numK_ << endl;
    cerr << "Size of location is " << ie - is + 1 << ", " << je - js + 1 <<
        ", " << ke - ks + 1 << endl;
    exit(0);
  }

  // in struct s is for index of sliced array, p is for index of parent array
  for (struct {int s; int p;} kk = {0, ks}; kk.s < array.numK_;
       kk.s++, kk.p++) {
    for (struct {int s; int p;} jj = {0, js}; jj.s < array.numJ_;
         jj.s++, jj.p++) {
      for (struct {int s; int p;} ii = {0, is}; ii.s < array.numI_;
           ii.s++, ii.p++) {
        (*this)(ii.p, jj.p, kk.p) = array(ii.s, jj.s, kk.s);
      }
    }
  }
}

template <typename T>
void multiArray3d<T>::GrowI() {
  multiArray3d<T> arr(numI_ + 1, numJ_, numK_);
  for (auto kk = 0; kk < arr.numK_; kk++) {
    for (auto jj = 0; jj < arr.numJ_; jj++) {
      for (auto ii = 0; ii < arr.numI_; ii++) {
        arr(ii, jj, kk) = (ii == arr.numI_ - 1) ? (*this)(ii - 1, jj, kk) :
            (*this)(ii, jj, kk);
      }
    }
  }
  *this = arr;
}

template <typename T>
void multiArray3d<T>::GrowJ() {
  multiArray3d<T> arr(numI_, numJ_ + 1, numK_);
  for (auto kk = 0; kk < arr.numK_; kk++) {
    for (auto jj = 0; jj < arr.numJ_; jj++) {
      for (auto ii = 0; ii < arr.numI_; ii++) {
        arr(ii, jj, kk) = (jj == arr.numJ_ - 1) ? (*this)(ii, jj - 1, kk) :
            (*this)(ii, jj, kk);
      }
    }
  }
  *this = arr;
}

template <typename T>
void multiArray3d<T>::GrowK() {
  multiArray3d<T> arr(numI_, numJ_, numK_ + 1);
  for (auto kk = 0; kk < arr.numK_; kk++) {
    for (auto jj = 0; jj < arr.numJ_; jj++) {
      for (auto ii = 0; ii < arr.numI_; ii++) {
        arr(ii, jj, kk) = (kk == arr.numK_ - 1) ? (*this)(ii, jj, kk - 1) :
            (*this)(ii, jj, kk);
      }
    }
  }
  *this = arr;
}

template <typename T>
void multiArray3d<T>::SameSizeResize(const int &ii, const int&jj,
                                     const int &kk) {
  if (this->Size() != (ii * jj * kk)) {
    cerr << "ERROR: Error in multiArray3d<T>::SameSizeResize. Attempting to "
         << "resize array of " << this->Size() << " cells to " <<
        ii * jj * kk << " cells." << endl;
    exit(0);
  }
    numI_ = ii;
    numJ_ = jj;
    numK_ = kk;
}

// member function to "zero out" the container with a supplied "zero"
template <typename T>
void multiArray3d<T>::Zero(const T &zero) {
  for (auto &val : data_) {
    val = zero;
  }
}

// operation overload for << - allows use of cout, cerr, etc.
template <typename T>
ostream &operator<<(ostream &os, const multiArray3d<T> &array) {
  os << "Size: " << array.NumI() << ", " << array.NumJ() << ", "
     << array.NumK() << endl;

  for (auto kk = 0; kk < array.NumK(); kk++) {
    for (auto jj = 0; jj < array.NumJ(); jj++) {
      for (auto ii = 0; ii < array.NumI(); ii++) {
        os << ii << ", " << jj << ", " << kk << ", " <<
            array(ii, jj, kk) << endl;
      }
    }
  }
  return os;
}

template <typename T>
void multiArray3d<T>::PutSlice(const multiArray3d<T> &array,
                               const interblock &inter, const int &d3,
                               const int &numG) {
  // array -- array to insert into *this
  // inter -- interblock data structure defining patches and orientation
  // d3 -- distance of directioni normal to patch to insert
  // numG -- number of ghost cells

  // check that number of cells to insert matches
  auto blkCell = (inter.Dir1EndFirst() - inter.Dir1StartFirst()) *
      (inter.Dir2EndFirst() - inter.Dir2StartFirst()) * d3;
  if (blkCell != array.Size()) {
    cerr << "ERROR: Error in multiArray3d<T>::PutSlice(). Number of cells "
            "being inserted does not match designated space to insert." << endl;
    cerr << "Direction 1, 2, 3 of multiArray3d<T> to insert into: "
         << inter.Dir1EndFirst() - inter.Dir1StartFirst() << ", "
         << inter.Dir2EndFirst() - inter.Dir2StartFirst() << ", " << d3 << endl;
    cerr << "Direction I, J, K of multiArray3d<T> to insert: " << array.NumI()
         << ", " << array.NumJ() << ", " << array.NumK() << endl;
    exit(0);
  }

  // adjust insertion indices if patch borders another interblock on the same
  // surface of the block
  auto adjS1 = (inter.Dir1StartInterBorderFirst()) ? numG : 0;
  auto adjE1 = (inter.Dir1EndInterBorderFirst()) ? numG : 0;
  auto adjS2 = (inter.Dir2StartInterBorderFirst()) ? numG : 0;
  auto adjE2 = (inter.Dir2EndInterBorderFirst()) ? numG : 0;

  // loop over cells to insert
  for (auto l3 = 0; l3 < d3; l3++) {
    for (auto l2 = adjS2;
         l2 < (inter.Dir2EndFirst() - inter.Dir2StartFirst() - adjE2); l2++) {
      for (auto l1 = adjS1;
           l1 < (inter.Dir1EndFirst() - inter.Dir1StartFirst() - adjE1); l1++) {
        // get acceptor and inserter indices
        auto indA = GetSwapLoc(l1, l2, l3, inter, true);
        auto indI = GetSwapLoc(l1, l2, l3, inter, false);

        // swap cell data
        (*this)(indA[0], indA[1], indA[2]) = array(indI[0], indI[1], indI[2]);
      }
    }
  }
}

/*Member function to pack an array into a buffer, swap it with its
 * interblock partner, and then unpack it into an array.*/
template <typename T>
void multiArray3d<T>::PackSwapUnpackMPI(const interblock &inter,
                                        const MPI_Datatype &MPI_cellData,
                                        const int &rank) {
  // inter -- interblock boundary for the swap
  // MPI_cellData -- MPI datatype to pass data type in array
  // rank -- processor rank

  // swap with mpi_send_recv_replace
  // pack data into buffer, but first get size
  auto bufSize = 0;
  auto tempSize = 0;
  MPI_Pack_size(this->Size(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for states
  bufSize += tempSize;
  // add size for 3 ints for multiArray3d dims
  MPI_Pack_size(3, MPI_INT, MPI_COMM_WORLD,
                &tempSize);
  bufSize += tempSize;

  auto *buffer = new char[bufSize];  // allocate buffer to pack data into

  // pack data into buffer
  auto numI = this->NumI();
  auto numJ = this->NumJ();
  auto numK = this->NumK();
  auto position = 0;
  MPI_Pack(&numI, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numJ, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numK, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&data_[0], this->Size(), MPI_cellData, buffer,
           bufSize, &position, MPI_COMM_WORLD);

  MPI_Status status;
  if (rank == inter.RankFirst()) {  // send/recv with second entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankSecond(), 1,
                         inter.RankSecond(), 1, MPI_COMM_WORLD, &status);
  } else {  // send/recv with first entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankFirst(), 1,
                         inter.RankFirst(), 1, MPI_COMM_WORLD, &status);
  }

  // put slice back into multiArray3d
  position = 0;
  MPI_Unpack(buffer, bufSize, &position, &numI, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numJ, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numK, 1, MPI_INT,
             MPI_COMM_WORLD);
  // resize slice
  this->SameSizeResize(numI, numJ, numK);

  MPI_Unpack(buffer, bufSize, &position, &data_[0],
             this->Size(), MPI_cellData, MPI_COMM_WORLD);

  delete[] buffer;
}

/* Function to swap slice using MPI. This is similar to the SwapSlice
function, but is called when the neighboring procBlocks are on different
processors.
*/
template <typename T>
void multiArray3d<T>::SwapSliceMPI(const interblock &inter, const int &rank,
                             const MPI_Datatype &MPI_cellData,
                             const int &numGhosts) {
  // inter -- interblock boundary information
  // rank -- processor rank
  // MPI_cellData -- MPI datatype for passing primVars, genArray
  // numGhosts -- number of ghost cells

  // Get indices for slice coming from block to swap
  auto is = 0;
  auto ie = 0;
  auto js = 0;
  auto je = 0;
  auto ks = 0;
  auto ke = 0;

  if (rank == inter.RankFirst()) {  // local block is first in interblock
    inter.FirstSliceIndices(is, ie, js, je, ks, ke, numGhosts);
  // local block is second in interblock
  } else if (rank == inter.RankSecond()) {
    inter.SecondSliceIndices(is, ie, js, je, ks, ke, numGhosts);
  } else {
    cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Processor rank does "
            "not match either of interblock ranks!" << endl;
    exit(0);
  }

  // get local state slice to swap
  auto slice = this->Slice(is, ie, js, je, ks, ke);

  // swap state slices with partner block
  slice.PackSwapUnpackMPI(inter, MPI_cellData, rank);

  // change interblocks to work with slice and ghosts
  auto interAdj = inter;

  // block to insert into is first in interblock
  if (rank == inter.RankFirst()) {
    interAdj.AdjustForSlice(true, numGhosts);
  } else {  // block to insert into is second in interblock, so pass swapped
            // version
    interAdj.AdjustForSlice(false, numGhosts);
  }

  // insert state slice into procBlock
  this->PutSlice(slice, interAdj, numGhosts, numGhosts);
}

/* Function to swap ghost cells between two blocks at an interblock
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
logic ensures that the ghost cells at the interblock boundary exactly match
their partner block as if there were no separation in the grid.
*/
template <typename T>
void multiArray3d<T>::SwapSlice(const interblock &inter,
                                multiArray3d<T> &array, const int &numG1,
                                const int &numG2) {
  // inter -- interblock boundary information
  // array -- second array involved in interblock boundary
  // numG1 -- number of ghost cells in first array
  // numG2 -- number of ghost cells in second array

  // Get indices for slice coming from first block to swap
  auto is1 = 0;
  auto ie1 = 0;
  auto js1 = 0;
  auto je1 = 0;
  auto ks1 = 0;
  auto ke1 = 0;

  inter.FirstSliceIndices(is1, ie1, js1, je1, ks1, ke1, numG1);

  // Get indices for slice coming from second block to swap
  auto is2 = 0;
  auto ie2 = 0;
  auto js2 = 0;
  auto je2 = 0;
  auto ks2 = 0;
  auto ke2 = 0;

  inter.SecondSliceIndices(is2, ie2, js2, je2, ks2, ke2, numG2);

  // get slices to swap
  auto slice1 = this->Slice(is1, ie1, js1, je1, ks1, ke1);
  auto slice2 = array.Slice(is2, ie2, js2, je2, ks2, ke2);

  // change interblocks to work with slice and ghosts
  interblock inter1 = inter;
  interblock inter2 = inter;
  inter1.AdjustForSlice(false, numG1);
  inter2.AdjustForSlice(true, numG2);

  // put slices in proper blocks
  this->PutSlice(slice2, inter2, numG2, numG2);
  array.PutSlice(slice1, inter1, numG1, numG1);
}


#endif

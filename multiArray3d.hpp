/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
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
#include "vector3d.hpp"

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
  multiArray3d() : data_(vector<T>(1)), numI_(1), numJ_(1), numK_(1) {}
  multiArray3d(const int&, const int&, const int&);
  multiArray3d(const int&, const int&, const int&, const T&);

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

  void GrowI();
  void GrowJ();
  void GrowK();

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

  T GetElem(const int &ii, const int &jj, const int &kk) const;

  multiArray3d<T> operator*(const T&) const;
  multiArray3d<T> operator/(const T&) const;
  multiArray3d<T> operator+(const T&) const;
  multiArray3d<T> operator-(const T&) const;

  template <typename TT>
  multiArray3d<T> operator*(const TT&) const;
  template <typename TT>
  multiArray3d<T> operator/(const TT&) const;
  template <typename TT>
  multiArray3d<T> operator+(const TT&) const;
  template <typename TT>
  multiArray3d<T> operator-(const TT&) const;

  template <typename TT>
  friend multiArray3d<TT> operator*(const TT&, const multiArray3d<TT>&);
  template <typename TT>
  friend multiArray3d<TT> operator/(const TT&, const multiArray3d<TT>&);
  template <typename TT>
  friend multiArray3d<TT> operator+(const TT&, const multiArray3d<TT>&);
  template <typename TT>
  friend multiArray3d<TT> operator-(const TT&, const multiArray3d<TT>&);

  template <typename T1, typename T2>
  friend multiArray3d<T1> operator*(const T2&, const multiArray3d<T1>&);
  template <typename T1, typename T2>
  friend multiArray3d<T1> operator/(const T2&, const multiArray3d<T1>&);
  template <typename T1, typename T2>
  friend multiArray3d<T1> operator+(const T2&, const multiArray3d<T1>&);
  template <typename T1, typename T2>
  friend multiArray3d<T1> operator-(const T2&, const multiArray3d<T1>&);

  multiArray3d<T> operator*(const multiArray3d<T> &) const;
  multiArray3d<T> operator/(const multiArray3d<T> &) const;
  multiArray3d<T> operator+(const multiArray3d<T> &) const;
  multiArray3d<T> operator-(const multiArray3d<T> &) const;

  void ClearResize(const int &ii, const int &jj, const int &kk) {
    *this = multiArray3d<T>(ii, jj, kk);
  }
  void ClearResize(const int &ii, const int &jj, const int &kk, const T &val) {
    *this = multiArray3d<T>(ii, jj, kk, val);
  }

  void SameSizeResize(const int &ii, const int&jj, const int &kk);

  template <typename TT>
  friend ostream &operator<<(ostream &, const multiArray3d<TT> &);

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

// constructor - default initialization
template <typename T>
multiArray3d<T>::multiArray3d(const int &ii, const int &jj, const int &kk) {
  // ii -- i-dimension
  // jj -- j-dimension
  // kk -- k-dimension

  numI_ = ii;
  numJ_ = jj;
  numK_ = kk;
  data_ = vector<T>(ii * jj * kk);
}

// constructor - uniform initialization
template <typename T>
multiArray3d<T>::multiArray3d(const int &ii, const int &jj, const int &kk,
                              const T &init) {
  // ii -- i-dimension
  // jj -- j-dimension
  // kk -- k-dimension
  // T -- data member to initialize entire class with

  numI_ = ii;
  numJ_ = jj;
  numK_ = kk;
  data_ = vector<T>(ii * jj * kk, init);
}

// member functions
// operator overload to multiply all data components by same factor
template <typename T>
multiArray3d<T> multiArray3d<T>::operator*(const T &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val *= factor;
  }
  return arr;
}

template <typename TT>
multiArray3d<TT> operator*(const TT &factor, const multiArray3d<TT> &arr) {
  return arr * factor;
}

template <typename T>
template <typename TT>
multiArray3d<T> multiArray3d<T>::operator*(const TT &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val = val * factor;
  }
  return arr;
}

template <typename T1, typename T2>
multiArray3d<T1> operator*(const T2 &factor, const multiArray3d<T1> &arr) {
  return arr * factor;
}


// operator overload to divide the same factor by all data components
template <typename T>
multiArray3d<T> multiArray3d<T>::operator/(const T &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val /= factor;
  }
  return arr;
}

template <typename TT>
multiArray3d<TT> operator/(const TT &factor, const multiArray3d<TT> &arr) {
  multiArray3d<TT> out = arr;
  for (TT& val : out.data_) {
    val = factor / val;
  }
  return out;
}

template <typename T>
template <typename TT>
multiArray3d<T> multiArray3d<T>::operator/(const TT &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val = val / factor;
  }
  return arr;
}

template <typename T1, typename T2>
multiArray3d<T1> operator/(const T2 &factor, const multiArray3d<T1> &arr) {
  multiArray3d<T1> out = arr;
  for (T1& val : out.data_) {
    val = factor / val;
  }
  return out;
}

// operator overload to add the same factor to all data components
template <typename T>
multiArray3d<T> multiArray3d<T>::operator+(const T &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val += factor;
  }
  return arr;
}

template <typename TT>
multiArray3d<TT> operator+(const TT &factor, const multiArray3d<TT> &arr) {
  return arr + factor;
}

template <typename T>
template <typename TT>
multiArray3d<T> multiArray3d<T>::operator+(const TT &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val = val + factor;
  }
  return arr;
}

template <typename T1, typename T2>
multiArray3d<T1> operator+(const T2 &factor, const multiArray3d<T1> &arr) {
  return arr + factor;
}

// operator overload to subtract the same factor from all data components
template <typename T>
multiArray3d<T> multiArray3d<T>::operator-(const T &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val -= factor;
  }
  return arr;
}

template <typename TT>
multiArray3d<TT> operator-(const TT &factor, const multiArray3d<TT> &arr) {
  multiArray3d<TT> out = arr;
  for (TT& val : out.data_) {
    val = factor - val;
  }
  return out;
}

template <typename T>
template <typename TT>
multiArray3d<T> multiArray3d<T>::operator-(const TT &factor) const {
  multiArray3d<T> arr = *this;
  for (T& val : arr.data_) {
    val = val - factor;
  }
  return arr;
}

template <typename T1, typename T2>
multiArray3d<T1> operator-(const T2 &factor, const multiArray3d<T1> &arr) {
  multiArray3d<T1> out = arr;
  for (T1& val : out.data_) {
    val = factor - val;
  }
  return out;
}

// operator overload for pointwise multiplication
template <typename T>
multiArray3d<T> multiArray3d<T>::operator*(const multiArray3d<T> &b) const {
  multiArray3d<T> arr = *this;
  for (unsigned int ii = 0; ii < data_.size(); ii++) {
    arr.data_[ii] = arr.data_[ii] * b.data_[ii];
  }
  return arr;
}

// operator overload to divide all data components by same factor
template <typename T>
multiArray3d<T> multiArray3d<T>::operator/(const multiArray3d<T> &b) const {
  multiArray3d<T> arr = *this;
  for (unsigned int ii = 0; ii < data_.size(); ii++) {
    arr.data_[ii] = arr.data_[ii] / b.data_[ii];
  }
  return arr;
}

// operator overload to add the same factor to all data components
template <typename T>
multiArray3d<T> multiArray3d<T>::operator+(const multiArray3d<T> &b) const {
  multiArray3d<T> arr = *this;
  for (unsigned int ii = 0; ii < data_.size(); ii++) {
    arr.data_[ii] = arr.data_[ii] + b.data_[ii];
  }
  return arr;
}

// operator overload to subtract the same factor from all data components
template <typename T>
multiArray3d<T> multiArray3d<T>::operator-(const multiArray3d<T> &b) const {
  multiArray3d<T> arr = *this;
  for (unsigned int ii = 0; ii < data_.size(); ii++) {
    arr.data_[ii] = arr.data_[ii] - b.data_[ii];
  }
  return arr;
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

// overload to print array contents to stream
template <typename TT>
ostream &operator<<(ostream &os, const multiArray3d<TT> &array) {
  os << "Size: " << array.numI_ << ", " << array.numJ_ << ", "
     << array.numK_ << endl;

  for (int kk = 0; kk < array.numK_; kk++) {
    for (int jj = 0; jj < array.numJ_; jj++) {
      for (int ii = 0; ii < array.numI_; ii++) {
        os << ii << ", " << jj << ", " << kk << ", " <<
            array(ii, jj, kk) << endl;
      }
    }
  }
  return os;
}

template <typename T>
void multiArray3d<T>::GrowI() {
  multiArray3d<T> arr(numI_ + 1, numJ_, numK_);
  for (int kk = 0; kk < arr.numK_; kk++) {
    for (int jj = 0; jj < arr.numJ_; jj++) {
      for (int ii = 0; ii < arr.numI_; ii++) {
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
  for (int kk = 0; kk < arr.numK_; kk++) {
    for (int jj = 0; jj < arr.numJ_; jj++) {
      for (int ii = 0; ii < arr.numI_; ii++) {
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
  for (int kk = 0; kk < arr.numK_; kk++) {
    for (int jj = 0; jj < arr.numJ_; jj++) {
      for (int ii = 0; ii < arr.numI_; ii++) {
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


#endif

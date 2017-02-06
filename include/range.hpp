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

#ifndef RANGEHEADERDEF  // only if the macro RANGEHEADERDEF
                        // is not defined execute these lines of code

#define RANGEHEADERDEF  // define the macro

/* This file contains the header and implementation for a range class that is
   used to index multiArray3d
 */

#include <iostream>  // ostream

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;

class range {
  int start_;
  int end_;

 public:
  // constructors
  range(const int &s, const int &e) : start_(s), end_(e) {}
  range(const int &ind) : range(ind, ind + 1) {}
  range(std::initializer_list<int> values) {
    if (values.size() > 2 || values.size() == 0) {
      cerr << "ERROR: Error in range::range(). Initializer list is invalid size!"
           << endl;
    } else if (values.size() == 2) {
      start_ = *values.begin();
      end_ = *(values.begin() + 1);
    } else {
      start_ = *values.begin();
      end_ = start_ + 1;
    }
  }

  // member functions
  int Start() const {return start_;}
  int End() const {return end_;}
  int Size() const {return end_ - start_;}
  bool IsInside(const range &r) const {
    return start_ >= r.start_ && start_ < r.end_ && end_ > r.start_ &&
        end_ <= r.end_;
  }
  bool IsValid() const {
    return end_ > start_;
  }
  void GrowStart() {start_++;}
  void GrowEnd() {end_++;}
  void ShrinkStart() {start_--;}
  void ShrinkEnd() {end_--;}

  // destructor
  ~range() noexcept {}
};

ostream &operator<<(ostream &os, const range &r);

#endif

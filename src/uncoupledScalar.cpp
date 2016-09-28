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

#include "uncoupledScalar.hpp"
#include "genArray.hpp"      // genArray

using std::cout;
using std::endl;
using std::cerr;

// member functions
// member function to multiply the scalars with a genArray
genArray uncoupledScalar::ArrayMult(genArray arr) const {
  arr[0] *= flowVar_;
  arr[1] *= flowVar_;
  arr[2] *= flowVar_;
  arr[3] *= flowVar_;
  arr[4] *= flowVar_;

  arr[5] *= turbVar_;
  arr[6] *= turbVar_;

  return arr;
}

// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, uncoupledScalar &scalar) {
  os << scalar.FlowVariable() << endl;
  os << scalar.TurbVariable() << endl;
  return os;
}

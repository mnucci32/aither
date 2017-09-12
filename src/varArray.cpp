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

#include <cstdlib>  // exit()
#include <iostream>  // cout
#include <cmath>
#include "varArray.hpp"

using std::cout;
using std::endl;
using std::cerr;

// ------------------------------------------------------------------
// functions for varArray class

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const varArray &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const primative &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const conserved &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const residual &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const source &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const flux &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}


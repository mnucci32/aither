/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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

#include <iostream>
#include "inviscidFlux.hpp"

using std::endl;
using std::ostream;

/* Member function to calculate the Roe flux, given the left and right
 * convective fluxes as well as the dissipation term.
*/
void inviscidFlux::RoeFlux(const inviscidFlux &right, const varArray &diss) {
  // right -- right convective flux
  // diss -- dissipation term

  (*this) += right - diss;
  (*this) *= 0.5;
}

// non-member functions
// -----------------------------------------------------------

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const inviscidFlux &flux) {
  for (auto rr = 0; rr < flux.Size(); rr++) {
    os << flux[rr] << endl;
  }
  return os;
}


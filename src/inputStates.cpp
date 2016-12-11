/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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

#include <iostream>     // cout
#include "inputStates.hpp"

using std::cout;
using std::endl;
using std::cerr;

ostream &operator<<(ostream &os, const icState &state) {
  os << "icState(tag=" << state.Tag() << "; pressure=" << state.Pressure()
     << "; density=" << state.Density() << "; velocity=[" << state.Velocity()
     << "]; turbulenceIntensity=" << state.TurbulenceIntensity()
     << "; eddyViscosityRatio=" << state.EddyViscosityRatio() << ")";
  return os;
}

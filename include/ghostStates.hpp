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

#ifndef GHOSTHEADERDEF
#define GHOSTHEADERDEF

#include <string>
#include "arrayView.hpp"
#include "tensor.hpp"
#include "primitive.hpp"

using std::string;

// forward class declarations
struct wallVars;
class input;

primitive GetGhostState(const primitiveView &interior, const string &bcType,
                        const vector3d<double> &areaVec, const double &wallDist,
                        const int &surf, const input &inputVars, const int &tag,
                        const physics &phys, wallVars &wVars, const int &layer,
                        const double &mu = 0.0,
                        const double &dt = 0.0, const primitive &stateN = {},
                        const vector3d<double> &pressGrad = {},
                        const tensor<double> &velGrad = {},
                        const double &avgMach = 0.0,
                        const double &maxMach = 0.0);

primitive ExtrapolateHoldMixture(const primitive &boundary,
                                 const double &factor,
                                 const primitiveView &interior);

#endif
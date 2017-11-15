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
#include <memory>
#include "arrayView.hpp"

using std::string;
using std::unique_ptr;

// forward class declarations
class turbModel;
struct wallVars;
class input;
class primitive;

primitive GetGhostState(
    const primitiveView &interior, const string &bcType, const vector3d<double> &areaVec,
    const double &wallDist, const int &surf, const input &inputVars,
    const int &tag, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, wallVars &wVars, const int &layer,
    const double &dt = 0.0, const primitive &stateN = {},
    const vector3d<double> &pressGrad = {}, const tensor<double> &velGrad = {},
    const double &avgMach = 0.0, const double &maxMach = 0.0);

#endif
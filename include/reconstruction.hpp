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

#ifndef RECONSTRUCTIONHEADERDEF
#define RECONSTRUCTIONHEADERDEF

/* This header contains functions to reconstruct the primitive variables from 
 * the cell centers to the face centers.
 */

#include <string>                  // string

using std::string;

// forward class declarations
class primitive;

// function definitions
// function to calculate reconstruction of state variables from cell
// center to cell face assuming value at cell center is constant over cell
// volume; zeroth order reconstruction results in first order accuracy
primitive FaceReconConst(const primitive &state) { return state; }

// function to calculate reconstruction of state variables from cell
// center to cell face this function uses muscle extrapolation resulting in
// higher order accuracy
primitive FaceReconMUSCL(const primitive &, const primitive &,
                         const primitive &, const double &, const string &,
                         const double &, const double &, const double &);

// calculate face reconstruction using 5th order weno scheme
primitive FaceReconWENO(const primitive &, const primitive &, const primitive &,
                        const primitive &, const primitive &, const double &,
                        const double &, const double &, const double &,
                        const double &, const bool &);

#endif

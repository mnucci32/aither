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
#include "arrayView.hpp"
#include "macros.hpp"

using std::string;

// forward class declarations
class primitive;

// function definitions
// function to calculate reconstruction of state variables from cell
// center to cell face assuming value at cell center is constant over cell
// volume; zeroth order reconstruction results in first order accuracy
primitive FaceReconConst(const primitive &state) { return state; }
primitive FaceReconConst(const primitiveView &state) {
  return state.CopyData();
}

// function to calculate reconstruction of state variables from cell
// center to cell face this function uses muscle extrapolation resulting in
// higher order accuracy
template <typename T>
primitive FaceReconMUSCL(const T &, const T &, const T &, const double &,
                         const string &, const double &, const double &,
                         const double &);

// calculate face reconstruction using 5th order weno scheme
template <typename T>
primitive FaceReconWENO(const T &, const T &, const T &, const T &, const T &,
                        const double &, const double &, const double &,
                        const double &, const double &, const bool &);

// function to reconstruct cell variables to the face using central
// differences
template <typename T>
auto FaceReconCentral(const T &varU, const T &varD,
                      const vector<double> &cellWidth) {
  // varU -- variable at the cell center of the upwind cell
  // varD -- variable at the cell center of the downwind cell
  // cellWidth -- width of cells in stencil
  MSG_ASSERT(cellWidth.size() == 2, "cell width size is unexpected");

  // get coefficients
  const auto coeffs = LagrangeCoeff(cellWidth, 1, 0, 0);

  // reconstruct with central difference
  return coeffs[0] * varD + coeffs[1] * varU;
}

// function to reconstruct cell variables to the face using central
// differences (4th order)
template <typename T>
auto FaceReconCentral4th(const T &varU2, const T &varU1, const T &varD1,
                         const T &varD2, const vector<double> &cellWidth) {
  // varU2 -- variable at the cell center of the second upwind cell
  // varU1 -- variable at the cell center of the first upwind cell
  // varD1 -- variable at the cell center of the first downwind cell
  // varD2 -- variable at the cell center of the second downwind cell
  // cellWidth -- width of cells in stencil
  MSG_ASSERT(cellWidth.size() == 4, "cell width size is unexpected");

  // get coefficients
  const auto coeffs = LagrangeCoeff(cellWidth, 3, 1, 1);

  // reconstruct with central difference
  return coeffs[0] * varU2 + coeffs[1] * varU1 + coeffs[2] * varD1 +
      coeffs[3] * varD2;
}


#endif

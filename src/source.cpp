/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (michael.nucci@gmail.com)

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
#include <cmath>
#include <memory>
#include "source.hpp"
#include "turbulence.hpp"
#include "transport.hpp"
#include "vector3d.hpp"
#include "tensor.hpp"
#include "matrix.hpp"
#include "arrayView.hpp"
#include "primitive.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const source &m) {
  for (auto rr = 0; rr < m.Size(); rr++) {
    os << m[rr] << endl;
  }
  return os;
}

// Member function to calculate the source terms for the species equations
squareMatrix source::CalcChemSrc(const unique_ptr<chemistry> &chem,
                                 const primitiveView &state,
                                 const double &temperature) {
  const auto src = chem->SourceTerms(state.RhoVec(), temperature);
  std::copy(std::begin(src), std::end(src), std::begin(*this));

  return squareMatrix();
}

// Member function to calculate the source terms for the turbulence equations
squareMatrix source::CalcTurbSrc(
    const unique_ptr<turbModel> &turb, const primitiveView &state,
    const tensor<double> &velGrad, const vector3d<double> &tGrad,
    const vector3d<double> &tkeGrad, const vector3d<double> &omegaGrad,
    const unique_ptr<transport> &trans, const double &vol, const double &mut,
    const double &f1, const double &f2, const double &phi) {
  // turb -- turbulence model
  // state -- primitive variables
  // velGrad -- velocity gradient
  // tGrad -- temperature gradient
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // f2 -- second blending coefficient
  // phi -- factor to reduce tke destruction by for des

  // calculate turbulent source terms
  vector<double> turbSrc(this->NumTurbulence(), 0.0);
  const auto srcJac =
      turb->CalcTurbSrc(state, velGrad, tkeGrad, omegaGrad, trans, vol, mut, f1,
                        f2, phi, turbSrc);

  // assign turbulent source terms
  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] = turbSrc[ii];
  }

  // return source jacobian
  return srcJac;
}

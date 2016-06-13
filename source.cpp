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

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include "source.hpp"
#include "turbulence.hpp"
#include "primVars.hpp"
#include "vector3d.hpp"
#include "tensor.hpp"
#include "matrix.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::unique_ptr;

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const source &src) {
  os << src.SrcMass() << endl;
  os << src.SrcMomX() << endl;
  os << src.SrcMomY() << endl;
  os << src.SrcMomZ() << endl;
  os << src.SrcEngy() << endl;
  os << src.SrcTke() << endl;
  os << src.SrcOmg() << endl;
  return os;
}

// Member function to calculate the source terms for the turbulence equations
squareMatrix source::CalcTurbSrc(const unique_ptr<turbModel> &turb,
                                 const primVars &state,
                                 const tensor<double> &velGrad,
                                 const vector3d<double> &tGrad,
                                 const vector3d<double> &tkeGrad,
                                 const vector3d<double> &omegaGrad,
                                 const sutherland &suth, const double &vol,
                                 const double &mut, const double &f1) {
  // turb -- turbulence model
  // state -- primative variables
  // velGrad -- velocity gradient
  // tGrad -- temperature gradient
  // tkeGrad -- tke gradient
  // omegaGrad -- omega gradient
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  // calculate turbulent source terms
  auto ksrc = 0.0;
  auto wsrc = 0.0;
  const auto srcJac = turb->CalcTurbSrc(state, velGrad, tkeGrad, omegaGrad,
                                        suth, vol, mut, f1, ksrc, wsrc);

  // assign turbulent source terms
  data_[5] = ksrc;
  data_[6] = wsrc;

  // return source jacobian
  return srcJac;
}

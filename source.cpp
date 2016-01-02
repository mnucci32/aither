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
#include "gradients.hpp"

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
double source::CalcTurbSrc(const unique_ptr<turbModel> &turb,
                         const primVars &state, const gradients &grads,
                         const sutherland &suth, const idealGas &eqnState,
                         const double &wallDist, const int &ii, const int &jj,
                         const int &kk) {
  // turb -- turbulence model
  // state -- primative variables
  // grads -- gradients
  // suth -- sutherland's law for viscosity
  // eqnState -- equation of state
  // wallDist -- distance to nearest viscous wall
  // ii -- cell i-location to calculate source terms at
  // jj -- cell j-location to calculate source terms at
  // kk -- cell k-location to calculate source terms at

  // get cell gradients
  const auto vGrad = grads.VelGradCell(ii, jj, kk);
  const auto kGrad = grads.TkeGradCell(ii, jj, kk);
  const auto wGrad = grads.OmegaGradCell(ii, jj, kk);

  // calculate turbulent source terms
  auto ksrc = 0.0;
  auto wsrc = 0.0;
  auto srcJacSpecRad = turb->CalcTurbSrc(state, vGrad, kGrad, wGrad, suth,
                                         eqnState, wallDist, ksrc, wsrc);

  // assign turbulent source terms
  data_[5] = ksrc;
  data_[6] = wsrc;

  // return spectral radius of source jacobian
  return srcJacSpecRad;
}

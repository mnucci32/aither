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

#include <cmath>  // sqrt
#include <string>
#include <memory>
#include <algorithm>  // max
#include "inviscidFlux.hpp"
#include "eos.hpp"
#include "thermodynamic.hpp"
#include "primitive.hpp"
#include "matrix.hpp"
#include "turbulence.hpp"
#include "utility.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::copysign;
using std::unique_ptr;



template <typename T1, typename T2>
void inviscidFlux::AUSMFlux(const T1 &left, const T2 &right,
                            const unique_ptr<eos> &eqnState,
                            const unique_ptr<thermodynamic> &thermo,
                            const vector3d<double> &area,
                            const double &sos, const double &mPlusLBar,
                            const double &mMinusRBar, const double &pPlus,
                            const double &pMinus) {
  static_assert(std::is_same<primitive, T1>::value ||
                    std::is_same<primitiveView, T1>::value,
                "T1 requires primitive or primativeView type");
  static_assert(std::is_same<primitive, T2>::value ||
                    std::is_same<primitiveView, T2>::value,
                "T2 requires primitive or primativeView type");

  // calculate left flux
  const auto vl = mPlusLBar * sos;
  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    (*this)[ii] = left.RhoN(ii) * vl;
  }
  (*this)[this->MomentumXIndex()] =
      left.Rho() * vl * left.U() + pPlus * left.P() * area.X();
  (*this)[this->MomentumYIndex()] =
      left.Rho() * vl * left.V() + pPlus * left.P() * area.Y();
  (*this)[this->MomentumZIndex()] =
      left.Rho() * vl * left.W() + pPlus * left.P() * area.Z();
  (*this)[this->EnergyIndex()] =
      left.Rho() * vl * left.Enthalpy(eqnState, thermo);
  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] =
        left.Rho() * vl * left.TurbulenceN(ii);
  }

  // calculate right flux (add contribution)
  const auto vr = mMinusRBar * sos;
  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    (*this)[ii] = right.RhoN(ii) * vr;
  }
  (*this)[this->MomentumXIndex()] +=
      right.Rho() * vr * right.U() + pMinus * right.P() * area.X();
  (*this)[this->MomentumYIndex()] +=
      right.Rho() * vr * right.V() + pMinus * right.P() * area.Y();
  (*this)[this->MomentumZIndex()] +=
      right.Rho() * vr * right.W() + pMinus * right.P() * area.Z();
  (*this)[this->EnergyIndex()] +=
      right.Rho() * vr * right.Enthalpy(eqnState, thermo);
  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] =
        right.Rho() * vr * right.TurbulenceN(ii);
  }
}


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


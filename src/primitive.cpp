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

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <cmath>
#include <algorithm>  // max
#include "primitive.hpp"
#include "input.hpp"               // input
#include "turbulence.hpp"          // turbModel
#include "utility.hpp"
#include "wallLaw.hpp"
#include "wallData.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::min;
using std::unique_ptr;

// member function to initialize a state with nondimensional values
void primitive::NondimensionalInitialize(const unique_ptr<eos> &eqnState,
                                        const input &inp,
                                        const unique_ptr<transport> &trans,
                                        const int &parBlock,
                                        const unique_ptr<turbModel> &turb) {
  // get initial condition state for parent block
  auto ic = inp.ICStateForBlock(parBlock);
  auto massFracs = ic.MassFractions();

  // DEBUG -- need to make sure species are going into correct index
  auto ii = 0;
  for (auto &mf : massFracs) {
    (*this)[ii] = mf.second * ic.Density();
    ii++;
  }

  (*this)[this->MomentumXIndex()] = ic.Velocity().X();
  (*this)[this->MomentumYIndex()] = ic.Velocity().Y();
  (*this)[this->MomentumZIndex()] = ic.Velocity().Z();

  (*this)[this->EnergyIndex()] = ic.Pressure();
  
  if (this->HasTurbulenceData()) {
    // Initialize turbulence quantities based off of specified turublence
    // intensity and eddy viscosity ratio. This is the default for
    // STAR-CCM+
    this->ApplyFarfieldTurbBC(this->Velocity(), ic.TurbulenceIntensity(),
                              ic.EddyViscosityRatio(), trans, eqnState, turb);
  }
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const primitive &prim) {
  for (auto rr = 0; rr < prim.Size(); rr++) {
    os << prim[rr] << endl;
  }
  return os;
}

primitive primitive::UpdateWithConsVars(
    const unique_ptr<eos> &eqnState, const unique_ptr<thermodynamic> &thermo,
    const varArrayView &du, const unique_ptr<turbModel> &turb) const {
  return UpdatePrimWithCons((*this), eqnState, thermo, du, turb);
}


// member function to apply farfield turbulence boundary conditions
// using the method in STAR-CCM+ involving turbulence intensity and
// eddy viscosity ratio
void primitive::ApplyFarfieldTurbBC(const vector3d<double> &vel,
                                   const double &turbInten,
                                   const double &viscRatio,
                                   const unique_ptr<transport> &trans,
                                   const unique_ptr<eos> &eqnState,
                                   const unique_ptr<turbModel> &turb) {
  // vel -- reference velocity (nondimensionalized)
  // turbInten -- turbulence intensity at farfield
  // viscRatio -- eddy viscosity ratio at farfield
  // trans -- viscous transport model
  // eqnState -- equation of state
  // turb --  turbulence model

  (*this)[this->TurbulenceIndex()] = 1.5 * pow(turbInten * vel.Mag(), 2.0);
  (*this)[this->TurbulenceIndex() + 1] =
      this->Rho() * this->Tke() /
      (viscRatio *
       trans->Viscosity(this->Temperature(eqnState), this->MassFractions()));
  this->LimitTurb(turb);
}

void primitive::LimitTurb(const unique_ptr<turbModel> &turb) {
  // Adjust turbulence variables to be above minimum if necessary
  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] =
        max((*this)[TurbulenceIndex() + ii], turb->TurbMinN(ii));
  }
}

// return element by element absolute value
primitive primitive::Abs() const {
  auto abs = *this;
  std::transform(std::begin(*this), std::end(*this), std::begin(abs),
                 [](const double &val) { return fabs(val); });
  return abs;
}

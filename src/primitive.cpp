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
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <cmath>
#include <algorithm>  // max
#include "primitive.hpp"
#include "input.hpp"               // input
#include "physicsModels.hpp"
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
void primitive::NondimensionalInitialize(const physics &phys, const input &inp,
                                         const int &parBlock) {
  // get initial condition state for parent block
  auto ic = inp.ICStateForBlock(parBlock);
  auto massFracs = ic.MassFractions();

  for (auto &mf : massFracs) {
    auto ind = inp.SpeciesIndex(mf.first);
    (*this)[ind] = mf.second * ic.Density();
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
                              ic.EddyViscosityRatio(), phys);
  }
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const primitive &prim) {
  for (auto rr = 0; rr < prim.Size(); rr++) {
    os << prim[rr] << endl;
  }
  return os;
}

primitive primitive::UpdateWithConsVars(const physics &phys,
                                        const varArrayView &du) const {
  return UpdatePrimWithCons((*this), phys, du);
}

// member function to apply farfield turbulence boundary conditions
// using the method in STAR-CCM+ involving turbulence intensity and
// eddy viscosity ratio
void primitive::ApplyFarfieldTurbBC(const vector3d<double> &vel,
                                    const double &turbInten,
                                    const double &viscRatio,
                                    const physics &phys) {
  // vel -- reference velocity (nondimensionalized)
  // turbInten -- turbulence intensity at farfield
  // viscRatio -- eddy viscosity ratio at farfield
  // phys -- physics models

  (*this)[this->TurbulenceIndex()] = 1.5 * pow(turbInten * vel.Mag(), 2.0);
  (*this)[this->TurbulenceIndex() + 1] =
      this->Rho() * this->Tke() /
      (viscRatio * phys.Transport()->Viscosity(this->Temperature(phys.EoS()),
                                               this->MassFractions()));
  this->LimitTurb(phys.Turbulence());
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

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

#include <vector>
#include <iostream>     // cout
#include <cstdlib>      // exit()
#include "transport.hpp"
#include "fluid.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;

// constructor for sutherland class
sutherland::sutherland(const vector<fluid> &fl, const double &tRef,
                       const double &rRef, const double &lRef,
                       const double &aRef) {
  auto numSpecies = fl.size();
  // size coefficient vectors
  viscC1_.reserve(numSpecies);
  viscS_.reserve(numSpecies);
  condC1_.reserve(numSpecies);
  condS_.reserve(numSpecies);
  muRef_.reserve(numSpecies);

  tRef_ = tRef;

  // get data from fluid class
  for (auto &f : fl) {
    auto viscCoeffs = f.ViscosityCoeffs();
    viscC1_.push_back(viscCoeffs[0]);
    viscS_.push_back(viscCoeffs[1]);

    auto condCoeffs = f.ConductivityCoeffs();
    condC1_.push_back(condCoeffs[0]);
    condS_.push_back(condCoeffs[1]);

    muRef_.push_back(viscCoeffs[0] * pow(tRef_, 1.5) / (tRef_ + viscCoeffs[1]));
  }

  // DEBUG
  // calculate reference viscosity for reference mixture and set scaling
  this->SetScaling(rRef, lRef, muRef_[0], aRef);
}

// Functions for sutherland class
// DEBUG -- use Wilke's method
double sutherland::Viscosity(const double &t) const {
  // Dimensionalize temperature
  const auto temp = t * tRef_;

  // Calculate viscosity
  const auto mu = (viscC1_[0] * pow(temp, 1.5)) / (temp + viscS_[0]);

  // Nondimensionalize viscosity
  return (mu / muRef_[0]);
}

double sutherland::EffectiveViscosity(const double &t) const {
  // Get viscosity and scale
  return this->Viscosity(t) * this->NondimScaling();
}

double sutherland::Lambda(const double &mu) const {
  // Calculate lambda (2nd coeff of viscosity)
  return bulkVisc_ - (2.0 / 3.0) * mu;
}

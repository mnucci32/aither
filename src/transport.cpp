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

#include <iostream>     // cout
#include <cstdlib>      // exit()
#include "transport.hpp"

using std::cout;
using std::endl;
using std::cerr;

// Functions for sutherland class
double sutherland::Viscosity(const double &t) const {
  // Dimensionalize temperature
  const auto temp = t * tRef_;

  // Calculate viscosity
  const auto mu = (cOne_ * pow(temp, 1.5)) / (temp + S_);

  // Nondimensionalize viscosity
  return (mu / muRef_);
}

double sutherland::EffectiveViscosity(const double &t) const {
  // Get viscosity and scale
  return this->Viscosity(t) * this->NondimScaling();
}

double sutherland::Lambda(const double &mu) const {
  // Calculate lambda (2nd coeff of viscosity)
  return bulkVisc_ - (2.0 / 3.0) * mu;
}

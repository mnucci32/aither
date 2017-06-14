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
#include "thermodynamic.hpp"
#include "utility.hpp"  // FindRoot

using std::cout;
using std::endl;
using std::cerr;

// Member functions for thermodynamic class

double caloricallyPerfect::TemperatureFromSpecEnergy(const double& e) const {
  const auto t = 1.0;  // cpg has constant Cv, so value of t is meaningless
  return e / this->Cv(t);
}

double thermallyPerfect::TemperatureFromSpecEnergy(const double& e) const {
  auto temperature = 0.0;
  auto func = [&](const double& t) {
    temperature = t;
    //cout << "t: " << temperature << endl;
    return e - this->SpecEnergy(t);
  };
  FindRoot(func, 1.0e-8, 1.0e4, 1.0e-8);
  
  return temperature;
}
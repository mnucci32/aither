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
#include "diffusion.hpp"
#include "fluid.hpp"
#include "macros.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;

// constructor for schmidt class
schmidt::schmidt(const vector<fluid> &fl) {
  diffCoeff_.reserve(fl.size());
  for (auto &f : fl) {
    diffCoeff_.push_back(f.SchmidtNumber());
  }
}

// member function to use Wilke's method to calculate mixture viscosity
double schmidt::WilkesDiff(const vector<double> &specDiff,
                           const vector<double> &mf) const {
  // specDiff -- vector of species viscosities
  // mf -- vector of species mass fractions
  MSG_ASSERT(mf.size() == specDiff.size(), "mismatch in species size");

  return 0.0;
}

// Functions for sutherland class
double schmidt::SpeciesDiffusion(const int &ii, const int &jj) const {
  MSG_ASSERT(ii < this->NumSpecies(), "Accessing index out of range");
  MSG_ASSERT(jj < this->NumSpecies(), "Accessing index out of range");
  // DEBUG
  return diffCoeff_[ii];
}

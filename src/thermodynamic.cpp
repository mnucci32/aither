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
#include "fluid.hpp"
#include "macros.hpp"

using std::cout;
using std::endl;
using std::cerr;

// constructor
caloricallyPerfect::caloricallyPerfect(const vector<fluid> &fl) {
  const auto numSpecies = fl.size();
  gamma_.reserve(numSpecies);
  for (auto &f : fl) {
    const auto n = f.N();
    gamma_.push_back((n + 1.0) / n);  // non-dim R = 1 / gamma
  }
}

thermallyPerfect::thermallyPerfect(const vector<fluid> &fl) {
  const auto numSpecies = fl.size();
  n_.reserve(numSpecies);
  vibTemp_.reserve(numSpecies);
  nonDimR_.reserve(numSpecies);
  for (auto &f : fl) {
    const auto n = f.N();
    n_.push_back(n);
    vibTemp_.push_back(f.VibrationalTemperature());
    nonDimR_.push_back(n / (n + 1.0));
  }
}

// ---------------------------------------------------------------------------
// Member functions for calorically perfect class
double caloricallyPerfect::Gamma(const double& t,
                                 const vector<double>& mf) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(mf.size()),
             "species size mismatch");
  auto gamma = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    gamma += mf[ss] * this->SpeciesGamma(t, ss);
  }
  return gamma;
}

double caloricallyPerfect::TemperatureFromSpecEnergy(
    const double& e, const vector<double>& mf) const {
  const auto t = 1.0;  // cpg has constant Cv, so value of t is meaningless
  return e / this->Cv(t, mf);
}


// ---------------------------------------------------------------------------
// thermally perfect functions
double thermallyPerfect::Cp(const double& t, const vector<double>& mf) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(mf.size()),
             "species size mismatch");
  auto cp = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    cp += mf[ss] * this->SpeciesCp(t, ss);
  }
  return cp;
}
double thermallyPerfect::Cv(const double& t, const vector<double>& mf) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(mf.size()),
             "species size mismatch");
  auto cv = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    cv += mf[ss] * this->SpeciesCv(t, ss);
  }
  return cv;
}

double thermallyPerfect::SpecEnergy(const double& t,
                                    const vector<double>& mf) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(mf.size()),
             "species size mismatch");
  auto e = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    e += mf[ss] * (nonDimR_[ss] * (n_[ss] * t + this->VibEqTerm(t, ss)));
  }
  return e;
}

double thermallyPerfect::SpecEnthalpy(const double& t,
                                    const vector<double>& mf) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(mf.size()),
             "species size mismatch");
  auto h = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    h += mf[ss] * (nonDimR_[ss] * ((n_[ss] + 1) * t + this->VibEqTerm(t, ss)));
  }
  return h;
}

double thermallyPerfect::TemperatureFromSpecEnergy(
    const double& e, const vector<double>& mf) const {
  auto temperature = 0.0;
  auto func = [&](const double& t) {
    temperature = t;
    return e - this->SpecEnergy(t, mf);
  };
  FindRoot(func, 1.0e-8, 1.0e4, 1.0e-8);
  
  return temperature;
}
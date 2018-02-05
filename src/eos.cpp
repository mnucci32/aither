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

#include <iostream>     // cout
#include <cstdlib>      // exit()
#include "eos.hpp"
#include "macros.hpp"
#include "fluid.hpp"

using std::cout;
using std::endl;
using std::cerr;

// constructor
idealGas::idealGas(const vector<fluid> &fl, const double &tRef,
                   const double &aRef) {
  const auto numSpecies = fl.size();
  gasConst_.reserve(numSpecies);
  for (auto &f : fl) {
    // nondimensionalize gas constant
    gasConst_.push_back(f.GasConstant() * tRef / (aRef * aRef));
  }
}

// Member functions for idealGas class
// These functions calculate values using the ideal gas equation of state
// P = rho * R * T
double idealGas::PressFromEnergy(const unique_ptr<thermodynamic> &thermo,
                                 const vector<double> &rho,
                                 const double &energy,
                                 const double &vel) const {
  const auto specEnergy = energy - 0.5 * vel * vel;
  const auto rhoSum = std::accumulate(rho.begin(), rho.end(), 0.0);
  vector<double> mf(rho.size());
  for (auto ii = 0U; ii < mf.size(); ++ii) {
    mf[ii] = rho[ii] / rhoSum;
  }
  const auto temperature = thermo->TemperatureFromSpecEnergy(specEnergy, mf);
  return this->PressureRT(rho, temperature);
}

double idealGas::PressureRT(const vector<double> &rho,
                            const double &temperature) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(rho.size()),
             "species size mismatch");
  auto p = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    p += rho[ss] * gasConst_[ss] * temperature;
  }
  return p;
}

double idealGas::SpecEnergy(const unique_ptr<thermodynamic> &thermo,
                            const double &t, const vector<double> &mf) const {
  return thermo->SpecEnergy(t, mf);
}

double idealGas::Energy(const double &specEn, const double &vel) const {
  return specEn + 0.5 * vel * vel;
}

double idealGas::SpeciesEnthalpy(const unique_ptr<thermodynamic> &thermo,
                                 const double &t, const double &vel,
                                 const int &ss) const {
  return thermo->SpeciesSpecEnthalpy(t, ss) + 0.5 * vel * vel;
}

double idealGas::Enthalpy(const unique_ptr<thermodynamic> &thermo,
                          const double &t, const double &vel,
                          const vector<double> &mf) const {
  return thermo->SpecEnthalpy(t, mf) + 0.5 * vel * vel;
}

double idealGas::SoS(const unique_ptr<thermodynamic> &thermo,
                     const double &pressure, const vector<double> &rho) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(rho.size()),
             "species size mismatch");
  const auto rhoSum = std::accumulate(rho.begin(), rho.end(), 0.0);
  vector<double> mf(rho.size());
  for (auto ii = 0U; ii < mf.size(); ++ii) {
    mf[ii] = rho[ii] / rhoSum;
  }
  const auto t = this->Temperature(pressure, rho);
  const auto gamma = thermo->Gamma(t, mf);
  return sqrt(gamma * pressure / rhoSum);
}

double idealGas::Temperature(const double &pressure,
                             const vector<double> &rho) const {
  MSG_ASSERT(this->NumSpecies() == static_cast<int>(rho.size()),
             "species size mismatch");
  auto rhoR = 0.0;
  for (auto ss = 0; ss < this->NumSpecies(); ++ss) {
    rhoR += rho[ss] * gasConst_[ss];
  }
  return pressure / rhoR;
}

double idealGas::DensityTP(const double &temp, const double &press,
                           const vector<double> &mf) const {
  auto R = this->MixtureGasConstant(mf);
  return press / (R * temp);
}

/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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
#include <cmath>
#include "transport.hpp"
#include "fluid.hpp"
#include "macros.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;

// constructor for sutherland class
sutherland::sutherland(const vector<fluid> &fl, const double &tRef,
                       const double &rRef, const double &lRef,
                       const double &aRef, const vector<double> &mixRef) {
  MSG_ASSERT(fl.size() == mixRef.size(), "mixture must be same size as fluids");

  auto numSpecies = fl.size();
  // size coefficient vectors
  viscC1_.reserve(numSpecies);
  viscS_.reserve(numSpecies);
  condC1_.reserve(numSpecies);
  condS_.reserve(numSpecies);
  molarMass_.reserve(numSpecies);

  tRef_ = tRef;

  vector<double> muSpecRef;
  muSpecRef.reserve(numSpecies);

  // get data from fluid class
  for (auto &f : fl) {
    auto viscCoeffs = f.ViscosityCoeffs();
    viscC1_.push_back(viscCoeffs[0]);
    viscS_.push_back(viscCoeffs[1]);

    auto condCoeffs = f.ConductivityCoeffs();
    condC1_.push_back(condCoeffs[0]);
    condS_.push_back(condCoeffs[1]);

    muSpecRef.push_back(viscCoeffs[0] * pow(tRef_, 1.5) /
                        (tRef_ + viscCoeffs[1]));
    molarMass_.push_back(f.MolarMass());
  }

  // calculate reference viscosity for reference mixture and set scaling
  muMixRef_ =
      numSpecies == 1 ? muSpecRef[0] : this->WilkesVisc(muSpecRef, mixRef);
  kNonDim_ = (aRef * aRef * muMixRef_) / tRef_;
  this->SetScaling(rRef, lRef, muMixRef_, aRef);
}

// member function to use Wilke's method to calculate mixture viscosity
double sutherland::WilkesVisc(const vector<double> &specVisc,
                              const vector<double> &mf) const {
  // specVisc -- vector of species viscosities
  // mf -- vector of species mass fractions
  MSG_ASSERT(mf.size() == specVisc.size(), "mismatch in species size");

  // calculate mole fractions
  auto moleFrac = this->MoleFractions(mf);

  auto mixtureVisc = 0.0;
  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    auto moleVisc = moleFrac[ii] * specVisc[ii];
    auto denom = 0.0;
    for (auto jj = 0; jj < this->NumSpecies(); ++jj) {
      denom += moleFrac[jj] / sqrt(1.0 + molarMass_[ii] / molarMass_[jj]) *
               pow(1.0 + sqrt(specVisc[ii] / specVisc[jj]) *
                             pow(molarMass_[jj] / molarMass_[ii], 0.25),
                   2.0);
    }
    mixtureVisc += moleVisc / denom;
  }
  return 4.0 / sqrt(2.0) * mixtureVisc;
}

// member function to use Wilke's method to calculate mixture conductivity
double sutherland::WilkesCond(const vector<double> &specCond,
                              const vector<double> &mf) const {
  // calculate mole fractions
  auto moleFrac = this->MoleFractions(mf);

  auto weightedAvg = 0.0;
  auto harmonicAvg = 0.0;
  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    weightedAvg += moleFrac[ii] * specCond[ii];
    harmonicAvg += moleFrac[ii] / specCond[ii];
  }
  harmonicAvg = 1.0 / harmonicAvg;
  return 0.5 * (weightedAvg + harmonicAvg);
}


// Functions for sutherland class
double sutherland::SpeciesViscosity(const double &t, const int &ii) const {
  MSG_ASSERT(ii < this->NumSpecies(), "Accessing index out of range");
  // Dimensionalize temperature
  const auto temp = t * tRef_;
  // Calculate viscosity
  const auto mu = (viscC1_[ii] * pow(temp, 1.5)) / (temp + viscS_[ii]);
  // Nondimensionalize viscosity
  return mu / muMixRef_;
}

double sutherland::SpeciesConductivity(const double &t, const int &ii) const {
  MSG_ASSERT(ii < this->NumSpecies(), "Accessing index out of range");
  // Dimensionalize temperature
  const auto temp = t * tRef_;
  // Calculate conductivity
  const auto k = (condC1_[ii] * pow(temp, 1.5)) / (temp + condS_[ii]);
  // Nondimensionalize conductivity
  return k / kNonDim_;
}


vector<double> sutherland::MoleFractions(const vector<double> &mf) const {
  MSG_ASSERT(static_cast<int>(mf.size()) == this->NumSpecies(),
             "mismatch in species size");

  vector<double> moleFrac(this->NumSpecies());
  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    moleFrac[ii] = mf[ii] / molarMass_[ii];
  }
  auto moleFracSum = std::accumulate(moleFrac.begin(), moleFrac.end(), 0.0);
  std::for_each(moleFrac.begin(), moleFrac.end(),
                [&moleFracSum](auto &val) { val /= moleFracSum; });
  return moleFrac;
}

// member function to use Wilke's method to calculate mixture viscosity
double sutherland::Viscosity(const double &t, const vector<double> &mf) const {
  MSG_ASSERT(static_cast<int>(mf.size()) == this->NumSpecies(),
             "mismatch in species size");

  if (this->NumSpecies() == 1) {
    return this->SpeciesViscosity(t, 0);
  } else {
    // calculate species viscosity and mole fractions
    vector<double> specVisc(this->NumSpecies());
    for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
      specVisc[ii] = this->SpeciesViscosity(t, ii);
    }
    return this->WilkesVisc(specVisc, mf);
  }
}

double sutherland::EffectiveViscosity(const double &t,
                                      const vector<double> &mf) const {
  // Get viscosity and scale
  return this->Viscosity(t, mf) * this->NondimScaling();
}

double sutherland::Lambda(const double &mu) const {
  // Calculate lambda (2nd coeff of viscosity)
  return bulkVisc_ - (2.0 / 3.0) * mu;
}

// member function to use Wilke's method to calculate mixture conductivity
double sutherland::Conductivity(const double &t,
                                const vector<double> &mf) const {
  if (this->NumSpecies() == 1) {
    return this->SpeciesConductivity(t, 0);
  } else {
    // calculate species conductivity and mole fractions
    vector<double> specCond(this->NumSpecies());
    for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
      specCond[ii] = this->SpeciesConductivity(t, ii);
    }
    return this->WilkesCond(specCond, mf);
  }
}

double sutherland::EffectiveConductivity(const double &t,
                                         const vector<double> &mf) const {
  // Get viscosity and scale
  return this->Conductivity(t, mf) * this->NondimScaling();
}

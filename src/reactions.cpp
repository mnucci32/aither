/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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

#include <string>
#include <iostream>
#include <memory>
#include "reactions.hpp"
#include "inputStates.hpp"
#include "input.hpp"
#include "macros.hpp"
#include "thermodynamic.hpp"

using std::string;
using std::cerr;
using std::cout;
using std::endl;
using std::unique_ptr;

// constructor
reaction::reaction(const string &str, const input &inp) {
  const auto ns = inp.NumSpecies();
  stoichReactants_.resize(ns, 0.0);
  stoichProducts_.resize(ns, 0.0);
  species_.reserve(ns);
  universalGasConst_ = inp.Fluid(0).UniversalGasConstant();

  // get species names
  for (auto ii = 0; ii < ns; ++ii) {
    species_.push_back(inp.Fluid(ii).Name());
  }

  // split into reaction and rate data
  const auto tokens = Tokenize(str, ":");
  MSG_ASSERT(tokens.size() == 2U, "problem with reaction data");

  // determine if forward only reaction
  auto rxEnd = tokens[0].find("<=>");
  const auto isEquilRx = rxEnd != string::npos;
  isForwardOnly_ = false;
  if (!isEquilRx) {
    rxEnd = tokens[0].find("=>");
    isForwardOnly_ = rxEnd != string::npos;
    modifyReactants_.resize(ns, 0.0);
  }
  MSG_ASSERT(isEquilRx || isForwardOnly_, "reaction type not recognized");

  // split reactants and products
  const auto prodStart = isEquilRx ? rxEnd + 3 : rxEnd + 2;
  const auto rxStr = tokens[0].substr(0, rxEnd);
  const auto reactants = Tokenize(rxStr, "+");
  const auto prodStr = tokens[0].substr(prodStart);
  const auto products = Tokenize(prodStr, "+");

  // split into coefficients and species
  for (const auto &pr : products) {
    const auto speciesStart = pr.find_first_not_of("0123456789.");
    const auto stoich =
        speciesStart == 0 ? 1.0 : std::stod(Trim(pr.substr(0, speciesStart)));
    const auto species = Trim(pr.substr(speciesStart));
    if (!inp.HaveSpecies(species)) {
      cerr << "species " << species << " is in reaction, but not in simulation"
           << endl;
      exit(EXIT_FAILURE);
    }
    stoichProducts_[inp.SpeciesIndex(species)] += stoich;
  }
  for (const auto &re : reactants) {
    const auto speciesStart = re.find_first_not_of("0123456789.");
    const auto stoich =
        speciesStart == 0 ? 1.0 : std::stod(Trim(re.substr(0, speciesStart)));
    const auto species = Trim(re.substr(speciesStart));
    if (!inp.HaveSpecies(species)) {
      cerr << "species " << species << " is in reaction, but not in simulation"
           << endl;
      exit(EXIT_FAILURE);
    }
    stoichReactants_[inp.SpeciesIndex(species)] += stoich;
  }

  // get rate data -----------------------------------------------------------
  const auto rateTokens = Tokenize(tokens[1], ";");
  MSG_ASSERT(rateTokens.size() >= 1, "missing rate data");
  for (const auto &rt : rateTokens) {
    const auto rateData = Tokenize(rt, "=", 1);
    MSG_ASSERT(rateData.size() == 2U, "rate data is wrong size");
    if (rateData[0] == "forwardRate") {
      if (rateData[1].find("arrhenius") == string::npos) {
        cerr << "ERROR: forwardRate must have arrhenius information" << endl;
        exit(EXIT_FAILURE);
      }
      const auto start = rateData[1].find("(") + 1;
      const auto end = rateData[1].find(")") + 1;
      const auto range = end - start + 1;  // +/-1 to ignore ()
      auto arrhenius = rateData[1].substr(start, range);
      const auto arrTokens = Tokenize(arrhenius, ",");
      MSG_ASSERT(arrTokens.size() == 3U, "arrhenius has wrong data size");
      for (const auto &arr : arrTokens) {
        const auto arrData = Tokenize(arr, "=");
        MSG_ASSERT(arrData.size() == 2U, "arrhenius data has wrong size");
        if (arrData[0] == "C") {
          arrheniusC_ = std::stod(arrData[1]);
        } else if (arrData[0] == "eta") {
          arrheniusEta_ = std::stod(arrData[1]);
        } else if (arrData[0] == "theta") {
          arrheniusTheta_ = std::stod(arrData[1]);
        } else {
          cerr << "ERROR: arrhenius data " << arr << " is not recognized"
               << endl;
          exit(EXIT_FAILURE);
        }
      }
    } else if (rateData[0] == "modifyReactants") {
      if (!isForwardOnly_) {
        cerr << "ERROR: modifyReactants is only supported for forward reactions"
             << endl;
        exit(EXIT_FAILURE);
      }
      const auto start = rateData[1].find("[") + 1;
      const auto end = rateData[1].find("]") - 1;
      const auto range = end - start + 1;  // +/-1 to ignore []
      auto mod = rateData[1].substr(start, range);
      auto modTokens = Tokenize(mod, ",");
      for (const auto &mt : modTokens) {
        // tokenize each token into name / mass fraction pairs
        auto species = Tokenize(mt, "=");
        if (species.size() != 2) {
          cerr << "ERROR. Problem with reading species. Substring is "
               << mod << endl;
          exit(EXIT_FAILURE);
        }
        if (!inp.HaveSpecies(species[0])) {
          cerr << "species " << species[0]
               << " is in reaction, but not in simulation" << endl;
          exit(EXIT_FAILURE);
        }
        modifyReactants_[inp.SpeciesIndex(species[0])] = std::stod(species[1]);
      }
    } else {
      cerr << "ERROR: can't identify rate data: " << rt << endl;
      exit(EXIT_FAILURE);
    }
  }
}

void reaction::Print(std::ostream &os) const {
  auto first = true;
  for (auto ii = 0U; ii < species_.size(); ++ii) {
    if (stoichReactants_[ii] > 0.0) {
      if (!first) {
        os << "+ ";
      }
      first = false;
      os << stoichReactants_[ii] << " " << species_[ii] << " ";
    }
  }
  isForwardOnly_ ? os << "=> " : os << "<=> ";
  first = true;
  for (auto ii = 0U; ii < species_.size(); ++ii) {
    if (stoichProducts_[ii] > 0.0) {
      if (!first) {
        os << "+ ";
      }
      first = false;
      os << stoichProducts_[ii] << " " << species_[ii] << " ";
    }
  }
  os << ": forwardRate=arrhenius(C=" << arrheniusC_
     << ", eta=" << arrheniusEta_ << ", theta=" << arrheniusTheta_ << ")";
  if (modifyReactants_.size() > 0) {
    os << "; modifyReactants=[";
    auto first = true;
    for (auto ii = 0U; ii < species_.size(); ++ii) {
      if (modifyReactants_[ii] > 0.0) {
        if (!first) {
          os << ", ";
        }
        first = false;
        os << species_[ii] << "=" << modifyReactants_[ii];
      }
    }
    os << "]";
  }
}

std::ostream &operator<<(std::ostream &os, const reaction &rx) {
  rx.Print(os);
  return os;
}

// member function to calculate equilibrium reaction rate
double reaction::EquilibriumRate(const double &t, const double &refP,
                                 const vector<double> &gibbsTerm) const {
  MSG_ASSERT(gibbsTerm.size() == stoichProducts_.size(),
             "species size mismatch");
  auto prodMinReac = 0.0;
  auto expTerm = 0.0;
  for (auto ss = 0U; ss < stoichProducts_.size(); ++ss) {
    auto specProdMinReac = stoichProducts_[ss] - stoichReactants_[ss];
    prodMinReac += specProdMinReac;
    expTerm += gibbsTerm[ss] * specProdMinReac;
  }
  const auto kp = std::exp(-expTerm);
  // reaction rate based on concentration
  return pow(refP / (universalGasConst_ * t), prodMinReac) * kp;
}
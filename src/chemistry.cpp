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

#include <algorithm>
#include <fstream>
#include <string>
#include <iostream>
#include <memory>
#include "chemistry.hpp"
#include "macros.hpp"
#include "utility.hpp"
#include "inputStates.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "primitive.hpp"
#include "physicsModels.hpp"

using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::unique_ptr;

squareMatrix chemistry::SourceJac(const primitive &state, const double &t,
                                  const vector<double> &gibbsTerm,
                                  const vector<double> &w, const physics &phys,
                                  const bool &isBlockMatrix) const {
  return isBlockMatrix ? squareMatrix(state.Size() - state.NumTurbulence())
                       : squareMatrix();
}

void reacting::ReadFromFile(const input &inp) {
  auto fname = inp.ChemistryMechanism() + ".mch";
  // open mechanism file -- first try run directory, then mechanism database
  ifstream mchFile(fname, std::ios::in);
  if (mchFile.fail()) {
    auto mechanismFile = GetEnvironmentVariable("AITHER_INSTALL_DIRECTORY") +
                         "/chemistryMechanisms/" + fname;
    mchFile.open(mechanismFile, std::ios::in);
    if (mchFile.fail()) {
      cerr << "ERROR: Error in reacting::ReadFromFile(). File " << fname
           << " did not open correctly!!!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  cout << "Reading reactions from " << fname << endl;
  string line = "";
  while (getline(mchFile, line)) {
    // remove leading and trailing whitespace and ignore comments
    line = Trim(line);

    if (line.length() > 0) {  // only proceed if line has data
      // read in reaction
      reaction rx(line, inp);
      cout << rx << endl;
      rx.Nondimensionalize(inp.TRef(), inp.LRef(), inp.ARef());
      reactions_.push_back(rx);
    }
  }
  cout << endl;

  // close mechanism file
  mchFile.close();
}

vector<double> reacting::SourceTerms(const vector<double> &rho, const double &t,
                                     const vector<double> &gibbsTerm,
                                     double &specRad) const {
  MSG_ASSERT(rho.size() == molarMass_.size(), "species size mismatch");

  vector<double> src(rho.size(), 0.0);
  if (t < freezingTemperature_) {  // no reactions
    return src;
  }

  // calculate mass fractions
  const auto rhoSum = std::accumulate(std::begin(rho), std::end(rho), 0.0);
  vector<double> mf;
  mf.reserve(rho.size());
  for (const auto &ri : rho) {
mf.push_back(ri / rhoSum);
  }

  vector<double> destruction(src.size(), 0.0);
  // loop over all species
  for (auto ss = 0U; ss < src.size(); ++ss) {
    // loop over all reactions
    for (const auto &rx : reactions_) {
      const auto prodMinReac = rx.StoichProduct(ss) - rx.StoichReactant(ss);
      const auto kf = rx.ForwardRate(t);
      auto fwdTerm = 1.0;
      for (auto ff = 0U; ff < src.size(); ++ff) {
        fwdTerm *= pow(rho[ff] / molarMass_[ff], rx.StoichReactant(ff));
      }

      const auto kb = rx.BackwardRate(t, refP_, gibbsTerm);
      auto bckTerm = 1.0;
      for (auto ff = 0U; ff < src.size(); ++ff) {
        bckTerm *= pow(rho[ff] / molarMass_[ff], rx.StoichProduct(ff));
      }
      destruction[ss] -= kb * bckTerm;
      src[ss] += prodMinReac * (kf * fwdTerm - kb * bckTerm);
    }
    destruction[ss] *= molarMass_[ss] / mf[ss];
    src[ss] *= molarMass_[ss];
  }

  specRad = *std::min_element(std::begin(destruction), std::end(destruction));
  return src;
}

squareMatrix reacting::SourceJac(const primitive &state, const double &t,
                                 const vector<double> &gibbsTerm,
                                 const vector<double> &w, const physics &phys,
                                 const bool &isBlockMatrix) const {
  MSG_ASSERT(gibbsTerm.size() == molarMass_.size(), "species size mismatch");
  MSG_ASSERT(w.size() == molarMass_.size(), "species size mismatch");
  if (!isBlockMatrix) {  // using spectral radius instead
    return squareMatrix();
  }

  auto chemJac = squareMatrix(state.Size() - state.NumTurbulence());
  if (t < freezingTemperature_) {  // no reactions
    return chemJac;
  }

  // calculate perturbation step for species equations
  constexpr auto epsilon = 1.0e-30;
  const auto rhoSum = state.Rho();
  const auto hRho = epsilon * rhoSum;
  auto specRad = 0.0;  // dummy value, not used

  for (auto cc = 0U; cc < w.size(); ++cc) {
    // preturb state - want derivative wrt to conservative variables, but for 
    // species terms conservative and primitive variables are the same - more
    // efficient to use primitive variables
    auto preturbed = state.RhoVec();
    preturbed[cc] += hRho;
    // compute chemistry source terms
    const auto wPreturbed = this->SourceTerms(preturbed, t, gibbsTerm, specRad);

    for (auto rr = 0U; rr < w.size(); ++rr) {
      chemJac(rr, cc) = (wPreturbed[rr] - w[rr]) / hRho;
    }
  }

  // calculate perturbation step for energy equation
  const auto maxEigV = state.Velocity().Mag() + state.SoS(phys);
  const auto hEnergy = epsilon * rhoSum * maxEigV * maxEigV;
  const auto ei = state.EnergyIndex();
  auto preturbedCons = state.ConsVars(phys);
  preturbedCons[ei] += hEnergy;
  const auto preturbed = primitive(preturbedCons, phys);
  const auto wPreturbed =
      this->SourceTerms(preturbed.RhoVec(), t, gibbsTerm, specRad);
  for (auto rr = 0U; rr < w.size(); ++rr) {
    chemJac(rr, ei) = (wPreturbed[rr] - w[rr]) / hEnergy;
  }

  return chemJac;
}

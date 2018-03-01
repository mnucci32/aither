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
#include "chemistry.hpp"
#include "macros.hpp"
#include "utility.hpp"
#include "inputStates.hpp"
#include "input.hpp"

using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

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
      reactions_.push_back(rx);
    }
  }
  cout << endl;

  // close mechanism file
  mchFile.close();
}

vector<double> reacting::SourceTerms(const vector<double>& rho,
                                     const double& t) const {
  MSG_ASSERT(rho.size() == molarMass_.size(), "species size mismatch");
  vector<double> src(rho.size(), 0.0);
  // loop over all species
  for (auto ss = 0U; ss < src.size(); ++ss) {
    // loop over all reactions
    for (const auto &rx : reactions_) {
      const auto prodMinReac = rx.StoichProduct(ss) - rx.StoichReactant(ss);
      const auto kf = rx.ForwardRate(t);
      auto fwdTerm = 0.0;
      for (auto ff = 0U; ff < src.size(); ++ff) {
        fwdTerm *= pow(rho[ff] / molarMass_[ff], rx.StoichReactant(ff));
      }

      const auto kb = rx.BackwardRate(t);
      auto bckTerm = 0.0;
      for (auto ff = 0U; ff < src.size(); ++ff) {
        bckTerm *= pow(rho[ff] / molarMass_[ff], rx.StoichProduct(ff));
      }
      src[ss] += prodMinReac * (kf * fwdTerm - kb * bckTerm);
    }
    src[ss] *= molarMass_[ss];
  }
  return src;
}

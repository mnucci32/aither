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

#include <iostream>     // cout
#include <vector>
#include <string>
#include <fstream>
#include "fluid.hpp"
#include "inputStates.hpp"
#include "utility.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

// construct initial condition from string
fluid::fluid(string &str, const string name) {
  const auto start = str.find("(") + 1;
  const auto end = str.find(")") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore ()
  auto state = str.substr(start, range);
  const auto id = str.substr(0, start - 1);
  if (id != name) {
    cerr << "ERROR. State condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // parameter counters
  auto nameCount = 0;
  auto mfCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with " << name << " parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "referenceMassFraction") {
      massFracRef_ = stod(RemoveTrailing(param[1], ","));
      mfCount++;
    } else if (param[0] == "name") {
      name_ = RemoveTrailing(param[1], ",");
      nameCount++;
    } else {
      cerr << "ERROR. " << name << " specifier " << param[0]
           << " is not recognized in fluid definition!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (nameCount != 1 || mfCount != 1) {
    cerr << "ERROR. For fluid 'name' and 'referenceMassFraction' must be "
            "specified"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (name != "fluid") {
    cerr << "ERROR. To specify fluid, properties must be enclosed in fluid()."
         << endl;
  }

  // get data from fluid database
  this->GetDatabaseProperties(name_);
}

void fluid::Nondimensionalize(const double &tRef, const double &rRef,
                              const double &aRef, const double &lRef) {
  if (!this->IsNondimensional()) {
    std::for_each(vibTemp_.begin(), vibTemp_.end(),
                  [&tRef](auto &val) { val /= tRef; });
    // converting hf & s from mol to kg values, then nondimensionalizing
    heatOfFormation_ /= molarMass_ * (aRef * aRef);
    refS_ /= molarMass_ / tRef * (aRef * aRef);
    molarMass_ /= rRef / pow(lRef, 3.0);
    refP_ /= rRef * aRef * aRef;
    refT_ /= tRef;
    universalGasConst_ /= aRef * aRef * rRef / (tRef * pow(lRef, 3.0));
    this->SetNondimensional(true);
  }
}

void fluid::GetDatabaseProperties(const string &name) {
  auto fname = name + ".dat";
  // open database file -- first try run directory, then fluid database
  ifstream datFile(fname, std::ios::in);
  if (datFile.fail()) {
    auto databaseFile = GetEnvironmentVariable("AITHER_INSTALL_DIRECTORY") +
                        "/fluidDatabase/" + fname;
    datFile.open(databaseFile, std::ios::in);
    if (datFile.fail()) {
      cerr << "ERROR: Error in fluid::GetDatabaseProperties(). File " << fname
           << " did not open correctly!!!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  string line = "";
  while (getline(datFile, line)) {
    // remove leading and trailing whitespace and ignore comments
    line = Trim(line);

    if (line.length() > 0) {  // only proceed if line has data
      // split line at variable separator
      auto tokens = Tokenize(line, ":", 2);
      // search to see if first token corresponds to any keywords
      auto key = tokens[0];

      if (key == "n") {
        n_ = std::stod(tokens[1]);
      } else if (key == "molarMass") {
        molarMass_ = std::stod(tokens[1]) / 1000.;  // convert to kg/mol
      } else if (key == "vibrationalTemperature") {
        vibTemp_ = ReadVectorXd(tokens[1]);
      } else if (key == "heatOfFormation") {
        heatOfFormation_ = std::stod(tokens[1]);
      } else if (key == "referencePressure") {
        refP_ = std::stod(tokens[1]);
      } else if (key == "referenceTemperature") {
        refT_ = std::stod(tokens[1]);
      } else if (key == "referenceEntropy") {
        refS_ = std::stod(tokens[1]);
      } else if (key == "sutherlandViscosityC1") {
        transportViscosity_[0] = std::stod(tokens[1]);
      } else if (key == "sutherlandViscosityS") {
        transportViscosity_[1] = std::stod(tokens[1]);
      } else if (key == "sutherlandConductivityC1") {
        transportConductivity_[0] = std::stod(tokens[1]);
      } else if (key == "sutherlandConductivityS") {
        transportConductivity_[1] = std::stod(tokens[1]);
      } else {
        cerr << "ERROR: Error in fluid::GetDatabaseProperties(). Property "
             << key << " is not recognized" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  // close database file
  datFile.close();
}

// function to read initial condition state from string
vector<fluid> ReadFluidList(ifstream &inFile, string &str) {
  vector<fluid> fluidList;
  auto openList = false;
  do {
    const auto start = openList ? 0 : str.find("<");
    const auto listOpened = str.find("<") == string::npos ? false : true;
    const auto end = str.find(">");
    openList = (end == string::npos) ? true : false;

    // test for fluid on current line
    // if < or > is alone on a line, should not look for fluid
    auto fluidPos = str.find("fluid");
    if (fluidPos != string::npos) {  // there is a fluid in current line
      string list;
      if (listOpened && openList) {  // list opened on current line, remains open
        list = str.substr(start + 1, string::npos);
      } else if (listOpened && !openList) {  // list opened/closed on current line
        const auto range = end - start - 1;
        list = str.substr(start + 1, range);  // +/- 1 to ignore <>
      } else if (!listOpened && openList) {  // list was open and remains open
        list = str.substr(start, string::npos);
      } else {  // list was open and is now closed
        const auto range = end - start;
        list = str.substr(start, range);
      }

      fluid fluidProps(list);
      fluidList.push_back(fluidProps);

      auto nextFluid = list.find("fluid");
      while (nextFluid != string::npos) {  // there are more fluids to read
        list.erase(0, nextFluid);  // remove commas separating states
        fluidProps = fluid(list);
        fluidList.push_back(fluidProps);
        nextFluid = list.find("fluid");
      }
    }

    if (openList) {
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);

  return fluidList;
}

ostream &operator<<(ostream &os, const fluid &fl) {
  auto vt = fl.VibrationalTemperature();
  os << "fluid(name=" << fl.Name()
     << "; referenceMassFraction=" << fl.MassFractionRef() << ")";
  return os;
}
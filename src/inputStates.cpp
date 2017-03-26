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
#include <vector>
#include <string>
#include "inputStates.hpp"

using std::cout;
using std::endl;
using std::cerr;

ostream &operator<<(ostream &os, const inputState &state) {
  state.Print(os);
  return os;
}

void icState::Print(ostream &os) const {
  os << "icState(tag=" << this->Tag() << "; pressure=" << this->Pressure()
     << "; density=" << this->Density() << "; velocity=[" << this->Velocity()
     << "]";
  if (this->SpecifiedTurbulence()) {
    os << "; turbulenceIntensity=" << this->TurbulenceIntensity()
       << "; eddyViscosityRatio=" << this->EddyViscosityRatio();
  }
  os << ")";
}

ostream &operator<<(ostream &os, const icState &state) {
  state.Print(os);
  return os;
}

void characteristic::Print(ostream &os) const {
  os << "characteristic(tag=" << this->Tag() << "; pressure=" << this->Pressure()
     << "; density=" << this->Density() << "; velocity=[" << this->Velocity()
     << "]";
  if (this->SpecifiedTurbulence()) {
    os << "; turbulenceIntensity=" << this->TurbulenceIntensity()
       << "; eddyViscosityRatio=" << this->EddyViscosityRatio();
  }
  os << ")";
}
ostream &operator<<(ostream &os, const characteristic &bc) {
  bc.Print(os);
  return os;
}

void stagnationInlet::Print(ostream &os) const {
  os << "stagnationInlet(tag=" << this->Tag() << "; p0="
     << this->StagnationPressure() << "; t0=" << this->StagnationTemperature()
     << "; direction=[" << this->Direction() << "]";
  if (this->SpecifiedTurbulence()) {
    os << "; turbulenceIntensity=" << this->TurbulenceIntensity()
       << "; eddyViscosityRatio=" << this->EddyViscosityRatio();
  }
  os << ")";
}

ostream &operator<<(ostream &os, const stagnationInlet &bc) {
  bc.Print(os);
  return os;
}

void pressureOutlet::Print(ostream &os) const {
  os << "pressureOutlet(tag=" << this->Tag() << "; pressure="
     << this->Pressure() << ")";
}

ostream &operator<<(ostream &os, const pressureOutlet &bc) {
  bc.Print(os);
  return os;
}

void supersonicInflow::Print(ostream &os) const {
  os << "supersonicInflow(tag=" << this->Tag() << "; pressure="
     << this->Pressure() << "; density=" << this->Density() << "; velocity=["
     << this->Velocity() << "]";
  if (this->SpecifiedTurbulence()) {
    os << "; turbulenceIntensity=" << this->TurbulenceIntensity()
       << "; eddyViscosityRatio=" << this->EddyViscosityRatio();
  }
  os << ")";
}

ostream &operator<<(ostream &os, const supersonicInflow &bc) {
  bc.Print(os);
  return os;
}

void subsonicOutflow::Print(ostream &os) const {
  os << "subsonicOutflow(tag=" << this->Tag() << "; pressure="
     << this->Pressure() << ")";
}

ostream &operator<<(ostream &os, const subsonicOutflow &bc) {
  bc.Print(os);
  return os;
}

void subsonicInflow::Print(ostream &os) const {
  os << "subsonicInflow(tag=" << this->Tag() << "; density=" << this->Density()
     << "; velocity=[" << this->Velocity() << "]";
  if (this->SpecifiedTurbulence()) {
    os << "; turbulenceIntensity=" << this->TurbulenceIntensity()
       << "; eddyViscosityRatio=" << this->EddyViscosityRatio();
  }
  os << ")";
}

ostream &operator<<(ostream &os, const subsonicInflow &bc) {
  bc.Print(os);
  return os;
}

void viscousWall::Print(ostream &os) const {
  os << "viscousWall(tag=" << this->Tag() << "; velocity=[" << this->Velocity()
     << "]";
  if (this->IsIsothermal()) {
    os << "; temperature=" << this->Temperature();
  } else if (specifiedHeatFlux_) {
    os << "; heatFlux=" << this->HeatFlux();
  }
  os << "; wallTreatment=" << wallTreatment_;
  if (this->IsWallLaw()) {
    os << "; vonKarmen=" << vonKarmen_ << "; wallConstant=" << wallConstant_;
  }
  os << ")";
}

ostream &operator<<(ostream &os, const viscousWall &bc) {
  bc.Print(os);
  return os;
}

void periodic::Print(ostream &os) const {
  os << "periodic(startTag=" << this->StartTag() << "; endTag=" << this->EndTag();
  if (this->IsTranslation()) {
    os << "; translation=[" << this->Translation() << "]";
  } else if (this->IsRotation()) {
    os << "; axis=[" << this->Axis() << "]; point=[" << this->Point()
       << "]; rotation=" << this->Rotation();
  }
  os << ")";
}

ostream &operator<<(ostream &os, const periodic &bc) {
  bc.Print(os);
  return os;
}


// function to trim leading and trailing whitespace from a string, and also
// remove data after a comment
string Trim(const string &s, const string &whitespace) {
  const string comment = "#";  // # is comment character for input file

  if (s.empty()) {
    return "";  // string is empty
  } else {
    // find index of first non whitespace character
    const auto sBegin = s.find_first_not_of(whitespace);
    if (sBegin == string::npos) {
      return "";  // string is only whitespace
    }

    // find index of last non whitespace character
    const auto sEnd = s.find_last_not_of(whitespace);
    const auto sRange = sEnd - sBegin + 1;  // range to trim string to
    auto temp = s.substr(sBegin, sRange);

    // find index of first comment character
    const auto tempComment = temp.find(comment);
    return temp.substr(0, tempComment);
  }
}

// function to tokenize a string based on a given character
vector<string> Tokenize(string str, const string &delimiter,
                        const unsigned int maxTokens) {
  // str -- string to tokenize
  // delimiter -- string to use as delimiter
  // maxTokens -- maximum number of tokens (if 0 (default), no max)

  vector<string> tokens;
  auto reachedMax = false;
  auto pos = str.find(delimiter);
  while (pos != string::npos && !reachedMax) {
    auto token = Trim(str.substr(0, pos));
    if (!token.empty()) {tokens.push_back(token);}
    // treat consecutive delimiters as single delimiter
    auto end = str.find_first_not_of(delimiter, pos);
    str.erase(0, end);
    if (maxTokens > 0 && maxTokens == tokens.size() - 1) {
      reachedMax = true;
    }
    pos = str.find(delimiter);
  }
  auto token = Trim(str);
  if (!token.empty()) {tokens.push_back(token);}
  return tokens;
}

// function to remove delimiter if it is last character
string RemoveTrailing(const string &str, const string &delimiter) {
  auto pos = str.rfind(delimiter);
  return (pos == str.length() - 1) ? str.substr(0, pos - 1) : str;
}

// function to read vector data from string
vector3d<double> ReadVector(const string &str) {
  const auto start = str.find("[") + 1;
  const auto end = str.find("]") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore []
  auto vec = str.substr(start, range);
  auto tokens = Tokenize(vec, ",");
  if (tokens.size() != 3) {
    cerr << "ERROR. Expected three components for vector, found "
         << tokens.size() << endl;
    cerr << "Vector string was " << vec << endl;
    exit(EXIT_FAILURE);
  }
  return {stod(tokens[0]), stod(tokens[1]), stod(tokens[2])};
}

// construct initial condition from string
icState::icState(string &str, const string name) {
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
  auto tagCount = 0;
  auto pressureCount = 0;
  auto densityCount = 0;
  auto velocityCount = 0;
  auto tiCount = 0;
  auto evrCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with " << name << " parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "pressure") {
      pressure_ = stod(RemoveTrailing(param[1], ","));
      pressureCount++;
    } else if (param[0] == "density") {
      density_ = stod(RemoveTrailing(param[1], ","));
      densityCount++;
    } else if (param[0] == "velocity") {
      velocity_ = ReadVector(RemoveTrailing(param[1], ","));
      velocityCount++;
    } else if (param[0] == "turbulenceIntensity") {
      this->SetSpecifiedTurbulence();
      turbIntensity_ = stod(RemoveTrailing(param[1], ","));
      tiCount++;
    } else if (param[0] == "eddyViscosityRatio") {
      eddyViscRatio_ = stod(RemoveTrailing(param[1], ","));
      evrCount++;
    } else if (param[0] == "tag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      tagCount++;
    } else {
      cerr << "ERROR. " << name << " specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (pressureCount != 1 || densityCount != 1 || velocityCount != 1) {
    cerr << "ERROR. For " << name << " pressure, density, and "
         << "velocity must be specified, and only specified once." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (tagCount > 1 || tiCount > 1 || evrCount > 1 || tiCount != evrCount) {
    cerr << "ERROR. For " << name << ", turbulenceIntensity, and "
         << "eddyViscosityRatio can only be specified once." << endl;
    cerr << "If either turbulenceIntensity or eddyViscosityRatio is specified "
         << "the other must be as well." << endl;
    exit(EXIT_FAILURE);
  }

  if (name != "icState" && tagCount != 1) {
    cerr << "ERROR. For " << name << ", tag must be specified." << endl;
    exit(EXIT_FAILURE);
  }
}

// construct stagnation inlet from string
stagnationInlet::stagnationInlet(string &str) {
  const auto start = str.find("(") + 1;
  const auto end = str.find(")") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore ()
  auto state = str.substr(start, range);
  const auto id = str.substr(0, start - 1);
  if (id != "stagnationInlet") {
    cerr << "ERROR. State condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // parameter counters
  auto tagCount = 0;
  auto p0Count = 0;
  auto t0Count = 0;
  auto directionCount = 0;
  auto tiCount = 0;
  auto evrCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with state condition parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "p0") {
      p0_ = stod(RemoveTrailing(param[1], ","));
      p0Count++;
    } else if (param[0] == "t0") {
      t0_ = stod(RemoveTrailing(param[1], ","));
      t0Count++;
    } else if (param[0] == "direction") {
      direction_ = ReadVector(RemoveTrailing(param[1], ",")).Normalize();
      directionCount++;
    } else if (param[0] == "turbulenceIntensity") {
      this->SetSpecifiedTurbulence();
      turbIntensity_ = stod(RemoveTrailing(param[1], ","));
      tiCount++;
    } else if (param[0] == "eddyViscosityRatio") {
      eddyViscRatio_ = stod(RemoveTrailing(param[1], ","));
      evrCount++;
    } else if (param[0] == "tag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      tagCount++;
    } else {
      cerr << "ERROR. stagnationInlet specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (p0Count != 1 || t0Count != 1 || directionCount != 1 || tagCount != 1) {
    cerr << "ERROR. For stagnationInlet p0, t0, tag, and "
         << "direction must be specified, and only specified once." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (tiCount > 1 || evrCount > 1 || tiCount != evrCount) {
    cerr << "ERROR. For stagnationInlet, turbulenceIntensity, and "
         << "eddyViscosityRatio can only be specified once." << endl;
    cerr << "If either turbulenceIntensity or eddyViscosityRatio is specified "
         << "the other must be as well." << endl;
    exit(EXIT_FAILURE);
  }
}

// construct pressureOutlet from string
pressureOutlet::pressureOutlet(string &str, const string name) {
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
  auto tagCount = 0;
  auto pressureCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with " << name << " parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "pressure") {
      pressure_ = stod(RemoveTrailing(param[1], ","));
      pressureCount++;
    } else if (param[0] == "tag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      tagCount++;
    } else {
      cerr << "ERROR. " << name << " specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (pressureCount != 1 || tagCount != 1) {
    cerr << "ERROR. For " << name << " pressure and tag must be specified, and "
         << "only specified once." << endl;
    exit(EXIT_FAILURE);
  }
}

// construct subsonic inflow from string
subsonicInflow::subsonicInflow(string &str) {
  const auto start = str.find("(") + 1;
  const auto end = str.find(")") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore ()
  auto state = str.substr(start, range);
  const auto id = str.substr(0, start - 1);
  if (id != "subsonicInflow") {
    cerr << "ERROR. State condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // parameter counters
  auto tagCount = 0;
  auto densityCount = 0;
  auto velocityCount = 0;
  auto tiCount = 0;
  auto evrCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with state condition parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "density") {
      density_ = stod(RemoveTrailing(param[1], ","));
      densityCount++;
    } else if (param[0] == "velocity") {
      velocity_ = ReadVector(RemoveTrailing(param[1], ","));
      velocityCount++;
    } else if (param[0] == "turbulenceIntensity") {
      this->SetSpecifiedTurbulence();
      turbIntensity_ = stod(RemoveTrailing(param[1], ","));
      tiCount++;
    } else if (param[0] == "eddyViscosityRatio") {
      eddyViscRatio_ = stod(RemoveTrailing(param[1], ","));
      evrCount++;
    } else if (param[0] == "tag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      tagCount++;
    } else {
      cerr << "ERROR. subsonicInlet specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (densityCount != 1 || velocityCount != 1 || tagCount != 1) {
    cerr << "ERROR. For subsonicInlet density, tag, and velocity "
         << "must be specified, and only specified once." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (tiCount > 1 || evrCount > 1 || tiCount != evrCount) {
    cerr << "ERROR. For subsonicInlet, turbulenceIntensity, and "
         << "eddyViscosityRatio can only be specified once." << endl;
    cerr << "If either turbulenceIntensity or eddyViscosityRatio is specified "
         << "the other must be as well." << endl;
    exit(EXIT_FAILURE);
  }
}

// construct viscousWall from string
viscousWall::viscousWall(string &str) {
  const auto start = str.find("(") + 1;
  const auto end = str.find(")") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore ()
  auto state = str.substr(start, range);
  const auto id = str.substr(0, start - 1);
  if (id != "viscousWall") {
    cerr << "ERROR. State condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // parameter counters
  auto tagCount = 0;
  auto velocityCount = 0;
  auto temperatureCount = 0;
  auto heatFluxCount = 0;
  auto wallTreatmentCount = 0;
  auto vonKarmenCount = 0;
  auto wallConstantCount = 0;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with state condition parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "velocity") {
      velocity_ = ReadVector(RemoveTrailing(param[1], ","));
      velocityCount++;
    } else if (param[0] == "temperature") {
      specifiedTemperature_ = true;
      temperature_ = stod(RemoveTrailing(param[1], ","));
      temperatureCount++;
    } else if (param[0] == "heatFlux") {
      specifiedHeatFlux_ = true;
      heatFlux_ = stod(RemoveTrailing(param[1], ","));
      heatFluxCount++;
    } else if (param[0] == "vonKarmen") {
      vonKarmen_ = stod(RemoveTrailing(param[1], ","));
      vonKarmenCount++;
    } else if (param[0] == "wallConstant") {
      wallConstant_ = stod(RemoveTrailing(param[1], ","));
      wallConstantCount++;
    } else if (param[0] == "wallTreatment") {
      wallTreatment_ = RemoveTrailing(param[1], ",");
      wallTreatmentCount++;
    } else if (param[0] == "tag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      tagCount++;
    } else {
      cerr << "ERROR. viscousWall specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (tagCount != 1) {
    cerr << "ERROR. For viscousWall tag must be specified, and only specified "
         << "once." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (velocityCount > 1 || temperatureCount > 1 || heatFluxCount > 1 ||
      wallTreatmentCount > 1 || vonKarmenCount > 1 || wallConstantCount > 1) {
    cerr << "ERROR. For viscousWall, velocity, heatFlux, temperature, "
         << "wallTreatment, vonKarmen, and wallConstant can only be specified "
         << "once." << endl;
    exit(EXIT_FAILURE);
  }
  if (specifiedHeatFlux_ && specifiedTemperature_) {
    cerr << "ERROR. For viscousWall can only specify temperature OR heatFlux."
         << endl;
    exit(EXIT_FAILURE);
  }
  if (wallTreatment_ != "lowRe" && wallTreatment_ != "wallLaw") {
    cerr << "ERROR. wallTreatment " << wallTreatment_ << " is not recognized!"
         << endl;
    cerr << "Choose lowRe or wallLaw" << endl;
    exit(EXIT_FAILURE);
  }
}

// construct periodic from string
periodic::periodic(string &str) {
  const auto start = str.find("(") + 1;
  const auto end = str.find(")") - 1;
  const auto range = end - start + 1;  // +/-1 to ignore ()
  auto state = str.substr(start, range);
  const auto id = str.substr(0, start - 1);
  if (id != "periodic") {
    cerr << "ERROR. State condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // parameter counters
  auto startTagCount = 0;
  auto endTagCount = 0;
  auto translationCount = 0;
  auto axisCount = 0;
  auto pointCount = 0;
  auto rotationCount = 0;
  auto specifiedTranslation = false;
  auto specifiedAxis = false;
  auto specifiedPoint = false;
  auto specifiedRotation = false;

  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with state condition parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "translation") {
      translation_ = ReadVector(RemoveTrailing(param[1], ","));
      specifiedTranslation = true;
      translationCount++;
    } else if (param[0] == "axis") {
      axis_ = ReadVector(RemoveTrailing(param[1], ",")).Normalize();
      specifiedAxis = true;
      axisCount++;
    } else if (param[0] == "point") {
      point_ = ReadVector(RemoveTrailing(param[1], ","));
      specifiedPoint = true;
      pointCount++;
    } else if (param[0] == "rotation") {
      rotation_ = stod(RemoveTrailing(param[1], ","));
      specifiedRotation = true;
      rotationCount++;
    } else if (param[0] == "startTag") {
      this->SetTag(stoi(RemoveTrailing(param[1], ",")));
      startTagCount++;
    } else if (param[0] == "endTag") {
      endTag_ = stoi(RemoveTrailing(param[1], ","));
      endTagCount++;
    } else {
      cerr << "ERROR. periodic specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (startTagCount != 1 || endTagCount != 1) {
    cerr << "ERROR. For periodic startTag and endTag must be specified, and "
         << "only specified once." << endl;
    exit(EXIT_FAILURE);
  }
  if (this->StartTag() == this->EndTag()) {
    cerr << "ERROR. For periodic startTag and endTag must be unique." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (translationCount > 1 || axisCount > 1 || rotationCount > 1 ||
      pointCount > 1) {
    cerr << "ERROR. For periodic, translation, axis, rotation, and point can "
         << "only be specified once." << endl;
    exit(EXIT_FAILURE);
  }
  if (specifiedTranslation &&
      (specifiedAxis || specifiedRotation || specifiedPoint)) {
    cerr << "ERROR. For periodic can only specify translation OR rotation."
         << endl;
    exit(EXIT_FAILURE);
  }
  if ((specifiedAxis && !(specifiedRotation && specifiedPoint)) ||
      (specifiedPoint && !(specifiedRotation && specifiedAxis)) ||
      (specifiedRotation && !(specifiedAxis && specifiedPoint))) {
    cerr << "ERROR. For periodic rotation must specify axis, point, & rotation."
         << endl;
    exit(EXIT_FAILURE);
  }
}


void CheckICTags(const vector<icState> &ics, const int &tag) {
  for (auto &ic : ics) {
    if (tag == ic.Tag()) {
      cerr << "ERROR: Initial condition tag " << tag << " is repeated!" << endl;
      exit(EXIT_FAILURE);
    }
  }
}

// function to read initial condition state from string
vector<icState> ReadICList(ifstream &inFile, string &str) {
  vector<icState> icList;
  auto openList = false;
  do {
    const auto start = openList ? 0 : str.find("<");
    const auto listOpened = str.find("<") == string::npos ? false : true;
    const auto end = str.find(">");
    openList = (end == string::npos) ? true : false;

    // test for state on current line
    // if < or > is alone on a line, should not look for icState
    auto statePos = str.find("icState");
    if (statePos != string::npos) {  // there is a state in current line
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

      icState ic(list);
      CheckICTags(icList, ic.Tag());
      icList.push_back(ic);

      auto nextState = list.find("icState");
      while (nextState != string::npos) {  // there are more states to read
        list.erase(0, nextState);  // remove commas separating states
        ic = icState(list);
        CheckICTags(icList, ic.Tag());
        icList.push_back(ic);
        nextState = list.find("icState");
      }
    }

    if (openList) {
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);

  return icList;
}

auto FindBCPosition(const string &str, const vector<string> &names,
                    string &type) {
  auto pos = string::npos;
  type = "none";
  for (auto &name : names) {
    auto currPos = str.find(name);
    if (currPos < pos) {
      pos = currPos;
      type = name;
    }
  }
  return pos;
}

void AddBCToList(const string &type, vector<unique_ptr<inputState>> &bcList,
                 string &list) {
  unique_ptr<inputState> bc(nullptr);
  if (type == "characteristic") {
    bc = unique_ptr<inputState>{std::make_unique<characteristic>(list)};
  } else if (type == "stagnationInlet") {
    bc = unique_ptr<inputState>{std::make_unique<stagnationInlet>(list)};
  } else if (type == "pressureOutlet") {
    bc = unique_ptr<inputState>{std::make_unique<pressureOutlet>(list)};
  } else if (type == "subsonicInflow") {
    bc = unique_ptr<inputState>{std::make_unique<subsonicInflow>(list)};
  } else if (type == "subsonicOutflow") {
    bc = unique_ptr<inputState>{std::make_unique<subsonicOutflow>(list)};
  } else if (type == "supersonicInflow") {
    bc = unique_ptr<inputState>{std::make_unique<supersonicInflow>(list)};
  } else if (type == "viscousWall") {
    bc = unique_ptr<inputState>{std::make_unique<viscousWall>(list)};
  } else if (type == "periodic") {
    bc = unique_ptr<inputState>{std::make_unique<periodic>(list)};
  } else {
    cerr << "ERROR. BC state " << type << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  // sanity check -- see if tag already exits
  for (auto &bcData : bcList) {
    if (bc->Tag() == bcData->Tag()) {
      cerr << "ERROR: Boundary state tag " << bc->Tag() << " is repeated!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  bcList.push_back(std::move(bc));
}

// function to read boundary condition list from string
vector<unique_ptr<inputState>> ReadBCList(ifstream &inFile, string &str) {
  vector<unique_ptr<inputState>> bcList;
  vector<string> bcNames {"characteristic", "stagnationInlet", "pressureOutlet",
        "subsonicInflow", "subsonicOutflow", "supersonicInflow", "viscousWall",
        "periodic"};
  auto openList = false;
  do {
    auto start = openList ? 0 : str.find("<");
    auto listOpened = str.find("<") == string::npos ? false : true;
    auto end = str.find(">");
    openList = (end == string::npos) ? true : false;

    // test for state on current line
    // if < or > is alone on a line, should not look for bc
    string type = "none";
    auto bcPos = FindBCPosition(str, bcNames, type);
    if (bcPos != string::npos) {  // there is a state in current line
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

      AddBCToList(type, bcList, list);

      auto nextState = FindBCPosition(list, bcNames, type);
      while (nextState != string::npos) {  // there are more states to read
        list.erase(0, nextState);  // remove commas separating states
        AddBCToList(type, bcList, list);
        nextState = FindBCPosition(list, bcNames, type);
      }
    }

    if (openList) {
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);

  return bcList;
}

// function to read a list of strings
// used for output variable specification
vector<string> ReadStringList(ifstream &inFile, string &str) {
  vector<string> strList;
  auto openList = false;
  do {
    auto start = openList ? 0 : str.find("<");
    auto listOpened = str.find("<") == string::npos ? false : true;
    auto end = str.find(">");
    openList = (end == string::npos) ? true : false;

    // test for argument on current line
    auto argPos = str.find_first_not_of(" \t\r\n,");
    if (argPos != string::npos) {  // there is an argument in current line
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

      // tokenize all arguments on current line and add to vector
      auto args = Tokenize(list, ",");
      for (auto &arg : args) {
        strList.push_back(arg);
      }
    }

    if (openList) {
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);

  return strList;
}


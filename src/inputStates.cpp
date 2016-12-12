/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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
     << "]; turbulenceIntensity=" << this->TurbulenceIntensity()
     << "; eddyViscosityRatio=" << this->EddyViscosityRatio() << ")";
}

ostream &operator<<(ostream &os, const icState &state) {
  state.Print(os);
  return os;
}

void characteristic::Print(ostream &os) const {
  os << "characteristic(tag=" << this->Tag() << "; pressure=" << this->Pressure()
     << "; density=" << this->Density() << "; velocity=[" << this->Velocity()
     << "]; turbulenceIntensity=" << this->TurbulenceIntensity()
     << "; eddyViscosityRatio=" << this->EddyViscosityRatio() << ")";
}
ostream &operator<<(ostream &os, const characteristic &bc) {
  bc.Print(os);
  return os;
}

void stagnationInlet::Print(ostream &os) const {
  os << "stagnationInlet(tag=" << this->Tag() << "; p0="
     << this->StagnationPressure() << "; t0=" << this->StagnationTemperature()
     << "; direction=[" << this->Direction() << "]; turbulenceIntensity="
     << this->TurbulenceIntensity() << "; eddyViscosityRatio="
     << this->EddyViscosityRatio() << ")";
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
     << this->Velocity() << "]; turbulenceIntensity="
     << this->TurbulenceIntensity() << "; eddyViscosityRatio="
     << this->EddyViscosityRatio() << ")";
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
     << "; velocity=[" << this->Velocity() << "]; turbulenceIntensity="
     << this->TurbulenceIntensity() << "; eddyViscosityRatio="
     << this->EddyViscosityRatio() << ")";
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
  }
  os << ")";
}

ostream &operator<<(ostream &os, const viscousWall &bc) {
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
    // find index of last non whitespace character
    const auto sEnd = s.find_last_not_of(whitespace);
    const auto sRange = sEnd - sBegin + 1;  // range to trim string to
    auto temp = s.substr(sBegin, sRange);

    // find index of first comment character
    const auto tempComment = temp.find(comment);
    const auto tempRange = tempComment - 0;  // find range of string to return

    return temp.substr(0, tempRange);
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
    auto token = str.substr(0, pos);
    tokens.push_back(Trim(token));
    // treat consecutive delimiters as single delimiter
    auto end = str.find_first_not_of(delimiter, pos);
    str.erase(0, end);
    if (maxTokens > 0 && maxTokens == tokens.size() - 1) {
      reachedMax = true;
    }
    pos = str.find(delimiter);
  }
  tokens.push_back(Trim(str));
  return tokens;
}

// function to remove delimiter if it is last character
string RemoveTrailing(const string &str, const string &delimiter) {
  auto pos = str.rfind(delimiter);
  return (pos == str.length() - 1) ? str.substr(0, pos - 1) : str;
}

// function to read vector data from string
vector3d<double> ReadVector(const string &str) {
  auto start = str.find("[");
  auto end = str.find("]");
  auto vec = str.substr(start + 1, end - 1);  // +/-1 to ignore []
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
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
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

  // set default values for optional parameters
  turbIntensity_ = 0.0;
  eddyViscRatio_ = 0.0;

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
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
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

  // set default values for optional parameters
  turbIntensity_ = 0.0;
  eddyViscRatio_ = 0.0;

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
      direction_ = ReadVector(RemoveTrailing(param[1], ","));
      directionCount++;
    } else if (param[0] == "turbulenceIntensity") {
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
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
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
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
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

  // set default values for optional parameters
  turbIntensity_ = 0.0;
  eddyViscRatio_ = 0.0;

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
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
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

  // set default values for optional parameters
  velocity_ = {0.0, 0.0, 0.0};
  temperature_ = 0.0;

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
      temperature_ = stoi(RemoveTrailing(param[1], ","));
      temperatureCount++;
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
  if (velocityCount > 1 || temperatureCount > 1) {
    cerr << "ERROR. For viscousWall, velocity and temperature can only be "
         << "specified once." << endl;
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
    auto start = openList ? 0 : str.find("<");
    auto listOpened = str.find("<") == string::npos ? false : true;
    auto end = str.find(">");
    openList = (end == string::npos) ? true : false;

    // test for state on current line
    // if < or > is alone on a line, should not look for icState
    auto statePos = str.find("icState");
    if (statePos != string::npos) {  // there is a state in current line
      string list;
      if (listOpened && openList) {  // list opened on current line, remains open
        list = str.substr(start + 1, end);
      } else if (listOpened && !openList) {  // list opened/closed on current line
        list = str.substr(start + 1, end - 1);  // +/- 1 to ignore <>
      } else if (!listOpened && openList) {  // list was open and remains open
        list = str.substr(start, end);
      } else {  // list was open and is now closed
        list = str.substr(start, end - 1);
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
        "subsonicInflow", "subsonicOutflow", "supersonicInflow", "viscousWall"};
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
        list = str.substr(start + 1, end);
      } else if (listOpened && !openList) {  // list opened/closed on current line
        list = str.substr(start + 1, end - 1);  // +/- 1 to ignore <>
      } else if (!listOpened && openList) {  // list was open and remains open
        list = str.substr(start, end);
      } else {  // list was open and is now closed
        list = str.substr(start, end - 1);
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
    // if < or > is alone on a line, should not look for icState
    auto argPos = str.find_first_not_of(" \t,");
    if (argPos != string::npos) {  // there is an argument in current line
      string list;
      if (listOpened && openList) {  // list opened on current line, remains open
        list = str.substr(start + 1, end);
      } else if (listOpened && !openList) {  // list opened/closed on current line
        list = str.substr(start + 1, end - 1);  // +/- 1 to ignore <>
      } else if (!listOpened && openList) {  // list was open and remains open
        list = str.substr(start, end);
      } else {  // list was open and is now closed
        list = str.substr(start, end - 1);
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


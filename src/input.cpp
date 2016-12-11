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

#include <chrono>     // timing capability
#include <iostream>   // cout
#include <iomanip>    // put_time
#include <fstream>    // ifstream
#include <cstdlib>    // exit()
#include <sstream>    // istringstream
#include <iterator>   // istring_iterator
#include <memory>     // make_unique
#include <algorithm>  // min
#include <string>
#include <vector>
#include "input.hpp"
#include "turbulence.hpp"
#include "inputStates.hpp"
#include "macros.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::cerr;
using std::istringstream;
using std::istream_iterator;
using std::unique_ptr;
using std::stoi;
using std::stod;

// constructor for input class
// initialize vector to have length of number of acceptable inputs to the code
input::input(const string &name) : simName_(name) {
  // default values for each variable
  gName_ = "";
  dt_ = -1.0;
  iterations_ = 1;
  pRef_ = -1.0;
  rRef_ = -1.0;
  lRef_ = 1.0;
  gamma_ = 1.4;
  gasConst_ = 287.058;
  velRef_ = {1.0, 0.0, 0.0};
  bc_ = vector<boundaryConditions>(1);
  timeIntegration_ = "explicitEuler";
  cfl_ = -1.0;
  kappa_ = -2.0;  // default to value outside of range to tell if higher
                  // order or constant method should be used
  limiter_ = "none";
  outputFrequency_ = 1;
  equationSet_ = "euler";
  tRef_ = -1.0;
  matrixSolver_ = "lusgs";
  matrixSweeps_ = 1;
  matrixRelaxation_ = 1.0;  // default is symmetric Gauss-Seidel
                            // with no overrelaxation
  timeIntTheta_ = 1.0;  // default results in implicit euler
  timeIntZeta_ = 0.0;  // default results in implicit euler
  nonlinearIterations_ = 1;  // default is 1 (steady)
  cflMax_ = 1.0;
  cflStep_ = 0.0;
  cflStart_ = 1.0;
  invFluxJac_ = "rusanov";  // default is approximate rusanov which is used
                            // with lusgs
  dualTimeCFL_ = -1.0;  // default value of -1; negative value means dual time
                       // stepping is not used
  inviscidFlux_ = "roe";  // default value is roe flux
  stagInletProps_[0] = 0.0;
  stagInletProps_[1] = 0.0;
  stagInletProps_[2] = 0.0;
  stagInletProps_[3] = 0.0;
  stagInletProps_[4] = 0.0;
  stagInletProps_[5] = 0.0;
  pressureOutlet_[0] = 0.0;
  pressureOutlet_[1] = 0.0;
  decompMethod_ = "cubic";  // default is cubic decomposition
  turbModel_ = "none";  // default turbulence model is none
  farfieldTurbInten_ = 0.01;
  farfieldEddyViscRatio_ = 10.0;

  // default to primative variables
  outputVariables_ = {"density", "vel_x", "vel_y", "vel_z", "pressure"};

  // keywords in the input file that the parser is looking for to define
  // variables
  vars_ = {"gridName",
           "timeStep",
           "iterations",
           "pressureRef",
           "densityRef",
           "lengthRef",
           "gamma",
           "gasConstant",
           "velocity",
           "timeIntegration",
           "faceReconstruction",
           "limiter",
           "outputFrequency",
           "equationSet",
           "temperatureRef",
           "matrixSolver",
           "matrixSweeps",
           "matrixRelaxation",
           "nonlinearIterations",
           "cflMax",
           "cflStep",
           "cflStart",
           "inviscidFluxJacobian",
           "dualTimeCFL",
           "inviscidFlux",
           "stagnationInlet",
           "pressureOutlet",
           "decompositionMethod",
           "turbulenceModel",
           "farfieldTurbulenceIntensity",
           "farfieldEddyViscosityRatio",
           "outputVariables",
           "initialConditions",
           "boundaryConditions"};
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

// function to remove delimiter if it is last character
string RemoveTrailing(const string &str, const string &delimiter) {
  auto pos = str.rfind(delimiter);
  return (pos == str.length() - 1) ? str.substr(0, pos - 1) : str;
}

// function to read initial condition state from string
icState ReadICState(string &str) {
  auto start = str.find("(");
  auto end = str.find(")");
  auto state = str.substr(start + 1, end - 1);  // +/-1 to ignore ()
  auto id = str.substr(0, start);
  if (id != "icState") {
    cerr << "ERROR. Initial condition specifier " << id << " is not recognized!"
         << endl;
    exit(EXIT_FAILURE);
  }
  auto tokens = Tokenize(state, ";");

  // erase portion used so multiple states in same string can easily be found
  str.erase(0, end);

  // paramter counters
  auto tagCount = 0;
  auto pressureCount = 0;
  auto densityCount = 0;
  auto velocityCount = 0;
  auto tiCount = 0;
  auto evrCount = 0;

  auto tag = 0;
  auto pressure = 0.0;
  auto density = 0.0;
  vector3d<double> velocity = {0.0, 0.0, 0.0};
  auto turbIntensity = 0.0;
  auto eddyViscRatio = 0.0;
  for (auto &token : tokens) {
    auto param = Tokenize(token, "=");
    if (param.size() != 2) {
      cerr << "ERROR. Problem with initial condition parameter " << token << endl;
      exit(EXIT_FAILURE);
    }

    if (param[0] == "pressure") {
      pressure = stod(RemoveTrailing(param[1], ","));
      pressureCount++;
    } else if (param[0] == "density") {
      density = stod(RemoveTrailing(param[1], ","));
      densityCount++;
    } else if (param[0] == "velocity") {
      velocity = ReadVector(RemoveTrailing(param[1], ","));
      velocityCount++;
    } else if (param[0] == "turbulenceIntensity") {
      turbIntensity = stod(RemoveTrailing(param[1], ","));
      tiCount++;
    } else if (param[0] == "eddyViscosityRatio") {
      eddyViscRatio = stod(RemoveTrailing(param[1], ","));
      evrCount++;
    } else if (param[0] == "tag") {
      tag = stoi(RemoveTrailing(param[1], ","));
      tagCount++;
    } else {
      cerr << "ERROR. Initial condition specifier " << param[0]
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // sanity checks
  // required variables
  if (pressureCount != 1 || densityCount != 1 || velocityCount != 1) {
    cerr << "ERROR. For initial condition state pressure, density, and "
         << "velocity must be specified, and only specified once." << endl;
    exit(EXIT_FAILURE);
  }
  // optional variables
  if (tagCount > 1 || tiCount > 1 || evrCount > 1 || tiCount != evrCount) {
    cerr << "ERROR. For initial condition state tag, turbulenceIntensity, and "
         << "eddyViscosityRatio can only be specified once." << endl;
    cerr << "If either turbulenceIntensity or eddyViscosityRatio is specified "
         << "the other must be as well." << endl;
    exit(EXIT_FAILURE);
  }

  if (tagCount == 1) {
    if (tiCount == 1) {
      return icState(velocity, density, pressure, turbIntensity, eddyViscRatio, tag);
    } else {
      return icState(velocity, density, pressure, tag);
    }
  } else {
    if (tiCount == 1) {
      return icState(velocity, density, pressure, turbIntensity, eddyViscRatio);
    } else {
      return icState(velocity, density, pressure);
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

      auto ic = ReadICState(list);
      icList.push_back(ic);

      auto nextState = list.find("icState");
      while (nextState != string::npos) {  // there are more states to read
        list.erase(0, nextState);  // remove commas separating states
        ic = ReadICState(list);
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

// function to print the time
void PrintTime() {
  auto now = std::chrono::system_clock::now();
  auto nowOut = std::chrono::system_clock::to_time_t(now);
  cout << std::put_time(std::localtime(&nowOut), "%c") << endl;
}

// function to read the input file and return the data as a member of the input
// class
// read the input file
void input::ReadInput(const int &rank) {
  // rank -- rank of processor

  if (rank == ROOTP) {
    cout << "##################################################################"
            "#########################################################" << endl;
    PrintTime();
    cout << endl;
    cout << "Parsing input file " << simName_ << endl << endl;
    cout << "Solver Inputs" << endl;
  }

  // open input file
  ifstream inFile(simName_, ios::in);
  if (inFile.fail()) {
    cerr << "ERROR: Error in input::ReadInput(). Input file " << simName_
         << " did not open correctly!!!" << endl;
    exit(EXIT_FAILURE);
  }

  string line = "";
  auto readingBCs = 0;
  auto blk = 0;
  auto surfCounter = 0;
  auto numBCBlks = 0;
  vector<boundaryConditions> tempBC(1);
  auto lEnd = 0;
  auto numSurf = 0;

  while (getline(inFile, line)) {  // while there are still lines in the input
    // file, execute loop

    // remove leading and trailing whitespace and ignore comments
    line = Trim(line);

    if (line.length() > 0) {  // only proceed if line has data
      // split line at variable separator
      auto tokens = Tokenize(line, ":", 2);

      // search to see if first token corresponds to any keywords
      auto key = tokens[0];

      // if first token matches a keyword or reading boundary condtions
      if (vars_.find(key) != vars_.end() || readingBCs > 0) {
        // if not yet reading BCs (readingBCs == 0), set variable in input
        // class to corresponding value and print assignment to std out
        if (key == "gridName") {
          gName_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->GridName() << endl;
          }
        } else if (key == "timeStep") {
          dt_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->Dt() << endl;
          }
        } else if (key == "iterations") {
          iterations_ = stoi(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->Iterations() << endl;
          }
        } else if (key == "pressureRef") {
          pRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->PRef() << endl;
          }
        } else if (key == "densityRef") {
          rRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->RRef() << endl;
          }
        } else if (key == "lengthRef") {
          lRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->LRef() << endl;
          }
        } else if (key == "gamma") {
          gamma_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->Gamma() << endl;
          }
        } else if (key == "gasConstant") {
          gasConst_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->R() << endl;
          }
        } else if (key == "velocity") {
          velRef_ = ReadVector(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->VelRef() << endl;
          }
        } else if (key == "timeIntegration") {
          timeIntegration_ = tokens[1];
          if (this->TimeIntegration() == "implicitEuler") {
            timeIntTheta_ = 1.0;
            timeIntZeta_ = 0.0;
          } else if (this->TimeIntegration() == "crankNicholson") {
            timeIntTheta_ = 0.5;
            timeIntZeta_ = 0.0;
          } else if (this->TimeIntegration() == "bdf2") {
            timeIntTheta_ = 1.0;
            timeIntZeta_ = 0.5;
          }

          if (rank == ROOTP) {
            cout << key << ": " << this->TimeIntegration() << endl;
          }
        } else if (key == "faceReconstruction") {
          if (tokens[1] == "upwind") {
            kappa_ = -1.0;
          } else if (tokens[1] == "fromm") {
            kappa_ = 0.0;
          } else if (tokens[1] == "quick") {
            kappa_ = 0.5;
          } else if (tokens[1] == "central") {
            kappa_ = 1.0;
          } else if (tokens[1] == "thirdOrder") {
            kappa_ = 1.0 / 3.0;
          } else if (tokens[1] == "constant") {
            kappa_ = -2.0;
          } else {
            // if string is not recognized, set kappa to number input
            kappa_ = stod(tokens[1]);
          }

          if (rank == ROOTP) {
            cout << key << ": " << tokens[1]
                 << " kappa = " << this->Kappa() << endl;
          }
          if ((this->Kappa() < -1.0) || (this->Kappa() > 1.0)) {
            cerr << "ERROR: Error in input::ReadInput(). Kappa value of "
                 << this->Kappa()
                 << " is not valid! Choose a value between -1.0 and 1.0."
                 << endl;
            exit(EXIT_FAILURE);
          }
        } else if (key == "limiter") {
          limiter_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->Limiter() << endl;
          }
        } else if (key == "outputFrequency") {
          outputFrequency_ = stoi(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->OutputFrequency() << endl;
          }
        } else if (key == "equationSet") {
          equationSet_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->EquationSet() << endl;
          }
        } else if (key == "temperatureRef") {
          tRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->TRef() << endl;
          }
        } else if (key == "matrixSolver") {
          matrixSolver_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->MatrixSolver() << endl;
          }
        } else if (key == "matrixSweeps") {
          matrixSweeps_ = stoi(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->MatrixSweeps() << endl;
          }
        } else if (key == "matrixRelaxation") {
          matrixRelaxation_ =
              stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->MatrixRelaxation() << endl;
          }
        } else if (key == "nonlinearIterations") {
          nonlinearIterations_ = stoi(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->NonlinearIterations() << endl;
          }
        } else if (key == "cflMax") {
          cflMax_ = stod(tokens[1]);  // double  (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->CFLMax() << endl;
          }
        } else if (key == "cflStep") {
          cflStep_ = stod(tokens[1]);  // double (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->CFLStep() << endl;
          }
        } else if (key == "cflStart") {
          cflStart_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->CFLStart() << endl;
          }
        } else if (key == "inviscidFluxJacobian") {
          invFluxJac_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->InvFluxJac() << endl;
          }
        } else if (key == "dualTimeCFL") {
          dualTimeCFL_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->DualTimeCFL() << endl;
          }
        } else if (key == "inviscidFlux") {
          inviscidFlux_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->InviscidFlux() << endl;
          }
        } else if (key == "stagnationInlet") {
          stagInletProps_[0] = stoi(tokens[1]);  // tag
          stagInletProps_[1] = stod(tokens[2]);  // stag pressure
          stagInletProps_[2] = stod(tokens[3]);  // stag temp
          stagInletProps_[3] = stod(tokens[4]);  // dir-x
          stagInletProps_[4] = stod(tokens[5]);  // dir-y
          stagInletProps_[5] = stod(tokens[6]);  // dir-z
          if (rank == ROOTP) {
            cout << key << ": " << this->StagInletTag() << " "
                 << this->StagInletP0() << " " << this->StagInletT0()
                 << " " << this->StagInletDx() << " "
                 << this->StagInletDy() << " " << this->StagInletDz() << endl;
          }
        } else if (key == "pressureOutlet") {
          pressureOutlet_[0] = stoi(tokens[1]);  // tag
          pressureOutlet_[1] = stod(tokens[2]);  // outlet pressure
          if (rank == ROOTP) {
            cout << key << ": " << this->PressureOutletTag()
                 << " " << this->PressureOutletP() << endl;
          }
        } else if (key == "decompositionMethod") {
          decompMethod_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->DecompMethod() << endl;
          }
        } else if (key == "turbulenceModel") {
          turbModel_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->TurbulenceModel() << endl;
          }
        } else if (key == "farfieldTurbulenceIntensity") {
          farfieldTurbInten_ = stod(tokens[1]);  // double (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->FarfieldTurbIntensity() << endl;
          }
        } else if (key == "farfieldEddyViscosityRatio") {
          farfieldEddyViscRatio_ = stod(tokens[1]);  // double (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->FarfieldEddyViscRatio() << endl;
          }
        } else if (key == "outputVariables") {
          // clear default variables from set
          outputVariables_.clear();
          auto specifiedVars = Tokenize(tokens[1], ",");
          for (auto &vars : specifiedVars) {
            outputVariables_.insert(vars);
          }
          if (rank == ROOTP) {
            cout << key << ": ";
            for (auto &vars : outputVariables_) {
              cout << vars << " ";
            }
            cout << endl;
          }

        } else if (key == "initialConditions") {
          ics_ = ReadICList(inFile, tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": <";
            for (auto ii = 0U; ii < ics_.size(); ++ii) {
              cout << ics_[ii];
              if (ii == ics_.size() - 1) {
                cout << ">" << endl;
              } else {
                cout << "," << endl << "                    ";
              }
            }
          }

          // reading BCs
          // -------------------------------------------------------------
        } else if (key == "boundaryConditions" || readingBCs > 0) {
          // read in boundary conditions and assign to boundaryConditions class
          if (readingBCs == 0) {  // variable read must be number of blocks if
            // first line of BCs
            numBCBlks = stoi(tokens[1]);
            tempBC.resize(numBCBlks);  // resize vector holding BCs to correct
            // number of blocks
            readingBCs++;
            lEnd++;
          } else if (readingBCs == lEnd) {  // variables must be number of i,
            // j, k surfaces for each block
            // set number of i, j, k surfaces and resize vectors
            // boundary conditions are space delimited
            tokens = Tokenize(line, " ");
            tempBC[blk].ResizeVecs(stoi(tokens[0]), stoi(tokens[1]),
                                   stoi(tokens[2]));

            // set total number of surfaces
            numSurf = stoi(tokens[0]) + stoi(tokens[1]) + stoi(tokens[2]);

            // set new target for when block data is finished
            lEnd += (numSurf + 1);
            readingBCs++;

          } else {  // assign BC block variables
            // boundary conditions are space delimited
            tokens = Tokenize(line, " ");
            tempBC[blk].AssignFromInput(surfCounter, tokens);

            surfCounter++;
            if (surfCounter == numSurf) {  // at end of block data, increment
              // block index, reset surface
              // counter
              blk++;
              surfCounter = 0;
            }

            readingBCs++;
          }

          // if block counter reaches number of blocks, BCs are finished (b/c
          // counter starts at 0), so assign BCs and write them out
          if (blk == numBCBlks) {
            bc_ = tempBC;
            if (rank == ROOTP) {
              cout << "boundaryConditions: " << this->NumBC() << endl;
              for (auto ll = 0; ll < this->NumBC(); ll++) {
                cout << "Block: " << ll << endl;
                cout << this->BC(ll) << endl;
              }
            }
          }
        }
      }
    }
  }


  // input file sanity checks
  this->CheckNonlinearIterations();
  this->CheckOutputVariables();
  this->CheckTurbulenceModel();

  if (rank == ROOTP) {
    cout << endl;
    cout << "Input file parse complete" << endl;
    cout << "##################################################################"
            "#########################################################" << endl
         << endl;
  }
}

// member function to calculate the cfl value for the step from the starting,
// ending, and step values
void input::CalcCFL(const int &ii) {
  cfl_ = std::min(cflStart_ + ii * cflStep_, cflMax_);
}

// member function to determine number of turbulence equations
int input::NumTurbEquations() const {
  auto numEqns = 0;
  if (this->IsTurbulent()) {
    numEqns = 2;
  }
  return numEqns;
}

// member function to determine number of equations to solver for
int input::NumEquations() const {
  auto numEqns = 0;
  if ((equationSet_ == "euler") ||
      (equationSet_ == "navierStokes")) {
    numEqns = this->NumFlowEquations();
  } else if (equationSet_ == "rans") {
    numEqns = this->NumFlowEquations() + this->NumTurbEquations();
  } else {
    cerr << "ERROR: Equations set is not recognized. Cannot determine number "
            "of equations!" << endl;
  }

  return numEqns;
}

// member function to determine of method is implicit or explicit
bool input::IsImplicit() const {
  if (timeIntegration_ == "implicitEuler" ||
      timeIntegration_ == "crankNicholson" ||
      timeIntegration_ == "bdf2") {
    return true;
  } else {
    return false;
  }
}

// member function to determine of method is vicous or inviscid
bool input::IsViscous() const {
  if (equationSet_ == "navierStokes" ||
      equationSet_ == "rans") {
    return true;
  } else {
    return false;
  }
}

// member function to determine of method is turbulent
bool input::IsTurbulent() const {
  if (equationSet_ == "rans") {
    return true;
  } else {
    return false;
  }
}

// member function to determine if solution should use a block matrix
bool input::IsBlockMatrix() const {
  if (this->IsImplicit() && (matrixSolver_ == "bdplur" ||
                             matrixSolver_ == "blusgs")) {
    return true;
  } else {
    return false;
  }
}

string input::OrderOfAccuracy() const {
  if (kappa_ == -2.0) {
    return "first";
  } else {
    return "second";
  }
}

unique_ptr<turbModel> input::AssignTurbulenceModel() const {
  // define turbulence model
  unique_ptr<turbModel> turb(nullptr);
  if (turbModel_ == "none") {
    turb = unique_ptr<turbModel>{std::make_unique<turbNone>()};
  } else if (turbModel_ == "kOmegaWilcox2006") {
    turb = unique_ptr<turbModel>{std::make_unique<turbKWWilcox>()};
  } else if (turbModel_ == "sst2003") {
    turb = unique_ptr<turbModel>{std::make_unique<turbKWSst>()};
  } else {
    cerr << "ERROR: Error in input::AssignTurbulenceModel(). Turbulence model "
         << turbModel_ << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
  return turb;
}

// member function to return the name of the simulation without the file
// extension i.e. "myInput.inp" would return "myInput"
string input::SimNameRoot() const {
  return simName_.substr(0, simName_.find("."));
}

// member function to check validity of nonlinear iterations
void input::CheckNonlinearIterations() {
  if (timeIntegration_ == "rk4" && nonlinearIterations_ != 4) {
    cerr << "WARNING: For RK4 method, nonlinear iterations should be set to "
         << 4 << " changing value from " << nonlinearIterations_ << " to "
         << 4 << endl;
    nonlinearIterations_ = 4;
  }

  if (timeIntegration_ == "euler" && nonlinearIterations_ != 1) {
    cerr << "WARNING: For euler method, nonlinear iterations should be set to "
         << 1 << " changing value from " << nonlinearIterations_ << " to " << 1
         << endl;
    nonlinearIterations_ = 1;
  }
}

// member function to check validity of the requested output variables
void input::CheckOutputVariables() {
  for (auto var : outputVariables_) {
    if (!this->IsTurbulent()) {  // can't have turbulent varibles output
      if (var == "tke" || var == "sdr" || var == "viscosityRatio" ||
          var.find("tkeGrad_") != string::npos ||
          var.find("sdrGrad_") != string::npos || var == "resid_tke" ||
          var == "resid_sdr") {
        cerr << "WARNING: Variable " << var <<
            " is not available for laminar simulations." << endl;
        outputVariables_.erase(var);
      }

      if (!this->IsViscous()) {  // can't have viscous variables output
        if (var.find("velGrad_") != string::npos
            || var.find("tempGrad_") != string::npos) {
          cerr << "WARNING: Variable " << var <<
              " is not available for inviscid simulations." << endl;
          outputVariables_.erase(var);
        }
      }
    }
  }
}

// member function to check that turbulence model makes sense with equation set
void input::CheckTurbulenceModel() const {
  if (equationSet_ == "rans" && turbModel_ == "none") {
    cerr << "ERROR: If solving RANS equations, must specify turbulence model!"
         << endl;
    exit(EXIT_FAILURE);
  }
  if (equationSet_ != "rans" && turbModel_ != "none") {
    cerr << "ERROR: Turbulence models are only valid for the RANS equation set!"
         << endl;
    exit(EXIT_FAILURE);
  }
}


// member function to calculate the coefficient used to scale the viscous
// spectral radius in the time step calculation
double input::ViscousCFLCoefficient() const {
  auto coeff = 1.0;
  if (this->Kappa() == 1.0) {  // central
    coeff = 4.0;
  } else if (this->Kappa() == -2.0) {  // first order upwind
    coeff = 2.0;
  }
  return coeff;
}

bool input::MatrixRequiresInitialization() const {
  // initialize matrix if using DPLUR / BDPLUR, or if using LUSGS / BLUSGS with
  // more than one sweep
  return (matrixSolver_ == "dplur" || matrixSolver_ == "bdplur" ||
          matrixSweeps_ > 1) ? true : false;
}

int input::NumberGhostLayers() const {
  if (this->OrderOfAccuracy() == "first") {
    return 1;
  } else {
    return 2;
  }
}

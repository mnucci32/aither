/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

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

#include <chrono>    // timing capability
#include <iostream>  // cout
#include <iomanip>   // put_time
#include <fstream>   // ifstream
#include <cstdlib>   // exit()
#include <sstream>   // istringstream
#include <iterator>  // istring_iterator
#include <memory>    // make_unique
#include <string>
#include <vector>
#include "input.hpp"
#include "turbulence.hpp"
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
input::input(const string &name) : simName_(name), vars_(35) {
  // default values for each variable
  gName_ = "";
  dt_ = -1.0;
  iterations_ = 1;
  pRef_ = -1.0;
  rRef_ = -1.0;
  lRef_ = 1.0;
  gamma_ = 1.4;
  gasConst_ = 287.058;
  velRef_.SetX(1.0);
  velRef_.SetY(0.0);
  velRef_.SetZ(0.0);
  vector<boundaryConditions> dummy(1);
  bc_ = dummy;
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

  // keywords in the input file that the parser is looking for to define
  // variables
  vars_[0] = "gridName:";
  vars_[1] = "timeStep:";
  vars_[2] = "iterations:";
  vars_[3] = "pressureRef:";
  vars_[4] = "densityRef:";
  vars_[5] = "lengthRef:";
  vars_[6] = "gamma:";
  vars_[7] = "gasConstant:";
  vars_[8] = "velocity:";
  vars_[9] = "timeIntegration:";
  vars_[10] = "cfl:";
  vars_[11] = "faceReconstruction:";
  vars_[12] = "limiter:";
  vars_[13] = "outputFrequency:";
  vars_[14] = "equationSet:";
  vars_[15] = "temperatureRef:";
  vars_[16] = "matrixSolver:";
  vars_[17] = "matrixSweeps:";
  vars_[18] = "matrixRelaxation:";
  vars_[19] = "timeIntTheta:";
  vars_[20] = "timeIntZeta:";
  vars_[21] = "nonlinearIterations:";
  vars_[22] = "cflMax:";
  vars_[23] = "cflStep:";
  vars_[24] = "cflStart:";
  vars_[25] = "inviscidFluxJacobian:";
  vars_[26] = "dualTimeCFL:";
  vars_[27] = "inviscidFlux:";
  vars_[28] = "stagnationInlet:";
  vars_[29] = "pressureOutlet:";
  vars_[30] = "decompositionMethod:";
  vars_[31] = "turbulenceModel:";
  vars_[32] = "farfieldTurbulenceIntensity:";
  vars_[33] = "farfieldEddyViscosityRatio:";

  // boundary conditions should be listed last
  vars_[vars_.size() - 1] = "boundaryConditions:";
}

// function to trim leading and trailing whitespace from a string, and also
// remove data after a comment
string trim(const string &s, const string &whitespace = " \t") {
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

// function to print the time
void PrintTime() {
  auto now = std::chrono::high_resolution_clock::now();
  auto nowOut = std::chrono::high_resolution_clock::to_time_t(now);
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
  vector3d<double> temp;
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
    line = trim(line);

    // split line into words
    istringstream buf(line);
    istream_iterator<string> beg(buf), end;
    vector<string> tokens(beg, end);

    // search to see if first token corresponds to any keywords
    // line must contain at least 2 tokens (keyword and value)
    if (tokens.size() >= 2) {
      for (auto ii = 0; ii < this->NumVars();
           ii++) {  // loop over all input variables and see if they match
                    // keywords

        if (tokens[0] == this->Vars(ii) ||
            readingBCs > 0) {  // if first token matches a keyword or reading
                               // boundary conditions

          // if not yet reading BCs (readingBCs == 0), set variable in input
          // class to corresponding value and print assignment to std out
          if (ii == 0 && readingBCs == 0) {
            gName_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->GridName() << endl;
            }
            continue;
          } else if (ii == 1 && readingBCs == 0) {
            dt_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->Dt() << endl;
            }
            continue;
          } else if (ii == 2 && readingBCs == 0) {
            iterations_ = stoi(tokens[1]);
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->Iterations() << endl;
            }
            continue;
          } else if (ii == 3 && readingBCs == 0) {
            pRef_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->PRef() << endl;
            }
            continue;
          } else if (ii == 4 && readingBCs == 0) {
            rRef_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->RRef() << endl;
            }
            continue;
          } else if (ii == 5 && readingBCs == 0) {
            lRef_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->LRef() << endl;
            }
            continue;
          } else if (ii == 6 && readingBCs == 0) {
            gamma_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->Gamma() << endl;
            }
            continue;
          } else if (ii == 7 && readingBCs == 0) {
            gasConst_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->R() << endl;
            }
            continue;
          } else if (ii == 8 && readingBCs == 0) {
            temp.SetX(stod(tokens[1]));  // double variable (stod)
            temp.SetY(stod(tokens[2]));
            temp.SetZ(stod(tokens[3]));
            velRef_ = temp;
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->VelRef() << endl;
            }
            continue;
          } else if (ii == 9 && readingBCs == 0) {
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
              cout << this->Vars(ii) << " " << this->TimeIntegration()
                   << endl;
            }
            continue;
          } else if (ii == 11 && readingBCs == 0) {
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
              cout << this->Vars(ii) << " " << tokens[1]
                   << " kappa = " << this->Kappa() << endl;
            }
            if ((this->Kappa() < -1.0) || (this->Kappa() > 1.0)) {
              cerr << "ERROR: Error in input::ReadInput(). Kappa value of "
                   << this->Kappa()
                   << " is not valid! Choose a value between -1.0 and 1.0."
                   << endl;
              exit(EXIT_FAILURE);
            }
            continue;
          } else if (ii == 12 && readingBCs == 0) {
            limiter_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->Limiter() << endl;
            }
            continue;
          } else if (ii == 13 && readingBCs == 0) {
            outputFrequency_ = stoi(tokens[1]);
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->OutputFrequency()
                   << endl;
            }
            continue;
          } else if (ii == 14 && readingBCs == 0) {
            equationSet_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->EquationSet() << endl;
            }
            continue;
          } else if (ii == 15 && readingBCs == 0) {
            tRef_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->TRef() << endl;
            }
            continue;
          } else if (ii == 16 && readingBCs == 0) {
            matrixSolver_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->MatrixSolver() << endl;
            }
            continue;
          } else if (ii == 17 && readingBCs == 0) {
            matrixSweeps_ = stoi(tokens[1]);
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->MatrixSweeps() << endl;
            }
            continue;
          } else if (ii == 18 && readingBCs == 0) {
            matrixRelaxation_ =
                stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->MatrixRelaxation()
                   << endl;
            }
            continue;
          } else if (ii == 21 && readingBCs == 0) {
            nonlinearIterations_ = stoi(tokens[1]);
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->NonlinearIterations()
                   << endl;
            }
            continue;
          } else if (ii == 22 && readingBCs == 0) {
            cflMax_ = stod(tokens[1]);  // double  (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->CFLMax() << endl;
            }
            continue;
          } else if (ii == 23 && readingBCs == 0) {
            cflStep_ = stod(tokens[1]);  // double (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->CFLStep() << endl;
            }
            continue;
          } else if (ii == 24 && readingBCs == 0) {
            cflStart_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->CFLStart() << endl;
            }
            continue;
          } else if (ii == 25 && readingBCs == 0) {
            invFluxJac_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->InvFluxJac() << endl;
            }
          } else if (ii == 26 && readingBCs == 0) {
            dualTimeCFL_ = stod(tokens[1]);  // double variable (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->DualTimeCFL() << endl;
            }
            continue;
          } else if (ii == 27 && readingBCs == 0) {
            inviscidFlux_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->InviscidFlux() << endl;
            }
          } else if (ii == 28 && readingBCs == 0) {
            stagInletProps_[0] = stoi(tokens[1]);  // tag
            stagInletProps_[1] = stod(tokens[2]);  // stag pressure
            stagInletProps_[2] = stod(tokens[3]);  // stag temp
            stagInletProps_[3] = stod(tokens[4]);  // dir-x
            stagInletProps_[4] = stod(tokens[5]);  // dir-y
            stagInletProps_[5] = stod(tokens[6]);  // dir-z
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->StagInletTag() << " "
                   << this->StagInletP0() << " " << this->StagInletT0()
                   << " " << this->StagInletDx() << " "
                   << this->StagInletDy() << " " << this->StagInletDz()
                   << endl;
            }
          } else if (ii == 29 && readingBCs == 0) {
            pressureOutlet_[0] = stoi(tokens[1]);  // tag
            pressureOutlet_[1] = stod(tokens[2]);  // outlet pressure
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->PressureOutletTag()
                   << " " << this->PressureOutletP() << endl;
            }
          } else if (ii == 30 && readingBCs == 0) {
            decompMethod_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->DecompMethod() << endl;
            }
            continue;
          } else if (ii == 31 && readingBCs == 0) {
            turbModel_ = tokens[1];
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->TurbulenceModel()
                   << endl;
            }
            continue;
          } else if (ii == 32 && readingBCs == 0) {
            farfieldTurbInten_ = stod(tokens[1]);  // double (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->FarfieldTurbIntensity()
                   << endl;
            }
            continue;
          } else if (ii == 33 && readingBCs == 0) {
            farfieldEddyViscRatio_ = stod(tokens[1]);  // double (stod)
            if (rank == ROOTP) {
              cout << this->Vars(ii) << " " << this->FarfieldEddyViscRatio()
                   << endl;
            }
            continue;

          // reading BCs
          // -------------------------------------------------------------
          } else if (ii == this->NumVars() - 1 || readingBCs > 0) {
            // read in boundary conditions and assign to boundaryConditions
            // class
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
              tempBC[blk].ResizeVecs(stoi(tokens[0]),
                                     stoi(tokens[1]),
                                     stoi(tokens[2]));

              // set total number of surfaces
              numSurf = stoi(tokens[0]) + stoi(tokens[1]) +
                        stoi(tokens[2]);

              lEnd += (numSurf +
                       1);  // set new target for when block data is finished
              readingBCs++;

            } else {  // assign BC block variables
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
                cout << this->Vars(this->NumVars() - 1) << " "
                     << this->NumBC() << endl << endl;
                for (auto ll = 0; ll < this->NumBC(); ll++) {
                  cout << "Block: " << ll << endl;
                  cout << this->BC(ll) << endl;
                }
              }
            }

            break;
          }
        }
      }
    } else {  // if there aren't 2 or more tokens just skip line
      continue;
    }
  }

  // input file sanity checks
  this->CheckNonlinearIterations();

  if (rank == ROOTP) {
    cout << endl;
    cout << "Input file parse complete" << endl;
    cout << "##################################################################"
            "#########################################################" << endl
         << endl;
  }
}

// member function to calculate the cfl_ value for the step from the starting,
// ending, and step values
void input::CalcCFL(const int &i) {
  if (i == 0) {  // first time step
    cfl_ = cflStart_;
  } else if (cflStart_ + i * cflStep_ >
             cflMax_) {  // if calculated value is higher than max, set to max
    cfl_ = cflMax_;
  } else {
    cfl_ = cflStart_ + i * cflStep_;
  }
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
  // initialize matrix if using DPLUR / BDPLUR, or if using LUSGS with more
  // than one sweep
  return (matrixSolver_ == "dplur" || matrixSolver_ == "bdplur" ||
          matrixSweeps_ > 1) ? true : false;
}

/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <time.h>    // strftime
#include <iostream>  // cout
#include <fstream>   // ifstream
#include <cstdlib>   // exit()
#include <sstream>   // istringstream
#include <iterator>  // istring_iterator
#include <string>
#include <vector>
#include "input.hpp"

#define ROOT 0

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::cerr;
using std::istringstream;
using std::istream_iterator;

// constructor for input class
// initialize vector to have length of number of acceptable inputs to the code
input::input() : vars_(33) {
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
  invFluxJac_ = "laxFriedrichs";  // default is approximate roe flux
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

  vars_[32] = "boundaryConditions:";  // boundary conditions should be listed
                                     // last
}

// function to trim leading and trailing whitespace from a string, and also
// remove data after a comment
string trim(const string &s, const string &whitespace = " \t") {
  const string comment = "#";  //# is comment character for input file

  if (s.empty()) {
    return "";  // string is empty
  } else {
    const unsigned int sBegin = s.find_first_not_of(
        whitespace);  // find index of first non whitespace character
    const int sEnd = s.find_last_not_of(
        whitespace);  // find index of last non whitespace character
    const int sRange = sEnd - sBegin + 1;  // range to trim string to
    string temp = s.substr(sBegin, sRange);

    const int tempComment =
        temp.find(comment);  // find index of first comment character
    const int tempRange = tempComment - 0;  // find range of string to return

    return temp.substr(0, tempRange);
  }
}

// function to print the time
void PrintTime() {
  time_t rawtime;
  struct tm timeinfo;

  time(&rawtime);
  localtime_r(&rawtime, &timeinfo);

  const int kBufLen = 100;  // hard coded max buffer length
  char buffer[kBufLen];
  strftime(buffer, kBufLen, "%c", &timeinfo);
  cout << buffer << endl;
}

// function to read the input file and return the data as a member of the input
// class
// read the input file
void input::ReadInput(const string &inputName, const int &rank) {
  // inputName -- name of input file
  // rank -- rank of processor

  if (rank == ROOT) {
    cout << "##################################################################"
            "#########################################################" << endl;
    PrintTime();
    cout << endl;
    cout << "Parsing input file " << inputName << endl << endl;
    cout << "Solver Inputs" << endl;
  }

  // open input file
  ifstream inFile;
  inFile.open(inputName.c_str(), ios::in);
  if (inFile.fail()) {
    cerr << "ERROR: Error in input::ReadInput(). Input file " << inputName
         << " did not open correctly!!!" << endl;
    exit(0);
  }

  string line;
  vector3d<double> temp;
  int readingBCs = 0;
  int blk = 0;
  int surfCounter = 0;
  int numBCBlks = 0;
  vector<boundaryConditions> tempBC(1);
  int lEnd = 0;
  int numSurf = 0;

  while (getline(inFile, line)) {  // while there are still lines in the input
                                   // file, execute loop

    line = trim(
        line);  // remove leading and trailing whitespace and ignore comments

    // split line into words
    istringstream buf(line);
    istream_iterator<string> beg(buf), end;
    vector<string> tokens(beg, end);

    // search to see if first token corresponds to any keywords
    // line must contain at least 2 tokens (keyword and value)
    if (tokens.size() >= 2) {
      for (int ii = 0; ii < (*this).NumVars();
           ii++) {  // loop over all input variables and see if they match
                    // keywords

        if (tokens[0] == (*this).Vars(ii) ||
            readingBCs > 0) {  // if first token matches a keyword or reading
                               // boundary conditions

          // if not yet reading BCs (readingBCs == 0), set variable in input
          // class to corresponding value and print assignment to std out
          if (ii == 0 && readingBCs == 0) {
            (*this).gName_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).GridName() << endl;
            }
            continue;
          } else if (ii == 1 && readingBCs == 0) {
            (*this).dt_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).Dt() << endl;
            }
            continue;
          } else if (ii == 2 && readingBCs == 0) {
            (*this).iterations_ = atoi(tokens[1].c_str());
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).Iterations() << endl;
            }
            continue;
          } else if (ii == 3 && readingBCs == 0) {
            (*this).pRef_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).PRef() << endl;
            }
            continue;
          } else if (ii == 4 && readingBCs == 0) {
            (*this).rRef_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).RRef() << endl;
            }
            continue;
          } else if (ii == 5 && readingBCs == 0) {
            (*this).lRef_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).LRef() << endl;
            }
            continue;
          } else if (ii == 6 && readingBCs == 0) {
            (*this).gamma_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).Gamma() << endl;
            }
            continue;
          } else if (ii == 7 && readingBCs == 0) {
            (*this).gasConst_ =
                atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).R() << endl;
            }
            continue;
          } else if (ii == 8 && readingBCs == 0) {
            temp.SetX(atof(tokens[1].c_str()));  // double variable (atof)
            temp.SetY(atof(tokens[2].c_str()));
            temp.SetZ(atof(tokens[3].c_str()));
            (*this).velRef_ = temp;
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).VelRef() << endl;
            }
            continue;
          } else if (ii == 9 && readingBCs == 0) {
            (*this).timeIntegration_ = tokens[1];
            if ((*this).TimeIntegration() == "implicitEuler") {
              (*this).timeIntTheta_ = 1.0;
              (*this).timeIntZeta_ = 0.0;
            } else if ((*this).TimeIntegration() == "crankNicholson") {
              (*this).timeIntTheta_ = 0.5;
              (*this).timeIntZeta_ = 0.0;
            } else if ((*this).TimeIntegration() == "bdf2") {
              (*this).timeIntTheta_ = 1.0;
              (*this).timeIntZeta_ = 0.5;
            }

            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).TimeIntegration()
                   << endl;
            }
            continue;
          } else if (ii == 11 && readingBCs == 0) {
            if (tokens[1] == "upwind") {
              (*this).kappa_ = -1.0;
            } else if (tokens[1] == "fromm") {
              (*this).kappa_ = 0.0;
            } else if (tokens[1] == "quick") {
              (*this).kappa_ = 0.5;
            } else if (tokens[1] == "central") {
              (*this).kappa_ = 1.0;
            } else if (tokens[1] == "thirdOrder") {
              (*this).kappa_ = 1.0 / 3.0;
            } else if (tokens[1] == "constant") {
              (*this).kappa_ = -2.0;
            } else {
              // if string is not recognized, set kappa to number input
              (*this).kappa_ = atof(tokens[1].c_str());
            }

            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << tokens[1]
                   << " kappa_ = " << (*this).Kappa() << endl;
            }
            if (((*this).Kappa() < -1.0) || ((*this).Kappa() > 1.0)) {
              cerr << "ERROR: Error in input::ReadInput(). Kappa value of "
                   << (*this).Kappa()
                   << " is not valid! Choose a value between -1.0 and 1.0."
                   << endl;
              exit(0);
            }
            continue;
          } else if (ii == 12 && readingBCs == 0) {
            (*this).limiter_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).Limiter() << endl;
            }
            continue;
          } else if (ii == 13 && readingBCs == 0) {
            (*this).outputFrequency_ = atoi(tokens[1].c_str());
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).OutputFrequency()
                   << endl;
            }
            continue;
          } else if (ii == 14 && readingBCs == 0) {
            (*this).equationSet_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).EquationSet() << endl;
            }
            continue;
          } else if (ii == 15 && readingBCs == 0) {
            (*this).tRef_ = atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).TRef() << endl;
            }
            continue;
          } else if (ii == 16 && readingBCs == 0) {
            (*this).matrixSolver_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).MatrixSolver() << endl;
            }
            continue;
          } else if (ii == 17 && readingBCs == 0) {
            (*this).matrixSweeps_ = atoi(tokens[1].c_str());
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).MatrixSweeps() << endl;
            }
            continue;
          } else if (ii == 18 && readingBCs == 0) {
            (*this).matrixRelaxation_ =
                atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).MatrixRelaxation()
                   << endl;
            }
            continue;
          } else if (ii == 21 && readingBCs == 0) {
            (*this).nonlinearIterations_ = atoi(tokens[1].c_str());
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).NonlinearIterations()
                   << endl;
            }
            continue;
          } else if (ii == 22 && readingBCs == 0) {
            (*this).cflMax_ = atof(tokens[1].c_str());  // double  (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).CFLMax() << endl;
            }
            continue;
          } else if (ii == 23 && readingBCs == 0) {
            (*this).cflStep_ = atof(tokens[1].c_str());  // double (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).CFLStep() << endl;
            }
            continue;
          } else if (ii == 24 && readingBCs == 0) {
            (*this).cflStart_ =
                atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).CFLStart() << endl;
            }
            continue;
          } else if (ii == 25 && readingBCs == 0) {
            (*this).invFluxJac_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).InvFluxJac() << endl;
            }
          } else if (ii == 26 && readingBCs == 0) {
            (*this).dualTimeCFL_ =
                atof(tokens[1].c_str());  // double variable (atof)
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).DualTimeCFL() << endl;
            }
            continue;
          } else if (ii == 27 && readingBCs == 0) {
            (*this).inviscidFlux_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).InviscidFlux() << endl;
            }
          } else if (ii == 28 && readingBCs == 0) {
            (*this).stagInletProps_[0] = atoi(tokens[1].c_str());  // tag
            (*this).stagInletProps_[1] =
                atof(tokens[2].c_str());  // stag pressure
            (*this).stagInletProps_[2] = atof(tokens[3].c_str());  // stag temp
            (*this).stagInletProps_[3] = atof(tokens[4].c_str());  // dir-x
            (*this).stagInletProps_[4] = atof(tokens[5].c_str());  // dir-y
            (*this).stagInletProps_[5] = atof(tokens[6].c_str());  // dir-z
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).StagInletTag() << " "
                   << (*this).StagInletP0() << " " << (*this).StagInletT0()
                   << " " << (*this).StagInletDx() << " "
                   << (*this).StagInletDy() << " " << (*this).StagInletDz()
                   << endl;
            }
          } else if (ii == 29 && readingBCs == 0) {
            (*this).pressureOutlet_[0] = atoi(tokens[1].c_str());  // tag
            (*this).pressureOutlet_[1] =
                atof(tokens[2].c_str());  // outlet pressure
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).PressureOutletTag()
                   << " " << (*this).PressureOutletP() << endl;
            }
          } else if (ii == 30 && readingBCs == 0) {
            (*this).decompMethod_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).DecompMethod() << endl;
            }
            continue;
          } else if (ii == 31 && readingBCs == 0) {
            (*this).turbModel_ = tokens[1];
            if (rank == ROOT) {
              cout << (*this).Vars(ii) << " " << (*this).TurbulenceModel()
                   << endl;
            }
            continue;
          // reading BCs
          // -------------------------------------------------------------
          } else if (ii == (*this).NumVars() - 1 || readingBCs > 0) {
            // read in boundary conditions and assign to boundaryConditions
            // class
            if (readingBCs == 0) {  // variable read must be number of blocks if
                                    // first line of BCs
              numBCBlks = atoi(tokens[1].c_str());
              tempBC.resize(numBCBlks);  // resize vector holding BCs to correct
                                         // number of blocks
              readingBCs++;
              lEnd++;
            } else if (readingBCs == lEnd) {  // variables must be number of i,
                                              // j, k surfaces for each block
              // set number of i, j, k surfaces and resize vectors
              tempBC[blk].ResizeVecs(atoi(tokens[0].c_str()),
                                     atoi(tokens[1].c_str()),
                                     atoi(tokens[2].c_str()));

              // set total number of surfaces
              numSurf = atoi(tokens[0].c_str()) + atoi(tokens[1].c_str()) +
                        atoi(tokens[2].c_str());

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
              (*this).bc_ = tempBC;
              if (rank == ROOT) {
                cout << (*this).Vars((*this).NumVars() - 1) << " "
                     << (*this).NumBC() << endl << endl;
                for (int ll = 0; ll < (*this).NumBC(); ll++) {
                  cout << "Block: " << ll << endl;
                  cout << (*this).BC(ll) << endl;
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

  if (rank == ROOT) {
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

// member function to determine number of equations to solver for
int input::NumEquations() const {
  int numEqns = 0;
  if (((*this).EquationSet() == "euler") ||
      ((*this).EquationSet() == "navierStokes")) {
    numEqns = 5;
  } else if ((*this).EquationSet() == "rans") {
    numEqns = 7;
  } else {
    cerr << "ERROR: Equations set is not recognized. Cannot determine number "
            "of equations!" << endl;
  }

  return numEqns;
}

// member function to determine of method is implicit or explicit
bool input::IsImplicit() const {
  if ((*this).TimeIntegration() == "implicitEuler" ||
      (*this).TimeIntegration() == "crankNicholson" ||
      (*this).TimeIntegration() == "bdf2") {
    return true;
  } else {
    return false;
  }
}

// member function to determine of method is vicous or inviscid
bool input::IsViscous() const {
  if ((*this).EquationSet() == "navierStokes" ||
      (*this).EquationSet() == "rans") {
    return true;
  } else {
    return false;
  }
}

// member function to determine of method is turbulent
bool input::IsTurbulent() const {
  if ((*this).EquationSet() == "rans") {
    return true;
  } else {
    return false;
  }
}

string input::OrderOfAccuracy() const {
  if ((*this).Kappa() == -2.0) {
    return "first";
  } else {
    return "second";
  }
}

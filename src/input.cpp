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
#include "eos.hpp"
#include "transport.hpp"
#include "thermodynamic.hpp"
#include "fluid.hpp"
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
input::input(const string &name, const string &resName) : simName_(name),
                                                          restartName_(resName) {
  // default values for each variable
  gName_ = "";
  dt_ = -1.0;
  iterations_ = 1;
  rRef_ = -1.0;
  tRef_ = -1.0;
  lRef_ = 1.0;
  aRef_ = 0.0;
  fluids_ = vector<fluid>(1);
  bc_ = vector<boundaryConditions>(1);
  timeIntegration_ = "explicitEuler";
  cfl_ = -1.0;
  faceReconstruction_ = "constant";
  viscousFaceReconstruction_ = "central";
  kappa_ = -2.0;  // default to value outside of range to tell if higher
                  // order or constant method should be used
  limiter_ = "none";
  outputFrequency_ = 1;
  equationSet_ = "euler";
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
  decompMethod_ = "cubic";  // default is cubic decomposition
  turbModel_ = "none";  // default turbulence model is none
  thermodynamicModel_ = "caloricallyPerfect";  // default to cpg
  equationOfState_ = "idealGas";  // default to ideal gas
  transportModel_ = "sutherland";  // default to sutherland
  restartFrequency_ = 0;  // default to not write restarts
  iterationStart_ = 0;  // default to start from iteration zero

  // default to primative variables
  outputVariables_ = {"density", "vel_x", "vel_y", "vel_z", "pressure"};
  wallOutputVariables_ = {};

  // keywords in the input file that the parser is looking for to define
  // variables
  vars_ = {"gridName",
           "timeStep",
           "iterations",
           "referenceDensity",
           "referenceTemperature",
           "referenceLength",
           "fluids",
           "timeIntegration",
           "faceReconstruction",
           "viscousFaceReconstruction",
           "limiter",
           "outputFrequency",
           "restartFrequency",
           "equationSet",
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
           "decompositionMethod",
           "turbulenceModel",
           "thermodynamicModel",
           "equationOfState",
           "transportModel",
           "outputVariables",
           "wallOutputVariables",
           "initialConditions",
           "boundaryStates",
           "boundaryConditions"};
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
        } else if (key == "referenceDensity") {
          rRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->RRef() << endl;
          }
        } else if (key == "referenceTemperature") {
          tRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->TRef() << endl;
          }
        } else if (key == "referenceLength") {
          lRef_ = stod(tokens[1]);  // double variable (stod)
          if (rank == ROOTP) {
            cout << key << ": " << this->LRef() << endl;
          }
        } else if (key == "fluids") {
          fluids_ = ReadFluidList(inFile, tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": <";
            for (auto ii = 0U; ii < fluids_.size(); ++ii) {
              cout << fluids_[ii];
              if (ii == fluids_.size() - 1) {
                cout << ">" << endl;
              } else {
                cout << "," << endl << "                    ";
              }
            }
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
            faceReconstruction_ = tokens[1];
            kappa_ = -1.0;
          } else if (tokens[1] == "fromm") {
            faceReconstruction_ = tokens[1];
            kappa_ = 0.0;
          } else if (tokens[1] == "quick") {
            faceReconstruction_ = tokens[1];
            kappa_ = 0.5;
          } else if (tokens[1] == "central") {
            faceReconstruction_ = tokens[1];
            kappa_ = 1.0;
          } else if (tokens[1] == "thirdOrder") {
            faceReconstruction_ = tokens[1];
            kappa_ = 1.0 / 3.0;
          } else if (tokens[1] == "constant" || tokens[1] == "weno" ||
                     tokens[1] == "wenoZ") {
            faceReconstruction_ = tokens[1];
          } else {
            // face reconstruction method not recognized
            cerr << "ERROR: Error in input::ReadInput(). Face reconstruction "
                 << "method " << tokens[1] << " is not recognized!";
                exit(EXIT_FAILURE);
          }

          if (rank == ROOTP) {
            cout << key << ": " << this->FaceReconstruction() << endl;
          }
        } else if (key == "viscousFaceReconstruction") {
          if (tokens[1] == "central") {
            viscousFaceReconstruction_ = tokens[1];
          } else if (tokens[1] == "centralFourth") {
            viscousFaceReconstruction_ = tokens[1];
          } else {
            // face reconstruction method not recognized
            cerr << "ERROR: Error in input::ReadInput(). Viscous face "
                 << "reconstruction method " << tokens[1] << " is not recognized!";
                exit(EXIT_FAILURE);
          }

          if (rank == ROOTP) {
            cout << key << ": " << this->ViscousFaceReconstruction() << endl;
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
        } else if (key == "restartFrequency") {
          restartFrequency_ = stoi(tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": " << this->RestartFrequency() << endl;
          }
        } else if (key == "equationSet") {
          equationSet_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->EquationSet() << endl;
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
        } else if (key == "thermodynamicModel") {
          thermodynamicModel_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->ThermodynamicModel() << endl;
          }
        } else if (key == "equationOfState") {
          equationOfState_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->EquationOfState() << endl;
          }
        } else if (key == "transportModel") {
          transportModel_ = tokens[1];
          if (rank == ROOTP) {
            cout << key << ": " << this->TransportModel() << endl;
          }
        } else if (key == "outputVariables") {
          // clear default variables from set
          outputVariables_.clear();
          auto specifiedVars = ReadStringList(inFile, tokens[1]);
          for (auto &vars : specifiedVars) {
            outputVariables_.insert(vars);
          }
          if (rank == ROOTP) {
            cout << key << ": <";
            auto count = 0U;
            auto numChars = 0U;
            for (auto &vars : outputVariables_) {
              if (count == outputVariables_.size() - 1) {
                cout << vars << ">" << endl;
              } else {
                cout << vars << ", ";
                numChars += vars.length();
                if (numChars >= 50) {  // if more than 50 chars, go to next line
                  cout << endl << "                  ";
                  numChars = 0U;
                }
              }
              count++;
            }
            cout << endl;
          }
        } else if (key == "wallOutputVariables") {
          // clear default variables from set
          wallOutputVariables_.clear();
          auto specifiedVars = ReadStringList(inFile, tokens[1]);
          for (auto &vars : specifiedVars) {
            wallOutputVariables_.insert(vars);
          }
          if (rank == ROOTP) {
            cout << key << ": <";
            auto count = 0U;
            auto numChars = 0U;
            for (auto &vars : wallOutputVariables_) {
              if (count == wallOutputVariables_.size() - 1) {
                cout << vars << ">" << endl;
              } else {
                cout << vars << ", ";
                numChars += vars.length();
                if (numChars >= 50) {  // if more than 50 chars, go to next line
                  cout << endl << "                      ";
                  numChars = 0U;
                }
              }
              count++;
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
        } else if (key == "boundaryStates") {
          bcStates_ = ReadBCList(inFile, tokens[1]);
          if (rank == ROOTP) {
            cout << key << ": <";
            for (auto ii = 0U; ii < bcStates_.size(); ++ii) {
              cout << *bcStates_[ii];
              if (ii == bcStates_.size() - 1) {
                cout << ">" << endl;
              } else {
                cout << "," << endl << "                 ";
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
  this->CheckWallOutputVariables();
  this->CheckTurbulenceModel();
  this->CheckSpecies();
  this->CheckNonreflecting();

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
  return (this->IsRANS()) ? 2 : 0;
}

// member function to determine number of equations to solver for
int input::NumEquations() const {
  auto numEqns = 0;
  if (equationSet_ == "euler" ||
      equationSet_ == "navierStokes" ||
      equationSet_ == "largeEddySimulation") {
    numEqns = this->NumFlowEquations();
  } else if (this->IsRANS()) {
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
  if (equationSet_ == "navierStokes" || this->IsTurbulent()) {
    return true;
  } else {
    return false;
  }
}

// member function to determine of method is turbulent
bool input::IsTurbulent() const {
  if (this->IsRANS() || this->IsLES()) {
    return true;
  } else {
    return false;
  }
}

// member function to determine if simulation is RANS
bool input::IsRANS() const {
  return (equationSet_ == "rans") ? true : false;
}

// member function to determine if simulation is LES
bool input::IsLES() const {
  return (equationSet_ == "largeEddySimulation") ? true : false;
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
  if (this->UsingConstantReconstruction()) {
    return "first";
  } else {
    return "second";
  }
}

bool input::UsingMUSCLReconstruction() const {
  auto muscl = false;
  if (faceReconstruction_ == "upwind" || faceReconstruction_ == "fromm" ||
      faceReconstruction_ == "quick" || faceReconstruction_ == "central" ||
      faceReconstruction_ == "thirdOrder") {
    muscl = true;
  }
  return muscl;
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
  } else if (turbModel_ == "sstdes") {
    turb = unique_ptr<turbModel>{std::make_unique<turbSstDes>()};
  } else if (turbModel_ == "wale") {
    turb = unique_ptr<turbModel>{std::make_unique<turbWale>()};
  } else {
    cerr << "ERROR: Error in input::AssignTurbulenceModel(). Turbulence model "
         << turbModel_ << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
  return turb;
}

// member function to get equation of state
unique_ptr<eos> input::AssignEquationOfState(
    const unique_ptr<thermodynamic> &thermo) {
  // get fluid
  auto fl = this->Fluid();
  // define equation of state
  unique_ptr<eos> eqnState(nullptr);
  if (equationOfState_ == "idealGas") {
    // nondimensional temperature is 1.0 (tRef_ / tRef_)
    eqnState = unique_ptr<eos>{
        std::make_unique<idealGas>(thermo, fl.GasConstant(), 1.0)};
  } else {
    cerr << "ERROR: Error in input::AssignEquationOfState(). Equation of state "
         << equationOfState_ << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
  // use equation of state to assign additional reference values
  const auto pRef = eqnState->PressureDim(rRef_, tRef_);
  aRef_ = eqnState->SoS(pRef, rRef_);
  return eqnState;
}

// member function to get transport model
unique_ptr<transport> input::AssignTransportModel() const {
  // define equation of state
  unique_ptr<transport> trans(nullptr);
  if (transportModel_ == "sutherland") {
    trans = unique_ptr<transport>{std::make_unique<sutherland>(
        tRef_, rRef_, lRef_, aRef_)};
  } else {
    cerr << "ERROR: Error in input::AssignTransportModel(). Transport model "
         << transportModel_ << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
  return trans;
}

// member function to get thermodynamic model
unique_ptr<thermodynamic> input::AssignThermodynamicModel() const {
  // get fluid
  auto fl = this->Fluid();
  // define equation of state
  unique_ptr<thermodynamic> thermo(nullptr);
  if (thermodynamicModel_ == "caloricallyPerfect") {
    thermo =
        unique_ptr<thermodynamic>{std::make_unique<caloricallyPerfect>(fl.N())};
  } else if (thermodynamicModel_ == "thermallyPerfect") {
    thermo = unique_ptr<thermodynamic>{std::make_unique<thermallyPerfect>(
        fl.N(), fl.VibrationalTemperature())};
  } else {
    cerr << "ERROR: Error in input::AssignThermodynamicModel(). Thermodynamic "
         << "model " << thermodynamicModel_ << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
  return thermo;
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

  if (timeIntegration_ == "explicitEuler" && nonlinearIterations_ != 1) {
    cerr << "WARNING: For euler method, nonlinear iterations should be set to "
         << 1 << " changing value from " << nonlinearIterations_ << " to " << 1
         << endl;
    nonlinearIterations_ = 1;
  }
}

// member function to check validity of the requested output variables
void input::CheckOutputVariables() {
  auto oVars = outputVariables_;
  for (auto &var : oVars) {
    if (!this->IsRANS()) {  // can't have RANS variables output
      if (var == "tke" || var == "sdr" || var.find("tkeGrad_") != string::npos ||
          var.find("sdrGrad_") != string::npos || var == "resid_tke" ||
          var == "resid_sdr" || var == "f1" || var == "f2") {
        cerr << "WARNING: Variable " << var <<
            " is not available for non-RANS simulations." << endl;
        outputVariables_.erase(var);
      }
    }

    if (!this->IsTurbulent()) {  // can't have turbulent variables output
      if (var == "viscosityRatio" || var == "turbulentViscosity") {
        cerr << "WARNING: Variable " << var <<
            " is not available for laminar simulations." << endl;
        outputVariables_.erase(var);
      }
    }

    if (!this->IsViscous()) {  // can't have viscous variables output
      if (var == "viscosity") {
        cerr << "WARNING: Variable " << var <<
            " is not available for inviscid simulations." << endl;
        outputVariables_.erase(var);
      }
    }
  }
}

// member function to check validity of the requested output variables
void input::CheckWallOutputVariables() {
  auto oVars = wallOutputVariables_;
  for (auto &var : oVars) {
    if (!this->IsViscous()) {  // can't have viscous variables output
      if (var == "yplus" || var == "heatFlux" || var == "shearStress" ||
          var == "frictionVelocity" || var == "viscosity") {
        cerr << "WARNING: Wall variable " << var
             << " is not available for inviscid simulations." << endl;
        outputVariables_.erase(var);
      }
    }
    if (!this->IsTurbulent()) {  // can't have turbulent variables output
      if (var == "viscosityRatio") {
        cerr << "WARNING: Wall variable " << var
             << " is not available for laminar simulations." << endl;
        outputVariables_.erase(var);
      }
    }
    if (!this->IsRANS()) {  // can't have RANS variables output
      if (var == "tke" || var == "sdr") {
        cerr << "WARNING: Wall variable " << var
             << " is not available for laminar simulations." << endl;
        outputVariables_.erase(var);
      }
    }
  }
}

// member function to check that turbulence model makes sense with equation
// set
void input::CheckTurbulenceModel() const {
  if (this->IsTurbulent() && turbModel_ == "none") {
    cerr << "ERROR: If solving RANS or LES equations, must specify turbulence "
         << "model!" << endl;
    exit(EXIT_FAILURE);
  }
  if (!this->IsTurbulent() && turbModel_ != "none") {
    cerr << "ERROR: Turbulence models are only valid for the RANS and LES "
         << "equation sets!" << endl;
    exit(EXIT_FAILURE);
  }
  if (this->IsRANS() && turbModel_ == "wale") {
    cerr << "ERROR: Equation set is RANS, but turbulence model is not!" << endl;
    exit(EXIT_FAILURE);
  }
  if (this->IsTurbulent() && !this->IsRANS() && turbModel_ != "wale") {
    cerr << "ERROR: Equation set is LES, but turbulence model is not!" << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to check that all species specified are defined
void input::CheckSpecies() const {
  // check all species in ICs are defined
  for (auto &ic : ics_) {
    if (!ic.IsFromFile()) {
      auto fracs = ic.MassFractions();  // get mass fractions for ICs
      // loop over all species and find if that species has been defined
      for (auto &species : fracs) {
        auto name = species.first;
        if (find_if(std::begin(fluids_), std::end(fluids_),
                    [&name](const fluid &fl) { return fl.Name() == name; }) ==
            std::end(fluids_)) {
          cerr << "ERROR: Species " << name << " is not defined!" << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // check all species in BCs are defined
  for (auto &bc : bcStates_) {
    auto fracs = bc->MassFractions();  // get mass fractions for BCs
    if (!fracs.empty()) {              // some BCs don't have mass fraction data
      for (auto &species : fracs) {
        auto name = species.first;
        if (find_if(std::begin(fluids_), std::end(fluids_),
                    [&name](const fluid &fl) { return fl.Name() == name; }) ==
            std::end(fluids_)) {
          cerr << "ERROR: Species " << name << " is not defined!" << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // check that there is at least one species defined
  if (fluids_.size() < 1) {
    cerr << "ERROR: At least one fluid species must be defined!" << endl;
    exit(EXIT_FAILURE);
  }
}

// check that nonreflecting BCs aren't used with explicit euler because
// solution at time N is not stored
void input::CheckNonreflecting() const {
  if (timeIntegration_ == "explicitEuler") {
    for (auto &bc : bcStates_) {
      if (bc->IsNonreflecting()) {
        cerr << "ERROR: Nonreflecting BCs cannot be used with explicitEuler "
                "time integration!"
             << endl;
        exit(EXIT_FAILURE);
      }
    }
  }
}

// member function to check that all species specified are defined
// vector of species comes from prescribed ic file
void input::CheckSpecies(const vector<string> &species) const {
  for (auto &name : species) {
    if (find_if(std::begin(fluids_), std::end(fluids_),
                [&name](const fluid &fl) { return fl.Name() == name; }) ==
        std::end(fluids_)) {
      cerr << "ERROR: Species " << name << " is not defined!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // check that there is at least one species from ic file
  if (species.size() < 1) {
    cerr << "ERROR: At least one fluid species must be defined in ic file!"
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
  auto layers = 0;
  if (this->UsingConstantReconstruction()) {
    layers = 1;
  } else if (this->UsingMUSCLReconstruction()) {
    layers = 2;
  } else if (this->UsingHigherOrderReconstruction()) {
    layers = 3;
  } else {
    cerr << "ERROR: Problem with face reconstruction. Not using any of the "
         << "supported methods!" << endl;
    exit(EXIT_FAILURE);
  }

  auto viscLayers = (this->ViscousFaceReconstruction() == "centralFourth") ? 2 : 1;
  return std::max(layers, viscLayers);
}

// member function to get the initial condition state for a given parent block
icState input::ICStateForBlock(const int &block) const {
  icState blockIC;
  auto foundExactMatch = false;
  auto foundDefaultMatch = false;

  // exact match trumps default
  for (auto &ic : ics_) {
    if (ic.Tag() == -1 && !foundExactMatch) {
      blockIC = ic;
      foundDefaultMatch = true;
    } else if (ic.Tag() == block) {
      blockIC = ic;
      foundExactMatch = true;
    }
  }

  // sanity check -- must specify default or exact match
  if (!foundDefaultMatch && !foundExactMatch) {
    cerr << "ERROR. Did not find default or specified initial condition state "
         << "for block " << block << endl;
    exit(EXIT_FAILURE);
  }

  return blockIC;
}

// member function to return boundary condition data index for a given tag
const shared_ptr<inputState> & input::BCData(const int &tag) const {
  for (auto &bcData : bcStates_) {
    // for non-periodic bcs tag = startTag = endTag
    if (bcData->Tag() == tag || bcData->EndTag() == tag) {
      return bcData;
    }
  }

  // if function doesn't return, tag not found
  cerr << "ERROR. Could not find data for boundary condition tag " << tag << endl;
  exit(EXIT_FAILURE);
}

void input::NondimensionalizeStateData(const unique_ptr<eos> &eqnState) {
  for (auto &state : bcStates_) {
    state->Nondimensionalize(rRef_, tRef_, lRef_, aRef_);
  }
  for (auto &ic : ics_) {
    ic.Nondimensionalize(rRef_, tRef_, lRef_, aRef_);
  }
}

void input::NondimensionalizeFluid() {
  for (auto &fl : fluids_) {
    fl.Nondimensionalize(tRef_);
  }
}

// default value is 0; code currently only supports single fluid flows
fluid input::Fluid(const int ind) const { return fluids_[ind]; }
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

#ifndef INPUTHEADERDEF  // only if the macro INPUTHEADERDEF is not defined
                        // execute these lines of code
#define INPUTHEADERDEF  // define the macro

#include <vector>                   // vector
#include <string>                   // string
#include <memory>                   // unique_ptr
#include "vector3d.hpp"
#include "boundaryConditions.hpp"

using std::vector;
using std::string;
using std::unique_ptr;

// forward class declaration
class turbModel;

class input {
  string simName_;  // simulation name
  string gName_;  // grid file name
  double dt_;  // time step
  int iterations_;  // number of iterations

  // variable names that are regcognized by the input file
  vector<string> vars_;
                        // parser
  double pRef_;  // reference pressure
  double rRef_;  // reference density
  double lRef_;  // reference length
  double gamma_;  // ratio of specific heats
  double gasConst_;  // gas constant of fluid
  vector3d<double> velRef_;  // reference velocity
  vector<boundaryConditions> bc_;  // vector of boundary conditions for each
                                  // block
  string timeIntegration_;  // time integration method
  double cfl_;  // cfl number for local time stepping
  double kappa_;  // kappa paramenter for MUSCL face reconstruction
  string limiter_;  // limiter to use in higher order calculations
  int outputFrequency_;  // how often to output results
  string equationSet_;  // which set of equations to solver Euler/Navier-Stokes
  double tRef_;  // reference temperature
  string matrixSolver_;  // matrix solver to solve Ax=b
  int matrixSweeps_;  // number of sweeps for matrix solver
  double matrixRelaxation_;  // relaxation parameter for matrix solver
  double timeIntTheta_;  // beam and warming time integration parameter
  double timeIntZeta_;  // beam and warming time integration parameter
  int nonlinearIterations_;  // number of nonlinear iterations for time accurate
                            // scheme
  double cflMax_;  // maximum cfl_ value
  double cflStep_;  // cfl_ step size for ramp
  double cflStart_;  // starting cfl_ number
  string invFluxJac_;  // inviscid flux jacobian
  double dualTimeCFL_;  // cfl_ number for dual time
  string inviscidFlux_;  // scheme for inviscid flux calculation
  double stagInletProps_[6];  // vector of stagnation inlet properties
  double pressureOutlet_[2];  // vector of pressure outlet properties
  string decompMethod_;  // method of decomposition for parallel problems
  string turbModel_;  // turbulence model
  double farfieldTurbInten_;  // turbulence intensity at farfield
  double farfieldEddyViscRatio_;  // eddy viscosity ratio at farfield

 public:
  // constructor
  explicit input(const string &);

  // move constructor and assignment operator
  input(input&&) noexcept = default;
  input& operator=(input&&) noexcept = default;

  // copy constructor and assignment operator
  input(const input&) = default;
  input& operator=(const input&) = default;

  // member functions
  string SimName() const {return simName_;}
  string SimNameRoot() const;
  string GridName() const {return gName_;}

  double Dt() const {return dt_;}

  int Iterations() const {return iterations_;}

  double PRef() const {return pRef_;}
  double RRef() const {return rRef_;}
  double LRef() const {return lRef_;}
  double TRef() const {return tRef_;}

  double Gamma() const {return gamma_;}
  double R() const {return gasConst_;}

  vector3d<double> VelRef() const {return velRef_;}

  boundaryConditions BC(const int &ind) const {return bc_[ind];}
  vector<boundaryConditions> AllBC() const {return bc_;}
  int NumBC() const {return bc_.size();}

  string TimeIntegration() const {return timeIntegration_;}

  double CFL() const {return cfl_;}
  void CalcCFL(const int &i);

  double Kappa() const {return kappa_;}

  string Limiter() const {return limiter_;}

  int OutputFrequency() const {return outputFrequency_;}

  string EquationSet() const {return equationSet_;}

  string MatrixSolver() const {return matrixSolver_;}
  int MatrixSweeps() const {return matrixSweeps_;}
  double MatrixRelaxation() const {return matrixRelaxation_;}

  double Theta() const {return timeIntTheta_;}
  double Zeta() const {return timeIntZeta_;}

  int NonlinearIterations() const {return nonlinearIterations_;}

  double CFLMax() const {return cflMax_;}
  double CFLStep() const {return cflStep_;}
  double CFLStart() const {return cflStart_;}

  string InvFluxJac() const {return invFluxJac_;}

  double DualTimeCFL() const {return dualTimeCFL_;}

  string InviscidFlux() const {return inviscidFlux_;}

  int StagInletTag() const {return static_cast<int> (stagInletProps_[0]);}
  double StagInletP0() const {return stagInletProps_[1];}
  double StagInletT0() const {return stagInletProps_[2];}
  double StagInletDx() const {return stagInletProps_[3];}
  double StagInletDy() const {return stagInletProps_[4];}
  double StagInletDz() const {return stagInletProps_[5];}

  int PressureOutletTag() const {
    return static_cast<int> (pressureOutlet_[0]);}
  double PressureOutletP() const {return pressureOutlet_[1];}

  string DecompMethod() const {return decompMethod_;}
  string TurbulenceModel() const {return turbModel_;}

  string Vars(const int &ind) const {return vars_[ind];}

  int NumVars() const {return vars_.size();}
  int NumEquations() const;

  void ReadInput(const int &);

  bool IsImplicit() const;
  bool IsViscous() const;
  bool IsTurbulent() const;

  string OrderOfAccuracy() const;

  double FarfieldTurbIntensity() const {return farfieldTurbInten_;}
  double FarfieldEddyViscRatio() const {return farfieldEddyViscRatio_;}

  unique_ptr<turbModel> AssignTurbulenceModel() const;

  void CheckNonlinearIterations();

  double ViscousCFLCoefficient() const;

  // destructor
  ~input() noexcept {}
};

// function declarations
void PrintTime();

#endif

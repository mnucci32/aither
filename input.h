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

#ifndef INPUTHEADERDEF             //only if the macro INPUTHEADERDEF is not defined execute these lines of code
#define INPUTHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h"
#include "boundaryConditions.h"

using std::vector;
using std::string;

class input {

  string gName;                         //grid file name
  double dt;                            //time step
  int iterations;                       //number of iterations
  vector<string> vars;                  //variable names that are regcognized by the input file parser
  double pRef;                          //reference pressure
  double rRef;                          //reference density
  double lRef;                          //reference length
  double gamma;                         //ratio of specific heats
  double gasConst;                      //gas constant of fluid
  vector3d<double> velRef;              //reference velocity
  vector<boundaryConditions> bc;        //vector of boundary conditions for each block
  string timeIntegration;               //time integration method
  double cfl;                           //cfl number for local time stepping
  double kappa;                         //kappa paramenter for MUSCL face reconstruction
  string limiter;                       //limiter to use in higher order calculations
  int outputFrequency;                  //how often to output results
  string equationSet;                   //which set of equations to solver Euler/Navier-Stokes
  double tRef;                          //reference temperature
  string matrixSolver;                  //matrix solver to solve Ax=b
  int matrixSweeps;                     //number of sweeps for matrix solver
  double matrixRelaxation;              //relaxation parameter for matrix solver
  double timeIntTheta;                  //beam and warming time integration parameter
  double timeIntZeta;                   //beam and warming time integration parameter
  int nonlinearIterations;              //number of nonlinear iterations for time accurate scheme
  double cflMax;                        //maximum cfl value
  double cflStep;                       //cfl step size for ramp
  double cflStart;                      //starting cfl number
  string invFluxJac;                    //inviscid flux jacobian
  double dualTimeCFL;                   //cfl number for dual time
  string inviscidFlux;                  //scheme for inviscid flux calculation
  double stagInletProps[6];             //vector of stagnation inlet properties
  double pressureOutlet[2];             //vector of pressure outlet properties
  string decompMethod;                  //method of decomposition for parallel problems
  string turbModel;                     //turbulence model

 public:
  //constructor
  input();

  //member functions
  string GridName()const{return gName;}

  double Dt()const{return dt;}

  int Iterations()const{return iterations;}

  double PRef()const{return pRef;}
  double RRef()const{return rRef;}
  double LRef()const{return lRef;}
  double TRef()const{return tRef;}

  double Gamma()const{return gamma;}
  double R()const{return gasConst;}

  vector3d<double> VelRef()const{return velRef;}

  boundaryConditions BC(const int &ind)const{return bc[ind];}
  vector<boundaryConditions> AllBC()const{return bc;}
  int NumBC()const{return bc.size();}

  string TimeIntegration()const{return timeIntegration;}

  double CFL()const{return cfl;}
  void CalcCFL(const int &i);

  double Kappa()const{return kappa;}

  string Limiter()const{return limiter;}

  int OutputFrequency()const{return outputFrequency;}

  string EquationSet()const{return equationSet;}

  string MatrixSolver()const{return matrixSolver;}
  int MatrixSweeps()const{return matrixSweeps;}
  double MatrixRelaxation()const{return matrixRelaxation;}

  double Theta()const{return timeIntTheta;}

  double Zeta()const{return timeIntZeta;}

  int NonlinearIterations()const{return nonlinearIterations;}

  double CFLMax()const{return cflMax;}
  double CFLStep()const{return cflStep;}
  double CFLStart()const{return cflStart;}

  string InvFluxJac()const{return invFluxJac;}

  double DualTimeCFL()const{return dualTimeCFL;}

  string InviscidFlux()const{return inviscidFlux;}

  int StagInletTag()const{return (int)stagInletProps[0];}
  double StagInletP0()const{return stagInletProps[1];}
  double StagInletT0()const{return stagInletProps[2];}
  double StagInletDx()const{return stagInletProps[3];}
  double StagInletDy()const{return stagInletProps[4];}
  double StagInletDz()const{return stagInletProps[5];}

  int PressureOutletTag()const{return (int)pressureOutlet[0];}
  double PressureOutletP()const{return pressureOutlet[1];}

  string DecompMethod()const{return decompMethod;}
  string TurbulenceModel()const{return turbModel;}

  string Vars(const int &ind)const{return vars[ind];}

  int NumVars()const{return vars.size();}
  int NumEquations()const;

  void ReadInput(const string &, const int &);

  //destructor
  ~input() {}

};

//function declarations
void PrintTime();


#endif

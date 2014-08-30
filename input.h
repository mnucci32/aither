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

 public:
  //constructor
  input();

  //member functions
  string GridName()const{return gName;}
  void SetGridName(const string &name){gName = name;}

  double Dt()const{return dt;}
  void SetDt(const double &a){dt = a;}

  int Iterations()const{return iterations;}
  void SetIterations(const int &i){iterations = i;}

  double PRef()const{return pRef;}
  void SetPRef(const double &a){pRef = a;}

  double RRef()const{return rRef;}
  void SetRRef(const double &a){rRef = a;}

  double LRef()const{return lRef;}
  void SetLRef(const double &a){lRef = a;}

  double TRef()const{return tRef;}
  void SetTRef(const double &a){tRef = a;}

  double Gamma()const{return gamma;}
  void SetGamma(const double &a){gamma = a;}

  double R()const{return gasConst;}
  void SetR(const double &a){gasConst = a;}

  vector3d<double> VelRef()const{return velRef;}
  void SetVelRef(const vector3d<double> &a){velRef = a;}

  boundaryConditions BC(const int &ind)const{return bc[ind];}
  int NumBC()const{return bc.size();}
  void SetBC(const vector<boundaryConditions> &a){bc = a;}
  void SetBCVec(const int &a);

  string TimeIntegration()const{return timeIntegration;}
  void SetTimeIntegration(const string &a){timeIntegration = a;}

  double CFL()const{return cfl;}
  void SetCFL(const double &a){cfl = a;}
  void CalcCFL(const int &i);

  double Kappa()const{return kappa;}
  void SetKappa(const double &a){kappa = a;}

  string Limiter()const{return limiter;}
  void SetLimiter(const string &a){limiter = a;}

  int OutputFrequency()const{return outputFrequency;}
  void SetOutputFrequency(const int &a){outputFrequency = a;}

  string EquationSet()const{return equationSet;}
  void SetEquationSet(const string &a){equationSet = a;}

  string MatrixSolver()const{return matrixSolver;}
  void SetMatrixSolver(const string &a){matrixSolver = a;}

  int MatrixSweeps()const{return matrixSweeps;}
  void SetMatrixSweeps(const int &i){matrixSweeps = i;}

  double MatrixRelaxation()const{return matrixRelaxation;}
  void SetMatrixRelaxation(const double &a){matrixRelaxation = a;}

  double Theta()const{return timeIntTheta;}
  void SetTheta(const double &a){timeIntTheta = a;}

  double Zeta()const{return timeIntZeta;}
  void SetZeta(const double &a){timeIntZeta = a;}

  int NonlinearIterations()const{return nonlinearIterations;}
  void SetNonlinearIterations(const int &i){nonlinearIterations = i;}

  double CFLMax()const{return cflMax;}
  void SetCFLMax(const double &a){cflMax = a;}

  double CFLStep()const{return cflStep;}
  void SetCFLStep(const double &a){cflStep = a;}

  double CFLStart()const{return cflStart;}
  void SetCFLStart(const double &a){cflStart = a;}

  string InvFluxJac()const{return invFluxJac;}
  void SetInvFluxJac(const string &a){invFluxJac = a;}

  string Vars(const int &ind)const{return vars[ind];}

  int NumVars()const{return vars.size();}

  //destructor
  ~input() {}

};

//function declarations
input ReadInput(const string &inputName);
void PrintTime();


#endif

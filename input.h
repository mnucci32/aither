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
  vector<double> stagInletProps;        //vector of stagnation inlet properties
  vector<double> pressureOutlet;        //vector of pressure outlet properties
  string decompMethod;                  //method of decomposition for parallel problems

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
  vector<boundaryConditions> AllBC()const{return bc;}
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

  double DualTimeCFL()const{return dualTimeCFL;}
  void SetDualTimeCFL(const double &a){dualTimeCFL = a;}

  string InviscidFlux()const{return inviscidFlux;}
  void SetInviscidFlux(const string &a){inviscidFlux = a;}

  int StagInletTag()const{return (int)stagInletProps[0];}
  void SetStagInletTag(const int &a){stagInletProps[0] = (double)a;}
  double StagInletP0()const{return stagInletProps[1];}
  void SetStagInletP0(const double &a){stagInletProps[1] = a;}
  double StagInletT0()const{return stagInletProps[2];}
  void SetStagInletT0(const double &a){stagInletProps[2] = a;}
  double StagInletDx()const{return stagInletProps[3];}
  void SetStagInletDx(const double &a){stagInletProps[3] = a;}
  double StagInletDy()const{return stagInletProps[4];}
  void SetStagInletDy(const double &a){stagInletProps[4] = a;}
  double StagInletDz()const{return stagInletProps[5];}
  void SetStagInletDz(const double &a){stagInletProps[5] = a;}

  int PressureOutletTag()const{return (int)pressureOutlet[0];}
  void SetPressureOutletTag(const int &a){pressureOutlet[0] = (double)a;}
  double PressureOutletP()const{return pressureOutlet[1];}
  void SetPressureOutletP(const double &a){pressureOutlet[1] = a;}

  string DecompMethod()const{return decompMethod;}
  void SetDecompMethod(const string &name){decompMethod = name;}

  string Vars(const int &ind)const{return vars[ind];}

  int NumVars()const{return vars.size();}

  //destructor
  ~input() {}

};

//function declarations
input ReadInput(const string &inputName, const int &rank);
void PrintTime();


#endif

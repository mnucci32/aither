#ifndef STATEVARSHEADERDEF             //only if the macro STATEVARSHEADERDEF is not defined execute these lines of code
#define STATEVARSHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "plot3d.h" //vector3d
#include "input.h" //inputVars
#include "eos.h"  // idealGas
#include "boundaryConditions.h" //boundaryConditions
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class stateVars {
  double rho;            //conserved variables at cell center
  double rhoU;
  double rhoV;
  double rhoW;
  double rhoE;

 public:
  //constructors
  stateVars() : rho(0.0), rhoU(0.0), rhoV(0.0), rhoW(0.0), rhoE(0.0) {}
  stateVars( double r, double e, vector3d<double> v) : rho(r), rhoU(r*v.X()), rhoV(r*v.Y()), rhoW(r*v.Z()), rhoE(r*e) {}

  //member functions
  void SetRho( const double &a){rho = a;}
  double Rho()const{return rho;}
  void SetRhoU( const double &a){rhoU = a;}
  double RhoU()const{return rhoU;}
  void SetRhoV( const double &a){rhoV = a;}
  double RhoV()const{return rhoV;}
  void SetRhoW( const double &a){rhoW = a;}
  double RhoW()const{return rhoW;}
  void SetRhoE( const double &a){rhoE = a;}
  double RhoE()const{return rhoE;}

  inline double Vx()const;
  inline double Vy()const;
  inline double Vz()const;
  inline double Energy()const;
  inline double Pressure(const idealGas&)const;
  inline double Enthalpy(const idealGas&)const;
  inline vector3d<double> Velocity()const;
  inline double Temperature(const idealGas&)const;

  //operator overloads for addition and subtraction of states
  stateVars operator + (const stateVars&)const;
  stateVars operator - (const stateVars&)const;
  stateVars operator * (const stateVars&)const;
  stateVars operator / (const stateVars&)const;

  friend stateVars operator + (const double&, const stateVars&);
  friend stateVars operator - (const double&, const stateVars&);
  friend stateVars operator * (const double&, const stateVars&);
  friend stateVars operator / (const double&, const stateVars&);
  friend ostream & operator<< (ostream &os, const stateVars&);

  //member function to calculate reconstruction of state variables from cell center to cell face 
  //assuming value at cell center is constant over cell volume; zeroth order reconstruction results in first order accuracy
  stateVars FaceReconConst()const{return *this;}     

  //member function to calculate reconstruction of state variables from cell center to cell face
  //this function uses muscle extrapolation resulting in higher order accuracy
  stateVars FaceReconMUSCL( const stateVars&, const stateVars&, const string&, const double&, const string&, const double=1.0, const double=1.0, const double=1.0)const;

  //member function to calculate Van Albada limiter function
  stateVars LimiterVanAlbada( const stateVars& )const;
  stateVars LimiterMinmod( const stateVars&, const stateVars&, const double )const;
  stateVars LimiterNone()const;

  //member function to return the state of the appropriate ghost cell
  stateVars GetGhostState( const string&, const vector3d<double>&, const string&, const input&, const idealGas&)const;

  /* //member function to advance the state using the explicit euler method (1st order accurate) */
  /* stateVars ExplicitEulerTimeAdvance( double, double, vector<double>& ); */

  /* //member function to advance the state using the 4th order Runge-Kutta method */
  /* stateVars RK4TimeAdvance( double, double, vector<double>&, int ); */

  //destructor
  ~stateVars() {}

};

//function definitions
//member function to calculate x velocity from conserved variables
double stateVars::Vx()const{
  return rhoU / rho ;
}

//member function to calculate x velocity from conserved variables
double stateVars::Vy()const{
  return rhoV / rho ;
}

//member function to calculate z velocity from conserved variables
double stateVars::Vz()const{
  return rhoW / rho ;
}

//member function to calculate total enthalpy from conserved variables
double stateVars::Energy()const{
  return rhoE / rho ;
}

//member function to calculate pressure from conserved variables and equation of state
double stateVars::Pressure(const idealGas &eqnState)const{

  vector3d<double> vel(rhoU/rho, rhoV/rho, rhoW/rho);
  double nonDimSpecE = rhoE/rho - 0.5 * vel.Mag() * vel.Mag();
  return eqnState.GetPressure(rho, nonDimSpecE);

}

//member function to calculate temperature from conserved variables and equation of state
double stateVars::Temperature(const idealGas &eqnState)const{

  double press = (*this).Pressure(eqnState);

  return eqnState.GetTemperature(press, rho);

}

//member function to calculate velocity from conserved variables
vector3d<double> stateVars::Velocity()const{
  vector3d<double> vel(rhoU/rho, rhoV/rho, rhoW/rho);
  return vel;
}



//member function to calculate enthalpy from conserved variables and equation of state
double stateVars::Enthalpy(const idealGas &eqnState)const{

  vector3d<double> vel(rhoU/rho, rhoV/rho, rhoW/rho);

  double nonDimSpecE = rhoE/rho - 0.5 * vel.Mag() * vel.Mag();
  double press = eqnState.GetPressure(rho, nonDimSpecE);

  return eqnState.GetEnthalpy(rhoE/rho, press, rho);

}

#endif

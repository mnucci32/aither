#ifndef PRIMVARSHEADERDEF             //only if the macro PRIMVARSHEADERDEF is not defined execute these lines of code
#define PRIMVARSHEADERDEF             //define the macro

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

class primVars {
  double rho;            //primative variables at cell center
  double u;
  double v;
  double w;
  double p;

 public:
  //constructors
  primVars() : rho(0.0), u(0.0), v(0.0), w(0.0), p(0.0) {}
  primVars( double r, double p, vector3d<double> v) : rho(r), u(v.X()), v(v.Y()), w(v.Z()), p(p) {}
  primVars( double a, double b, double c, double d, double e) : rho(a), u(b), v(c), w(d), p(e) {}

  //member functions
  void SetRho( const double &a){rho = a;}
  double Rho()const{return rho;}
  void SetU( const double &a){u = a;}
  double U()const{return u;}
  void SetV( const double &a){v = a;}
  double V()const{return v;}
  void SetW( const double &a){w = a;}
  double W()const{return w;}
  void SetP( const double &a){p = a;}
  double P()const{return p;}

  inline vector3d<double> Velocity()const;

  inline double Energy(const idealGas&)const;
  inline double Enthalpy(const idealGas&)const;
  inline double Temperature(const idealGas&)const;
  inline double SoS(const idealGas&)const;

  inline vector<double> ConsVars(const idealGas&)const;

  //operator overloads for addition and subtraction of states
  primVars operator + (const primVars&)const;
  primVars operator - (const primVars&)const;
  primVars operator * (const primVars&)const;
  primVars operator / (const primVars&)const;

  friend primVars operator + (const double&, const primVars&);
  friend primVars operator - (const double&, const primVars&);
  friend primVars operator * (const double&, const primVars&);
  friend primVars operator / (const double&, const primVars&);
  friend ostream & operator<< (ostream &os, const primVars&);

  //member function to calculate reconstruction of state variables from cell center to cell face 
  //assuming value at cell center is constant over cell volume; zeroth order reconstruction results in first order accuracy
  primVars FaceReconConst()const{return *this;}     

  //member function to calculate reconstruction of state variables from cell center to cell face
  //this function uses muscle extrapolation resulting in higher order accuracy
  primVars FaceReconMUSCL( const primVars&, const primVars&, const string&, const double&, const string&, const double=1.0, const double=1.0, const double=1.0)const;

  //member function to calculate Van Albada limiter function
  primVars LimiterVanAlbada( const primVars& )const;
  primVars LimiterMinmod( const primVars&, const primVars&, const double )const;
  primVars LimiterNone()const;

  //member function to return the state of the appropriate ghost cell
  primVars GetGhostState( const string&, const vector3d<double>&, const string&, const input&, const idealGas&)const;

  //destructor
  ~primVars() {}

};

//function definitions

//member function to calculate temperature from conserved variables and equation of state
double primVars::Temperature(const idealGas &eqnState)const{
  return eqnState.GetTemperature(p, rho);
}

//member function to calculate velocity from conserved variables
vector3d<double> primVars::Velocity()const{
  vector3d<double> vel(u, v, w);
  return vel;
}

//member function to calculate total enthalpy from conserved variables
double primVars::Energy(const idealGas &eqnState)const{
  return eqnState.GetEnergy( eqnState.GetSpecEnergy(p, rho), (*this).Velocity().Mag() );
}

//member function to calculate speed of sound from primative varialbes
double primVars::SoS(const idealGas &eqnState)const{
  return sqrt(eqnState.Gamma() * p / rho);
}

//member function to calculate enthalpy from conserved variables and equation of state
double primVars::Enthalpy(const idealGas &eqnState)const{
  return eqnState.GetEnthalpy((*this).Energy(eqnState), p, rho);
}

//member function to calculate conserved variables from primative variables
vector<double> primVars::ConsVars(const idealGas &eqnState)const{
  vector<double> cv (5);
  cv[0] = rho;
  cv[1] = rho * u;
  cv[2] = rho * v;
  cv[3] = rho * w;
  cv[4] = rho * (*this).Energy(eqnState);
  return cv;
}


#endif

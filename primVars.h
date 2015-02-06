#ifndef PRIMVARSHEADERDEF             //only if the macro PRIMVARSHEADERDEF is not defined execute these lines of code
#define PRIMVARSHEADERDEF             //define the macro

/* This header contains the primVars class.

The primVars class stores the primative variables for the Euler and Navier-Stokes equations [rho, u, v, w, P]. It contains several member functions
to manipulate these primative varibles. It also contains member functions to extrapolate the primative variables from the cell center to the cell
face using constant and MUSCL reconstruction. It is has a member function to supply a ghost state given a boundary condition and boundary cell.  */

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "plot3d.h" //vector3d
#include "input.h" //inputVars
#include "eos.h"  // idealGas
#include "matrix.h" //colMatrix
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
  double data[5];        //primative variables at cell center

 public:
  //constructors
 primVars() : data{0.0, 0.0, 0.0, 0.0, 0.0} {}
 primVars( double r, double p, vector3d<double> v) : data{r, v.X(), v.Y(), v.Z(), p} {}
 primVars( double a, double b, double c, double d, double e) : data{a, b, c, d, e} {}

  //member functions
  void SetRho( const double &a){data[0] = a;}
  double Rho()const{return data[0];}
  void SetU( const double &a){data[1] = a;}
  double U()const{return data[1];}
  void SetV( const double &a){data[2] = a;}
  double V()const{return data[2];}
  void SetW( const double &a){data[3] = a;}
  double W()const{return data[3];}
  void SetP( const double &a){data[4] = a;}
  double P()const{return data[4];}

  void NondimensionalInitialize(const idealGas &, const vector3d<double>&);

  inline vector3d<double> Velocity()const;

  inline double Energy(const idealGas&)const;
  inline double Enthalpy(const idealGas&)const;
  inline double Temperature(const idealGas&)const;
  inline double SoS(const idealGas&)const;

  inline colMatrix ConsVars(const idealGas&)const;
  primVars UpdateWithConsVars(const idealGas&, const colMatrix&)const;

  //operator overloads for addition and subtraction of states
  primVars operator + (const primVars&)const;
  primVars operator - (const primVars&)const;
  primVars operator * (const primVars&)const;
  primVars operator / (const primVars&)const;

  primVars operator + (const double&);
  primVars operator - (const double&);
  primVars operator * (const double&);
  primVars operator / (const double&);

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
  primVars FaceReconMUSCL( const primVars&, const primVars&, const double&, const string&, const double=1.0, const double=1.0, const double=1.0)const;

  //member function to calculate Van Albada limiter function
  primVars LimiterVanAlbada( const primVars& )const;
  primVars LimiterMinmod( const primVars&, const primVars&, const double )const;
  primVars LimiterNone()const;

  //member function to return the state of the appropriate ghost cell
  primVars GetGhostState( const string&, const vector3d<double>&, const string&, const input&, const idealGas&, const int=1)const;

  //destructor
  ~primVars() {}

};

//function definitions

//member function to calculate temperature from conserved variables and equation of state
double primVars::Temperature(const idealGas &eqnState)const{
  return eqnState.GetTemperature(data[4], data[0]);
}

//member function to calculate velocity from conserved variables
vector3d<double> primVars::Velocity()const{
  vector3d<double> vel(data[1], data[2], data[3]);
  return vel;
}

//member function to calculate total enthalpy from conserved variables
double primVars::Energy(const idealGas &eqnState)const{
  return eqnState.GetEnergy( eqnState.GetSpecEnergy(data[4], data[0]), (*this).Velocity().Mag() );
}

//member function to calculate speed of sound from primative varialbes
double primVars::SoS(const idealGas &eqnState)const{
  return sqrt(eqnState.Gamma() * data[4] / data[0]);
}

//member function to calculate enthalpy from conserved variables and equation of state
double primVars::Enthalpy(const idealGas &eqnState)const{
  return eqnState.GetEnthalpy((*this).Energy(eqnState), data[4], data[0]);
}

//member function to calculate conserved variables from primative variables
colMatrix primVars::ConsVars(const idealGas &eqnState)const{
  colMatrix cv(5);
  cv.SetData(0, data[0]);
  cv.SetData(1, data[0] * data[1]);
  cv.SetData(2, data[0] * data[2]);
  cv.SetData(3, data[0] * data[3]);
  cv.SetData(4, data[0] * (*this).Energy(eqnState));
  return cv;
}


#endif

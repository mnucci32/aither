#ifndef EOSHEADERDEF             //only if the macro EOSHEADERDEF is not defined execute these lines of code
#define EOSHEADERDEF             //define the macro

/* This header file contains the idealGas and sutherland classes. 

The ideal gas class stores the ratio of specific heats (gamma),
and the gas constant R. It contains several member functions to calculate state variables using the equation of state
P = rho * (gamma - 1) * specificEnergy for the euler equations and P = rho * R * temperature for the Navier-Stokes equations.

The sutherland class stores a reference temperature and viscosity, as well as the sutherland coefficients. It is used for calculating
a temperature dependent viscosity for the Navier-Stokes equations. */

#include <vector>  //vector
#include <string>  //string
#include <math.h>       //sqrt

using std::vector;
using std::string;

class idealGas {
  const double gamma;
  const double gasConst;


 public:
  //constructor
  idealGas() : gamma(1.4), gasConst(287.058) {}
  idealGas( const double a, const double b ) : gamma(a), gasConst(b) {}

  //member functions
  double GetPressure(const double &rho, const double &specEn)const;
  double GetPressFromEnergy(const double &rho, const double &energy, const double &vel)const;
  double GetDensity(const double &pressure, const double &specEn)const;
  double GetSpecEnergy(const double &pressure, const double &rho)const;
  double GetEnergy(const double &specEn, const double &vel)const;
  double GetEnthalpy(const double &energy, const double &pressure, const double &rho)const;
  double GetSoS(const double &pressure, const double &rho)const;
  double Gamma()const{return gamma;}
  double GasConst()const{return gasConst;}
  double GetPrandtl()const{return (4.0*gamma)/(9.0*gamma-5.0);}

  double GetTemperature(const double &pressure, const double &rho)const;
  double GetConductivity(const double &mu)const{return mu/( (*this).GetPrandtl()*(gamma-1.0) );} //nondimensional version (R=1/gamma)
  double GetDensityTP(const double &temp, const double &press)const{return press*gamma/temp;}

  //destructor
  ~idealGas() {}

};

//function declarations


class sutherland {

  double cOne;
  double S;
  double tRef;
  double muRef;
  double bulkVisc;

 public:
  //constructors
  //Stoke's hypothesis -- bulk viscosity = 0
  //Sutherland's Law -- mu = muref * (C1 * Tref^1.5) / (T + S)
 sutherland() : cOne(1.458e-6), S(110.4), tRef(288.15), muRef(cOne * pow(tRef,1.5)/(tRef+S)), bulkVisc(0.0) {}
 sutherland(const double a, const double b, const double c) : cOne(a), S(b), tRef(c), muRef(cOne * pow(c,1.5)/(c+S)), bulkVisc(0.0) {}
  sutherland(const double t) : cOne(1.458e-6), S(110.4), tRef(t), muRef(cOne * pow(t,1.5)/(t+S)), bulkVisc(0.0) {}

  //member functions
  double GetViscosity(const double &)const;
  double GetLambda(const double &)const;
  double ConstC1()const{return cOne;}
  double ConstS()const{return S;}
  double TRef()const{return tRef;}
  double MuRef()const{return muRef;}

  void InitializeTRef(const double t){
    tRef = t;
    muRef = cOne * pow(t,1.5)/(t+S);
  }
};

#endif

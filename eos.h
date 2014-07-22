#ifndef EOSHEADERDEF             //only if the macro EOSHEADERDEF is not defined execute these lines of code
#define EOSHEADERDEF             //define the macro

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
 sutherland() : cOne(1.458e-6), S(110.4), tRef(288.15), muRef(cOne * pow(tRef,1.5)/(tRef+S)), bulkVisc(0.0) {}
 sutherland(const double a, const double b, const double c) : cOne(a), S(b), tRef(c), muRef(cOne * pow(tRef,1.5)/(tRef+S)), bulkVisc(0.0) {}
  sutherland(const double t) : cOne(1.458e-6), S(110.4), tRef(t), muRef(cOne * pow(tRef,1.5)/(tRef+S)), bulkVisc(0.0) {}

  //member functions
  double GetViscosity(const double &)const;
  double GetLambda(const double &)const;
  double ConstC1()const{return cOne;}
  double ConstS()const{return S;}
  double TRef()const{return tRef;}
  double MuRef()const{return muRef;}
};

#endif

#include <cstdlib>      //exit()
#include "eos.h"
#include <iostream>     //cout

using std::cout;
using std::endl;
using std::cerr;

//member functions for idealGas class
double idealGas::GetPressure(const double &rho, const double &specEn)const{
  return (gamma-1.0) * rho * specEn;
}

double idealGas::GetPressFromEnergy(const double &rho, const double &energy, const double &vel)const{
  return (gamma-1.0) * rho * (energy - 0.5 * vel * vel);
}

double idealGas::GetDensity(const double &pressure, const double &specEn)const{
  return pressure / ((gamma-1.0) * specEn);
}

double idealGas::GetSpecEnergy(const double &pressure, const double &rho)const{
  return pressure / ((gamma-1.0) * rho);
}

double idealGas::GetEnergy(const double &specEn, const double &vel)const{
  return specEn + 0.5 * vel * vel;
}

double idealGas::GetEnthalpy(const double &energy, const double &pressure, const double &rho)const{
  return energy + pressure / rho;
}

double idealGas::GetSoS(const double &pressure, const double &rho)const{
  return sqrt(gamma * pressure / rho);
}


double idealGas::GetTemperature(const double &pressure, const double &rho)const{
  return pressure * gamma / rho;
}


//functions for sutherland class
double sutherland::GetViscosity(const double &t)const{

  //dimensionalize temperature
  double temp = t * tRef;

  //calculate viscosity
  double mu = (cOne * pow(temp, 1.5)) / (temp + S);

  //nondimensionalize viscosity
  return mu / muRef;
}

double sutherland::GetLambda(const double &m)const{
  //calculate lambda
  return bulkVisc - (2.0/3.0) * m;
}

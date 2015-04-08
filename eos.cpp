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

#include <cstdlib>      //exit()
#include "eos.h"
#include <iostream>     //cout

using std::cout;
using std::endl;
using std::cerr;

//member functions for idealGas class
//these functions calculate values using the ideal gas equation of state P = rho * R * T (for Navier-Stokes) or P = (g-1) * rho * e (for Euler)
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
  //calculate lambda (2nd coeff of viscosity) 
  return bulkVisc - (2.0/3.0) * m;
}

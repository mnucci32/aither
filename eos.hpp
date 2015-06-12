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

// Only if the macro EOSHEADERDEF is not defined execute these lines of code
#ifndef EOSHEADERDEF
#define EOSHEADERDEF  // define the macro

/* This header file contains the idealGas and sutherland classes. 

The ideal gas class stores the ratio of specific heats (gamma),
and the gas constant R. It contains several member functions to calculate state variables using the equation of state
P = rho * (gamma - 1) * specificEnergy for the euler equations and P = rho * R * temperature for the Navier-Stokes equations.

The sutherland class stores a reference temperature and viscosity, as well as the sutherland coefficients. It is used for calculating
a temperature dependent viscosity for the Navier-Stokes equations. */

#include <math.h>  // sqrt
#include <vector>  // vector
#include <string>  // string

using std::vector;
using std::string;

class idealGas {
  const double gamma_;
  const double gasConst_;

 public:
  // Constructor
  idealGas() : gamma_(1.4), gasConst_(287.058) {}
  idealGas(const double &a, const double &b) : gamma_(a), gasConst_(b) {}

  // Member functions
  double GetPressure(const double &rho, const double &specEn) const;
  double GetPressFromEnergy(const double &rho, const double &energy,
                            const double &vel) const;
  double GetDensity(const double &pressure, const double &specEn) const;
  double GetSpecEnergy(const double &pressure, const double &rho) const;
  double GetEnergy(const double &specEn, const double &vel) const;
  double GetEnthalpy(const double &energy, const double &pressure,
                     const double &rho) const;
  double GetSoS(const double &pressure, const double &rho) const;
  double Gamma() const {return gamma_;}
  double GasConst() const {return gasConst_;}
  double GetPrandtl() const {return (4.0*gamma_)/(9.0*gamma_-5.0);}

  double GetTemperature(const double &pressure, const double &rho) const;

  // nondimensional version (R=1/gamma_)
  double GetConductivity(const double &mu) const {
    return mu/( (*this).GetPrandtl()*(gamma_-1.0) );}
  // Nondimensional version (R=1/gamma_)
  double GetTurbConductivity(const double &eddyVisc, const double &prt) const {
    return eddyVisc/( prt * (gamma_-1.0) );}
  double GetDensityTP(const double &temp, const double &press) const {
    return press*gamma_/temp;}

  // Destructor
  ~idealGas() {}
};

// Function declarations

class sutherland {
  double cOne_;
  double S_;
  double tRef_;
  double muRef_;
  double bulkVisc_;

 public:
  // Constructors
  // Stoke's hypothesis -- bulk viscosity = 0
  // Sutherland's Law -- mu = muref * (C1 * Tref^1.5) / (T + S_)
  sutherland() : cOne_(1.458e-6), S_(110.4), tRef_(288.15),
                 muRef_(cOne_ * pow(tRef_, 1.5)/(tRef_+S_)), bulkVisc_(0.0) {}

  sutherland(const double a, const double b, const double c) :
     cOne_(a), S_(b), tRef_(c), muRef_(cOne_ * pow(c, 1.5)/(c+S_)),
     bulkVisc_(0.0) {}

  explicit sutherland(const double t) : cOne_(1.458e-6), S_(110.4),
                                        tRef_(t),
                                        muRef_(cOne_ * pow(t, 1.5)/(t+S_)),
                                        bulkVisc_(0.0) {}

  // Member functions
  double GetViscosity(const double &) const;
  double GetLambda(const double &) const;
  double ConstC1() const {return cOne_;}
  double ConstS() const {return S_;}
  double TRef() const {return tRef_;}
  double MuRef() const {return muRef_;}
};

#endif

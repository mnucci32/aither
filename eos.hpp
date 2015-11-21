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
and the gas constant R. It contains several member functions to calculate
state variables using the equation of state
P = rho * (gamma - 1) * specificEnergy for the euler equations and
P = rho * R * temperature for the Navier-Stokes equations.

The sutherland class stores a reference temperature and viscosity,
as well as the sutherland coefficients. It is used for calculating
a temperature dependent viscosity for the Navier-Stokes equations. */

#include <math.h>  // sqrt
#include <vector>  // vector
#include <string>  // string
#include "vector3d.hpp"

using std::vector;
using std::string;

class idealGas {
  const double gamma_;
  const double gasConst_;

 public:
  // Constructor
  idealGas(const double &a, const double &b) : gamma_(a), gasConst_(b) {}
  idealGas() : idealGas(1.4, 287.058) {}

  // move constructor and assignment operator
  idealGas(idealGas&&) noexcept = default;
  idealGas& operator=(idealGas&&) noexcept = default;

  // copy constructor and assignment operator
  idealGas(const idealGas&) = default;
  idealGas& operator=(const idealGas&) = default;

  // Member functions
  double Pressure(const double &rho, const double &specEn) const;
  double PressFromEnergy(const double &rho, const double &energy,
                         const double &vel) const;
  double Density(const double &pressure, const double &specEn) const;
  double SpecEnergy(const double &pressure, const double &rho) const;
  double Energy(const double &specEn, const double &vel) const;
  double Enthalpy(const double &energy, const double &pressure,
                  const double &rho) const;
  double SoS(const double &pressure, const double &rho) const;
  double Gamma() const {return gamma_;}
  double GasConst() const {return gasConst_;}
  double Prandtl() const {return (4.0 * gamma_) / (9.0 * gamma_ - 5.0);}

  double Temperature(const double &pressure, const double &rho) const;

  // nondimensional version (R=1/gamma_)
  double Conductivity(const double &mu) const {
    return mu / (this->Prandtl() * (gamma_ - 1.0) );}
  // Nondimensional version (R=1/gamma_)
  double TurbConductivity(const double &eddyVisc, const double &prt) const {
    return eddyVisc / ( prt * (gamma_ - 1.0) );}
  double DensityTP(const double &temp, const double &press) const {
    return press * gamma_ / temp;}

  // Destructor
  ~idealGas() noexcept {}
};

// Function declarations

class sutherland {
  const double cOne_;
  const double S_;
  const double tRef_;
  const double muRef_;
  const double bulkVisc_;
  const double reRef_;
  const double mRef_;
  const double scaling_;
  const double invScaling_;

 public:
  // Constructors
  // Stoke's hypothesis -- bulk viscosity = 0
  // Sutherland's Law -- mu = muref * (C1 * Tref^1.5) / (T + S_)
  sutherland(const double &c, const double &s, const double &t,
             const double &r, const double &p, const double &l,
             const vector3d<double> &vel, const idealGas &eos) :
      cOne_(c), S_(s), tRef_(t), muRef_(cOne_ * pow(tRef_, 1.5) / (tRef_ + S_)),
      bulkVisc_(0.0), reRef_(r * vel.Mag() * l / muRef_),
      mRef_(vel.Mag() / eos.SoS(p, r)), scaling_(mRef_ / reRef_),
      invScaling_(reRef_ / mRef_) {}
  sutherland(const double &t, const double &r, const double &l, const double &p,
             const vector3d<double> &vel, const idealGas &eos) :
      sutherland(1.458e-6, 110.4, t, r, p, l, vel, eos) {}

  explicit sutherland(const double &t) : cOne_(1.458e-6), S_(110.4),
                                         tRef_(t),
                                         muRef_(cOne_ * pow(t, 1.5)/(t+S_)),
                                         bulkVisc_(0.0), reRef_(0.0),
                                         mRef_(0.0), scaling_(0.0),
                                         invScaling_(0.0) {}
  sutherland() : sutherland(288.15) {}

  // move constructor and assignment operator
  sutherland(sutherland&&) noexcept = default;
  sutherland& operator=(sutherland&&) noexcept = default;

  // copy constructor and assignment operator
  sutherland(const sutherland&) = default;
  sutherland& operator=(const sutherland&) = default;

  // Member functions
  double Viscosity(const double&) const;
  double EffectiveViscosity(const double&) const;
  double Lambda(const double&) const;
  double ConstC1() const {return cOne_;}
  double ConstS() const {return S_;}
  double TRef() const {return tRef_;}
  double MuRef() const {return muRef_;}
  double ReRef() const {return reRef_;}
  double MRef() const {return mRef_;}
  double NondimScaling() const {return scaling_;}
  double InvNondimScaling() const {return invScaling_;}

  // Destructor
  ~sutherland() noexcept {}
};

#endif

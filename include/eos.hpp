/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
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
#include "vector3d.hpp"

// abstract base class for equation of state
class eos {

 public:
  // Constructor
  eos() {}

  // move constructor and assignment operator
  eos(eos&&) noexcept = default;
  eos& operator=(eos&&) noexcept = default;

  // copy constructor and assignment operator
  eos(const eos&) = default;
  eos& operator=(const eos&) = default;

  // Member functions for abstract base class
  virtual double PressFromEnergy(const double &rho, const double &energy,
                                 const double &vel) const = 0;
  virtual double PressureRT(const double &rho,
                            const double &temperature) const = 0;
  virtual double SpecEnergy(const double &pressure,
                            const double &rho) const = 0;
  virtual double Energy(const double &specEn, const double &vel) const = 0;
  virtual double Enthalpy(const double &energy, const double &pressure,
                          const double &rho) const = 0;
  virtual double SoS(const double &pressure, const double &rho) const = 0;
  virtual double Temperature(const double &pressure,
                             const double &rho) const = 0;
  virtual double DensityTP(const double &temp, const double &press) const = 0;

  // delete this later once moved to thermo model
  virtual double Gamma() const = 0;
  virtual double GasConst() const = 0;
  virtual double Prandtl() const = 0;
  virtual double SpecificHeat() const = 0;
  virtual double Conductivity(const double &) const = 0;
  virtual double TurbConductivity(const double &, const double &) const = 0;

  // Destructor
  virtual ~eos() noexcept {}
};


class idealGas : public eos {
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
  double PressFromEnergy(const double &rho, const double &energy,
                         const double &vel) const override;
  double PressureRT(const double &rho, const double &temperature) const override;
  double SpecEnergy(const double &pressure, const double &rho) const override;
  double Energy(const double &specEn, const double &vel) const override;
  double Enthalpy(const double &energy, const double &pressure,
                  const double &rho) const override;
  double SoS(const double &pressure, const double &rho) const override;
  double Temperature(const double &pressure, const double &rho) const override;
  // nondimensional version (R=1/gamma)
  double DensityTP(const double &temp, const double &press) const override {
    return press * gamma_ / temp;
  }

  double Gamma() const override {return gamma_;}
  double GasConst() const override {return gasConst_;}
  double Prandtl() const override {return (4.0 * gamma_) / (9.0 * gamma_ - 5.0);}
  double SpecificHeat() const override {return 1.0 / (gamma_ - 1.0);}
  double Conductivity(const double &mu) const override {
    return mu * this->SpecificHeat() / this->Prandtl();
  }
  double TurbConductivity(const double &eddyVisc,
                          const double &prt) const override {
    return eddyVisc * this->SpecificHeat() / prt;
  }

  // Destructor
  ~idealGas() noexcept {}
};


#endif

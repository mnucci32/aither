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

// This header file contains the equation of state classes

#include <memory>
#include "vector3d.hpp"
#include "thermodynamic.hpp"

using std::unique_ptr;

// forward class declaration
class fluid;

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
  virtual int NumSpecies() const = 0;
  virtual double GasConstant(const int &ii) const = 0;
  virtual const vector<double> &GasConstants() const = 0;
  virtual double PressFromEnergy(const unique_ptr<thermodynamic> &thermo,
                                 const vector<double> &rho,
                                 const double &energy,
                                 const double &vel) const = 0;
  virtual double PressureRT(const vector<double> &rho,
                            const double &temperature) const = 0;
  virtual double SpecEnergy(const unique_ptr<thermodynamic> &thermo,
                            const double &t,
                            const vector<double> &mf) const = 0;
  virtual double Energy(const double &specEn, const double &vel) const = 0;
  virtual double Enthalpy(const unique_ptr<thermodynamic> &thermo,
                          const double &t, const double &vel,
                          const vector<double> &mf) const = 0;
  virtual double SoS(const unique_ptr<thermodynamic> &thermo,
                     const double &pressure,
                     const vector<double> &rho) const = 0;
  virtual double Temperature(const double &pressure,
                             const vector<double> &rho) const = 0;
  virtual double DensityTP(const double &temp, const double &press) const = 0;

  // Destructor
  virtual ~eos() noexcept {}
};


// The idealGas class uses the ideal gas law to calculate state variables
// The ideal gas equation of state is P = rho * R * T. In 
// nondimensional from it is P = rho * T / gammaRef

class idealGas : public eos {
  vector<double> gasConst_;

 public:
  // Constructor
  idealGas(const vector<fluid> &, const double &, const double &);

  // Member functions
  int NumSpecies() const override { return gasConst_.size(); }
  double GasConstant(const int &ii) const override { return gasConst_[ii]; }
  const vector<double> &GasConstants() const override { return gasConst_; }
  double PressFromEnergy(const unique_ptr<thermodynamic> &thermo,
                         const vector<double> &rho, const double &energy,
                         const double &vel) const override;
  double PressureRT(const vector<double> &rho,
                    const double &temperature) const override;
  double SpecEnergy(const unique_ptr<thermodynamic> &thermo, const double &t,
                    const vector<double> &mf) const override;
  double Energy(const double &specEn, const double &vel) const override;
  double Enthalpy(const unique_ptr<thermodynamic> &thermo, const double &t,
                  const double &vel, const vector<double> &mf) const override;
  double SoS(const unique_ptr<thermodynamic> &thermo, const double &pressure,
             const vector<double> &rho) const override;
  double Temperature(const double &pressure,
                     const vector<double> &rho) const override;
  double DensityTP(const double &temp, const double &press) const override;

  // Destructor
  ~idealGas() noexcept {}
};


#endif

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

#ifndef THERMOHEADERDEF
#define THERMOHEADERDEF

// This header file contains the thermodynamic model classes
#include <iostream>
#include <cmath>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

// abstract base class for thermodynamic model
class thermodynamic {

 public:
  // Constructor
  thermodynamic() {}

  // move constructor and assignment operator
  thermodynamic(thermodynamic&&) noexcept = default;
  thermodynamic& operator=(thermodynamic&&) noexcept = default;

  // copy constructor and assignment operator
  thermodynamic(const thermodynamic&) = default;
  thermodynamic& operator=(const thermodynamic&) = default;

  // Member functions for abstract base class
  virtual double Gamma(const double &t) const = 0;
  virtual double Prandtl(const double &t) const = 0;
  virtual double Cp(const double &t) const = 0;
  virtual double Cv(const double &t) const = 0;
  virtual double SpecEnergy(const double& t) const = 0;
  virtual double SpecEnthalpy(const double& t) const = 0;
  virtual double TemperatureFromSpecEnergy(const double& e) const = 0;

  // Destructor
  virtual ~thermodynamic() noexcept {}
};

// thermodynamic model for calorically perfect gas.
// Cp and Cv are constants
class caloricallyPerfect : public thermodynamic {
  const double gamma_;

 public:
  // Constructor
  explicit caloricallyPerfect(const double &n) : gamma_(1.0 / n + 1.0) {}
  caloricallyPerfect() : caloricallyPerfect(1.4) {}

  // move constructor and assignment operator
  caloricallyPerfect(caloricallyPerfect&&) noexcept = default;
  caloricallyPerfect& operator=(caloricallyPerfect&&) noexcept = default;

  // copy constructor and assignment operator
  caloricallyPerfect(const caloricallyPerfect&) = default;
  caloricallyPerfect& operator=(const caloricallyPerfect&) = default;

  // Member functions
  double Gamma(const double &t) const override {return gamma_;}
  double Prandtl(const double& t) const override {
    return (4.0 * gamma_) / (9.0 * gamma_ - 5.0);
  }
  double Cp(const double& t) const override { return 1.0 / (gamma_ - 1.0); }
  double Cv(const double& t) const override {
    return 1.0 / (gamma_ * (gamma_ - 1.0));
  }

  double SpecEnergy(const double& t) const override {return this->Cv(t) * t;}
  double SpecEnthalpy(const double& t) const override {return this->Cp(t) * t;}
  double TemperatureFromSpecEnergy(const double& e) const override;

  // Destructor
  ~caloricallyPerfect() noexcept {}
};

// thermodynamic model for thermally perfect gas
// Cp and Cv are functions of T
class thermallyPerfect : public thermodynamic {
  const double n_;
  const vector<double> vibTemp_;
  const double nonDimR_;

  // private member functions
  double ThetaV(const double& t, const int& ii) const {
    return vibTemp_[ii] / (2.0 * t);
  }

  double VibEqCpCvTerm(const double &t) const {
    auto vibEq = 0.0;
    for (auto ii = 0U; ii < vibTemp_.size(); ++ii) {
      const auto tv = this->ThetaV(t, ii);
      vibEq += pow(tv / sinh(tv), 2.0);
    }
    return vibEq;
  }

  double VibEqTerm(const double &t) const {
    auto vibEq = 0.0;
    for (auto &vt : vibTemp_) {
      vibEq += vt / (exp(vt / t) - 1.0);
    }
    return vibEq;
  }


 public:
  // Constructor
  thermallyPerfect(const double& n, const vector<double>& vt)
      : n_(n), vibTemp_(vt), nonDimR_(n / (n + 1.0)) {}

  // Member functions
  double Gamma(const double& t) const override {
    return this->Cp(t) / this->Cv(t);
  }
  double Prandtl(const double& t) const override {
    return (4.0 * this->Gamma(t)) / (9.0 * this->Gamma(t) - 5.0);
  }
  double Cp(const double& t) const override {
    return nonDimR_ * ((n_ + 1.0) + this->VibEqCpCvTerm(t));
  }
  double Cv(const double& t) const override {
    return nonDimR_ * (n_ + this->VibEqCpCvTerm(t));
  }

  // DEBUG -- fix this for multiple vib temps
  double SpecEnergy(const double& t) const override {
    return nonDimR_ * (n_ * t + this->VibEqTerm(t));
  }
  // DEBUG -- fix this for multiple vib temps
  double SpecEnthalpy(const double& t) const override {
    return nonDimR_ * ((n_ + 1) * t + this->VibEqTerm(t));
  }
  
  double TemperatureFromSpecEnergy(const double& e) const override;

  // Destructor
  ~thermallyPerfect() noexcept {}
};

#endif

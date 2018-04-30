/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (michael.nucci@gmail.com)

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

// forward class declarations
class fluid;

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
  virtual int NumSpecies() const = 0;
  virtual const double& N(const int& ss) const = 0;
  virtual const double& R(const int& ss) const = 0;
  virtual const double& Hf(const int& ss) const = 0;
  virtual const double& S0(const int& ss) const = 0;

  double SpeciesGamma(const double& t, const int& ss) const {
    return this->SpeciesCp(t, ss) / this->SpeciesCv(t, ss);
  }
  double Gamma(const double& t, const vector<double>& mf) const {
    return this->Cp(t, mf) / this->Cv(t, mf);
  }
  double Prandtl(const double& t, const vector<double>& mf) const {
    const auto gamma = this->Gamma(t, mf);
    return (4.0 * gamma) / (9.0 * gamma - 5.0);
  }
  double Cp(const double& t, const vector<double>& mf) const;
  double Cv(const double& t, const vector<double>& mf) const;
  virtual double SpeciesSpecEnergy(const double& t, const int& ss) const = 0;
  double SpecEnergy(const double& t, const vector<double>& mf) const;
  virtual double SpeciesSpecEnthalpy(const double& t, const int& ss) const = 0;
  double SpecEnthalpy(const double& t, const vector<double>& mf) const;
  virtual double TemperatureFromSpecEnergy(const double& e,
                                           const vector<double>& mf) const = 0;
  virtual double SpeciesCp(const double& t, const int& ss) const = 0;
  virtual double SpeciesCv(const double& t, const int& ss) const = 0;
  virtual double SpeciesGibbsMinStdState(const double& t,
                                         const int& ss) const = 0;
  virtual vector<double> GibbsMinimization(const double& t) const = 0;

  // Destructor
  virtual ~thermodynamic() noexcept {}
};

// thermodynamic model for calorically perfect gas.
// Cp and Cv are constants
class caloricallyPerfect : public thermodynamic {
  vector<double> n_;
  vector<double> gasConst_;
  vector<double> hf_;
  vector<double> s0_;

 public:
  // Constructor
  caloricallyPerfect(const vector<fluid>&, const double&, const double&);

  // Member functions
  int NumSpecies() const override { return n_.size(); }
  const double& N(const int &ss) const override { return n_[ss]; }
  const double& R(const int &ss) const override { return gasConst_[ss]; }
  const double& Hf(const int& ss) const override { return hf_[ss]; }
  const double& S0(const int& ss) const override { return s0_[ss]; }
  void SubtractS0(const int& ss, const double& val) { s0_[ss] -= val; }
  double SpeciesSpecEnergy(const double& t, const int& ss) const override {
    return hf_[ss] + this->SpeciesCv(t, ss) * t;
  }
  double SpeciesSpecEnthalpy(const double& t, const int& ss) const override {
    return hf_[ss] + this->SpeciesCp(t, ss) * t;
  }
  double TemperatureFromSpecEnergy(const double& e,
                                   const vector<double>& mf) const override;
  double SpeciesCp(const double& t, const int& ss) const override {
    return this->R(ss) * (this->N(ss) + 1.0);
  }
  double SpeciesCv(const double& t, const int& ss) const override {
    return this->R(ss) * this->N(ss);
  }
  double SpeciesGibbsMinStdState(const double& t, const int& ss) const override;
  vector<double> GibbsMinimization(const double& t) const override;

  // Destructor
  ~caloricallyPerfect() noexcept {}
};

// thermodynamic model for thermally perfect gas
// Cp and Cv are functions of T
class thermallyPerfect : public caloricallyPerfect {
  vector<vector<double>> vibTemp_;

  // private member functions
  double ThetaV(const double& t, const int &ss, const int& ii) const {
    return vibTemp_[ss][ii] / (2.0 * t);
  }

  double VibEqCpCvTerm(const double &t, const int &ss) const {
    auto vibEq = 0.0;
    for (auto ii = 0U; ii < vibTemp_[ss].size(); ++ii) {
      const auto tv = this->ThetaV(t, ss, ii);
      vibEq += pow(tv / sinh(tv), 2.0);
    }
    return vibEq;
  }

  double VibEqTerm(const double &t, const int &ss) const {
    auto vibEq = 0.0;
    for (auto &vt : vibTemp_[ss]) {
      vibEq += vt / (exp(vt / t) - 1.0);
    }
    return vibEq;
  }

  double GibbsVibTerm(const double &t, const int &ss) const {
    auto gibbsVib = 0.0;
    for (auto &vt : vibTemp_[ss]) {
      gibbsVib += log(1.0 - exp(-vt / t));
    }
    return gibbsVib * t;
  }

 public:
  // Constructor
  thermallyPerfect(const vector<fluid>& fl, const double& tRef,
                   const double& aRef);

  // Member functions
  double SpeciesSpecEnergy(const double& t, const int& ss) const override {
    return this->Hf(ss) + caloricallyPerfect::SpeciesCv(t, ss) * t + 
           this->R(ss) * this->VibEqTerm(t, ss);
  }
  double SpeciesSpecEnthalpy(const double& t, const int& ss) const override {
    return this->Hf(ss) + caloricallyPerfect::SpeciesCp(t, ss) * t + 
           this->R(ss) * this->VibEqTerm(t, ss);
  }
  double TemperatureFromSpecEnergy(const double& e,
                                   const vector<double>& mf) const override;
  double SpeciesCp(const double& t, const int& ss) const override {
    return caloricallyPerfect::SpeciesCp(t, ss) +
           this->R(ss) * this->VibEqCpCvTerm(t, ss);
  }
  double SpeciesCv(const double& t, const int& ss) const override {
    return caloricallyPerfect::SpeciesCv(t, ss) +
           this->R(ss) * this->VibEqCpCvTerm(t, ss);
  }
  double SpeciesGibbsMinStdState(const double& t,
                                 const int& ss) const override {
    return caloricallyPerfect::SpeciesGibbsMinStdState(t, ss) +
           this->R(ss) * this->GibbsVibTerm(t, ss);
  }

  // Destructor
  ~thermallyPerfect() noexcept {}
};

#endif

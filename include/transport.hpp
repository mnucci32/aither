/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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

#ifndef TRANSPORTHEADERDEF
#define TRANSPORTHEADERDEF

// This header file contains the classes for all viscous transport models

#include <math.h>  // sqrt
#include <memory>
#include <vector>
#include "vector3d.hpp"
#include "thermodynamic.hpp"

using std::vector;
using std::unique_ptr;

// forward class declarations
class fluid;

// abstract base class for transport models
class transport {
  double scaling_;
  double invScaling_;

 protected:
  void SetScaling(const double &rho, const double &l, const double &mu,
                  const double &a) {
    scaling_ = mu / (rho * a * l);
    invScaling_ = 1.0 / scaling_;
  }

 public:
  // Constructor
  transport() : scaling_(0), invScaling_(0) {}
  transport(const double &rho, const double &l, const double &mu,
            const double &a) {
    this->SetScaling(rho, l, mu, a);
  }

  // move constructor and assignment operator
  transport(transport&&) noexcept = default;
  transport& operator=(transport&&) noexcept = default;

  // copy constructor and assignment operator
  transport(const transport&) = default;
  transport& operator=(const transport&) = default;

  // Member functions for abstract base class
  virtual int NumSpecies() const = 0;
  virtual double SpeciesViscosity(const double &, const int &) const = 0;
  virtual double SpeciesConductivity(const double &, const int &) const = 0;
  virtual double Viscosity(const double &, const vector<double> &) const = 0;
  virtual double EffectiveViscosity(const double &,
                                    const vector<double> &) const = 0;
  virtual double Lambda(const double &) const = 0;
  virtual double TRef() const = 0;
  virtual double MuRef() const = 0;
  virtual double Conductivity(const double &, const vector<double> &) const = 0;
  virtual double EffectiveConductivity(const double &,
                                       const vector<double> &) const = 0;
  virtual double TurbConductivity(const double &, const double &,
                                  const double &,
                                  const unique_ptr<thermodynamic> &,
                                  const vector<double> &mf) const = 0;
  double NondimScaling() const { return scaling_; }
  double InvNondimScaling() const {return invScaling_;}
  virtual vector<double> MoleFractions(const vector<double> &) const = 0;

  // Destructor
  virtual ~transport() noexcept {}
};


// this class models viscous transport using Sutherland's law
class sutherland : public transport {
  vector<double> viscC1_;
  vector<double> viscS_;
  vector<double> condC1_;
  vector<double> condS_;
  vector<double> molarMass_;
  double tRef_;
  double muMixRef_;
  double kNonDim_;
  double bulkVisc_ = 0.0;

  // private member functions
  double WilkesVisc(const vector<double> &, const vector<double> &) const;
  double WilkesCond(const vector<double> &, const vector<double> &) const;

 public:
  // Constructors
  // Stoke's hypothesis -- bulk viscosity = 0
  // Sutherland's Law -- mu = muref * (C1 * Tref^1.5) / (T + S_)
  sutherland(const vector<fluid> &, const double &, const double &,
             const double &, const double &, const vector<double> &);

  // move constructor and assignment operator
  sutherland(sutherland&&) noexcept = default;
  sutherland& operator=(sutherland&&) noexcept = default;

  // copy constructor and assignment operator
  sutherland(const sutherland&) = default;
  sutherland& operator=(const sutherland&) = default;

  // Member functions
  int NumSpecies() const override { return molarMass_.size(); }
  double SpeciesViscosity(const double &, const int &) const override;
  double SpeciesConductivity(const double &, const int &) const override;
  double Viscosity(const double &, const vector<double> &) const override;
  double EffectiveViscosity(const double &,
                            const vector<double> &) const override;
  double Lambda(const double &) const override;
  double TRef() const override {return tRef_;}
  double MuRef() const override {return muMixRef_;}
  double Conductivity(const double &, const vector<double> &) const override;
  double EffectiveConductivity(const double &,
                               const vector<double> &) const override;
  double TurbConductivity(const double &eddyVisc, const double &prt,
                          const double &t,
                          const unique_ptr<thermodynamic> &thermo,
                          const vector<double> &mf) const override {
    return eddyVisc * thermo->Cp(t, mf) / prt;
  }
  vector<double> MoleFractions(const vector<double> &) const override;

  // Destructor
  ~sutherland() noexcept {}
};


// --------------------------------------------------------------------------
// function declarations



#endif

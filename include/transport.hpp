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

#ifndef TRANSPORTHEADERDEF
#define TRANSPORTHEADERDEF

// This header file contains the classes for all viscous transport models

#include <math.h>  // sqrt
#include <memory>
#include "vector3d.hpp"
#include "eos.hpp"
#include "thermodynamic.hpp"

using std::unique_ptr;

// abstract base class for transport models
class transport {
  const double reRef_;
  const double mRef_;
  const double scaling_;
  const double invScaling_;

 public:
  // Constructor
  transport(const double &re, const double &m)
      : reRef_(re), mRef_(m), scaling_(m / re), invScaling_(re / m) {}

  // move constructor and assignment operator
  transport(transport&&) noexcept = default;
  transport& operator=(transport&&) noexcept = default;

  // copy constructor and assignment operator
  transport(const transport&) = default;
  transport& operator=(const transport&) = default;

  // Member functions for abstract base class
  virtual double Viscosity(const double&) const = 0;
  virtual double EffectiveViscosity(const double&) const = 0;
  virtual double Lambda(const double&) const = 0;
  virtual double ConstC1() const = 0;
  virtual double ConstS() const = 0;
  virtual double TRef() const = 0;
  virtual double MuRef() const = 0;
  virtual double Conductivity(const double &,
                              const unique_ptr<thermodynamic> &) const = 0;
  virtual double TurbConductivity(const double &, const double &,
                                  const unique_ptr<thermodynamic> &) const = 0;

  double ReRef() const {return reRef_;}
  double MRef() const {return mRef_;}
  double NondimScaling() const {return scaling_;}
  double InvNondimScaling() const {return invScaling_;}

  // Destructor
  virtual ~transport() noexcept {}
};


// this class models viscous transport using Sutherland's law
class sutherland : public transport {
  const double cOne_;
  const double S_;
  const double tRef_;
  const double muRef_;
  const double bulkVisc_;

 public:
  // Constructors
  // Stoke's hypothesis -- bulk viscosity = 0
  // Sutherland's Law -- mu = muref * (C1 * Tref^1.5) / (T + S_)
  sutherland(const double &c, const double &s, const double &t, const double &r,
             const double &p, const double &l, const vector3d<double> &vel,
             const unique_ptr<eos> &eqnState)
      : transport(r * vel.Mag() * l / (c * pow(t, 1.5) / (t + s)),
                  vel.Mag() / eqnState->SoS(p, r)),
        cOne_(c),
        S_(s),
        tRef_(t),
        muRef_(cOne_ * pow(tRef_, 1.5) / (tRef_ + S_)),
        bulkVisc_(0.0) {}
  sutherland(const double &t, const double &r, const double &l, const double &p,
             const vector3d<double> &vel, const unique_ptr<eos> &eqnState)
      : sutherland(1.458e-6, 110.4, t, r, p, l, vel, eqnState) {}

  // move constructor and assignment operator
  sutherland(sutherland&&) noexcept = default;
  sutherland& operator=(sutherland&&) noexcept = default;

  // copy constructor and assignment operator
  sutherland(const sutherland&) = default;
  sutherland& operator=(const sutherland&) = default;

  // Member functions
  double Viscosity(const double&) const override;
  double EffectiveViscosity(const double&) const override;
  double Lambda(const double&) const override;
  double ConstC1() const override {return cOne_;}
  double ConstS() const override {return S_;}
  double TRef() const override {return tRef_;}
  double MuRef() const override {return muRef_;}
  double Conductivity(const double &mu,
                      const unique_ptr<thermodynamic> &thermo) const override {
    return mu * thermo->SpecificHeat() / thermo->Prandtl();
  }
  double TurbConductivity(
      const double &eddyVisc, const double &prt,
      const unique_ptr<thermodynamic> &thermo) const override {
    return eddyVisc * thermo->SpecificHeat() / prt;
  }

  // Destructor
  ~sutherland() noexcept {}
};

#endif

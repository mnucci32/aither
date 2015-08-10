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

#ifndef TURBHEADERDEF  // only if the macro TURBHEADERDEF is not defined execute
                       // these lines of code
#define TURBHEADERDEF  // define the macro

/* This defines the classes for the turbulence  models implemented in the code.
 */

#include <math.h>  // sqrt
#include <vector>  // vector
#include <string>  // string
#include "vector3d.hpp"  // vector3d
#include "tensor.hpp"  // tensor

using std::vector;
using std::string;

// forward class declaration
class primVars;
class sutherland;

class turbModel {
  string eddyViscMethod_;

 public:
  // constructor
  turbModel() : eddyViscMethod_("boussinesq") {}
  explicit turbModel(const string &meth) : eddyViscMethod_(meth) {}

  // member functions
  string EddyViscMethod() const {return eddyViscMethod_;}

  virtual double EddyVisc(const primVars &state,
                          const tensor<double> &vGrad,
                          const sutherland &suth) const = 0;
  virtual double EddyViscNoLim(const primVars &state) const = 0;
  virtual double TurbPrandtlNumber() const = 0;
  virtual double Production1(const primVars &state,
                             const tensor<double> &velGrad,
                             const sutherland &suth) const = 0;
  virtual double Production2(const primVars &state,
                             const tensor<double> &velGrad,
                             const sutherland &suth) const = 0;
  virtual double Destruction1(const primVars &state) const = 0;
  virtual double Destruction2(const primVars &state,
                              const tensor<double> &velGrad,
                              const sutherland &suth) const = 0;
  virtual double CrossDiff2(const primVars &state,
                            const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const = 0;
  virtual double MolecDiff1Coeff() const = 0;
  virtual double MolecDiff2Coeff() const = 0;

  virtual void Print() const = 0;

  tensor<double> BoussinesqReynoldsStress(const primVars &state,
                                          const tensor<double> &velGrad,
                                          const sutherland &suth) const;

  // destructor
  virtual ~turbModel() {}
};

class turbNone : public turbModel {
 public:
  // constructor
  turbNone() : turbModel() {}
  explicit turbNone(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars &state,
                  const tensor<double> &vGrad,
                  const sutherland &suth) const override {return 0.0;}
  double EddyViscNoLim(const primVars &state) const override {return 0.0;}
  double TurbPrandtlNumber() const override {return 0.9;}
  double Production1(const primVars &state, const tensor<double> &velGrad,
                     const sutherland &suth) const override {return 0.0;}
  double Production2(const primVars &state, const tensor<double> &velGrad,
                     const sutherland &suth) const override {return 0.0;}
  double Destruction1(const primVars &state) const override {return 0.0;}
  double Destruction2(const primVars &state,
                      const tensor<double> &velGrad,
                      const sutherland &suth) const override {
    return 0.0;
  }
  double CrossDiff2(const primVars &state, const vector3d<double> &kGrad,
                    const vector3d<double> &wGrad) const override {
    return 0.0;
  }
  double MolecDiff1Coeff() const override {return 0.0;}
  double MolecDiff2Coeff() const override {return 0.0;}

  void Print() const override;

  // destructor
  ~turbNone() {}
};

class turbKWWilcox : public turbModel {
  const double alpha_ = 0.52;
  const double betaStar_ = 0.09;
  const double sigma_ = 0.5;
  const double sigmaStar_ = 0.6;
  const double sigmaD0_ = 0.125;
  const double beta0_ = 0.0708;
  const double clim_ = 0.875;
  const double prt_ = 8.0 / 9.0;
  // These aren't used yet
  const double ks_ = 1.0e-6;
  const double minK_ = 1.0e-15;
  const double minW_ = 1.0e-15;

 public:
  // constructor
  turbKWWilcox() : turbModel() {}
  explicit turbKWWilcox(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars&, const tensor<double> &,
                  const sutherland &) const override;
  double EddyViscNoLim(const primVars&) const override;
  double TurbPrandtlNumber() const override {return prt_;}
  double Production1(const primVars &, const tensor<double> &,
                     const sutherland &) const override;
  double Production2(const primVars &, const tensor<double> &,
                     const sutherland &) const override;
  double Destruction1(const primVars &) const override;
  double Destruction2(const primVars &, const tensor<double> &,
                      const sutherland &) const override;
  double CrossDiff2(const primVars &, const vector3d<double> &,
                    const vector3d<double> &) const override;
  double MolecDiff1Coeff() const override {return sigmaStar_;}
  double MolecDiff2Coeff() const override {return sigma_;}

  double Alpha() const {return alpha_;}
  double BetaStar() const {return betaStar_;}
  double Sigma() const {return sigma_;}
  double SigmaStar() const {return sigmaStar_;}
  double SigmaD0() const {return sigmaD0_;}
  double Beta0() const {return beta0_;}
  double CLim() const {return clim_;}

  double SigmaD(const vector3d<double>&, const vector3d<double>&) const;
  double Xw(const primVars&, const tensor<double>&, const sutherland&) const;
  double FBeta(const primVars&, const tensor<double>&, const sutherland&) const;
  double Beta(const primVars&, const tensor<double>&, const sutherland&) const;
  tensor<double> StrainKI(const tensor<double>&) const;
  double OmegaTilda(const primVars&, const tensor<double>&,
                    const sutherland&) const;

  void Print() const override;

  // destructor
  ~turbKWWilcox() {}
};

// function declarations

#endif

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

class turbModel {
  string eddyViscMethod_;

 public:
  // constructor
  turbModel() : eddyViscMethod_("boussinesq") {}
  explicit turbModel(const string &meth) : eddyViscMethod_(meth) {}

  // member functions
  string EddyViscMethod() const {return eddyViscMethod_;}

  virtual double EddyVisc(const primVars &state,
                          const tensor<double> &vGrad) const = 0;
  virtual double EddyViscNoLim(const primVars &state) const = 0;
  virtual double TurbPrandtlNumber() const = 0;
  virtual double Eqn1ProductionCoeff() const = 0;
  virtual double Eqn2ProductionCoeff() const = 0;
  virtual double Eqn1DissipationCoeff() const = 0;
  virtual double Eqn2DissipationCoeff(const primVars &state,
                                      const tensor<double> &velGrad) const = 0;
  virtual double Eqn2CrossDiffCoeff(const vector3d<double> &kGrad,
                                    const vector3d<double> &wGrad) const = 0;
  virtual double Eqn1MolecDiffCoeff() const = 0;
  virtual double Eqn2MolecDiffCoeff() const = 0;

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
                  const tensor<double> &vGrad) const override {return 0.0;}
  double EddyViscNoLim(const primVars &state) const override {return 0.0;}
  double TurbPrandtlNumber() const override {return 0.9;}
  double Eqn1ProductionCoeff() const override {return 0.0;}
  double Eqn2ProductionCoeff() const override {return 0.0;}
  double Eqn1DissipationCoeff() const override {return 0.0;}
  double Eqn2DissipationCoeff(const primVars &state,
                              const tensor<double> &velGrad) const override {
    return 0.0;
  }
  double Eqn2CrossDiffCoeff(const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const override {
    return 0.0;
  }
  double Eqn1MolecDiffCoeff() const override {return 0.0;}
  double Eqn2MolecDiffCoeff() const override {return 0.0;}


  // destructor
  ~turbNone() {}
};

class turbKWWilcox : public turbModel {
 public:
  // constructor
  turbKWWilcox() : turbModel() {}
  explicit turbKWWilcox(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars&, const tensor<double> &) const override;
  double EddyViscNoLim(const primVars&) const override;
  double TurbPrandtlNumber() const override {return 8.0 / 9.0;}
  double Eqn1ProductionCoeff() const override {return 1.0;}
  double Eqn2ProductionCoeff() const override {return (*this).Alpha();}
  double Eqn1DissipationCoeff() const override {return (*this).BetaStar();}
  double Eqn2DissipationCoeff(const primVars &state,
                              const tensor<double> &velGrad) const override {
    return (*this).Beta(state, velGrad);
  }
  double Eqn2CrossDiffCoeff(const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const override {
    return (*this).SigmaD(kGrad, wGrad);
  }
  double Eqn1MolecDiffCoeff() const override {return (*this).SigmaStar();}
  double Eqn2MolecDiffCoeff() const override {return (*this).Sigma();}

  double Alpha() const {return 0.52;}
  double BetaStar() const {return 0.09;}
  double Sigma() const {return 0.5;}
  double SigmaStar() const {return 0.6;}
  double SigmaD0() const {return 0.125;}
  double Beta0() const {return 0.0708;}
  double CLim() const {return 0.875;}

  double SigmaD(const vector3d<double>&, const vector3d<double>&) const;
  double Xw(const primVars&, const tensor<double>&) const;
  double FBeta(const primVars&, const tensor<double>&) const;
  double Beta(const primVars&, const tensor<double>&) const;
  tensor<double> StrainKI(const tensor<double>&) const;
  double OmegaTilda(const primVars&, const tensor<double>&) const;

  // destructor
  ~turbKWWilcox() {}
};

// function declarations

#endif

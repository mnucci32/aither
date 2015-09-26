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
class idealGas;

class turbModel {
  string eddyViscMethod_;

 public:
  // constructor
  turbModel() : eddyViscMethod_("boussinesq") {}
  explicit turbModel(const string &meth) : eddyViscMethod_(meth) {}

  // member functions
  string EddyViscMethod() const {return eddyViscMethod_;}
  virtual double EddyViscNoLim(const primVars &state) const;
  virtual double TurbPrandtlNumber() const {return 0.9;}
  tensor<double> BoussinesqReynoldsStress(const primVars &state,
                                          const tensor<double> &velGrad,
                                          const sutherland &suth,
                                          const idealGas &eos,
                                          const double &wallDist) const;
  double ReynoldsStressDDotVelGrad(const primVars &state,
                                   const tensor<double> &velGrad,
                                   const sutherland &suth,
                                   const idealGas &eos,
                                   const double &wallDist) const;
  double TkeDestruction(const primVars &state) const;
  double OmegaDestruction(const primVars &state) const;
  double CrossDiffusion(const primVars &state,
                        const vector3d<double> &kGrad,
                        const vector3d<double> &wGrad) const;

  // abstract functions
  virtual double EddyVisc(const primVars &state,
                          const tensor<double> &vGrad,
                          const sutherland &suth, const idealGas &eos,
                          const double &wallDist) const = 0;
  virtual double Production1(const primVars &state,
                             const tensor<double> &velGrad,
                             const sutherland &suth, const idealGas &eos,
                             const double &wallDist) const = 0;
  virtual double Destruction1(const primVars &state) const = 0;
  virtual double Production2(const primVars &state,
                             const tensor<double> &velGrad,
                             const vector3d<double> &kGrad,
                             const vector3d<double> &wGrad,
                             const sutherland &suth, const idealGas &eos,
                             const double &wallDist) const = 0;
  virtual double Destruction2(const primVars &state,
                              const tensor<double> &velGrad,
                              const vector3d<double> &kGrad,
                              const vector3d<double> &wGrad,
                              const sutherland &suth, const idealGas &eos,
                              const double &wallDist) const = 0;
  virtual double CrossDiff2(const primVars &state,
                            const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad,
                            const sutherland &suth, const idealGas &eos,
                            const double &walldist) const = 0;
  virtual double MolecDiff1Coeff(const primVars &state,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad,
                                 const sutherland &suth, const idealGas &eos,
                                 const double &walldist) const = 0;
  virtual double MolecDiff2Coeff(const primVars &state,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad,
                                 const sutherland &suth, const idealGas &eos,
                                 const double &walldist) const = 0;

  virtual void Print() const = 0;

  // destructor
  virtual ~turbModel() {}
};

class turbNone : public turbModel {
 public:
  // constructor
  turbNone() : turbModel() {}
  explicit turbNone(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars &state,  const tensor<double> &vGrad,
                  const sutherland &suth, const idealGas &eos,
                  const double &wallDist) const override {return 0.0;}
  double EddyViscNoLim(const primVars &state) const override {return 0.0;}
  double Production1(const primVars &state, const tensor<double> &velGrad,
                     const sutherland &suth, const idealGas &eos,
                     const double &wallDist) const override {return 0.0;}
  double Destruction1(const primVars &state) const override {return 0.0;}
  double Production2(const primVars &state, const tensor<double> &velGrad,
                     const vector3d<double> &kGrad,
                     const vector3d<double> &wGrad,
                     const sutherland &suth, const idealGas &eos,
                     const double &wallDist) const override {return 0.0;}
  double Destruction2(const primVars &state, const tensor<double> &velGrad,
                      const vector3d<double> &kGrad,
                      const vector3d<double> &wGrad, const sutherland &suth,
                      const idealGas &eos, const double &wallDist)
      const override {return 0.0;}
  double CrossDiff2(const primVars &state, const vector3d<double> &kGrad,
                    const vector3d<double> &wGrad, const sutherland &suth,
                    const idealGas &eos, const double &wallDist)
      const override {return 0.0;}
  double MolecDiff1Coeff(const primVars &state, const vector3d<double> &kGrad,
                         const vector3d<double> &wGrad,
                         const sutherland &suth, const idealGas &eos,
                         const double &walldist) const override {return 0.0;}
  double MolecDiff2Coeff(const primVars &state, const vector3d<double> &kGrad,
                         const vector3d<double> &wGrad,
                         const sutherland &suth, const idealGas &eos,
                         const double &walldist) const override {return 0.0;}

  void Print() const override;

  // destructor
  ~turbNone() {}
};

class turbKWWilcox : public turbModel {
  const double gamma_ = 0.52;
  const double betaStar_ = 0.09;
  const double sigma_ = 0.5;
  const double sigmaStar_ = 0.6;
  const double sigmaD0_ = 0.125;
  const double beta0_ = 0.0708;
  const double clim_ = 0.875;
  const double prt_ = 8.0 / 9.0;

  // private member functions
  double SigmaD(const vector3d<double>&, const vector3d<double>&) const;
  double Xw(const primVars&, const tensor<double>&, const sutherland&) const;
  double FBeta(const primVars&, const tensor<double>&, const sutherland&) const;
  double Beta(const primVars&, const tensor<double>&, const sutherland&) const;
  tensor<double> StrainKI(const tensor<double>&) const;
  double OmegaTilda(const primVars&, const tensor<double>&,
                    const sutherland&) const;

 public:
  // constructor
  turbKWWilcox() : turbModel() {}
  explicit turbKWWilcox(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars&, const tensor<double> &,
                  const sutherland &, const idealGas &,
                  const double &) const override;
  double TurbPrandtlNumber() const override {return prt_;}
  double Production1(const primVars &state, const tensor<double> &velGrad,
                     const sutherland &suth, const idealGas &eos,
                     const double &wallDist) const override {
    return this->ReynoldsStressDDotVelGrad(state, velGrad, suth, eos, wallDist);
  }
  double Destruction1(const primVars &state) const override {
    return betaStar_ * this->TkeDestruction(state);
  }
  double Production2(const primVars &, const tensor<double> &,
                     const vector3d<double> &, const vector3d<double> &,
                     const sutherland &, const idealGas &,
                     const double &) const override;
  double Destruction2(const primVars &, const tensor<double> &,
                      const vector3d<double> &, const vector3d<double> &,
                      const sutherland &, const idealGas &,
                      const double &) const override;
  double CrossDiff2(const primVars &, const vector3d<double> &,
                    const vector3d<double> &, const sutherland &,
                    const idealGas &, const double &) const override;
  double MolecDiff1Coeff(const primVars &state, const vector3d<double> &kGrad,
                         const vector3d<double> &wGrad,
                         const sutherland &suth, const idealGas &eos,
                         const double &walldist) const override {
    return sigmaStar_;
  }
  double MolecDiff2Coeff(const primVars &state, const vector3d<double> &kGrad,
                         const vector3d<double> &wGrad,
                         const sutherland &suth, const idealGas &eos,
                         const double &walldist) const override {
    return sigma_;
  }

  double Gamma() const {return gamma_;}
  double BetaStar() const {return betaStar_;}
  double Sigma() const {return sigma_;}
  double SigmaStar() const {return sigmaStar_;}
  double SigmaD0() const {return sigmaD0_;}
  double Beta0() const {return beta0_;}
  double CLim() const {return clim_;}

  void Print() const override;

  // destructor
  ~turbKWWilcox() {}
};

class turbKWSst : public turbModel {
  const double betaStar_ = 0.09;
  const double sigmaK1_ = 0.85;
  const double sigmaK2_ = 1.0;
  const double sigmaW1_ = 0.5;
  const double sigmaW2_ = 0.856;
  const double beta1_ = 0.075;
  const double beta2_ = 0.0828;
  const double gamma1_ = 5.0 / 9.0;
  const double gamma2_ = 0.44;
  const double a1_ = 0.31;
  const double prt_ = 0.9;
  const double kProd2Dest_ = 10.0;

  // private member functions
  double CDkw(const primVars &, const vector3d<double> &,
              const vector3d<double> &) const;
  double F1(const primVars &, const vector3d<double> &,
            const vector3d<double> &, const sutherland &, const idealGas&,
            const double &) const;
  double F2(const primVars &, const sutherland &, const idealGas &,
            const double &) const;
  double BlendedCoeff(const double &, const double &, const primVars &,
                      const vector3d<double> &, const vector3d<double> &,
                      const sutherland &, const idealGas &,
                      const double &) const;

 public:
  // constructor
  turbKWSst() : turbModel() {}
  explicit turbKWSst(const string &meth) : turbModel(meth) {}

  // member functions
  double EddyVisc(const primVars &, const tensor<double> &,
                  const sutherland &, const idealGas &,
                  const double &) const override;
  double Production1(const primVars &, const tensor<double> &,
                     const sutherland &, const idealGas &,
                     const double &) const override;
  double Destruction1(const primVars &state) const override {
    return betaStar_ * this->TkeDestruction(state);
  }
  double Production2(const primVars &, const tensor<double> &,
                     const vector3d<double> &, const vector3d<double> &,
                     const sutherland &, const idealGas &,
                     const double &) const override;
  double Destruction2(const primVars &, const tensor<double> &,
                      const vector3d<double> &, const vector3d<double> &,
                      const sutherland &, const idealGas &,
                      const double &) const override;

  double CrossDiff2(const primVars &, const vector3d<double> &,
                    const vector3d<double> &, const sutherland &,
                    const idealGas &, const double &) const override;

  double MolecDiff1Coeff(const primVars &, const vector3d<double> &,
                         const vector3d<double> &, const sutherland &,
                         const idealGas &, const double &) const override;
  double MolecDiff2Coeff(const primVars &, const vector3d<double> &,
                         const vector3d<double> &, const sutherland &,
                         const idealGas &, const double &) const override;

  double Gamma1() const {return gamma1_;}
  double Beta1() const {return beta1_;}
  double SigmaK1() const {return sigmaK1_;}
  double SigmaW1() const {return sigmaW1_;}

  double Gamma2() const {return gamma2_;}
  double Beta2() const {return beta2_;}
  double SigmaK2() const {return sigmaK2_;}
  double SigmaW2() const {return sigmaW2_;}

  double A1() const {return a1_;}
  double BetaStar() const {return betaStar_;}

  void Print() const override;

  // destructor
  ~turbKWSst() {}
};


// function declarations

#endif

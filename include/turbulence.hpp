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

#ifndef TURBHEADERDEF  // only if the macro TURBHEADERDEF is not defined execute
                       // these lines of code
#define TURBHEADERDEF  // define the macro

/* This defines the classes for the turbulence  models implemented in the code.
 */

#include <math.h>  // sqrt
#include <string>  // string
#include <memory>
#include <algorithm>
#include "vector3d.hpp"  // vector3d
#include "tensor.hpp"  // tensor
#include "arrayView.hpp"  // primitiveView

using std::string;
using std::unique_ptr;

// forward class declaration
class primitive;
class transport;
class squareMatrix;

class turbModel {
  const string eddyViscMethod_;

 public:
  // constructor
  turbModel() : eddyViscMethod_("boussinesq") {}
  explicit turbModel(const string &meth) : eddyViscMethod_(meth) {}

  // move constructor and assignment operator
  turbModel(turbModel&&) = default;
  turbModel& operator=(turbModel&&) = default;

  // copy constructor and assignment operator
  turbModel(const turbModel&) = default;
  turbModel& operator=(const turbModel&) = default;

  // member functions
  string EddyViscMethod() const {return eddyViscMethod_;}
  tensor<double> MeanStrainRate(const tensor<double> &) const;
  virtual double EddyViscNoLim(const primitive &state) const;
  virtual double EddyViscNoLim(const primitiveView &state) const;
  virtual double TurbPrandtlNumber() const {return 0.9;}
  virtual double TkeMin() const {return 1.0e-20;}
  virtual double OmegaMin() const {return 1.0e-20;}
  virtual double TurbMinN(const int &ii) const { 
    MSG_ASSERT(ii == 0 || ii ==1, "turbulence index out of range");
    return ii == 0 ? this->TkeMin() : this->OmegaMin();
  }
  virtual double SigmaK(const double &f1) const {return 0.0;}
  virtual double SigmaW(const double &f1) const {return 0.0;}
  virtual double WallSigmaK() const {return 0.0;}
  virtual double WallSigmaW() const {return 0.0;}
  virtual bool UseUnlimitedEddyVisc() const {return false;}
  virtual bool UsePhi() const {return false;}
  virtual double EddyVisc(const primitive &state,
                          const tensor<double> &vGrad,
                          const unique_ptr<transport> &trans,
                          const double &f2,
                          const double &length) const {return 0.0;}
  virtual double WallBeta() const {return 1.0;}
  virtual double BetaStar() const {return 0.0;}
  virtual double SrcSpecRad(const primitive &state,
                            const unique_ptr<transport> &trans, const double &vol,
                            const double &phi = 1.0) const {return 0.0;}
  virtual squareMatrix InviscidJacobian(const primitive &state,
                                        const unitVec3dMag<double> &fArea,
                                        const bool &positive) const;
  virtual squareMatrix InviscidConvJacobian(
      const primitive &state, const unitVec3dMag<double> &fArea) const;
  virtual squareMatrix InviscidDissJacobian(
      const primitive &state, const unitVec3dMag<double> &fArea) const;
  virtual double InviscidCellSpecRad(const primitive &state,
                                     const unitVec3dMag<double> &fAreaL,
                                     const unitVec3dMag<double> &fAreaR) const;
  virtual double InviscidFaceSpecRad(const primitive &state,
                                     const unitVec3dMag<double> &fArea,
                                     const bool &positive) const;
  virtual squareMatrix ViscousJacobian(const primitive &state,
                                       const unitVec3dMag<double> &fArea,
                                       const double &mu,
                                       const unique_ptr<transport> &trans,
                                       const double &dist, const double &mut,
                                       const double &f1) const;
  virtual double ViscCellSpecRad(const primitive &state,
                                 const unitVec3dMag<double> &fAreaL,
                                 const unitVec3dMag<double> &fAreaR,
                                 const double &mu,
                                 const unique_ptr<transport> &trans,
                                 const double &vol, const double &mut,
                                 const double &f1) const {
    return 0.0;
  }
  virtual double ViscFaceSpecRad(const primitive &state,
                                 const unitVec3dMag<double> &fArea,
                                 const double &mu,
                                 const unique_ptr<transport> &trans,
                                 const double &dist, const double &mut,
                                 const double &f1) const {
    return 0.0;
  }
  tensor<double> BoussinesqReynoldsStress(const primitive &state,
                                          const tensor<double> &velGrad,
                                          const unique_ptr<transport> &trans,
                                          const double &mut) const;
  double ReynoldsStressDDotVelGrad(const primitive &state,
                                   const tensor<double> &velGrad,
                                   const unique_ptr<transport> &trans,
                                   const double &mut) const;
  double TkeDestruction(const primitive &state, const double &phi = 1.0) const;
  double OmegaDestruction(const primitive &state) const;
  double CrossDiffusion(const primitive &state,
                        const vector3d<double> &kGrad,
                        const vector3d<double> &wGrad) const;

  double CellSpectralRadius(const primitive &state,
                            const unitVec3dMag<double> &fAreaL,
                            const unitVec3dMag<double> &fAreaR,
                            const double &mu,
                            const unique_ptr<transport> &trans,
                            const double &vol, const double &mut,
                            const double &f1, const bool &addSrc) const;
  double FaceSpectralRadius(const primitive &state,
                            const unitVec3dMag<double> &fArea, const double &mu,
                            const unique_ptr<transport> &trans,
                            const double &dist, const double &mut,
                            const double &f1, const bool &positive) const;
  virtual squareMatrix CalcTurbSrc(const primitive &state,
                                   const tensor<double> &velGrad,
                                   const vector3d<double> &kGrad,
                                   const vector3d<double> &wGrad,
                                   const unique_ptr<transport> &trans,
                                   const double &vol,
                                   const double &mut, const double &f1,
                                   const double &f2, const double &width,
                                   vector<double> &turbSrc) const;
  virtual squareMatrix TurbSrcJac(const primitive &state,
                                  const double &beta,
                                  const unique_ptr<transport> &trans,
                                  const double &vol,
                                  const double &phi = 1.0) const;

  // abstract functions (need one for abstract base class)
  virtual void EddyViscAndBlending(const primitive &state,
                                   const tensor<double> &vGrad,
                                   const vector3d<double> &kGrad,
                                   const vector3d<double> &wGrad,
                                   const double &mu, const double &wallDist,
                                   const unique_ptr<transport> &trans,
                                   const double &length, double &mut,
                                   double &f1, double &f2) const = 0;
  virtual void Print() const = 0;

  // destructor
  virtual ~turbModel() noexcept {}
};

class turbNone : public turbModel {
 public:
  // constructor
  turbNone() : turbModel() {}
  explicit turbNone(const string &meth) : turbModel(meth) {}

  // move constructor and assignment operator
  turbNone(turbNone &&model) noexcept : turbModel(std::move(model)) {}
  turbNone& operator=(turbNone&&) = default;

  // copy constructor and assignment operator
  turbNone(const turbNone &model) : turbModel(model) {}
  turbNone& operator=(const turbNone&) = default;

  // member functions
  void EddyViscAndBlending(const primitive &state, const tensor<double> &vGrad,
                           const vector3d<double> &kGrad,
                           const vector3d<double> &wGrad, const double &mu,
                           const double &wallDist,
                           const unique_ptr<transport> &trans,
                           const double &length, double &mut, double &f1,
                           double &f2) const override {}
  double EddyViscNoLim(const primitive &state) const override { return 0.0; }
  double EddyViscNoLim(const primitiveView &state) const override {
    return 0.0;
  }
  double InviscidCellSpecRad(
      const primitive &state, const unitVec3dMag<double> &fAreaL,
      const unitVec3dMag<double> &fAreaR) const override {
    return 0.0;
  }
  double InviscidFaceSpecRad(const primitive &state,
                             const unitVec3dMag<double> &fArea,
                             const bool &postive) const override {
    return 0.0;
  }

  squareMatrix InviscidJacobian(const primitive &state,
                                const unitVec3dMag<double> &fArea,
                                const bool &positive) const override;
  squareMatrix InviscidConvJacobian(
      const primitive &state, const unitVec3dMag<double> &fArea) const override;
  squareMatrix InviscidDissJacobian(
      const primitive &state, const unitVec3dMag<double> &fArea) const override;

  double TkeMin() const override {return 0.0;}
  double OmegaMin() const override {return 0.0;}

  void Print() const override;

  // destructor
  ~turbNone() noexcept {}
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
  double Xw(const primitive &, const tensor<double> &,
            const unique_ptr<transport> &) const;
  double FBeta(const primitive &, const tensor<double> &,
               const unique_ptr<transport> &) const;
  double Beta(const primitive &, const tensor<double> &,
              const unique_ptr<transport> &) const;
  tensor<double> StrainKI(const tensor<double> &) const;
  double OmegaTilda(const primitive&, const tensor<double>&,
                    const unique_ptr<transport>&) const;

 public:
  // constructor
  turbKWWilcox() : turbModel() {}
  explicit turbKWWilcox(const string &meth) : turbModel(meth) {}

  // move constructor and assignment operator
  turbKWWilcox(turbKWWilcox &&model) noexcept : turbModel(std::move(model)) {}
  turbKWWilcox& operator=(turbKWWilcox&&) = default;

  // copy constructor and assignment operator
  turbKWWilcox(const turbKWWilcox &model) : turbModel(model) {}
  turbKWWilcox& operator=(const turbKWWilcox&) = default;

  // member functions
  squareMatrix CalcTurbSrc(const primitive &, const tensor<double> &,
                           const vector3d<double> &, const vector3d<double> &,
                           const unique_ptr<transport> &, const double &,
                           const double &, const double &, const double &,
                           const double &, vector<double> &) const override;
  double EddyVisc(const primitive&, const tensor<double> &,
                  const unique_ptr<transport> &, const double &,
                  const double &) const override;
  void EddyViscAndBlending(const primitive &, const tensor<double> &,
                           const vector3d<double> &, const vector3d<double> &,
                           const double &, const double &,
                           const unique_ptr<transport> &, const double &,
                           double &, double &, double &) const override;
  bool UseUnlimitedEddyVisc() const override { return true; }
  double SrcSpecRad(const primitive &, const unique_ptr<transport> &,
                    const double &, const double & = 1.0) const override;
  squareMatrix ViscousJacobian(const primitive &,
                               const unitVec3dMag<double> &,
                               const double &, const unique_ptr<transport> &,
                               const double &, const double &,
                               const double &) const override;
  double ViscCellSpecRad(const primitive &,
                         const unitVec3dMag<double> &,
                         const unitVec3dMag<double> &,
                         const double &, const unique_ptr<transport> &,
                         const double &, const double &,
                         const double &) const override;
  double ViscFaceSpecRad(const primitive &,
                         const unitVec3dMag<double> &,
                         const double &, const unique_ptr<transport> &,
                         const double &, const double &,
                         const double &) const override;

  squareMatrix TurbSrcJac(const primitive &, const double &,
                          const unique_ptr<transport> &, const double &,
                          const double & = 1.0) const override;

  double TurbPrandtlNumber() const override {return prt_;}
  double WallBeta() const override {return beta0_;}

  double TurbLengthScale(const primitive &state,
                         const unique_ptr<transport> &) const;

  double Gamma() const {return gamma_;}
  double BetaStar() const override {return betaStar_;}
  double Sigma() const {return sigma_;}
  double SigmaStar() const {return sigmaStar_;}
  double SigmaD0() const {return sigmaD0_;}
  double Beta0() const {return beta0_;}
  double CLim() const {return clim_;}

  double SigmaK(const double &f1) const override {return this->SigmaStar();}
  double SigmaW(const double &f1) const override {return this->Sigma();}
  double WallSigmaK() const override {return this->SigmaStar();}
  double WallSigmaW() const override {return this->Sigma();}

  void Print() const override;

  // destructor
  ~turbKWWilcox() noexcept {}
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
  double F1(const double &, const double &, const double &) const;
  double F2(const double &, const double &) const;
  double Alpha1(const primitive &, const unique_ptr<transport> &,
                const double &) const;
  double Alpha2(const primitive &, const unique_ptr<transport> &, const double &,
                const double &) const;
  double Alpha3(const primitive &, const double &, const double &) const;

 public:
  // constructor
  turbKWSst() : turbModel() {}
  explicit turbKWSst(const string &meth) : turbModel(meth) {}

  // move constructor and assignment operator
  turbKWSst(turbKWSst &&model) noexcept : turbModel(std::move(model)) {}
  turbKWSst& operator=(turbKWSst&&) = default;

  // copy constructor and assignment operator
  turbKWSst(const turbKWSst &model) : turbModel(model) {}
  turbKWSst& operator=(const turbKWSst&) = default;

  // member functions
  double BlendedCoeff(const double &, const double &, const double &) const;
  double CDkw(const primitive &, const vector3d<double> &,
              const vector3d<double> &) const;
  virtual squareMatrix CalcTurbSrc(
      const primitive &, const tensor<double> &, const vector3d<double> &,
      const vector3d<double> &, const unique_ptr<transport> &, const double &,
      const double &, const double &, const double &, const double &, 
      vector<double> &) const override;
  double EddyVisc(const primitive &, const tensor<double> &,
                  const unique_ptr<transport> &, const double &,
                  const double &) const override;
  void EddyViscAndBlending(const primitive &, const tensor<double> &,
                           const vector3d<double> &, const vector3d<double> &,
                           const double &, const double &,
                           const unique_ptr<transport> &, const double &,
                           double &, double &, double &) const override;

  virtual double SrcSpecRad(const primitive &, const unique_ptr<transport> &,
                            const double &,
                            const double & = 1.0) const override;
  squareMatrix ViscousJacobian(const primitive &, const unitVec3dMag<double> &,
                               const double &, const unique_ptr<transport> &,
                               const double &, const double &,
                               const double &) const override;
  double ViscCellSpecRad(const primitive &,
                         const unitVec3dMag<double> &,
                         const unitVec3dMag<double> &,
                         const double &, const unique_ptr<transport> &,
                         const double &, const double &,
                         const double &) const override;
  double ViscFaceSpecRad(const primitive &,
                         const unitVec3dMag<double> &,
                         const double &, const unique_ptr<transport> &,
                         const double &, const double &,
                         const double &) const override;

  virtual squareMatrix TurbSrcJac(const primitive &, const double &,
                                  const unique_ptr<transport> &,
                                  const double &,
                                  const double & = 1.0) const override;

  double WallBeta() const override {return beta1_;}
  double TurbPrandtlNumber() const override {return prt_;}

  double Gamma1() const {return gamma1_;}
  double Beta1() const {return beta1_;}
  double SigmaK1() const {return sigmaK1_;}
  double SigmaW1() const {return sigmaW1_;}

  double Gamma2() const {return gamma2_;}
  double Beta2() const {return beta2_;}
  double SigmaK2() const {return sigmaK2_;}
  double SigmaW2() const {return sigmaW2_;}

  double A1() const {return a1_;}
  double BetaStar() const override {return betaStar_;}
  double TkeProd2DestRatio() const {return kProd2Dest_;}
  double TurbLengthScale(const primitive &state,
                         const unique_ptr<transport> &) const;

  // use coefficients from 1 because they are smaller
  // this is used for TSL flux jacobian, so smaller will help increase
  // diagonal dominance of implicit operator
  double SigmaK(const double &f1) const override {
    return this->BlendedCoeff(sigmaK1_, sigmaK2_, f1);
  }
  double SigmaW(const double &f1) const override {
    return this->BlendedCoeff(sigmaW1_, sigmaW2_, f1);
  }
  double WallSigmaK() const override {return sigmaK1_;}
  double WallSigmaW() const override {return sigmaW1_;}

  void Print() const override;

  // destructor
  virtual ~turbKWSst() noexcept {}
};


class turbSstDes : public turbKWSst {
  const double cdes1_ = 0.78;
  const double cdes2_ = 0.61;

  // private member functions
  double Phi(const primitive &state, const double &cdes, const double &width,
             const double &f2, const unique_ptr<transport> &trans) const;

 public:
  // constructor
  turbSstDes() : turbKWSst() {}
  explicit turbSstDes(const string &meth) : turbKWSst(meth) {}

  // move constructor and assignment operator
  turbSstDes(turbSstDes &&model) noexcept : turbKWSst(std::move(model)) {}
  turbSstDes& operator=(turbSstDes&&) = default;

  // copy constructor and assignment operator
  turbSstDes(const turbSstDes &model) : turbKWSst(model) {}
  turbSstDes& operator=(const turbSstDes&) = default;


  squareMatrix CalcTurbSrc(const primitive &, const tensor<double> &,
                           const vector3d<double> &, const vector3d<double> &,
                           const unique_ptr<transport> &, const double &,
                           const double &, const double &, const double &,
                           const double &, vector<double> &) const override;

  double SrcSpecRad(const primitive &, const unique_ptr<transport> &,
                    const double &, const double &) const override;

  double CDes1() const {return cdes1_;}
  double CDes2() const {return cdes2_;}
  double CDes(const double &f1) const {
    return this->BlendedCoeff(cdes1_, cdes2_, f1);
  }
  bool UsePhi() const override {return true;}

  void Print() const override;

  // destructor
  ~turbSstDes() noexcept {}
};


class turbWale : public turbModel {
  const double cw_ = 0.544;

  // private member functions

 public:
  // constructor
  turbWale() : turbModel() {}
  explicit turbWale(const string &meth) : turbModel(meth) {}

  // move constructor and assignment operator
  turbWale(turbWale &&model) noexcept : turbModel(std::move(model)) {}
  turbWale& operator=(turbWale&&) = default;

  // copy constructor and assignment operator
  turbWale(const turbWale &model) : turbModel(model) {}
  turbWale& operator=(const turbWale&) = default;

  double EddyVisc(const primitive &state, const tensor<double> &vGrad,
                  const unique_ptr<transport> &trans, const double &f2,
                  const double &length) const override;

  void EddyViscAndBlending(const primitive &state, const tensor<double> &vGrad,
                           const vector3d<double> &kGrad,
                           const vector3d<double> &wGrad, const double &mu,
                           const double &wallDist,
                           const unique_ptr<transport> &trans,
                           const double &length, double &mut, double &f1,
                           double &f2) const override {
    f1 = 1.0;
    f2 = 0.0;
    mut = this->EddyVisc(state, vGrad, trans, f2, length);
  }

  double Cw() const {return cw_;}
  tensor<double> SigmaD(const tensor<double> &vGrad) const;
  void Print() const override;

  // destructor
  ~turbWale() noexcept {}
};



// function declarations
template <typename T>
double EddyViscUnlimited(const T &state);

#endif

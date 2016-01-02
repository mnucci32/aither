/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

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

#include <cstdlib>    // exit()
#include <iostream>   // cout
#include <algorithm>  // max
#include "turbulence.hpp"
#include "primVars.hpp"  // primVars
#include "eos.hpp"       // sutherland

using std::cout;
using std::endl;
using std::cerr;
using std::max;
using std::min;

// -------------------------------------------------------------------------
// member functions for turbModel class, all turbulence models inherit from
// this class

// member function to return the eddy viscosity calculated without the stress
// limiter
double turbModel::EddyViscNoLim(const primVars &state) const {
  return state.Rho() * state.Tke() / state.Omega();
}

tensor<double> turbModel::BoussinesqReynoldsStress(
    const primVars &state, const tensor<double> &velGrad,
    const sutherland &suth, const double &mut) const {
  const auto lambda = suth.Lambda(mut);

  tensor<double> I;
  I.Identity();
  const auto tau =
      lambda * velGrad.Trace() * I + mut * (velGrad + velGrad.Transpose())
      - 2.0 / 3.0 * state.Rho() * state.Tke() * I;
  return tau;
}

// Member function for calculation of tke production
double turbModel::ReynoldsStressDDotVelGrad(const primVars &state,
                                            const tensor<double> &velGrad,
                                            const sutherland &suth,
                                            const double &mut) const {
  const auto tau = this->BoussinesqReynoldsStress(state, velGrad, suth, mut);
  return tau.DoubleDotTrans(velGrad);
}

// member function for destruction of tke
double turbModel::TkeDestruction(const primVars &state) const {
  return state.Rho() * state.Tke() * state.Omega();
}

// member function for destruction of omega
double turbModel::OmegaDestruction(const primVars &state) const {
  return state.Rho() * state.Omega() * state.Omega();
}

// member function for cross diffusion term without coefficient
double turbModel::CrossDiffusion(const primVars &state,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad) const {
  return state.Rho() / state.Omega() * kGrad.DotProd(wGrad);
}

// -------------------------------------------------------------------------
// member functions for the turbNone class

// member function to print out turbulence variables
void turbNone::Print() const {
  cout << "No Turbulence Model" << endl;
}

// member function to calculate turbulence source terms
double turbNone::CalcTurbSrc(const primVars &state,
                             const tensor<double> &velGrad,
                             const vector3d<double> &kGrad,
                             const vector3d<double> &wGrad,
                             const sutherland &suth, const idealGas &eos,
                             const double &wallDist, double &ksrc,
                             double &wsrc) const {
  // set k and omega source terms to zero
  ksrc = 0.0;
  wsrc = 0.0;

  // return source jacobian spectral radius
  return this->SpecRad(state, suth);
}

// ---------------------------------------------------------------------
// K-Omega Wilcox member functions

// member function to return the eddy viscosity calculated with the stress
// limiter
double turbKWWilcox::EddyVisc(const primVars &state,
                              const tensor<double> &vGrad,
                              const sutherland &suth, const double &f2) const {
  return state.Rho() * state.Tke() / this->OmegaTilda(state, vGrad, suth);
}

double turbKWWilcox::SigmaD(const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const {
  // kGrad -- tke gradient
  // wGrad -- omega gradient

  if (kGrad.DotProd(wGrad) <= 0.0) {
    return 0.0;
  } else {
    return sigmaD0_;
  }
}

double turbKWWilcox::Beta(const primVars &state,
                          const tensor<double> &velGrad,
                          const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  return beta0_ * this->FBeta(state, velGrad, suth);
}

double turbKWWilcox::FBeta(const primVars &state,
                           const tensor<double> &velGrad,
                           const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  auto xw = this->Xw(state, velGrad, suth);
  return (1.0 + 85.0 * xw) / (1.0 + 100.0 * xw);
}

double turbKWWilcox::Xw(const primVars &state,
                        const tensor<double> &velGrad,
                        const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity

  const auto vorticity = 0.5 * (velGrad - velGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  // vorticity is asymmetric but vorticity * vorticity is symmetric
  return fabs( (vorticity * vorticity).DoubleDotTrans(this->StrainKI(velGrad))
               / pow(betaStar_ * state.Omega(), 3.0) )
      * pow(suth.NondimScaling(), 3.0);
}

tensor<double> turbKWWilcox::StrainKI(const tensor<double> &velGrad) const {
  // velGrad -- velocity gradient
  tensor<double> I;
  I.Identity();
  return 0.5 * (velGrad + velGrad.Transpose() - velGrad.Trace() * I);
}

double turbKWWilcox::OmegaTilda(const primVars &state,
                                const tensor<double> &velGrad,
                                const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity

  tensor<double> I;
  I.Identity();
  const auto sHat = 0.5 * (velGrad + velGrad.Transpose()) - 1.0 / 3.0 *
      velGrad.Trace() * I;

  // using DoubleDotTrans instead of DoubleDot for speed
  // since tensors are symmetric, result is the same
  return std::max(state.Omega(), suth.NondimScaling() * clim_ *
                  sqrt(2.0 * sHat.DoubleDotTrans(sHat) / betaStar_));
}

// member function to calculate turbulence source terms and return source
// spectral radius
double turbKWWilcox::CalcTurbSrc(const primVars &state,
                                 const tensor<double> &velGrad,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad,
                                 const sutherland &suth, const idealGas &eos,
                                 const double &wallDist, double &ksrc,
                                 double &wsrc) const {
  // calculate tke destruction
  const auto tkeDest = suth.InvNondimScaling() * betaStar_ *
      this->TkeDestruction(state);

  // calculate omega destruction
  const auto omgDest = suth.InvNondimScaling() *
      this->Beta(state, velGrad, suth) * this->OmegaDestruction(state);

  // calculate tke production
  const auto mut = this->EddyVisc(state, velGrad, suth, 0.0);
  const auto tkeProd = suth.NondimScaling() *
      this->ReynoldsStressDDotVelGrad(state, velGrad, suth, mut);

  // calculate omega production
  const auto omgProd = gamma_ * state.Omega() / state.Tke() * tkeProd;

  // calculate omega cross diffusion
  const auto omgCd = suth.NondimScaling() * this->SigmaD(kGrad, wGrad) *
      this->CrossDiffusion(state, kGrad, wGrad);

  // assign source term values
  ksrc = tkeProd - tkeDest;
  wsrc = omgProd - omgDest + omgCd;

  // return spectral radius of source jacobian
  return this->SpecRad(state, suth);
}

// member function to calculate the eddy viscosity, and the molecular diffusion
// coefficients. The eddy viscosity is used in the viscous flux calculation
// for the Navier-Stokes equations, and the molecular diffusion coefficients are
// used in the viscous flux calculation for the turbulence equations
double turbKWWilcox::EddyViscAndMolecDiffCoeff(const primVars &state,
                                               const tensor<double> &velGrad,
                                               const vector3d<double> &kGrad,
                                               const vector3d<double> &wGrad,
                                               const sutherland &suth,
                                               const idealGas &eos,
                                               const double &wallDist,
                                               double &sigmaK,
                                               double &sigmaW) const {
  // calculate limited eddy (effective) viscosity
  const auto mut = this->EddyViscNoLim(state) * suth.NondimScaling();

  // calculate blended coefficients
  sigmaK = sigmaStar_ * mut;
  sigmaW = sigma_ * mut;

  return mut;
}

// member function to calculate the spectral radius of the source jacobian
double turbKWWilcox::SpecRad(const primVars &state,
                             const sutherland &suth) const {
  return -2.0 * betaStar_ * state.Omega() * suth.InvNondimScaling();
}

// member function to print out turbulence variables
void turbKWWilcox::Print() const {
  cout << "Eddy Viscosity Method: " << this->EddyViscMethod() << endl;
  cout << "Gamma: " << gamma_ << endl;
  cout << "Beta*: " << betaStar_ << endl;
  cout << "Sigma: " << sigma_ << endl;
  cout << "Sigma*: " << sigmaStar_ << endl;
  cout << "SigmaD0: " << sigmaD0_ << endl;
  cout << "Beta0: " << beta0_ << endl;
  cout << "Clim: " << clim_ << endl;
  cout << "Turbulent Prandtl Number: " << prt_ << endl;
}

// ----------------------------------------------------------------------
// K-Omega SST member functions

// member function for eddy viscosity with limiter
double turbKWSst::EddyVisc(const primVars &state, const tensor<double> &vGrad,
                           const sutherland &suth, const double &f2) const {
  const auto strainRate = 0.5 * (vGrad + vGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  const auto meanStrainRate = sqrt(2.0 * strainRate.DoubleDotTrans(strainRate));
  return state.Rho() * a1_ * state.Tke() /
      max(a1_ * state.Omega(), suth.NondimScaling() * meanStrainRate * f2);
}

double turbKWSst::CDkw(const primVars &state, const vector3d<double> &kGrad,
                       const vector3d<double> &wGrad) const {
  return max(2.0 * state.Rho() * sigmaW2_ / state.Omega() *
             kGrad.DotProd(wGrad), 1.0e-10);
}

double turbKWSst::F1(const double &alpha1, const double &alpha2,
                     const double &alpha3) const {
  const auto arg1 = min(max(alpha1, alpha2), alpha3);
  return tanh(pow(arg1, 4.0));
}

double turbKWSst::F2(const double &alpha1, const double &alpha2) const {
  const auto arg2 = max(2.0 * alpha1, alpha2);
  return tanh(arg2 * arg2);
}

double turbKWSst::BlendedCoeff(const double &coeff1, const double &coeff2,
                               const double &f1) const {
  return f1 * coeff1 + (1.0 - f1) * coeff2;
}

double turbKWSst::Alpha1(const primVars &state, const sutherland &suth,
                         const double &wallDist) const {
  return suth.NondimScaling() * sqrt(state.Tke()) /
      (betaStar_ * state.Omega() * wallDist);
}

double turbKWSst::Alpha2(const primVars &state, const sutherland &suth,
                         const idealGas &eos, const double &wallDist) const {
  return suth.NondimScaling() * suth.NondimScaling() *
      500.0 * suth.Viscosity(state.Temperature(eos)) /
      (wallDist * wallDist * state.Rho() * state.Omega());
}

double turbKWSst::Alpha3(const primVars &state, const double &wallDist,
                         const double &cdkw) const {
  return 4.0 * state.Rho() * sigmaW2_ * state.Tke() /
      (cdkw * wallDist * wallDist);
}

// member function to calculate turbulence source terms and source spectral
// radius
double turbKWSst::CalcTurbSrc(const primVars &state,
                              const tensor<double> &velGrad,
                              const vector3d<double> &kGrad,
                              const vector3d<double> &wGrad,
                              const sutherland &suth, const idealGas &eos,
                              const double &wallDist, double &ksrc,
                              double &wsrc) const {
  // calculate blending functions
  const auto alpha1 = this->Alpha1(state, suth, wallDist);
  const auto alpha2 = this->Alpha2(state, suth, eos, wallDist);
  const auto cdkw = this->CDkw(state, kGrad, wGrad);
  const auto alpha3 = this->Alpha3(state, wallDist, cdkw);
  const auto f1 = this->F1(alpha1, alpha2, alpha3);
  const auto f2 = this->F2(alpha1, alpha2);

  // calculate blended coefficients
  const auto gamma = this->BlendedCoeff(gamma1_, gamma2_, f1);
  const auto beta = this->BlendedCoeff(beta1_, beta2_, f1);

  // calculate tke destruction
  const auto tkeDest = suth.InvNondimScaling() * betaStar_ *
      this->TkeDestruction(state);

  // calculate omega destruction
  const auto omgDest = suth.InvNondimScaling() * beta *
      this->OmegaDestruction(state);

  // calculate tke production
  const auto mut = this->EddyVisc(state, velGrad, suth, f2);
  const auto tkeProd =
      min(suth.NondimScaling() *
          this->ReynoldsStressDDotVelGrad(state, velGrad, suth, mut),
          kProd2Dest_ * tkeDest);

  // calculate omega production
  const auto omgProd = gamma * state.Rho() / mut * tkeProd;

  // calculate omega cross diffusion
  // Using CDkw instead of whole cross diffusion term
  // Both Loci/CHEM and SU2 use this implementation
  const auto omgCd = suth.NondimScaling() * (1.0 - f1) * cdkw;

  // assign source term values
  ksrc = tkeProd - tkeDest;
  wsrc = omgProd - omgDest + omgCd;

  // return spectral radius of source jacobian
  return this->SpecRad(state, suth);
}

// member function to calculate the eddy viscosity, and the molecular diffusion
// coefficients. The eddy viscosity is used in the viscous flux calculation
// for the Navier-Stokes equations, and the molecular diffusion coefficients are
// used in the viscous flux calculation for the turbulence equations
double turbKWSst::EddyViscAndMolecDiffCoeff(const primVars &state,
                                            const tensor<double> &velGrad,
                                            const vector3d<double> &kGrad,
                                            const vector3d<double> &wGrad,
                                            const sutherland &suth,
                                            const idealGas &eos,
                                            const double &wallDist,
                                            double &sigmaK,
                                            double &sigmaW) const {
  // calculate blending functions
  const auto alpha1 = this->Alpha1(state, suth, wallDist);
  const auto alpha2 = this->Alpha2(state, suth, eos, wallDist);
  const auto cdkw = this->CDkw(state, kGrad, wGrad);
  const auto alpha3 = this->Alpha3(state, wallDist, cdkw);
  const auto f1 = this->F1(alpha1, alpha2, alpha3);
  const auto f2 = this->F2(alpha1, alpha2);

  // calculate limited eddy (effective) viscosity
  const auto mut = this->EddyVisc(state, velGrad, suth, f2) *
      suth.NondimScaling();

  // calculate blended coefficients
  sigmaK = this->BlendedCoeff(sigmaK1_, sigmaK2_, f1) * mut;
  sigmaW = this->BlendedCoeff(sigmaW1_, sigmaW2_, f1) * mut;

  return mut;
}

// member function to calculate the spectral radius of the source jacobian
double turbKWSst::SpecRad(const primVars &state, const sutherland &suth) const {
  // DEBUG
  return -2.0 * betaStar_ * state.Omega() * suth.InvNondimScaling() * 1.0e-6;
}

// member function to print out turbulence variables
void turbKWSst::Print() const {
  cout << "Eddy Viscosity Method: " << this->EddyViscMethod() << endl;
  cout << "A1: " << a1_ << endl;
  cout << "Beta*: " << betaStar_ << endl;
  cout << "Turbulent Prandtl Number: " << prt_ << endl;
  cout << "Production to Destruction Ratio: " << kProd2Dest_ << endl;

  cout << "Sigma K1: " << sigmaK1_ << endl;
  cout << "Sigma W1: " << sigmaW1_ << endl;
  cout << "Beta1: " << beta1_ << endl;
  cout << "Gamma1: " << gamma1_ << endl;

  cout << "Sigma K2: " << sigmaK2_ << endl;
  cout << "Sigma W2: " << sigmaW2_ << endl;
  cout << "Beta2: " << beta2_ << endl;
  cout << "Gamma2: " << gamma2_ << endl;
}

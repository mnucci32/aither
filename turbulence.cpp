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

// member funtion to print out turbulence variables
void turbNone::Print() const {
  cout << "No Turbulence Model" << endl;
}

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
    const sutherland &suth, const idealGas &eos, const double &wallDist) const {
  double mut = this->EddyVisc(state, velGrad, suth, eos, wallDist);
  double lambda = suth.Lambda(mut);

  tensor<double> I;
  I.Identity();
  tensor<double> tau =
      lambda * velGrad.Trace() * I + mut * (velGrad + velGrad.Transpose())
      - 2.0 / 3.0 * state.Rho() * state.Tke() * I;
  return tau;
}

// Member function for calculation of tke production
double turbModel::ReynoldsStressDDotVelGrad(const primVars &state,
                                            const tensor<double> &velGrad,
                                            const sutherland &suth,
                                            const idealGas &eos,
                                            const double &wallDist) const {
  tensor<double> tau = this->BoussinesqReynoldsStress(state, velGrad, suth, eos,
                                                      wallDist);
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


// ---------------------------------------------------------------------
// K-Omega Wilcox member functions

// member function to return the eddy viscosity calculated with the stress
// limiter
double turbKWWilcox::EddyVisc(const primVars &state,
                              const tensor<double> &vGrad,
                              const sutherland &suth, const idealGas &eos,
                              const double &wallDist) const {
  return state.Rho() * state.Tke() / this->OmegaTilda(state, vGrad, suth);
}

// member function for the production of tke
double turbKWWilcox::Production1(const primVars &state,
                                 const tensor<double> &velGrad,
                                 const sutherland &suth, const idealGas &eos,
                                 const double &wallDist) const {
  return suth.NondimScaling() *
      this->ReynoldsStressDDotVelGrad(state, velGrad, suth, eos, wallDist);
}

// member function for the destruction of tke
double turbKWWilcox::Destruction1(const primVars &state,
                                  const sutherland &suth) const {
  return suth.InvNondimScaling() * betaStar_ * this->TkeDestruction(state);
}

// member function for production of omega
double turbKWWilcox::Production2(const primVars &state,
                                 const tensor<double> &velGrad,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad,
                                 const sutherland &suth, const idealGas &eos,
                                 const double &wallDist) const {
  return gamma_ * state.Omega() / state.Tke() *
      this->Production1(state, velGrad, suth, eos, wallDist);
}

// member function for destruction of omega
double turbKWWilcox::Destruction2(const primVars &state,
                                  const tensor<double> &velGrad,
                                  const vector3d<double> &kGrad,
                                  const vector3d<double> &wGrad,
                                  const sutherland &suth, const idealGas &eos,
                                  const double &wallDist) const {
  return suth.InvNondimScaling() * this->Beta(state, velGrad, suth) *
      this->OmegaDestruction(state);
}

// member funtion for cross diffusion term in omega equation
double turbKWWilcox::CrossDiff2(const primVars &state,
                                const vector3d<double> &kGrad,
                                const vector3d<double> &wGrad,
                                const sutherland &suth, const idealGas &eos,
                                const double &wallDist) const {
  return suth.NondimScaling() * this->SigmaD(kGrad, wGrad) *
      this->CrossDiffusion(state, kGrad, wGrad);
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
  double xw = this->Xw(state, velGrad, suth);
  return (1.0 + 85.0 * xw) / (1.0 + 100.0 * xw);
}

double turbKWWilcox::Xw(const primVars &state,
                        const tensor<double> &velGrad,
                        const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity

  tensor<double> vorticity = 0.5 * (velGrad - velGrad.Transpose());

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
  tensor<double> sHat = 0.5 * (velGrad + velGrad.Transpose()) - 1.0 / 3.0 *
      velGrad.Trace() * I;

  // using DoubleDotTrans instead of DoubleDot for speed
  // since tensors are symmetric, result is the same
  return std::max(state.Omega(), suth.NondimScaling() * clim_ *
                  sqrt(2.0 * sHat.DoubleDotTrans(sHat) / betaStar_));
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
                           const sutherland &suth, const idealGas &eos,
                           const double &wallDist) const {
  tensor<double> strainRate = 0.5 * (vGrad + vGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  double meanStrainRate = sqrt(2.0 * strainRate.DoubleDotTrans(strainRate));
  return state.Rho() * a1_ * state.Tke() /
      max(a1_ * state.Omega(), suth.NondimScaling() * meanStrainRate *
          this->F2(state, suth, eos, wallDist));
}

// member function for production of tke
double turbKWSst::Production1(const primVars &state,
                              const tensor<double> &velGrad,
                              const sutherland &suth, const idealGas &eos,
                              const double &wallDist) const {
  return min(suth.NondimScaling() *
             this->ReynoldsStressDDotVelGrad(state, velGrad, suth, eos,
                                             wallDist),
             kProd2Dest_ * this->Destruction1(state, suth));
}

// member function for the destruction of tke
double turbKWSst::Destruction1(const primVars &state,
                               const sutherland &suth) const {
  return suth.InvNondimScaling() * betaStar_ * this->TkeDestruction(state);
}

// member function for production of omega
double turbKWSst::Production2(const primVars &state,
                              const tensor<double> &velGrad,
                              const vector3d<double> &kGrad,
                              const vector3d<double> &wGrad,
                              const sutherland &suth, const idealGas &eos,
                              const double &wallDist) const {
  return this->BlendedCoeff(gamma1_, gamma2_, state, kGrad, wGrad, suth, eos,
                            wallDist) * state.Omega() / state.Tke() *
      this->Production1(state, velGrad, suth, eos, wallDist);
}

// member function for destruction of omega
double turbKWSst::Destruction2(const primVars &state,
                               const tensor<double> &velGrad,
                               const vector3d<double> &kGrad,
                               const vector3d<double> &wGrad,
                               const sutherland &suth, const idealGas &eos,
                               const double &wallDist) const {
  return suth.InvNondimScaling() *
      this->BlendedCoeff(beta1_, beta2_, state, kGrad, wGrad, suth, eos,
                            wallDist) * this->OmegaDestruction(state);
}

// member funtion for cross diffusion term in omega equation
double turbKWSst::CrossDiff2(const primVars &state,
                             const vector3d<double> &kGrad,
                             const vector3d<double> &wGrad,
                             const sutherland &suth, const idealGas &eos,
                             const double &wallDist) const {
  // DEBUG - using CDkw instead of CrossDiffusion (keeps value positive)
  return suth.NondimScaling() *
      (1.0 - this->F1(state, kGrad, wGrad, suth, eos, wallDist)) *
      this->CDkw(state, kGrad, wGrad);

  // return suth.NondimScaling() * 2.0 *
  //     (1.0 - this->F1(state, kGrad, wGrad, suth, eos, wallDist)) *
  //     sigmaW2_ * this->CrossDiffusion(state, kGrad, wGrad);
}

double turbKWSst::CDkw(const primVars &state, const vector3d<double> &kGrad,
                       const vector3d<double> &wGrad) const {
  return max(2.0 * state.Rho() * sigmaW2_ / state.Omega() *
             kGrad.DotProd(wGrad), 10e-10);
}

double turbKWSst::F1(const primVars &state, const vector3d<double> &kGrad,
                     const vector3d<double> &wGrad, const sutherland &suth,
                     const idealGas &eos, const double &wallDist) const {
  double arg1 = min(max(this->Alpha1(state, suth, wallDist),
                        this->Alpha2(state, suth, eos, wallDist)),
                    this->Alpha3(state, kGrad, wGrad, wallDist));
  return tanh(pow(arg1, 4.0));
}

double turbKWSst::F2(const primVars &state, const sutherland &suth,
                     const idealGas &eos, const double &wallDist) const {
  double arg2 = max(2.0 * this->Alpha1(state, suth, wallDist),
                    this->Alpha2(state, suth, eos, wallDist));
  return tanh(arg2 * arg2);
}

double turbKWSst::BlendedCoeff(const double &coeff1, const double &coeff2,
                               const primVars &state,
                               const vector3d<double> &kGrad,
                               const vector3d<double> &wGrad,
                               const sutherland &suth, const idealGas &eos,
                               const double &wallDist) const {
  double f1 = this->F1(state, kGrad, wGrad, suth, eos, wallDist);
  return f1 * coeff1 + (1.0 - f1) * coeff2;
}

double turbKWSst::MolecDiff1Coeff(const primVars &state,
                                  const vector3d<double> &kGrad,
                                  const vector3d<double> &wGrad,
                                  const sutherland &suth, const idealGas &eos,
                                  const double &wallDist) const {
  return this->BlendedCoeff(sigmaK1_, sigmaK2_, state, kGrad, wGrad, suth, eos,
                            wallDist);
}

double turbKWSst::MolecDiff2Coeff(const primVars &state,
                                  const vector3d<double> &kGrad,
                                  const vector3d<double> &wGrad,
                                  const sutherland &suth, const idealGas &eos,
                                  const double &wallDist) const {
  return this->BlendedCoeff(sigmaW1_, sigmaW2_, state, kGrad, wGrad, suth, eos,
                            wallDist);
}

double turbKWSst::Alpha1(const primVars &state, const sutherland &suth,
                         const double &wallDist) const {
  // DEBUG - set to 0 if k is negative (ghost cells)
  return (state.Tke() <= 0.0) ? 0.0 :
      suth.NondimScaling() * sqrt(state.Tke()) /
      (betaStar_ * state.Omega() * wallDist);
}

double turbKWSst::Alpha2(const primVars &state, const sutherland &suth,
                         const idealGas &eos, const double &wallDist) const {
  return suth.NondimScaling() * suth.NondimScaling() *
      500.0 * suth.Viscosity(state.Temperature(eos)) /
      (wallDist * wallDist * state.Rho() * state.Omega());
}

double turbKWSst::Alpha3(const primVars &state, const vector3d<double> &kGrad,
                         const vector3d<double> &wGrad,
                         const double &wallDist) const {
  return 4.0 * state.Rho() * sigmaW2_ * state.Tke() /
      (this->CDkw(state, kGrad, wGrad) * wallDist * wallDist);
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

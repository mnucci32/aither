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
  // state -- primative variables
  return state.Rho() * state.Tke() / state.Omega();
}

/* member function to calculate the Reynolds stress tensor using the Boussinesq
   approximation.

   tau = lambda * trace(velGrad) * I + mut * (velGrad + transpose(velGrad))
         - 2/3 * rho * k * I

   In the above equation tau is the reynolds stress tensor, mut is the turbulent
   eddy viscosity, velGrad is the velocity gradient, and trace and transpose are
   operators on a tensor that take the trace and transpose respectively.
*/
tensor<double> turbModel::BoussinesqReynoldsStress(
    const primVars &state, const tensor<double> &velGrad,
    const sutherland &suth, const double &mut) const {
  // state -- primative variables
  // velGrad -- velocity gradient tensor
  // suth -- sutherland's law for viscosity
  // mut -- turbulent eddy viscosity

  const auto lambda = suth.Lambda(mut);  // 2nd eddy viscosity

  tensor<double> I;
  I.Identity();
  return lambda * velGrad.Trace() * I + mut * (velGrad + velGrad.Transpose())
      - 2.0 / 3.0 * state.Rho() * state.Tke() * I;
}

// Member function for calculation of tke production
// P = tau_ij * velGrad_ij
double turbModel::ReynoldsStressDDotVelGrad(const primVars &state,
                                            const tensor<double> &velGrad,
                                            const sutherland &suth,
                                            const double &mut) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for eddy viscosity
  // mut -- turbulent eddy viscosity

  const auto tau = this->BoussinesqReynoldsStress(state, velGrad, suth, mut);
  return tau.DoubleDotTrans(velGrad);
}

// member function for destruction of tke
// Dk = rho * k * w
double turbModel::TkeDestruction(const primVars &state) const {
  // state -- primative variables
  return state.Rho() * state.Tke() * state.Omega();
}

// member function for destruction of omega
// Dw = rho * w * w
double turbModel::OmegaDestruction(const primVars &state) const {
  // state -- primative variables
  return state.Rho() * state.Omega() * state.Omega();
}

// member function for cross diffusion term without coefficient
// CD = rho / w * kGrad (dot) wGrad
// In the above equation kGrad and wGrad are the tke and omega gradients.
double turbModel::CrossDiffusion(const primVars &state,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad) const {
  // state -- primative variables
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  return state.Rho() / state.Omega() * kGrad.DotProd(wGrad);
}

// member function to calculate the spectral radius of the turbulence equations
double turbModel::SpectralRadius(const primVars &state,
                                 const unitVec3dMag<double> &fAreaL,
                                 const unitVec3dMag<double> &fAreaR,
                                 const double &mu, const sutherland &suth,
                                 const double &vol, const double &mut,
                                 const double &f1, const bool &addSrc) const {
  // state -- primative variables
  // fAreaL -- area at left face
  // fAreaR -- area at right face
  // mu -- laminar viscosity
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // addSrc -- flag to determine if source jacobian spectral radius should be
  //           included

  auto specRad = this->InviscidSpecRad(state, fAreaL, fAreaR);
  // factor of 2 because viscous spectral radius is not halved (Blazek 6.53)
  specRad += 2.0 * this->ViscSpecRad(state, fAreaL, fAreaR, mu, suth, vol,
                                     mut, f1);
  if (addSrc) {
    // minus sign because source terms are on RHS
    specRad -= this->SrcSpecRad(state, suth, vol);
  }

  return specRad;
}

// member function to calculate inviscid spectral radius
// df_dq = [vel (dot) area   0
//                0          vel (dot) area]
double turbModel::InviscidSpecRad(const primVars &state,
                                  const unitVec3dMag<double> &fAreaL,
                                  const unitVec3dMag<double> &fAreaR) const {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face

  auto normAvg = (0.5 * (fAreaL.UnitVector() +
                         fAreaR.UnitVector())).Normalize();
  auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  return state.Velocity().DotProd(normAvg) * fMag;
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
                             const sutherland &suth, const double &vol,
                             const double &turbVisc, const double &f1,
                             double &ksrc, double &wsrc) const {
  // set k and omega source terms to zero
  ksrc = 0.0;
  wsrc = 0.0;

  // return source jacobian spectral radius
  return this->SrcSpecRad(state, suth, vol);
}

// ---------------------------------------------------------------------
// K-Omega Wilcox member functions

// member function to return the eddy viscosity calculated with the stress
// limiter
double turbKWWilcox::EddyVisc(const primVars &state,
                              const tensor<double> &vGrad,
                              const sutherland &suth, const double &f2) const {
  // state -- primative variables
  // vGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  // f2 -- SST blending coefficient (not used in Wilcox K-W)
  return state.Rho() * state.Tke() / this->OmegaTilda(state, vGrad, suth);
}

// member function to calculate the cross diffusion coefficient term
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

// member function to calculate coefficient for omega destruction
double turbKWWilcox::Beta(const primVars &state,
                          const tensor<double> &velGrad,
                          const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  return beta0_ * this->FBeta(state, velGrad, suth);
}

// member function to calculate coefficient used in beta calculation
double turbKWWilcox::FBeta(const primVars &state,
                           const tensor<double> &velGrad,
                           const sutherland &suth) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  auto xw = this->Xw(state, velGrad, suth);
  return (1.0 + 85.0 * xw) / (1.0 + 100.0 * xw);
}

// member function to calculate vortex stretching coefficient
// used in fbeta calculation
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
  return fabs(
      (vorticity.MatMult(vorticity)).DoubleDotTrans(this->StrainKI(velGrad))
      / pow(betaStar_ * state.Omega(), 3.0) )
      * pow(suth.NondimScaling(), 3.0);
}

// member function to calculate term used in calculation of Xw
tensor<double> turbKWWilcox::StrainKI(const tensor<double> &velGrad) const {
  // velGrad -- velocity gradient
  tensor<double> I;
  I.Identity();
  return 0.5 * (velGrad + velGrad.Transpose() - velGrad.Trace() * I);
}

// member function to calculate adjusted omega used in eddy viscosity limiter
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
                                 const sutherland &suth, const double &vol,
                                 const double &mut, const double &f1,
                                 double &ksrc, double &wsrc) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // suth -- sutherland's law
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // ksrc -- source term for tke equation
  // wsrc -- source term for omega equation

  // calculate tke destruction
  const auto tkeDest = suth.InvNondimScaling() * betaStar_ *
      this->TkeDestruction(state);

  // calculate omega destruction
  const auto omgDest = suth.InvNondimScaling() *
      this->Beta(state, velGrad, suth) * this->OmegaDestruction(state);

  // calculate tke production
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
  return this->SrcSpecRad(state, suth, vol);
}

// member function to calculate the eddy viscosity, and the blending
// coefficients.
void turbKWWilcox::EddyViscAndBlending(const primVars &state,
                                       const tensor<double> &velGrad,
                                       const vector3d<double> &kGrad,
                                       const vector3d<double> &wGrad,
                                       const double &mu,
                                       const double &wallDist,
                                       const sutherland &suth,
                                       double &mut, double &f1,
                                       double &f2) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // mu -- laminar viscosity
  // walldist -- distance to nearest viscous wall
  // suth -- sutherland's law for viscosity
  // mut -- turbulent eddy viscosity
  // f1 -- first blending coefficient
  // f2 -- second blending coefficient

  f1 = 1.0;
  f2 = 0.0;

  // return eddy viscosity, scaled for nondimensional equations
  mut = this->EddyVisc(state, velGrad, suth, f2);
}


/* member function to calculate the spectral radius of the source jacobian
   Source jacobian from Wilcox used.
   
   ds_dq = [ -2 * betaStar * w   0
                   0            -2 * beta * w]
   
   This is a diagonal matrix so eigenvalues are trivial. Since betaStar is
   always larger than beta, this eigenvalue is used
*/
double turbKWWilcox::SrcSpecRad(const primVars &state,
                                const sutherland &suth,
                                const double &vol) const {
  // state -- primative variables
  // suth -- sutherland's law for viscosity
  // vol -- cell volume

  // return spectral radius scaled for nondimensional equations
  return -2.0 * betaStar_ * state.Omega() * vol * suth.InvNondimScaling();
}

// member function to calculate viscous spectral radius
// dfv_dq = [ (area / vol) * (nu + sigmaStar * nut)    0
//                           0                (area / vol) * (nu + sigma * nut)]
double turbKWWilcox::ViscSpecRad(const primVars &state,
                                 const unitVec3dMag<double> &fAreaL,
                                 const unitVec3dMag<double> &fAreaR,
                                 const double &mu, const sutherland &suth,
                                 const double &vol, const double &mut,
                                 const double &f1) const {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // mu -- laminar viscosity
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());

  // Wilcox method uses unlimited eddy viscosity
  return suth.NondimScaling() * fMag * fMag / (vol * state.Rho()) *
      (mu + this->SigmaK(f1) * this->EddyViscNoLim(state));
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
  // state -- primative variables
  // vGrad -- velocity gradient
  // suth -- sutherland's law for viscosity
  // f2 -- SST blending coefficient

  const auto strainRate = 0.5 * (vGrad + vGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  const auto meanStrainRate = sqrt(2.0 * strainRate.DoubleDotTrans(strainRate));
  return state.Rho() * a1_ * state.Tke() /
      max(a1_ * state.Omega(), suth.NondimScaling() * meanStrainRate * f2);
}

// member function to calculate cross diffusion term
double turbKWSst::CDkw(const primVars &state, const vector3d<double> &kGrad,
                       const vector3d<double> &wGrad) const {
  return max(2.0 * state.Rho() * sigmaW2_ / state.Omega() *
             kGrad.DotProd(wGrad), 1.0e-10);
}

// member function to calculate blending function
double turbKWSst::F1(const double &alpha1, const double &alpha2,
                     const double &alpha3) const {
  const auto arg1 = min(max(alpha1, alpha2), alpha3);
  return tanh(pow(arg1, 4.0));
}

// member function to calculate blending function
double turbKWSst::F2(const double &alpha1, const double &alpha2) const {
  const auto arg2 = max(2.0 * alpha1, alpha2);
  return tanh(arg2 * arg2);
}

// member function to calculate a blended coefficient
double turbKWSst::BlendedCoeff(const double &coeff1, const double &coeff2,
                               const double &f1) const {
  // coeff1 -- coefficient from set 1
  // coeff2 -- coefficient from set 2
  // f1 -- blending term

  return f1 * coeff1 + (1.0 - f1) * coeff2;
}

// member function to calculate blending term
double turbKWSst::Alpha1(const primVars &state, const sutherland &suth,
                         const double &wallDist) const {
  return suth.NondimScaling() * sqrt(state.Tke()) /
      (betaStar_ * state.Omega() * wallDist);
}

// member function to calculate blending term
double turbKWSst::Alpha2(const primVars &state, const sutherland &suth,
                         const double &mu, const double &wallDist) const {
  return suth.NondimScaling() * suth.NondimScaling() *
      500.0 * mu / (wallDist * wallDist * state.Rho() * state.Omega());
}

// member function to calculate blending term
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
                              const sutherland &suth, const double &vol,
                              const double &mut, const double &f1,
                              double &ksrc, double &wsrc) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // suth -- sutherland's law
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // ksrc -- source term for tke equation
  // wsrc -- source term for omega equation

  // calculate cross diffusion coefficient
  const auto cdkw = this->CDkw(state, kGrad, wGrad);

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
  return this->SrcSpecRad(state, suth, vol);
}

// member function to calculate the eddy viscosity, and the blending
// coefficients.
void turbKWSst::EddyViscAndBlending(const primVars &state,
                                    const tensor<double> &velGrad,
                                    const vector3d<double> &kGrad,
                                    const vector3d<double> &wGrad,
                                    const double &mu,
                                    const double &wallDist,
                                    const sutherland &suth,
                                    double &mut, double &f1,
                                    double &f2) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // mu -- laminar viscosity
  // walldist -- distance to nearest viscous wall
  // suth -- sutherland's law for viscosity
  // mut -- turbulent eddy viscosity
  // f1 -- first blending coefficient
  // f2 -- second blending coefficient

  // calculate blending functions
  const auto alpha1 = this->Alpha1(state, suth, wallDist);
  const auto alpha2 = this->Alpha2(state, suth, mu, wallDist);
  const auto cdkw = this->CDkw(state, kGrad, wGrad);
  const auto alpha3 = this->Alpha3(state, wallDist, cdkw);
  f1 = this->F1(alpha1, alpha2, alpha3);
  f2 = this->F2(alpha1, alpha2);

  // return limited eddy viscosity
  mut = this->EddyVisc(state, velGrad, suth, f2);
}


/* member function to calculate the spectral radius of the source jacobian
   Source jacobian from Wilcox used.
   
   ds_dq = [ -2 * betaStar * w   0
                   0            -2 * beta * w]
   
   This is a diagonal matrix so eigenvalues are trivial. Since betaStar is
   always larger than beta, this eigenvalue is used
*/
double turbKWSst::SrcSpecRad(const primVars &state,
                             const sutherland &suth,
                             const double &vol) const {
  // state -- primative variables
  // suth -- sutherland's law for viscosity

  // return spectral radius scaled for nondimensional equations
  return -2.0 * betaStar_ * state.Omega() * vol * suth.InvNondimScaling();
}


// member function to calculate viscous spectral radius
// dfv_dq = [ (area / vol) * (nu + sigmaStar * nut)    0
//                           0                (area / vol) * (nu + sigma * nut)]
double turbKWSst::ViscSpecRad(const primVars &state,
                              const unitVec3dMag<double> &fAreaL,
                              const unitVec3dMag<double> &fAreaR,
                              const double &mu, const sutherland &suth,
                              const double &vol, const double &mut,
                              const double &f1) const {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // mu -- laminar viscosity
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  const auto sigmaK = this->SigmaK(f1);

  return suth.NondimScaling() * fMag * fMag / (vol * state.Rho()) *
      (mu + sigmaK * mut);
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

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

#include <cstdlib>    // exit()
#include <iostream>   // cout
#include <algorithm>  // max
#include "turbulence.hpp"
#include "primitive.hpp"   // primitive
#include "arrayView.hpp"   // primitiveView
#include "transport.hpp"  // transport model
#include "matrix.hpp"     // squareMatrix

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
double turbModel::EddyViscosityNoLim(const primitiveView &state) const {
  return state.Rho() * state.Tke() / state.Omega();
}

// member function to return mean strain rate
tensor<double> turbModel::MeanStrainRate(const tensor<double> &vGrad) const {
  return 0.5 * (vGrad + vGrad.Transpose());
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
    const primitiveView &state, const tensor<double> &velGrad,
    const unique_ptr<transport> &trans, const double &mut) const {
  // state -- primitive variables
  // velGrad -- velocity gradient tensor
  // trans -- viscous transport model
  // mut -- turbulent eddy viscosity

  const auto lambda = trans->Lambda(mut);  // 2nd eddy viscosity
  tensor<double> I;
  I.Identity();
  return lambda * velGrad.Trace() * I + mut * (velGrad + velGrad.Transpose())
      - 2.0 / 3.0 * state.Rho() * state.Tke() * I;
}

// Member function for calculation of tke production
// P = tau_ij * velGrad_ij
double turbModel::ReynoldsStressDDotVelGrad(const primitiveView &state,
                                            const tensor<double> &velGrad,
                                            const unique_ptr<transport> &trans,
                                            const double &mut) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // trans -- viscous transport model
  // mut -- turbulent eddy viscosity
  const auto tau = this->BoussinesqReynoldsStress(state, velGrad, trans, mut);
  return tau.DoubleDotTrans(velGrad);
}

// member function for destruction of tke
// Dk = rho * k * w
double turbModel::TkeDestruction(const primitiveView &state,
                                 const double &phi) const {
  // state -- primitive variables
  return state.Rho() * state.Tke() * state.Omega() * phi;
}

// member function for destruction of omega
// Dw = rho * w * w
double turbModel::OmegaDestruction(const primitiveView &state) const {
  // state -- primitive variables
  return state.Rho() * state.Omega() * state.Omega();
}

// member function for cross diffusion term without coefficient
// CD = rho / w * kGrad (dot) wGrad
// In the above equation kGrad and wGrad are the tke and omega gradients.
double turbModel::CrossDiffusion(const primitiveView &state,
                                 const vector3d<double> &kGrad,
                                 const vector3d<double> &wGrad) const {
  // state -- primitive variables
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  return state.Rho() / state.Omega() * kGrad.DotProd(wGrad);
}

// member function to calculate inviscid flux jacobian
// v = vel (dot) area
// df_dq = [v +/- |v|     0     ]
//         [   0       v +/- |v|]
squareMatrix turbModel::InvJac(const primitiveView &state,
                               const unitVec3dMag<double> &fArea,
                               const bool &positive) const {
  // state -- primitive variables at face
  // fArea -- face area
  // positive -- flag to determine whether to add/subtract dissipation

  return positive ? 0.5 * (this->InviscidConvJacobian(state, fArea) +
                           this->InviscidDissJacobian(state, fArea)) :
      0.5 * (this->InviscidConvJacobian(state, fArea) -
             this->InviscidDissJacobian(state, fArea));
}

// function to calculate convective portion of inviscid flux jacobian
squareMatrix turbModel::InviscidConvectiveJacobian(
    const primitiveView &state, const unitVec3dMag<double> &fArea) const {
  // state -- primitive variables at face
  // fArea -- face area

  const auto velNorm = state.Velocity().DotProd(fArea.UnitVector());
  const auto diag = velNorm * fArea.Mag();
  squareMatrix jacobian(2);
  jacobian(0, 0) = diag;
  jacobian(1, 1) = diag;
  return jacobian;
}

// function to calculate dissipative portion of inviscid flux jacobian
squareMatrix turbModel::InviscidDissipationJacobian(
    const primitiveView &state, const unitVec3dMag<double> &fArea) const {
  // state -- primitive variables at face
  // fArea -- face area

  const auto velNorm = state.Velocity().DotProd(fArea.UnitVector());
  const auto diag = fabs(velNorm) * fArea.Mag();
  squareMatrix jacobian(2);
  jacobian(0, 0) = diag;
  jacobian(1, 1) = diag;
  return jacobian;
}


// member function to calculate inviscid spectral radius
// df_dq = [vel (dot) area   0
//                0          vel (dot) area]
double turbModel::InviscidCellSpectralRadius(
    const primitiveView &state, const unitVec3dMag<double> &fAreaL,
    const unitVec3dMag<double> &fAreaR) const {
  // state -- primitive variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  auto normAvg =
      (0.5 * (fAreaL.UnitVector() + fAreaR.UnitVector())).Normalize();
  auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  return fabs(state.Velocity().DotProd(normAvg)) * fMag;
}

double turbModel::InviscidFaceSpectralRadius(const primitiveView &state,
                                             const unitVec3dMag<double> &fArea,
                                             const bool &positive) const {
  // state -- primitive variables
  // fArea -- face area
  // positive -- add or subtract dissipation term

  const auto velNorm = state.Velocity().DotProd(fArea.UnitVector());
  // returning absolute value because it is expected that spec rad is positive
  return positive ? 0.5 * fArea.Mag() * fabs(velNorm + fabs(velNorm)) :
      0.5 * fArea.Mag() * fabs(velNorm - fabs(velNorm));
  return InviscidFaceSpectralRadius(state, fArea, positive);
}

// member function to calculate viscous flux jacobian for models with no
// viscous contribution
squareMatrix turbModel::ViscJac(const primitiveView &state,
                                const unitVec3dMag<double> &fArea,
                                const double &mu,
                                const unique_ptr<transport> &trans,
                                const double &dist, const double &mut,
                                const double &f1) const {
  // state -- primitive variables
  // fAreaL -- face area for left face
  // mu -- laminar viscosity
  // trans -- unique_ptr<transport>'s law for viscosity
  // dist -- distance from cell center to cell center across face
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  return squareMatrix();
}

// member function to calculate turbulence source terms
squareMatrix turbModel::CalcTurbSrc(
    const primitiveView &state, const tensor<double> &velGrad,
    const vector3d<double> &kGrad, const vector3d<double> &wGrad,
    const unique_ptr<transport> &trans, const double &vol,
    const double &turbVisc, const double &f1, const double &f2,
    const double &width, vector<double> &turbSrc) const {
  // set source terms to zero
  for (auto &val : turbSrc) {
    val = 0.0;
  }

  // return source jacobian spectral radius
  return this->TurbSrcJac(state, 0.0, trans, vol);
}

squareMatrix turbModel::TurbSrcJac(const primitiveView &state, 
                                   const double &beta,
                                   const unique_ptr<transport> &trans,
                                   const double &vol,
                                   const double &phi) const {
  return squareMatrix();
}

// -------------------------------------------------------------------------
// member functions for the turbNone class

// member function to print out turbulence variables
void turbNone::Print() const {
  cout << "No Turbulence Model" << endl;
}

squareMatrix turbNone::InvJac(const primitiveView &state,
                              const unitVec3dMag<double> &fArea,
                              const bool &positive) const {
  // state -- primitive variables at face
  // fArea -- face area
  // positive -- flag to determine whether to add/subtract spectral radius

  return squareMatrix();
}

squareMatrix turbNone::InviscidConvectiveJacobian(
    const primitiveView &state, const unitVec3dMag<double> &fArea) const {
  // state -- primitive variables at face
  // fArea -- face area

  return squareMatrix();
}

squareMatrix turbNone::InviscidDissipationJacobian(
    const primitiveView &state, const unitVec3dMag<double> &fArea) const {
  // state -- primitive variables at face
  // fArea -- face area

  return squareMatrix();
}

// ---------------------------------------------------------------------
// K-Omega Wilcox member functions

// member function to return the eddy viscosity calculated with the stress
// limiter
double turbKWWilcox::EddyVisc(const primitive &state,
                              const tensor<double> &vGrad,
                              const unique_ptr<transport> &trans,
                              const double &f2, const double &length) const {
  // state -- primitive variables
  // vGrad -- velocity gradient
  // trans -- viscous transport model
  // f2 -- SST blending coefficient (not used in Wilcox K-W)
  // length -- length scale (not used in Wilcox K-W)
  return state.Rho() * state.Tke() / this->OmegaTilda(state, vGrad, trans);
}

// member function to calculate the cross diffusion coefficient term
double turbKWWilcox::SigmaD(const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const {
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  return kGrad.DotProd(wGrad) <= 0.0 ? 0.0 : sigmaD0_;
}

// member function to calculate coefficient for omega destruction
double turbKWWilcox::Beta(const primitiveView &state,
                          const tensor<double> &velGrad,
                          const unique_ptr<transport> &trans) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // trans -- viscous transport model
  return beta0_ * this->FBeta(state, velGrad, trans);
}

// member function to calculate coefficient used in beta calculation
double turbKWWilcox::FBeta(const primitiveView &state,
                           const tensor<double> &velGrad,
                           const unique_ptr<transport> &trans) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // trans -- viscous transport model
  auto xw = this->Xw(state, velGrad, trans);
  return (1.0 + 85.0 * xw) / (1.0 + 100.0 * xw);
}

// member function to calculate vortex stretching coefficient
// used in fbeta calculation
double turbKWWilcox::Xw(const primitiveView &state,
                        const tensor<double> &velGrad,
                        const unique_ptr<transport> &trans) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // trans -- viscous transport model
  const auto vorticity = 0.5 * (velGrad - velGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  // vorticity is asymmetric but vorticity * vorticity is symmetric
  return fabs(
      (vorticity.MatMult(vorticity)).DoubleDotTrans(this->StrainKI(velGrad))
      / pow(betaStar_ * state.Omega(), 3.0) )
      * pow(trans->NondimScaling(), 3.0);
}

// member function to calculate term used in calculation of Xw
tensor<double> turbKWWilcox::StrainKI(const tensor<double> &velGrad) const {
  // velGrad -- velocity gradient
  tensor<double> I;
  I.Identity();
  return 0.5 * (velGrad + velGrad.Transpose() - velGrad.Trace() * I);
}

// member function to calculate adjusted omega used in eddy viscosity limiter
double turbKWWilcox::OmegaTilda(const primitive &state,
                                const tensor<double> &velGrad,
                                const unique_ptr<transport> &trans) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // trans -- viscous transport model

  tensor<double> I;
  I.Identity();
  const auto sHat = 0.5 * (velGrad + velGrad.Transpose()) - 1.0 / 3.0 *
      velGrad.Trace() * I;

  // using DoubleDotTrans instead of DoubleDot for speed
  // since tensors are symmetric, result is the same
  return std::max(state.Omega(), trans->NondimScaling() * clim_ *
                  sqrt(2.0 * sHat.DoubleDotTrans(sHat) / betaStar_));
}

// member function to calculate turbulence source terms and return source
// jacobian
squareMatrix turbKWWilcox::CalcTurbSrc(const primitiveView &state,
                                       const tensor<double> &velGrad,
                                       const vector3d<double> &kGrad,
                                       const vector3d<double> &wGrad,
                                       const unique_ptr<transport> &trans,
                                       const double &vol,
                                       const double &mut, const double &f1,
                                       const double &f2, const double &width,
                                       vector<double> &turbSrc) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // turbSrc -- source terms

  MSG_ASSERT(turbSrc.size() == 2, "turbulence source vector is wrong size");

  // calculate tke destruction
  const auto tkeDest = trans->InvNondimScaling() * betaStar_ *
      this->TkeDestruction(state);

  // calculate omega destruction
  const auto beta = this->Beta(state, velGrad, trans);
  const auto omgDest = trans->InvNondimScaling() * beta *
      this->OmegaDestruction(state);

  // calculate tke production
  auto tkeProd = trans->NondimScaling() *
      this->ReynoldsStressDDotVelGrad(state, velGrad, trans, mut);
  tkeProd = max(tkeProd, 0.0);

  // calculate omega production
  auto omgProd = gamma_ * state.Omega() / state.Tke() * tkeProd;
  omgProd = max(omgProd, 0.0);

  // calculate omega cross diffusion
  const auto omgCd = trans->NondimScaling() * this->SigmaD(kGrad, wGrad) *
      this->CrossDiffusion(state, kGrad, wGrad);

  // assign source term values
  turbSrc[0] = tkeProd - tkeDest;
  turbSrc[1] = omgProd - omgDest + omgCd;

  // return source jacobian
  return this->TurbSrcJac(state, beta, trans, vol);
}

// member function to calculate the eddy viscosity, and the blending
// coefficients.
void turbKWWilcox::EddyViscAndBlending(const primitive &state,
                                       const tensor<double> &velGrad,
                                       const vector3d<double> &kGrad,
                                       const vector3d<double> &wGrad,
                                       const double &mu,
                                       const double &wallDist,
                                       const unique_ptr<transport> &trans,
                                       const double &length,
                                       double &mut, double &f1,
                                       double &f2) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // mu -- laminar viscosity
  // walldist -- distance to nearest viscous wall
  // trans -- viscous transport model
  // length -- cell length scale
  // mut -- turbulent eddy viscosity
  // f1 -- first blending coefficient
  // f2 -- second blending coefficient

  f1 = 1.0;
  f2 = 0.0;

  // return eddy viscosity
  mut = this->EddyVisc(state, velGrad, trans, f2, length);
}


/* member function to calculate the spectral radius of the source jacobian
   Source jacobian from Wilcox used.
   
   ds_dq = [ -2 * betaStar * w   0
                   0            -2 * beta * w]
   
   This is a diagonal matrix so eigenvalues are trivial. Since betaStar is
   always larger than beta, this eigenvalue is used
*/
double turbKWWilcox::SrcSpecRad(const primitiveView &state,
                                const unique_ptr<transport> &trans,
                                const double &vol, const double &phi) const {
  // state -- primitive variables
  // trans -- viscous transport model
  // vol -- cell volume

  // return spectral radius scaled for nondimensional equations
  return -2.0 * betaStar_ * state.Omega() * vol * trans->InvNondimScaling();
}

squareMatrix turbKWWilcox::TurbSrcJac(const primitiveView &state,
                                      const double &beta,
                                      const unique_ptr<transport> &trans,
                                      const double &vol,
                                      const double &phi) const {
  // state -- primitive variables
  // beta -- destruction coefficient for omega equation
  // trans -- viscous transport model
  // vol -- cell volume
  // phi -- factor to reduce tke destruction for des

  squareMatrix jac(2);
  jac(0, 0) =
      -2.0 * betaStar_ * state.Omega() * vol * trans->InvNondimScaling();
  jac(1, 1) = -2.0 * beta * state.Omega() * vol * trans->InvNondimScaling();

  // return jacobian scaled for nondimensional equations
  return jac;
}

// member function to calculate viscous flux jacobian
// dfv_dq = [ (nu + sigmaStar * nut) / dist           0               ]
//          [             0                  (nu + sigma * nut) / dist]
squareMatrix turbKWWilcox::ViscJac(const primitiveView &state,
                                   const unitVec3dMag<double> &fArea,
                                   const double &mu,
                                   const unique_ptr<transport> &trans,
                                   const double &dist, const double &mut,
                                   const double &f1) const {
  // state -- primitive variables
  // fAreaL -- face area for left face
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center across face
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  auto length = fArea.Mag() / dist;
  squareMatrix jacobian(2);
  // Wilcox method uses unlimited eddy viscosity
  jacobian(0, 0) = trans->NondimScaling() * length / state.Rho() *
                   (mu + this->SigmaK(f1) * this->EddyViscosityNoLim(state));
  jacobian(1, 1) = trans->NondimScaling() * length / state.Rho() *
                   (mu + this->SigmaW(f1) * this->EddyViscosityNoLim(state));
  return jacobian;
}

// member function to calculate viscous spectral radius
// dfv_dq = [ (area / vol) * (nu + sigmaStar * nut)    0
//                           0                (area / vol) * (nu + sigma * nut)]
double turbKWWilcox::ViscousCellSpectralRadius(
    const primitiveView &state, const unitVec3dMag<double> &fAreaL,
    const unitVec3dMag<double> &fAreaR, const double &mu,
    const unique_ptr<transport> &trans, const double &vol, const double &mut,
    const double &f1) const {
  // state -- primitive variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  const auto length = fMag * fMag / vol;
  // Wilcox method uses unlimited eddy viscosity
  return trans->NondimScaling() * length / state.Rho() *
         (mu + this->SigmaK(f1) * this->EddyViscosityNoLim(state));
}

double turbKWWilcox::ViscousFaceSpectralRadius(
    const primitiveView &state, const unitVec3dMag<double> &fArea,
    const double &mu, const unique_ptr<transport> &trans, const double &dist,
    const double &mut, const double &f1) const {
  // state -- primitive variables
  // fArea -- face area
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  const auto length = fArea.Mag() / dist;
  // Wilcox method uses unlimited eddy viscosity
  return trans->NondimScaling() * length / state.Rho() *
         (mu + this->SigmaK(f1) * this->EddyViscosityNoLim(state));
}

double turbKWWilcox::TurbLengthScale(const primitiveView &state,
                                     const unique_ptr<transport> &trans) const {
  return sqrt(state.Tke()) / (betaStar_ * state.Omega()) *
         trans->NondimScaling();
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
double turbKWSst::EddyVisc(const primitive &state, const tensor<double> &vGrad,
                           const unique_ptr<transport> &trans, const double &f2,
                           const double &length) const {
  // state -- primitive variables
  // vGrad -- velocity gradient
  // trans -- viscous transport model
  // f2 -- SST blending coefficient
  // length -- length scale (not used for SST)

  const auto strainRate = this->MeanStrainRate(vGrad);

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  const auto meanStrainRate = sqrt(2.0 * strainRate.DoubleDotTrans(strainRate));
  return state.Rho() * a1_ * state.Tke() /
      max(a1_ * state.Omega(), trans->NondimScaling() * meanStrainRate * f2);
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
double turbKWSst::Alpha1(const primitive &state,
                         const unique_ptr<transport> &trans,
                         const double &wallDist) const {
  return trans->NondimScaling() * sqrt(state.Tke()) /
         (betaStar_ * state.Omega() * (wallDist + EPS));
}

// member function to calculate blending term
double turbKWSst::Alpha2(const primitive &state,
                         const unique_ptr<transport> &trans, const double &mu,
                         const double &wallDist) const {
  return trans->NondimScaling() * trans->NondimScaling() * 500.0 * mu /
         ((wallDist + EPS) * (wallDist + EPS) * state.Rho() * state.Omega());
}

// member function to calculate blending term
double turbKWSst::Alpha3(const primitive &state, const double &wallDist,
                         const double &cdkw) const {
  return 4.0 * state.Rho() * sigmaW2_ * state.Tke() /
         (cdkw * (wallDist + EPS) * (wallDist + EPS));
}

// member function to calculate turbulence source terms and source jacobian
squareMatrix turbKWSst::CalcTurbSrc(
    const primitiveView &state, const tensor<double> &velGrad,
    const vector3d<double> &kGrad, const vector3d<double> &wGrad,
    const unique_ptr<transport> &trans, const double &vol, const double &mut,
    const double &f1, const double &f2, const double &width, 
    vector<double> &turbSrc) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // turbSrc -- source terms

  MSG_ASSERT(turbSrc.size() == 2, "turbulence source vector is wrong size");

  // calculate cross diffusion coefficient
  const auto cdkw = this->CDkw(state, kGrad, wGrad);

  // calculate blended coefficients
  const auto gamma = this->BlendedCoeff(gamma1_, gamma2_, f1);
  const auto beta = this->BlendedCoeff(beta1_, beta2_, f1);

  // calculate tke destruction
  const auto tkeDest = trans->InvNondimScaling() * betaStar_ *
      this->TkeDestruction(state);

  // calculate omega destruction
  const auto omgDest = trans->InvNondimScaling() * beta *
      this->OmegaDestruction(state);

  // calculate tke production
  auto tkeProd = min(trans->NondimScaling() *
                     this->ReynoldsStressDDotVelGrad(state, velGrad, trans, mut),
                     kProd2Dest_ * tkeDest);
  tkeProd = max(tkeProd, 0.0);

  // calculate omega production
  auto omgProd = gamma * state.Rho() / mut * tkeProd;
  omgProd = max(omgProd, 0.0);

  // calculate omega cross diffusion
  // Using CDkw instead of whole cross diffusion term
  // Both Loci/CHEM and SU2 use this implementation
  const auto omgCd = trans->NondimScaling() * (1.0 - f1) * cdkw;

  // assign source term values
  turbSrc[0] = tkeProd - tkeDest;
  turbSrc[1] = omgProd - omgDest + omgCd;

  // return spectral radius of source jacobian
  return this->TurbSrcJac(state, beta, trans, vol);
}

// member function to calculate the eddy viscosity, and the blending
// coefficients.
void turbKWSst::EddyViscAndBlending(const primitive &state,
                                    const tensor<double> &velGrad,
                                    const vector3d<double> &kGrad,
                                    const vector3d<double> &wGrad,
                                    const double &mu,
                                    const double &wallDist,
                                    const unique_ptr<transport> &trans,
                                    const double &length,
                                    double &mut, double &f1,
                                    double &f2) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // mu -- laminar viscosity
  // walldist -- distance to nearest viscous wall
  // trans -- viscous transport model
  // length -- length scale
  // mut -- turbulent eddy viscosity
  // f1 -- first blending coefficient
  // f2 -- second blending coefficient

  // calculate blending functions
  const auto alpha1 = this->Alpha1(state, trans, wallDist);
  const auto alpha2 = this->Alpha2(state, trans, mu, wallDist);
  const auto cdkw = this->CDkw(state, kGrad, wGrad);
  const auto alpha3 = this->Alpha3(state, wallDist, cdkw);
  f1 = this->F1(alpha1, alpha2, alpha3);
  f2 = this->F2(alpha1, alpha2);

  // return limited eddy viscosity
  mut = this->EddyVisc(state, velGrad, trans, f2, length);
}


/* member function to calculate the spectral radius of the source jacobian
   Source jacobian from Wilcox used.
   
   ds_dq = [ -2 * betaStar * w   0
                   0            -2 * beta * w]
   
   This is a diagonal matrix so eigenvalues are trivial. Since betaStar is
   always larger than beta, this eigenvalue is used
*/
double turbKWSst::SrcSpecRad(const primitiveView &state,
                             const unique_ptr<transport> &trans,
                             const double &vol, const double &phi) const {
  // state -- primitive variables
  // trans -- unique_ptr<transport>'s law for viscosity

  // return spectral radius scaled for nondimensional equations
  return -2.0 * betaStar_ * state.Omega() * vol * trans->InvNondimScaling();
}

squareMatrix turbKWSst::TurbSrcJac(const primitiveView &state,
                                   const double &beta,
                                   const unique_ptr<transport> &trans,
                                   const double &vol,
                                   const double &phi) const {
  // state -- primitive variables
  // beta -- destruction coefficient for omega equation
  // trans -- viscous transport model
  // vol -- cell volume
  // phi -- factor to reduce tke destruction for des

  squareMatrix jac(2);
  jac(0, 0) = -2.0 * betaStar_ * state.Omega() * phi * vol *
      trans->InvNondimScaling();
  jac(1, 1) = -2.0 * beta * state.Omega() * vol * trans->InvNondimScaling();

  // return jacobian scaled for nondimensional equations
  return jac;
}

// member function to calculate viscous flux jacobian
// dfv_dq = [ (nu + sigmaStar * nut) / dist           0               ]
//          [             0                  (nu + sigma * nut) / dist]
squareMatrix turbKWSst::ViscJac(const primitiveView &state,
                                const unitVec3dMag<double> &fArea,
                                const double &mu,
                                const unique_ptr<transport> &trans,
                                const double &dist, const double &mut,
                                const double &f1) const {
  // state -- primitive variables
  // fArea -- face area
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center across face
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  auto length = fArea.Mag() / dist;
  squareMatrix jacobian(2);
  jacobian(0, 0) = trans->NondimScaling() * length / state.Rho() *
                   (mu + this->SigmaK(f1) * mut);
  jacobian(1, 1) = trans->NondimScaling() * length / state.Rho() *
                   (mu + this->SigmaW(f1) * mut);
  return jacobian;
}

// member function to calculate viscous spectral radius
// dfv_dq = [ (area / vol) * (nu + sigmaStar * nut)    0
//                           0                (area / vol) * (nu + sigma * nut)]
double turbKWSst::ViscousCellSpectralRadius(
    const primitiveView &state, const unitVec3dMag<double> &fAreaL,
    const unitVec3dMag<double> &fAreaR, const double &mu,
    const unique_ptr<transport> &trans, const double &vol, const double &mut,
    const double &f1) const {
  // state -- primitive variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient

  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  const auto length = fMag * fMag / vol;
  return trans->NondimScaling() * length / state.Rho() *
         (mu + this->SigmaK(f1) * mut);
}

double turbKWSst::ViscousFaceSpectralRadius(
    const primitiveView &state, const unitVec3dMag<double> &fArea,
    const double &mu, const unique_ptr<transport> &trans, const double &dist,
    const double &mut, const double &f1) const {
  // state -- primitive variables
  // fArea -- face area
  // mu -- laminar viscosity
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  const auto length = fArea.Mag() / dist;
  return trans->NondimScaling() * length / state.Rho() *
         (mu + this->SigmaK(f1) * mut);
}

double turbKWSst::TurbLengthScale(const primitiveView &state,
                                  const unique_ptr<transport> &trans) const {
  return sqrt(state.Tke()) / (betaStar_ * state.Omega()) *
         trans->NondimScaling();
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

double turbSstDes::Phi(const primitiveView &state, const double &cdes,
                       const double &width, const double &f2,
                       const unique_ptr<transport> &trans) const {
  return std::max(
      (1.0 - f2) * this->TurbLengthScale(state, trans) / (cdes * width), 1.0);
}

// member function to calculate turbulence source terms and source jacobian
squareMatrix turbSstDes::CalcTurbSrc(
    const primitiveView &state, const tensor<double> &velGrad,
    const vector3d<double> &kGrad, const vector3d<double> &wGrad,
    const unique_ptr<transport> &trans, const double &vol, const double &mut,
    const double &f1, const double &f2, const double &width, 
    vector<double> &turbSrc) const {
  // state -- primitive variables
  // velGrad -- velocity gradient
  // kGrad -- tke gradient
  // wGrad -- omega gradient
  // trans -- viscous transport model
  // vol -- cell volume
  // mut -- turbulent viscosity
  // f1 -- first blending coefficient
  // turbSrc -- source terms

  MSG_ASSERT(turbSrc.size() == 2, "turbulence source vector is wrong size");

  // calculate cross diffusion coefficient
  const auto cdkw = this->CDkw(state, kGrad, wGrad);

  // calculate blended coefficients
  const auto gamma = this->BlendedCoeff(this->Gamma1(), this->Gamma2(), f1);
  const auto beta = this->BlendedCoeff(this->Beta1(), this->Beta2(), f1);
  const auto cdes = this->BlendedCoeff(cdes1_, cdes2_, f1);

  // calculate tke destruction
  const auto phi = this->Phi(state, cdes, width, f2, trans);
  const auto tkeDest =
      trans->InvNondimScaling() * this->TkeDestruction(state, phi);

  // calculate omega destruction
  const auto omgDest = trans->InvNondimScaling() * beta *
      this->OmegaDestruction(state);

  // calculate tke production
  auto tkeProd = min(trans->NondimScaling() *
                     this->ReynoldsStressDDotVelGrad(state, velGrad, trans, mut),
                     this->TkeProd2DestRatio() * tkeDest);
  tkeProd = max(tkeProd, 0.0);

  // calculate omega production
  auto omgProd = gamma * state.Rho() / mut * tkeProd;
  omgProd = max(omgProd, 0.0);

  // calculate omega cross diffusion
  // Using CDkw instead of whole cross diffusion term
  // Both Loci/CHEM and SU2 use this implementation
  const auto omgCd = trans->NondimScaling() * (1.0 - f1) * cdkw;

  // assign source term values
  turbSrc[0] = tkeProd - tkeDest;
  turbSrc[1] = omgProd - omgDest + omgCd;

  // return spectral radius of source jacobian
  return this->TurbSrcJac(state, beta, trans, vol, phi);
}


double turbSstDes::SrcSpecRad(const primitiveView &state,
                              const unique_ptr<transport> &trans,
                              const double &vol, const double &phi) const {
  // state -- primitive variables
  // trans -- viscous transport model

  // return spectral radius scaled for nondimensional equations
  auto beta = this->Beta2();  // using beta2 b/c it is larger than beta1
  auto jac = this->TurbSrcJac(state, beta, trans, vol, phi);
  return -1.0 * jac.MaxAbsValOnDiagonal();
}


// member function to print out turbulence variables
void turbSstDes::Print() const {
  cout << "Eddy Viscosity Method: " << this->EddyViscMethod() << endl;
  cout << "A1: " << this->A1() << endl;
  cout << "Beta*: " << this->BetaStar() << endl;
  cout << "Turbulent Prandtl Number: " << this->TurbPrandtlNumber() << endl;
  cout << "Production to Destruction Ratio: " << this->TkeProd2DestRatio()
       << endl;

  cout << "Sigma K1: " << this->SigmaK1() << endl;
  cout << "Sigma W1: " << this->SigmaW1() << endl;
  cout << "Beta1: " << this->Beta1() << endl;
  cout << "Gamma1: " << this->Gamma1() << endl;
  cout << "CDES1: " << this->CDes1() << endl;

  cout << "Sigma K2: " << this->SigmaK2() << endl;
  cout << "Sigma W2: " << this->SigmaW2() << endl;
  cout << "Beta2: " << this->Beta2() << endl;
  cout << "Gamma2: " << this->Gamma2() << endl;
  cout << "CDES2: " << this->CDes2() << endl;
}

// member function to print out turbulence variables
void turbWale::Print() const {
  cout << "Eddy Viscosity Method: " << this->EddyViscMethod() << endl;
  cout << "Cw: " << this->Cw() << endl;
}

// member function to determine eddy viscosity
double turbWale::EddyVisc(const primitive &state, const tensor<double> &vGrad,
                          const unique_ptr<transport> &trans, const double &f2,
                          const double &length) const {
  // state -- primitive variables
  // vGrad -- velocity gradient
  // trans -- viscous transport model
  // f2 -- SST blending coefficient (not used in WALE)
  // length -- cell length

  const auto strainRate = this->MeanStrainRate(vGrad);
  const auto sigmaD = this->SigmaD(vGrad);
  const auto sigmaDDoubleDot = sigmaD.DoubleDotTrans(sigmaD);

  const auto velGradTerm = pow(sigmaDDoubleDot, 1.5) /
      (pow(strainRate.DoubleDotTrans(strainRate), 2.5) +
       pow(sigmaDDoubleDot, 1.25) + EPS);

  const auto lengthScale = pow(cw_ * length, 2.0);

  return lengthScale * velGradTerm;
}

tensor<double> turbWale::SigmaD(const tensor<double> &vGrad) const {
  const auto vGradSq = vGrad.MatMult(vGrad);

  tensor<double> I;
  I.Identity();

  return 0.5 * (vGradSq + vGradSq.Transpose()) - 1.0 / 3.0 * I * vGradSq.Trace();
}

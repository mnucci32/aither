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

#include <cstdlib>    // exit()
#include <iostream>   // cout
#include <memory>
#include "wallLaw.hpp"
#include "primVars.hpp"  // primVars
#include "eos.hpp"       // idealGas
#include "utility.hpp"   // TauShear
#include "turbulence.hpp"
#include "wallData.hpp"

using std::cout;
using std::endl;
using std::cerr;

// -------------------------------------------------------------------------
wallVars wallLaw::AdiabaticBCs(const vector3d<double> &area,
                               const vector3d<double> &velWall,
                               const idealGas &eos, const sutherland &suth,
                               const unique_ptr<turbModel> &turb) {
  // initialize wallVars
  wallVars wVars;
  wVars.heatFlux_ = 0.0;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  this->CalcRecoveryFactor(eos);
  // this is correct form of crocco-busemann, typo in Nichols & Nelson (2004)
  auto tW = state_.Temperature(eos) +
            0.5 * recoveryFactor_ * velTanMag * velTanMag / eos.SpecificHeat();
  // set wall properties
  this->SetWallVars(tW, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate constants
    this->UpdateGamma(eos);
    this->UpdateConstants(wVars.heatFlux_);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    wVars.yplus_ = yplus;
    return this->CalcYplusRoot(yplus);
  };

  // iteratively solve for y+
  FindRoot(func, 1.0e1, 1.0e4, 1.0e-8);

  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.temperature_ = tW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  return wVars;
}

wallVars wallLaw::HeatFluxBCs(const vector3d<double> &area,
                              const vector3d<double> &velWall,
                              const idealGas &eos, const sutherland &suth,
                              const unique_ptr<turbModel> &turb,
                              const double &heatFluxW) {
  // initialize wallVars
  wallVars wVars;
  wVars.heatFlux_ = heatFluxW;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  this->CalcRecoveryFactor(eos);
  
  // set wall properties - guess wall temperature equals interior temperature
  wVars.temperature_ = state_.Temperature(eos);
  this->SetWallVars(wVars.temperature_, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall temperature from croco-busemann
    wVars.temperature_ = this->CalcWallTemperature(eos, wVars.heatFlux_);
    this->SetWallVars(wVars.temperature_, eos, suth);
    this->UpdateGamma(eos);
    this->UpdateConstants(wVars.heatFlux_);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    wVars.yplus_ = yplus;
    return this->CalcYplusRoot(yplus);
  };
  
  // iteratively solve for y+
  FindRoot(func, 1.0e1, 1.0e4, 1.0e-8);
    
  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  return wVars;
}

wallVars wallLaw::IsothermalBCs(const vector3d<double> &area,
                                const vector3d<double> &velWall,
                                const idealGas &eos, const sutherland &suth,
                                const unique_ptr<turbModel> &turb,
                                const double &tW) {
  // initialize wallVars
  wallVars wVars;
  wVars.temperature_ = tW;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall properties
  this->CalcRecoveryFactor(eos);
  this->SetWallVars(wVars.temperature_, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall heat flux from croco-busemann equation
    this->UpdateGamma(eos);
    wVars.heatFlux_ = this->CalcHeatFlux(eos);
    // calculate constants
    this->UpdateConstants(wVars.heatFlux_);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    wVars.yplus_ = yplus;
    return this->CalcYplusRoot(yplus);
  };

  // iteratively solve for y+
  FindRoot(func, 1.0e1, 1.0e4, 1.0e-8);
  
  // calculate turbulent eddy viscosity from yplus
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  return wVars;
}

void wallLaw::UpdateGamma(const idealGas &eos) {
  // calculate constants
  gamma_ = recoveryFactor_ * uStar_ * uStar_ / (2.0 * eos.SpecificHeat() * tW_);
}

void wallLaw::UpdateConstants(const double &heatFluxW) {
  // calculate constants
  beta_ = heatFluxW * muW_ / (rhoW_ * tW_ * kW_ * uStar_);  // 0 for adiabatic
  q_ = sqrt(beta_ * beta_ + 4.0 * gamma_);
  phi_ = std::asin(-beta_ / q_);
}

void wallLaw::CalcYplusWhite() {
  yplusWhite_ =
      std::exp((vonKarmen_ / sqrt(gamma_)) *
               (std::asin((2.0 * gamma_ * uplus_ - beta_) / q_) - phi_)) *
      yplus0_;
}

double wallLaw::CalcHeatFlux(const idealGas &eos) const {
  // calculate wall heat flux from croco-busemann equation
  auto tmp =
      (state_.Temperature(eos) / tW_ - 1.0 + gamma_ * uplus_ * uplus_) / uplus_;
  return tmp * (rhoW_ * tW_ * kW_ * uStar_) / muW_;
}

double wallLaw::CalcWallTemperature(const idealGas &eos,
                                    const double &heatFluxW) const {
  return state_.Temperature(eos) +
         recoveryFactor_ * uStar_ * uStar_ * uplus_ * uplus_ /
             (2.0 * eos.SpecificHeat() +
              heatFluxW * muW_ / (rhoW_ * kW_ * uStar_));
}

void wallLaw::SetWallVars(const double &tW, const idealGas &eos,
                          const sutherland &suth) {
  tW_ = tW;
  rhoW_ = eos.DensityTP(tW_, state_.P());

  // get wall viscosity, conductivity from wall temperature
  muW_ = suth.EffectiveViscosity(tW_);
  kW_ = eos.Conductivity(muW_);
}

double wallLaw::CalcYplusRoot(const double &yplus) const {
  constexpr auto sixth = 1.0 / 6.0;
  const auto ku = vonKarmen_ * uplus_;
  return yplus - (uplus_ + yplusWhite_ -
                  yplus0_ * (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0)));
}

void wallLaw::EddyVisc(const idealGas &eos, const sutherland &suth) {
  const auto dYplusWhite =
      2.0 * yplusWhite_ * vonKarmen_ * sqrt(gamma_) / q_ *
      sqrt(std::max(1.0 - pow(2.0 * gamma_ * uplus_ - beta_, 2.0) / (q_ * q_),
                    0.0));
  const auto ku = vonKarmen_ * uplus_;
  mutW_ = muW_ * (1.0 + dYplusWhite -
                  vonKarmen_ * yplus0_ * (1.0 + ku + 0.5 * ku * ku)) -
          suth.EffectiveViscosity(state_.Temperature(eos));
  mutW_ = std::max(mutW_, 0.0);
}

void wallLaw::CalcVelocities(const double &yplus, const double &u) {
  // calculate u* and u+ from y+
  uplus_ = (wallDist_ * rhoW_ * u) / (muW_ * yplus);
  uStar_ = u / uplus_;
}

void wallLaw::CalcTurbVars(const unique_ptr<turbModel> &turb,
                           const idealGas &eos, const sutherland &suth,
                           double &kWall, double &wWall) {
  this->EddyVisc(eos, suth);
  // calculate turbulence variables from eddy viscosity
  auto wi = 6.0 * muW_ / (turb->WallBeta() * rhoW_ * wallDist_ * wallDist_);
  wi *= suth.NondimScaling();
  auto wo = uStar_ / (sqrt(turb->BetaStar()) * vonKarmen_ * wallDist_);
  wo *= suth.NondimScaling();
  wWall = sqrt(wi * wi + wo * wo);
  kWall = wWall * mutW_ / state_.Rho() * suth.InvNondimScaling();
}

void wallLaw::CalcRecoveryFactor(const idealGas &eos) {
  recoveryFactor_ = pow(eos.Prandtl(), 1.0 / 3.0);
}
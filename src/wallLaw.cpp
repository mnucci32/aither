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

using std::cout;
using std::endl;
using std::cerr;

// -------------------------------------------------------------------------
vector3d<double> wallLaw::AdiabaticBCs(const vector3d<double> &area,
                                       const vector3d<double> &velWall,
                                       const idealGas &eos,
                                       const sutherland &suth,
                                       const unique_ptr<turbModel> &turb,
                                       double &kWall, double &wWall) {
  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  constexpr auto heatFluxW = 0.0;

  // get wall temperature from crocco-busemann equation
  this->CalcRecoveryFactor(eos);
  auto tW =
      state_.Temperature(eos) /
      (1.0 +
       0.5 * (eos.Gamma() - 1.0) * recoveryFactor_ * velTanMag * velTanMag);
  // set wall properties
  this->SetWallVars(tW, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate constants
    this->UpdateGamma(eos);
    this->UpdateConstants(heatFluxW);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    return this->CalcYplusRoot(yplus);
  };
  
  // iteratively solve for y+
  FindRoot(func, 1.0e-5, 1.0e5, 1.0e-8);
    
  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, kWall, wWall);
  }

  return this->GhostVelocity(area);
}

vector3d<double> wallLaw::HeatFluxBCs(const vector3d<double> &area,
                                      const vector3d<double> &velWall,
                                      const idealGas &eos,
                                      const sutherland &suth,
                                      const unique_ptr<turbModel> &turb,
                                      const double &heatFluxW, double &tWall,
                                      double &kWall, double &wWall) {
  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  this->CalcRecoveryFactor(eos);
  
  // set wall properties - guess wall temperature equals interior temperature
  tWall = state_.Temperature(eos);
  this->SetWallVars(tWall, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall temperature from croco-busemann
    tWall = this->CalcWallTemperature(eos, heatFluxW);
    this->SetWallVars(tWall, eos, suth);
    this->UpdateGamma(eos);
    this->UpdateConstants(heatFluxW);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    return this->CalcYplusRoot(yplus);
  };
  
  // iteratively solve for y+
  FindRoot(func, 1.0e-5, 1.0e5, 1.0e-8);
    
  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, kWall, wWall);
  }

  return this->GhostVelocity(area);
}


vector3d<double> wallLaw::IsothermalBCs(const vector3d<double> &area,
                                        const vector3d<double> &velWall,
                                        const idealGas &eos,
                                        const sutherland &suth,
                                        const unique_ptr<turbModel> &turb,
                                        const double &tW, double &heatFluxW,
                                        double &kWall, double &wWall,
                                        double &mutWall) {
  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall properties
  this->CalcRecoveryFactor(eos);
  this->SetWallVars(tW, eos, suth);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall heat flux from croco-busemann equation
    this->UpdateGamma(eos);
    heatFluxW = this->CalcHeatFlux(eos);
    // calculate constants
    this->UpdateConstants(heatFluxW);
    // calculate y+ from White & Christoph
    this->CalcYplusWhite();
    // calculate root of y+ equation
    return this->CalcYplusRoot(yplus);
  };
  
  // iteratively solve for y+
  FindRoot(func, 1.0e-5, 1.0e5, 1.0e-8);
  
  // calculate turbulent eddy viscosity from yplus
  // use compressible form of equation (Nichols & Nelson 2004)
  // calculate turbulence variables from eddy viscosity
  if (isRANS_) {
    this->CalcTurbVars(turb, eos, suth, kWall, wWall);
    mutWall = mutW_;
  }

  return this->GhostVelocity(area);
}

vector3d<double> wallLaw::GhostVelocity(const vector3d<double> &area) const {
  const auto shearStressMag = this->ShearStressMag();
  const auto velInt = state_.Velocity();
  const auto velIntNormal = velInt.DotProd(area) * area;
  const auto velIntTan = velInt - velIntNormal;
  // use 2x wall distance as gradient length
  const auto velGhostTan =
      shearStressMag / (muW_ + mutW_) * 2.0 * wallDist_ + velIntTan;
  return velGhostTan - velIntNormal;
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
  uStar_ = yplus * muW_ / (rhoW_ * wallDist_);
  uplus_ = u / uStar_;
}

void wallLaw::CalcTurbVars(const unique_ptr<turbModel> &turb,
                           const idealGas &eos, const sutherland &suth,
                           double &kWall, double &wWall) {
  this->EddyVisc(eos, suth);
  // calculate turbulence variables from eddy viscosity
  // DEBUG -- check nondimensionalization
  auto wi = 6.0 * muW_ / (turb->WallBeta() * rhoW_ * wallDist_ * wallDist_);
  wi *= suth.NondimScaling();
  auto wo = uStar_ / (sqrt(turb->BetaStar()) * vonKarmen_ * wallDist_);
  wo *= suth.NondimScaling();
  wWall = sqrt(wi * wi + wo * wo);
  kWall = wWall * mutW_ / state_.Rho();
}

void wallLaw::CalcRecoveryFactor(const idealGas &eos) {
  recoveryFactor_ = pow(eos.Prandtl(), 1.0 / 3.0);
}
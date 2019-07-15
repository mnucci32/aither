/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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
#include "primitive.hpp"       // primitive
#include "physicsModels.hpp"
#include "utility.hpp"        // find root
#include "wallData.hpp"

using std::cout;
using std::endl;
using std::cerr;

// -------------------------------------------------------------------------
wallVars wallLaw::AdiabaticBCs(const vector3d<double> &area,
                               const vector3d<double> &velWall,
                               const vector<double> &wallMf,
                               const physics &phys, const bool &isLower) {
  // initialize wallVars
  wallVars wVars;
  wVars.heatFlux_ = 0.0;
  wVars.mf_ = wallMf;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  auto t = state_.Temperature(phys.EoS());
  this->CalcRecoveryFactor(phys.Thermodynamic(), t);
  // this is correct form of crocco-busemann, typo in Nichols & Nelson (2004)
  auto tW = t + 0.5 * recoveryFactor_ * velTanMag * velTanMag /
                    phys.Thermodynamic()->Cp(t, state_.MassFractions());
  // set wall properties
  this->SetWallVars(tW, phys);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate constants
    this->UpdateGamma(phys.Thermodynamic(), tW);
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
    this->CalcTurbVars(phys, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.temperature_ = tW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  if (!isLower) {
    wVars.shearStress_ *= -1.0;
  }
  return wVars;
}

wallVars wallLaw::HeatFluxBCs(const vector3d<double> &area,
                              const vector3d<double> &velWall,
                              const vector<double> &wallMf, const physics &phys,
                              const double &heatFluxW, const bool &isLower) {
  // initialize wallVars
  wallVars wVars;
  wVars.heatFlux_ = heatFluxW;
  wVars.mf_ = wallMf;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  wVars.temperature_ = state_.Temperature(phys.EoS());
  this->CalcRecoveryFactor(phys.Thermodynamic(), wVars.temperature_);
  
  // set wall properties - guess wall temperature equals interior temperature
  this->SetWallVars(wVars.temperature_, phys);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall temperature from croco-busemann
    wVars.temperature_ =
        this->CalcWallTemperature(phys, wVars.heatFlux_);
    this->SetWallVars(wVars.temperature_, phys);
    this->UpdateGamma(phys.Thermodynamic(), wVars.temperature_);
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
    this->CalcTurbVars(phys, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  if (!isLower) {
    wVars.shearStress_ *= -1.0;
  }
  return wVars;
}

wallVars wallLaw::IsothermalBCs(const vector3d<double> &area,
                                const vector3d<double> &velWall,
                                const vector<double> &wallMf,
                                const physics &phys, const double &tW,
                                const bool &isLower) {
  // initialize wallVars
  wallVars wVars;
  wVars.temperature_ = tW;
  wVars.mf_ = wallMf;

  // get tangential velocity
  const auto vel = state_.Velocity() - velWall;
  const auto velTan = vel - vel.DotProd(area) * area;
  const auto velTanMag = velTan.Mag();

  // get wall properties
  this->CalcRecoveryFactor(phys.Thermodynamic(), tW);
  this->SetWallVars(wVars.temperature_, phys);

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    this->CalcVelocities(yplus, velTanMag);
    // calculate wall heat flux from croco-busemann equation
    this->UpdateGamma(phys.Thermodynamic(), tW);
    wVars.heatFlux_ = this->CalcHeatFlux(phys.EoS());
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
    this->CalcTurbVars(phys, wVars.tke_, wVars.sdr_);
  }

  wVars.density_ = rhoW_;
  wVars.viscosity_ = muW_;
  wVars.turbEddyVisc_ = mutW_;
  wVars.frictionVelocity_ = uStar_;
  wVars.shearStress_ = this->ShearStressMag() * velTan / velTanMag;
  if (!isLower) {
    wVars.shearStress_ *= -1.0;
  }
  return wVars;
}

void wallLaw::UpdateGamma(const unique_ptr<thermodynamic> &thermo,
                          const double &t) {
  // calculate constants
  gamma_ = recoveryFactor_ * uStar_ * uStar_ /
           (2.0 * thermo->Cp(t, state_.MassFractions()) * tW_);
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

double wallLaw::CalcHeatFlux(const unique_ptr<eos> &eqnState) const {
  // calculate wall heat flux from croco-busemann equation
  auto tmp =
      (state_.Temperature(eqnState) / tW_ - 1.0 + gamma_ * uplus_ * uplus_) /
      uplus_;
  return tmp * (rhoW_ * tW_ * kW_ * uStar_) / muW_;
}

double wallLaw::CalcWallTemperature(const physics &phys,
                                    const double &heatFluxW) const {
  const auto t = state_.Temperature(phys.EoS());
  return t + recoveryFactor_ * uStar_ * uStar_ * uplus_ * uplus_ /
                 (2.0 * phys.Thermodynamic()->Cp(t, state_.MassFractions()) +
                  heatFluxW * muW_ / (rhoW_ * kW_ * uStar_));
}

void wallLaw::SetWallVars(const double &tW, const physics &phys) {
  tW_ = tW;
  rhoW_ = phys.EoS()->DensityTP(tW_, state_.P(), state_.MassFractions());

  // get wall viscosity, conductivity from wall temperature
  muW_ = phys.Transport()->EffectiveViscosity(tW_, state_.MassFractions());
  kW_ = phys.Transport()->EffectiveConductivity(tW_, state_.MassFractions());
}

double wallLaw::CalcYplusRoot(const double &yplus) const {
  constexpr auto sixth = 1.0 / 6.0;
  const auto ku = vonKarmen_ * uplus_;
  return yplus - (uplus_ + yplusWhite_ -
                  yplus0_ * (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0)));
}

void wallLaw::EddyVisc(const physics &phys) {
  const auto dYplusWhite =
      2.0 * yplusWhite_ * vonKarmen_ * sqrt(gamma_) / q_ *
      sqrt(std::max(1.0 - pow(2.0 * gamma_ * uplus_ - beta_, 2.0) / (q_ * q_),
                    0.0));
  const auto ku = vonKarmen_ * uplus_;
  mutW_ = muW_ * (1.0 + dYplusWhite -
                  vonKarmen_ * yplus0_ * (1.0 + ku + 0.5 * ku * ku)) -
          phys.Transport()->EffectiveViscosity(state_.Temperature(phys.EoS()),
                                               state_.MassFractions());
  mutW_ = std::max(mutW_, 0.0);
}

void wallLaw::CalcVelocities(const double &yplus, const double &u) {
  // calculate u* and u+ from y+
  uplus_ = (wallDist_ * rhoW_ * u) / (muW_ * yplus);
  uStar_ = u / uplus_;
}

void wallLaw::CalcTurbVars(const physics &phys, double &kWall, double &wWall) {
  this->EddyVisc(phys);
  // calculate turbulence variables from eddy viscosity
  auto wi = 6.0 * muW_ /
            (phys.Turbulence()->WallBeta() * rhoW_ * wallDist_ * wallDist_);
  wi *= phys.Transport()->NondimScaling();
  auto wo =
      uStar_ / (sqrt(phys.Turbulence()->BetaStar()) * vonKarmen_ * wallDist_);
  wo *= phys.Transport()->NondimScaling();
  wWall = sqrt(wi * wi + wo * wo);
  kWall = wWall * mutW_ / state_.Rho() * phys.Transport()->InvNondimScaling();
}

void wallLaw::CalcRecoveryFactor(const unique_ptr<thermodynamic> &thermo,
                                 const double &t) {
  recoveryFactor_ = pow(thermo->Prandtl(t, state_.MassFractions()), 1.0 / 3.0);
}
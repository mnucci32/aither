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
vector3d<double> wallLaw::WallShearStress(
    const primVars &state, const vector3d<double> &area, const idealGas &eos,
    const sutherland &suth, const unique_ptr<turbModel> &turb,
    const double &wallDist, double &kWall, double &wWall) const {
  // get tangential velocity
  auto vel = state.Velocity();
  auto velTan = vel - vel.DotProd(area) * area;
  auto velTanMag = velTan.Mag();

  // get wall temperature from crocco-busemann equation
  auto recoveryFactor = pow(eos.Prandtl(), 1.0 / 3.0);
  auto tW =
      state.Temperature(eos) /
      (1.0 +
       0.5 * (eos.Gamma() - 1.0) * recoveryFactor * velTanMag * velTanMag);
  auto muW = suth.EffectiveViscosity(tW);

  // get wall density from equation of state
  auto rhoW = eos.DensityTP(tW, state.P());

  // declare variables needed after shear stress loop
  auto yplusWhite = 0.0, beta = 0.0, phi = 0.0, q = 0.0, gamma = 0.0;
  auto yplus0 = std::exp(-vonKarmen_ * wallConst_);
  auto uStar = 0.0, uplus = 0.0;

  auto func = [&](const double &yplus) {
    // calculate u* and u+ from y+
    uStar = yplus * muW / (rhoW * wallDist);
    uplus = velTanMag / uStar;

    // calculate constants
    gamma = recoveryFactor * uStar * uStar / (2.0 * eos.SpecificHeat() * tW);
    beta = 0.0;  // due to zero heat flux
    q = sqrt(beta * beta + 4.0 * gamma);
    phi = std::asin(-beta / q);

    // calculate y+ from White & Christoph
    yplusWhite = std::exp((vonKarmen_ / sqrt(gamma)) *
                          (std::asin((2.0 * gamma * uplus - beta) / q) - phi)) *
                 yplus0;

    // calculate y+
    constexpr auto sixth = 1.0 / 6.0;
    auto ku = vonKarmen_ * uplus;
    return yplus - (uplus + yplusWhite -
                    yplus0 * (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0)));
  };
  
  FindRoot(func, 1.0e-5, 1.0e5, 1.0e-8);
  auto shearStressMag = uStar * uStar * rhoW;
  
  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  auto dYplusWhite =
      2.0 * yplusWhite * vonKarmen_ * sqrt(gamma) / q *
      sqrt(std::max(1.0 - pow(2.0 * gamma * uplus - beta, 2.0) / (q * q), 0.0));
  auto ku = vonKarmen_ * uplus;
  auto mutW = muW * (1.0 + dYplusWhite -
                     vonKarmen_ * yplus0 * (1.0 + ku + 0.5 * ku * ku)) -
              suth.EffectiveViscosity(state.Temperature(eos));
  mutW = std::max(mutW, 0.0);

  // calculate turbulence variables from eddy viscosity
  // DEBUG -- check nondimensionalization
  auto wi = 6.0 * muW / (turb->WallBeta() * rhoW * wallDist * wallDist);
  wi *= suth.NondimScaling() * suth.NondimScaling();
  auto wo = uStar / (sqrt(turb->BetaStar()) * vonKarmen_ * wallDist);
  wo *= suth.NondimScaling();
  wWall = sqrt(wi * wi + wo * wo);
  kWall = wWall * mutW / state.Rho();

  return shearStressMag * velTan / velTanMag;
}

vector3d<double> wallLaw::IsothermalWallShearStress(
    const primVars &state, const vector3d<double> &area,
    const tensor<double> &velGradW, const vector3d<double> &tGradW,
    const idealGas &eos, const sutherland &suth,
    const unique_ptr<turbModel> &turb, const double &wallDist,
    const double &tW, double &heatFlux, double &kWall, double &wWall) const {
  // get tangential velocity
  auto vel = state.Velocity();
  auto velTan = vel - vel.DotProd(area) * area;
  auto velTanMag = velTan.Mag();

  // get wall density from equation of state
  auto rhoW = eos.DensityTP(tW, state.P());

  // get wall viscosity from wall temperature
  auto muW = suth.Viscosity(tW);

  // get wall recovery factor
  auto recoveryFactor = pow(eos.Prandtl(), 1.0 / 3.0);

  // iteratively solve for wall shear stress and heat flux
  auto shearStress = TauShear(velGradW, area, muW, 0.0, suth);
  auto shearStressMag = shearStress.Mag();

  auto k = eos.Conductivity(muW);
  heatFlux = k * tGradW.DotProd(area);

  // declare variables that are needed after shear stress loop
  auto yplusWhite = 0.0, beta = 0.0, phi = 0.0, q = 0.0, gamma = 0.0;
  auto yplus0 = std::exp(-vonKarmen_ * wallConst_);
  auto uStar = 0.0, uplus = 0.0;
  
  auto diff = 1.0;
  auto counter = 0;
  while (diff > 1.0e-10 && counter < 100) {
    // assign old shear stress value and heat flux
    auto shearStressMagOld = shearStressMag;

    // calculate friction velocity and u+
    uStar = sqrt(shearStressMag / rhoW);
    uplus = velTanMag / uStar;

    // calculate constants
    gamma = recoveryFactor * uStar * uStar / (2.0 * eos.SpecificHeat() * tW);
    beta = heatFlux * muW / (rhoW * tW * k * uStar);
    q = sqrt(beta * beta + 4.0 * gamma);
    phi = std::asin(-beta / q);
    
    // calculate y+ from White & Christoph
    yplusWhite = std::exp((vonKarmen_ / sqrt(gamma)) *
                              (std::asin((2.0 * gamma * uplus - beta) / q)) -
                          phi) *
                 yplus0;

    // calculate y+
    constexpr auto sixth = 1.0 / 6.0;
    auto ku = vonKarmen_ * uplus;
    auto yplus = uplus + yplusWhite -
                 yplus0 * (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0));

    // calculate new wall shear stress
    shearStressMag = pow(yplus * muW / (rhoW * wallDist), 2.0) * rhoW;

    // calculate new heat flux from crocco-busemann equation
    auto tmp =
        (state.Temperature(eos) / tW - 1.0 + gamma * uStar * uStar) / uplus;
    heatFlux = tmp * (rhoW * tW * k * uStar) / muW;

    diff = fabs(shearStressMag - shearStressMagOld);
    counter++;
  }

  // calculate turbulent eddy viscosity from wall shear stress
  // use compressible form of equation (Nichols & Nelson 2004)
  auto dYplusWhite = 2.0 * yplusWhite * vonKarmen_ * sqrt(gamma) / q *
                     sqrt(1.0 - pow(2.0 * gamma * uplus - beta, 2.0) / (q * q));
  auto ku = vonKarmen_ * uplus;
  auto mutW = muW * (1.0 + dYplusWhite -
                     vonKarmen_ * yplus0 * (1.0 + ku + 0.5 * ku * ku)) -
              suth.Viscosity(state.Temperature(eos));

  // calculate turbulence variables from eddy viscosity
  auto wi = 6.0 * muW / (turb->WallBeta() * rhoW * wallDist * wallDist);
  auto wo = uStar / (sqrt(turb->BetaStar()) * vonKarmen_ * wallDist);
  wWall = sqrt(wi * wi + wo * wo);
  kWall = wWall * mutW / state.Rho();

  return shearStressMag * velTan / velTanMag;
}
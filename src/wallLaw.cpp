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
double wallLaw::WallShearStress(const primVars &state,
                                const vector3d<double> &area,
                                const tensor<double> &velGradW,
                                const idealGas &eos, const sutherland &suth,
                                const double &wallDist,
                                const double &muW,
                                const double &mutW) const {
    // get tangential velocity
    auto vel = state.Velocity();
    auto velTan = vel - vel.DotProd(area) * area;
    auto velTanMag = velTan.Mag();

    // get wall temperature from crocco-busemann equation
    auto recoveryFactor = pow(eos.Prandtl(), 1.0 / 3.0);
    auto tW = state.Temperature(eos) /
        (1.0 + 0.5 * (eos.Gamma() - 1.0) * recoveryFactor * velTanMag *
         velTanMag);

    // get wall density from equation of state
    auto rhoW = eos.DensityTP(tW, state.P());

    // iteratively solve for wall shear stress
    auto shearStress = TauShear(velGradW, area, muW, mutW, suth);
    auto shearStressMag = shearStress.Mag();

    auto diff = 1.0;
    auto counter = 0;
    while (diff > 1.0e-10 && counter < 100) {
      // assign old shear stress value
      auto shearStressMagOld = shearStressMag;

      // calculate friction velocity and u+
      auto uStar = sqrt(shearStressMag / rhoW);
      auto uplus = velTanMag / uStar;

      // calculate constants
      auto gamma = recoveryFactor * uStar * uStar / 
          (2.0 * eos.SpecificHeat() * tW);
      constexpr auto beta = 0.0;  // due to zero heat flux
      auto q = sqrt(4.0 * gamma);
      auto phi = std::asin(0.0);
      auto yplus0 = std::exp(-vonKarmen_ * wallConst_);

      // calculate y+ from White & Christoph
      auto yplusWhite = std::exp((vonKarmen_ / sqrt(gamma)) *
          (std::asin((2.0 * gamma * uplus - beta) / q)) - phi) * yplus0;

      // calculate y+
      constexpr auto sixth = 1.0 / 6.0;
      auto ku = vonKarmen_ * uplus;
      auto yplus = uplus + yplusWhite - yplus0 *
          (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0));

      // calculate new wall shear stress
      shearStressMag = pow(yplus * muW / (rhoW * wallDist), 2.0) * rhoW;

      diff = fabs(shearStressMag - shearStressMagOld);
      counter++;
  }

  return shearStressMag;
}

double wallLaw::IsothermalWallShearStress(const primVars &state,
                                          const vector3d<double> &area,
                                          const tensor<double> &velGradW,
                                          const vector3d<double> &tGradW,
                                          const idealGas &eos, 
                                          const sutherland &suth,
                                          const unique_ptr<turbModel> &turb,
                                          const double &wallDist,
                                          const double &muW, const double &mutW,
                                          const double &tW, 
                                          double &heatFlux) const {
    // get tangential velocity
    auto vel = state.Velocity();
    auto velTan = vel - vel.DotProd(area) * area;
    auto velTanMag = velTan.Mag();

    // get wall density from equation of state
    auto rhoW = eos.DensityTP(tW, state.P());

    // get wall recovery factor
    auto recoveryFactor = pow(eos.Prandtl(), 1.0 / 3.0);

    // iteratively solve for wall shear stress and heat flux
    auto shearStress = TauShear(velGradW, area, muW, mutW, suth);
    auto shearStressMag = shearStress.Mag();
    
    auto k = eos.Conductivity(muW);
    auto kt = eos.TurbConductivity(mutW, turb->TurbPrandtlNumber());
    heatFlux = (k + kt) * tGradW.DotProd(area);

    auto diff = 1.0;
    auto counter = 0;
    while (diff > 1.0e-10 && counter < 100) {
      // assign old shear stress value and heat flux
      auto shearStressMagOld = shearStressMag;
      auto heatFluxOld = heatFlux;

      // calculate friction velocity and u+
      auto uStar = sqrt(shearStressMag / rhoW);
      auto uplus = velTanMag / uStar;

      // calculate constants
      auto gamma = recoveryFactor * uStar * uStar / 
          (2.0 * eos.SpecificHeat() * tW);
      auto beta = heatFlux * muW / (rhoW * tW * k * uStar);
      auto q = sqrt(beta * beta + 4.0 * gamma);
      auto phi = std::asin(-beta / q);
      auto yplus0 = std::exp(-vonKarmen_ * wallConst_);

      // calculate y+ from White & Christoph
      auto yplusWhite = std::exp((vonKarmen_ / sqrt(gamma)) *
          (std::asin((2.0 * gamma * uplus - beta) / q)) - phi) * yplus0;

      // calculate y+
      constexpr auto sixth = 1.0 / 6.0;
      auto ku = vonKarmen_ * uplus;
      auto yplus = uplus + yplusWhite - yplus0 *
          (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0));

      // calculate new wall shear stress
      shearStressMag = pow(yplus * muW / (rhoW * wallDist), 2.0) * rhoW;

      // calculate new heat flux from crocco-busemann equation
      auto tmp = (state.Temperature(eos) / tW - 1.0 + gamma * uStar * uStar) /
          uplus;
      heatFlux = tmp * (rhoW * tW * k * uStar) / muW;

      diff = fabs(shearStressMag - shearStressMagOld);
      diff += fabs(heatFlux - heatFluxOld);
      counter++;
  }

  return shearStressMag;
}
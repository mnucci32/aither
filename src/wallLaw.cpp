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
#include "wallLaw.hpp"
#include "primVars.hpp"  // primVars
#include "eos.hpp"       // idealGas
#include "utility.hpp"   // TauShear

using std::cout;
using std::endl;
using std::cerr;

// -------------------------------------------------------------------------
double wallLaw::Gamma(const idealGas &eos, const double &uStar, 
                      const double &tW) const {
    auto recoveryFactor = pow(eos.Prandtl(), 1.0 / 3.0);
    return recoveryFactor * uStar * uStar / (2.0 * eos.SpecificHeat() * tW);
}

double wallLaw::WallShearStress(const primVars &state, 
                                const vector3d<double> &area,
                                const tensor<double> &velGradW, 
                                const idealGas &eos, 
                                const double &wallDist,
                                const double &muW,
                                const double &mutW,
                                const double &tW,
                                const double &rhoW,
                                const sutherland &suth) const {
    auto shearStress = TauShear(velGradW, area, muW, mutW, suth);
    auto shearStressMag = shearStress.Mag();
    auto uStar = sqrt(shearStressMag / rhoW);

    auto gamma = this->Gamma(eos, uStar, tW);
    auto beta = this->Beta();
    auto q = std::sqrt(beta * beta + 4.0 * gamma);
    auto phi = std::asin(-beta / q);
    auto yplus0 = std::exp(-vonKarmen_ * wallConst_);

    auto vel = state.Velocity();
    auto velTan = vel - vel.DotProd(area) * area;
    auto velTanMag = velTan.Mag();

    auto uplus = velTanMag / uStar;

    auto yplusWhite = std::exp((vonKarmen_ / std::sqrt(gamma)) *
        (std::asin((2.0 * gamma * uplus - beta) / q)) - phi) * yplus0;

    constexpr auto sixth = 1.0 / 6.0;
    auto ku = vonKarmen_ * uplus;
    auto yplus = uplus + yplusWhite - yplus0 *
        (1.0 + ku + 0.5 * ku * ku + sixth * pow(ku, 3.0));

    // calculate new wall shear stress
    shearStressMag = pow(yplus * muW / (rhoW * wallDist), 2.0) * rhoW;

    return shearStressMag;
}
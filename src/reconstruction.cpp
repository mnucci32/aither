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

#include <string>
#include "reconstruction.hpp"
#include "primitive.hpp"
#include "limiter.hpp"

using std::string;

// member function to calculate reconstruction of primitive variables from cell
// center to cell face this function uses MUSCL extrapolation resulting in
// higher order accuracy
/*

____________________|_____________________
|         |       Ul|Ur        |         |
|         |         |          |         |
|   Ui-1  |   Ui    |   Ui+1   |   Ui+2  |
|         |         |          |         |
|         |         |          |         |
|_________|_________|__________|_________|
                    |
     face to reconstruct state at
<--UW2---><--UW1---><---DW---->

The diagram above shows the stencil used to reconstruct the left and right
states for a given face. For each reconstruction (right and left) only 3 points
are needed, two upwind points and one downwind point. For the left
reconstruction, Ui is the first upwind point, Ui-1 is the second upwind point,
and Ui+1 is the downwind point. For the right reconstruction Ui+1 is the first
upwind point, Ui+2 is the second upwind point, and Ui is the downwind point.

For left reconstruction the MUSCL scheme goes as follows:
Ui+1/2 = Ui + 0.25 * ( (1-K) * (Ui - Ui-1) + (1+K) * (Ui+1 - Ui) )

The above equation assumes there is no limiter, and that the grid spacing is
uniform. In the above equation K refers to the parameter kappa which can be
varied to produce different reconstructions. Acceptable values of K are [-1,1]
(inclusive)

K = -1 -- fully upwind reconstruction - linear - results in 2nd order accuracy
K = 0 -- fromm scheme - linear - results in 2nd order accuracy
K = 0.5 -- QUICK scheme - parabolic - results in 2nd order accuracy
K = 1 -- central scheme - linear - results in 2nd order accuracy but is unstable
without dissipation added
K = 1/3 -- third order - parabolic - results 2nd order accuracy with lowest
error
(all order of accuracy estimates assume constant fluxes on cell faces which
limits order of accuracy to 2)

With limiters the equation looks as follows:
Ui+1/2 = Ui + 0.25 * (Ui - Ui-1) * ( (1-K) * L  + (1+K) * R * Linv )

L represents the limiter function which ranges between 0 and 1. A value of 1
implies there is no limiter and the full accuracy of the scheme is achieved. A
value of 0 implies that the solution has been limited to first order accuracy.
An in between value results in an order of accuracy between first and second.
Linv is the inverse of the limiter.

R represents the divided difference that the limiter is a function of.
R = (Ui - Ui-1) / (Ui+1 - Ui)

The MUSCL scheme can be extended to nonuniform grids by adding in terms
representing the difference in size between the cells. In the above diagram the
values UW2, UW1, and DW represent the length of the second upwind, first upwind,
and downwind cells respectively. dP and dM represent the factors due to the
change in cell size between the first upwind to downwind and first upwind to
second upwind cells.

dP = (UW1 + DW) / (2.0 * UW)
dM = (UW + UW2) / (2.0 * UW)
R = ((Ui - Ui-1) / dP) / ((Ui+1 - Ui) / dM)

Ui+1/2 = Ui + 0.25 * ((Ui - Ui-1) / dM) * ( (1-K) * L  + (1+K) * R * Linv )

*/
template <typename T>
primitive FaceReconMUSCL(const T &upwind2, const T &upwind1, const T &downwind1,
                         const double &kappa, const string &lim,
                         const double &uw, const double &uw2,
                         const double &dw) {
  // upwind2 -- upwind cell furthest from the face at which the primitive is
  //            being reconstructed.
  // upwind1 -- upwind cell nearest to the face at which the primitive is
  //            being reconstructed.
  // downwind1 -- downwind cell.
  // kappa -- parameter that determines which scheme is implemented
  // uw -- length of upwind cell
  // uw2 -- length of furthest upwind cell
  // dw -- length of downwind cell

  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "FaceReconMUSCL requires primitive or primativeView type");

  const auto dPlus = (uw + uw) / (uw + dw);
  const auto dMinus = (uw + uw) / (uw + uw2);

  // divided differences to base limiter on; eps must be listed to left of
  // primitive
  const auto r = (EPS + (downwind1 - upwind1) * dPlus) /
      (EPS + (upwind1 - upwind2) * dMinus);

  primitive limiter(upwind2.Size(), upwind2.NumSpecies());
  primitive invLimiter(upwind2.Size(), upwind2.NumSpecies());
  if (lim == "none") {
    limiter = LimiterNone(limiter.Size(), limiter.NumSpecies());
    invLimiter = limiter;
  } else if (lim == "vanAlbada") {
    limiter = LimiterVanAlbada(r);
    invLimiter = LimiterVanAlbada(1.0 / r);
  } else if (lim == "minmod") {
    limiter = LimiterMinmod(upwind1 - upwind2, downwind1 - upwind1, kappa);
    invLimiter = limiter / r;
  } else {
    cerr << "ERROR: Limiter " << lim << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  // calculate reconstructed state at face using MUSCL method with limiter
  return upwind1 + 0.25 * ((upwind1 - upwind2) * dMinus) *
    ((1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter);
}

// member function for higher order reconstruction via weno
template <typename T>
primitive FaceReconWENO(const T &upwind3, const T &upwind2, const T &upwind1,
                        const T &downwind1, const T &downwind2,
                        const double &uw3, const double &uw2, const double &uw1,
                        const double &dw1, const double &dw2,
                        const bool &isWenoZ) {
  // make sure template type is correct
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "FaceReconWENO requires primitive or primativeView type");

  // get candidate smaller stencils
  const vector<double> cellWidth = {uw3, uw2, uw1, dw1, dw2};

  constexpr auto degree = 2;
  constexpr auto up1Loc = 2;
  const auto coeffs0 = LagrangeCoeff(cellWidth, degree, 2, up1Loc);
  const auto stencil0 = coeffs0[0] * upwind3 + coeffs0[1] * upwind2 +
      coeffs0[2] * upwind1;

  const auto coeffs1 = LagrangeCoeff(cellWidth, degree, 1, up1Loc);
  const auto stencil1 = coeffs1[0] * upwind2 + coeffs1[1] * upwind1 +
      coeffs1[2] * downwind1;

  const auto coeffs2 = LagrangeCoeff(cellWidth, degree, 0, up1Loc);
  const auto stencil2 = coeffs2[0] * upwind1 + coeffs2[1] * downwind1 +
      coeffs2[2] * downwind2;

  // get coefficients for large stencil
  const auto fullCoeffs = LagrangeCoeff(cellWidth, 4, 2, up1Loc);

  // linear weights
  const auto lw0 = fullCoeffs[0] / coeffs0[0];
  const auto lw1 = fullCoeffs[4] / coeffs2[2];
  const auto lw2 = 1.0 - lw0 - lw1;

  const auto beta0 = Beta0(uw3, uw2, uw1, upwind3, upwind2, upwind1);
  const auto beta1 = Beta1(uw2, uw1, dw1, upwind2, upwind1, downwind1);
  const auto beta2 = Beta2(uw1, dw1, dw2, upwind1, downwind1, downwind2);

  // calculate nonlinear weights
  primitive nlw0(upwind3.Size(), upwind3.NumSpecies());
  primitive nlw1(upwind3.Size(), upwind3.NumSpecies());
  primitive nlw2(upwind3.Size(), upwind3.NumSpecies());
  if (isWenoZ) {
    // using weno-z weights with q = 2
    const auto tau5 = (beta0 - beta2).Abs();
    constexpr auto eps = 1.0e-40;
    nlw0 = lw0 * (1.0 + (tau5 / (eps + beta0)).Squared());
    nlw1 = lw1 * (1.0 + (tau5 / (eps + beta1)).Squared());
    nlw2 = lw2 * (1.0 + (tau5 / (eps + beta2)).Squared());
  } else {  // standard WENO
    // calculate nonlinear weights
    constexpr auto eps = 1.0e-6;
    nlw0 = lw0 / (eps + beta0).Squared();
    nlw1 = lw1 / (eps + beta1).Squared();
    nlw2 = lw2 / (eps + beta2).Squared();
  }

  // normalize weights
  const auto sum_nlw = nlw0 + nlw1 + nlw2;
  nlw0 /= sum_nlw;
  nlw1 /= sum_nlw;
  nlw2 /= sum_nlw;

  // return weighted contribution of each stencil
  return nlw0 * stencil0 + nlw1 * stencil1 + nlw2 * stencil2;
}

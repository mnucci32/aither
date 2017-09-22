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

#include <memory>
#include <cmath>
#include <algorithm>  // max
#include "spectralRadius.hpp"
#include "primitive.hpp"
#include "eos.hpp"
#include "transport.hpp"
#include "thermodynamic.hpp"

using std::max;
using std::min;
using std::unique_ptr;

/*Function to return the inviscid spectral radius for one direction (i, j, or k)
given a cell state, equation of state, and 2 face area vectors

L = 0.5 * (A1 + A2) * (|Vn| + SoS)

In the above equation L is the spectral radius in either the i, j, or k
direction. A1 and A2 are the two face areas in that direction. Vn is the
cell velocity normal to that direction. SoS is the speed of sound at the cell
 */
double InvCellSpectralRadius(const primitive &state,
                             const unitVec3dMag<double> &fAreaL,
                             const unitVec3dMag<double> &fAreaR,
                             const unique_ptr<thermodynamic> &thermo,
                             const unique_ptr<eos> &eqnState) {
  // state -- primitive state variables
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // thermo -- thermodynamic model
  // eqnState -- equation of state

  // normalize face areas
  const auto normAvg = (0.5 * (fAreaL.UnitVector() +
                               fAreaR.UnitVector())).Normalize();
  // average area magnitude
  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());

  // return spectral radius
  return (fabs(state.Velocity().DotProd(normAvg)) +
          state.SoS(thermo, eqnState)) *
         fMag;
}

double InvFaceSpectralRadius(const primitive &state,
                             const unitVec3dMag<double> &fArea,
                             const unique_ptr<thermodynamic> &thermo,
                             const unique_ptr<eos> &eqnState) {
  // state -- primitive state variables
  // fArea -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state

  // return spectral radius
  return 0.5 * fArea.Mag() *
         (fabs(state.Velocity().DotProd(fArea.UnitVector())) +
          state.SoS(thermo, eqnState));
}

/*Function to calculate the viscous spectral radius for one direction (i, j, or
k).

L = max(4/(3*rho), g/rho) * mu/Pr * A^2 / V

In the above equation L is the viscous spectral radius for a given direction (i,
j, or k). Rho is the density at the cell center. G is gamma, mu is viscosity,
and Pr is the Prandtl number (all at the cell center). A is the average face area
of the given direction (i, j, k), and V is the cell volume. This implementation
comes from Blazek.
 */
double ViscCellSpectralRadius(
    const primitive &state, const unitVec3dMag<double> &fAreaL,
    const unitVec3dMag<double> &fAreaR, const unique_ptr<thermodynamic> &thermo,
    const unique_ptr<eos> &eqnState, const unique_ptr<transport> &trans,
    const double &vol, const double &mu, const double &mut,
    const unique_ptr<turbModel> &turb) {
  // state -- primitive state variables
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // vol -- cell volume
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model

  // average area magnitude
  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  const auto t = state.Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * state.Rho()), thermo->Gamma(t) / state.Rho());
  // viscous term
  const auto viscTerm = trans->NondimScaling() *
      (mu / thermo->Prandtl(t) +  mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return maxTerm * viscTerm * fMag * fMag / vol;
}

double ViscFaceSpectralRadius(
    const primitive &state, const unitVec3dMag<double> &fArea,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<eos> &eqnState,
    const unique_ptr<transport> &trans, const double &dist, const double &mu,
    const double &mut, const unique_ptr<turbModel> &turb) {
  // state -- primitive state variables
  // fArea -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // dist -- distacne from cell center to cell center
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model

  const auto t = state.Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * state.Rho()), thermo->Gamma(t) / state.Rho());
  // viscous term
  const auto viscTerm = trans->NondimScaling() *
      (mu / thermo->Prandtl(t) +  mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return fArea.Mag() / dist * maxTerm * viscTerm;
}

double CellSpectralRadius(
    const primitive &state, const unitVec3dMag<double> &fAreaL,
    const unitVec3dMag<double> &fAreaR, const unique_ptr<thermodynamic> &thermo,
    const unique_ptr<eos> &eqnState, const unique_ptr<transport> &trans,
    const double &vol, const double &mu, const double &mut,
    const unique_ptr<turbModel> &turb, const bool &isViscous) {
  // state -- primitive state variables
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // vol -- cell volume
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model
  // isViscous -- flag that is true if simulation is viscous

  auto specRad = InvCellSpectralRadius(state, fAreaL, fAreaR, thermo, eqnState);

  if (isViscous) {
    // factor 2 because viscous spectral radius is not halved (Blazek 6.53)
    specRad +=
        2.0 * ViscCellSpectralRadius(state, fAreaL, fAreaR, thermo, eqnState,
                                     trans, vol, mu, mut, turb);
  }
  return specRad;
}

double FaceSpectralRadius(const primitive &state,
                          const unitVec3dMag<double> &fArea,
                          const unique_ptr<thermodynamic> &thermo,
                          const unique_ptr<eos> &eqnState,
                          const unique_ptr<transport> &trans,
                          const double &dist, const double &mu,
                          const double &mut, const unique_ptr<turbModel> &turb,
                          const bool &isViscous) {
  // state -- primitive state variables
  // fAreaL -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model
  // isViscous -- flag that is true if simulation is viscous

  auto specRad = InvFaceSpectralRadius(state, fArea, thermo, eqnState);

  if (isViscous) {
    specRad += ViscFaceSpectralRadius(state, fArea, thermo, eqnState, trans,
                                      dist, mu, mut, turb);
  }
  return specRad;
}


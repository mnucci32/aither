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

#ifndef SPECTRALRADIUSHEADERDEF
#define SPECTRALRADIUSHEADERDEF

/* This header contains the functions to calculate the inviscid and viscous
 * spectral radii
 */

#include <string>                  // string
#include <memory>                  // unique_ptr
#include <algorithm>               // max
#include "vector3d.hpp"            // vector3d
#include "primitive.hpp"
#include "arrayView.hpp"
#include "eos.hpp"
#include "transport.hpp"
#include "thermodynamic.hpp"
#include "turbulence.hpp"

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
template <typename T>
double InvCellSpectralRadius(const T &state,
                             const unitVec3dMag<double> &fAreaL,
                             const unitVec3dMag<double> &fAreaR,
                             const unique_ptr<thermodynamic> &thermo,
                             const unique_ptr<eos> &eqnState) {
  // state -- primitive state variables
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

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

template <typename T>
double InvFaceSpectralRadius(const T &state,
                             const unitVec3dMag<double> &fArea,
                             const unique_ptr<thermodynamic> &thermo,
                             const unique_ptr<eos> &eqnState) {
  // state -- primitive state variables
  // fArea -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

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
template <typename T>
double ViscCellSpectralRadius(
    const T &state, const unitVec3dMag<double> &fAreaL,
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
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  // average area magnitude
  const auto fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  const auto t = state.Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * state.Rho()),
          thermo->Gamma(t, state.MassFractions()) / state.Rho());
  // viscous term
  const auto viscTerm =
      trans->NondimScaling() * (mu / thermo->Prandtl(t, state.MassFractions()) +
                                mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return maxTerm * viscTerm * fMag * fMag / vol;
}

template <typename T>
double ViscFaceSpectralRadius(
    const T &state, const unitVec3dMag<double> &fArea,
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
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  const auto t = state.Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * state.Rho()),
          thermo->Gamma(t, state.MassFractions()) / state.Rho());
  // viscous term
  const auto viscTerm =
      trans->NondimScaling() * (mu / thermo->Prandtl(t, state.MassFractions()) +
                                mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return fArea.Mag() / dist * maxTerm * viscTerm;
}

template <typename T>
double CellSpectralRadius(
    const T &state, const unitVec3dMag<double> &fAreaL,
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
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  auto specRad = InvCellSpectralRadius(state, fAreaL, fAreaR, thermo, eqnState);

  if (isViscous) {
    // factor 2 because viscous spectral radius is not halved (Blazek 6.53)
    specRad +=
        2.0 * ViscCellSpectralRadius(state, fAreaL, fAreaR, thermo, eqnState,
                                     trans, vol, mu, mut, turb);
  }
  return specRad;
}

template <typename T>
double FaceSpectralRadius(const T &state,
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
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  auto specRad = InvFaceSpectralRadius(state, fArea, thermo, eqnState);

  if (isViscous) {
    specRad += ViscFaceSpectralRadius(state, fArea, thermo, eqnState, trans,
                                      dist, mu, mut, turb);
  }
  return specRad;
}

#endif

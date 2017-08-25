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

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <cmath>
#include <algorithm>  // max
#include "primVars.hpp"
#include "input.hpp"               // input
#include "turbulence.hpp"          // turbModel
#include "utility.hpp"
#include "wallLaw.hpp"
#include "wallData.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::min;
using std::unique_ptr;

primVars::primVars(const genArray &a, const bool &prim,
                   const unique_ptr<eos> &eqnState,
                   const unique_ptr<thermodynamic> &thermo,
                   const unique_ptr<turbModel> &turb) {
  // a -- array of conservative or primative variables
  // prim -- flag that is true if variable a is primative variables
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // turb -- turbulence model

  if (prim) {  // genArray is primative variables
    for (auto ii = 0; ii < NUMVARS; ii++) {
      data_[ii] = a[ii];
    }
  } else {  // genArray is conserved variables
    data_[0] = a[0];
    data_[1] = a[1] / a[0];
    data_[2] = a[2] / a[0];
    data_[3] = a[3] / a[0];
    const auto energy = a[4] / a[0];
    data_[4] = eqnState->PressFromEnergy(thermo, data_[0], energy,
                                         this->Velocity().Mag());
    data_[5] = a[5] / a[0];
    data_[6] = a[6] / a[0];
  }

  // Adjust turbulence variables to be above minimum if necessary
  this->LimitTurb(turb);
}

// member function to initialize a state with nondimensional values
void primVars::NondimensionalInitialize(const unique_ptr<eos> &eqnState,
                                        const input &inp,
                                        const unique_ptr<transport> &trans,
                                        const int &parBlock,
                                        const unique_ptr<turbModel> &turb) {
  // get initial condition state for parent block
  auto ic = inp.ICStateForBlock(parBlock);

  data_[0] = ic.Density();
  data_[1] = ic.Velocity().X();
  data_[2] = ic.Velocity().Y();
  data_[3] = ic.Velocity().Z();
  data_[4] = ic.Pressure();
  
  if (inp.IsRANS()) {
    // Initialize turbulence quantities based off of specified turublence
    // intensity and eddy viscosity ratio. This is the default for
    // STAR-CCM+
    this->ApplyFarfieldTurbBC(this->Velocity(), ic.TurbulenceIntensity(),
                              ic.EddyViscosityRatio(), trans, eqnState, turb);
  } else {
    data_[5] = 0.0;
    data_[6] = 0.0;
  }
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const primVars &prim) {
  os << prim.Rho() << endl;
  os << prim.U() << endl;
  os << prim.V() << endl;
  os << prim.W() << endl;
  os << prim.P() << endl;
  os << prim.Tke() << endl;
  os << prim.Omega() << endl;
  return os;
}

// member function to calculate reconstruction of primative variables from cell
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
primVars primVars::FaceReconMUSCL(const primVars &primUW2,
                                  const primVars &primDW1, const double &kappa,
                                  const string &lim, const double &uw,
                                  const double &uw2, const double &dw) const {
  // primUW2 -- upwind cell furthest from the face at which the primative is
  //            being reconstructed.
  // primUW1 -- upwind cell nearest to the face at which the primative is
  //            being reconstructed.
  // primDW1 -- downwind cell.
  // kappa -- parameter that determines which scheme is implemented
  // uw -- length of upwind cell
  // uw2 -- length of furthest upwind cell
  // dw -- length of downwind cell

  const auto primUW1 = *this;

  const auto dPlus = (uw + uw) / (uw + dw);
  const auto dMinus = (uw + uw) / (uw + uw2);

  // divided differences to base limiter on; eps must be listed to left of
  // primVars
  const auto r = (EPS + (primDW1 - primUW1) * dPlus) /
      (EPS + (primUW1 - primUW2) * dMinus);

  primVars limiter;
  primVars invLimiter;
  if (lim == "none") {
    limiter = LimiterNone();
    invLimiter = limiter;
  } else if (lim == "vanAlbada") {
    limiter = LimiterVanAlbada(r);
    invLimiter = LimiterVanAlbada(1.0 / r);
  } else if (lim == "minmod") {
    limiter = LimiterMinmod(primUW1 - primUW2, primDW1 - primUW1, kappa);
    invLimiter = limiter / r;
  } else {
    cerr << "ERROR: Limiter " << lim << " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }

  // calculate reconstructed state at face using MUSCL method with limiter
  return primUW1 + 0.25 * ((primUW1 - primUW2) * dMinus) *
    ((1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter);
}

// member function for higher order reconstruction via weno
primVars primVars::FaceReconWENO(const primVars &upwind2,
                                 const primVars &upwind3,
                                 const primVars &downwind1,
                                 const primVars &downwind2, const double &uw1,
                                 const double &uw2, const double &uw3,
                                 const double &dw1, const double &dw2,
                                 const bool &isWenoZ) const {
  // get candidate smaller stencils
  const vector<double> cellWidth = {uw3, uw2, uw1, dw1, dw2};

  constexpr auto degree = 2;
  constexpr auto up1Loc = 2;
  const auto coeffs0 = LagrangeCoeff(cellWidth, degree, 2, up1Loc);
  const auto stencil0 = coeffs0[0] * upwind3 + coeffs0[1] * upwind2 +
      coeffs0[2] * (*this);

  const auto coeffs1 = LagrangeCoeff(cellWidth, degree, 1, up1Loc);
  const auto stencil1 = coeffs1[0] * upwind2 + coeffs1[1] * (*this) +
      coeffs1[2] * downwind1;

  const auto coeffs2 = LagrangeCoeff(cellWidth, degree, 0, up1Loc);
  const auto stencil2 = coeffs2[0] * (*this) + coeffs2[1] * downwind1 +
      coeffs2[2] * downwind2;

  // get coefficients for large stencil
  const auto fullCoeffs = LagrangeCoeff(cellWidth, 4, 2, up1Loc);

  // linear weights
  const auto lw0 = fullCoeffs[0] / coeffs0[0];
  const auto lw1 = fullCoeffs[4] / coeffs2[2];
  const auto lw2 = 1.0 - lw0 - lw1;

  const auto beta0 = Beta0(uw3, uw2, uw1, upwind3, upwind2, (*this));
  const auto beta1 = Beta1(uw2, uw1, dw1, upwind2, (*this), downwind1);
  const auto beta2 = Beta2(uw1, dw1, dw2, (*this), downwind1, downwind2);

  // calculate nonlinear weights
  primVars nlw0, nlw1, nlw2;
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


// member function to calculate minmod limiter
primVars primVars::LimiterMinmod(const primVars &upwind,
                                 const primVars &downwind,
                                 const double &kap) const {
  // upwind -- upwind state (primative)
  // downwind -- downwind state (primative)
  // kap -- MUSCL parameter kappa

  primVars limiter;

  // calculate minmod parameter beta
  const auto beta = (3.0 - kap) / (1.0 - kap);

  // calculate minmod limiter
  for (auto ii = 0; ii < NUMVARS; ii++) {
    auto sign = 0.0;
    if (upwind.data_[ii] > 0.0) {
      sign = 1.0;
    } else if (upwind.data_[ii] < 0.0) {
      sign = -1.0;
    }
    limiter.data_[ii] = sign *
        max(0.0, min(fabs(upwind.data_[ii]),
                     sign * downwind.data_[ii] * beta));
  }

  return limiter;
}

// member function to calculate Van Albada limiter
primVars primVars::LimiterVanAlbada(const primVars &r) const {
  // r -- ratio of divided differences

  auto limiter = (r + r * r) / (1 + r * r);
  // if value is negative, return zero
  for (auto ii = 0; ii < NUMVARS; ii++) {
    limiter.data_[ii] = max(0.0, limiter.data_[ii]);
  }
  return limiter;
}

// member function to return no limiter
primVars primVars::LimiterNone() const {
  // for no limiter return all 1s
  primVars limiter(1.0);
  return limiter;
}

// member function to return the state of the appropriate ghost cell
/*

In the diagram below Ui represents the interior state. This should be represnted
by (*this). Ug1 and Ug2 are the first and second layers of ghost cells
respectively. They are the output of the function. The input "layer" determines
whether Ug1 or Ug2 is output for BCs that use extrapolation. Reflected boundaries
(slipWall, viscousWall) do not use the layer input. Instead Ui+1 must be
specified as (*this) to output Ug2.

____________________|___________
|         |         |          |
|         |         |          |
|   Ug2   |   Ug1   |   Ui     |
|         |         |          |
|         |         |          |
|_________|_________|__________|
                    |
             boundary face

Currently the following boundary conditions are supported: slipWall,
viscousWall, characteristic, stagnationInlet, pressureOutlet, subsonicInflow,
subsonicOutflow, supersonicInflow, supersonicOutflow
*/
primVars primVars::GetGhostState(
    const string &bcType, const vector3d<double> &areaVec,
    const double &wallDist, const int &surf, const input &inputVars,
    const int &tag, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, wallVars &wVars, const int &layer,
    const double &dt, const primVars &stateN, const vector3d<double> &pressGrad,
    const tensor<double> &velGrad, const double &avgMach,
    const double &maxMach) const {
  // bcType -- type of boundary condition to supply ghost cell for
  // areaVec -- unit area vector of boundary face
  // surf -- surface type [1-6]
  // wallDist -- distance from cell center to nearest wall boundary
  // inputVar -- all input variables
  // tag -- boundary condition tag
  // eqnState -- equation of state
  // trans -- unique_ptr<transport> model for viscosity
  // turb -- turbulence model
  // layer -- layer of ghost cell to return (1st (closest) or 2nd (farthest))
  // dt -- cell time step nearest to wall boundary
  // stateN -- solution at boundary adjacent cell at time n
  // pressGrad -- pressure gradient in adjcent cell
  // velGrad -- velocity gradient in adjacent cell
  // avgMach -- average mach number on surface patch
  // maxMach -- maximum mach number on surface patch

  // the instance of primVars being acted upon should be the interior cell
  // bordering the boundary

  // set ghost state equal to boundary state to start
  auto ghostState = (*this);

  // face area vector (should always point out of domain)
  // at lower surface normal should point out of domain for ghost cell calc
  const auto isLower = surf % 2 == 1;
  const auto normArea = isLower? -1.0 * areaVec : areaVec;

  // slip wall boundary condition
  // ----------------------------------------------------------------------
  // slipWall is implemented as a reflection of the interior states so that
  // there is no mass flow across the boundary.
  if (bcType == "slipWall") {  // for slip wall state should be reflected across
                               // boundary face, density and pressure stay equal
                               // to the boundary cell
    const auto stateVel = this->Velocity();
    const auto normVelCellCenter = stateVel.DotProd(normArea);

    // for a slip wall the velocity of the boundary cell center is reflected
    // across the boundary face to get the velocity at the ghost cell center
    const auto ghostVel = stateVel - 2.0 * normArea * normVelCellCenter;

    ghostState.data_[1] = ghostVel.X();
    ghostState.data_[2] = ghostVel.Y();
    ghostState.data_[3] = ghostVel.Z();

    // numerical BCs for rho and pressure, same as boundary state
    // numerical BCs for turbulence variables

  // viscous wall boundary condition
  // -------------------------------------------------------------------------
  // viscous wall uses the interior density and pressure, but flips the sign on
  // the velocity so that the velocity at the boundary is 0.
  } else if (bcType == "viscousWall") {  // for viscous wall velocity at face
                                       // should be 0.0, density and pressure
                                       // stay equal to the boundary cell
    const auto & bcData = inputVars.BCData(tag);

    // ghost cell velocity at cell center is set to opposite of velocity at
    // boundary cell center so that velocity at face will be zero
    // only true for low-Re wall treatment
    const auto velWall = bcData->Velocity();
    const auto ghostVel = 2.0 * velWall - this->Velocity();
    ghostState.data_[1] = ghostVel.X();
    ghostState.data_[2] = ghostVel.Y();
    ghostState.data_[3] = ghostVel.Z();

    if (bcData->IsIsothermal()) {  //-----------------------------------------
      const auto tWall = bcData->Temperature();
      // for wall law ghost velocity and turbulence variables calculated
      // simultaneously
      if (bcData->IsWallLaw()) {
        wallLaw wl(bcData->VonKarmen(), bcData->WallConstant(), *this, wallDist,
                   inputVars.IsRANS());
        wVars = wl.IsothermalBCs(normArea, velWall, eqnState, thermo, trans,
                                 turb, tWall, isLower);

        if (wVars.SwitchToLowRe()) {
          const auto tGhost = 2.0 * tWall - this->Temperature(eqnState);
          ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());
        } else {
          // use wall law heat flux to get ghost cell density
          // need turbulent contribution because eddy viscosity is not 0 at wall
          const auto kappa = trans->Conductivity(wVars.viscosity_,
                                                 wVars.temperature_, thermo) +
                             trans->TurbConductivity(
                                 wVars.turbEddyVisc_, turb->TurbPrandtlNumber(),
                                 wVars.temperature_, thermo);
          // 2x wall distance as gradient length
          const auto tGhost = tWall - wVars.heatFlux_ / kappa * 2.0 * wallDist;
          ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());
        }

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState.data_[5] = 2.0 * wVars.tke_ - this->Tke();
          ghostState.data_[6] = 2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState.data_[5] = layer * ghostState.data_[5] - wVars.tke_;
            ghostState.data_[6] = layer * ghostState.data_[6] - wVars.sdr_;
          }
        }
      } else {  // low-Re wall treatment
        const auto tGhost = 2.0 * tWall - this->Temperature(eqnState);
        ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());
      }
    } else if (bcData->IsConstantHeatFlux()) {  //-----------------------------
      // must nondimensionalize heat flux
      const auto qWall = bcData->HeatFlux();
      if (bcData->IsWallLaw()) {
        wallLaw wl(bcData->VonKarmen(), bcData->WallConstant(), *this, wallDist,
                   inputVars.IsRANS());
        wVars = wl.HeatFluxBCs(normArea, velWall, eqnState, thermo, trans, turb,
                               qWall, isLower);

        if (wVars.SwitchToLowRe()) {
          // don't need turbulent contribution b/c eddy viscosity is 0 at wall
          const auto t = this->Temperature(eqnState);
          const auto mu = trans->EffectiveViscosity(t);
          const auto kappa = trans->Conductivity(mu, t, thermo);
          // 2x wall distance as gradient length
          const auto tGhost =
              this->Temperature(eqnState) - qWall / kappa * 2.0 * wallDist;
          ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());

        } else {
          // use wall law wall temperature to get ghost cell density
          const auto tGhost =
              2.0 * wVars.temperature_ - this->Temperature(eqnState);
          ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());
        }

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState.data_[5] = 2.0 * wVars.tke_ - this->Tke();
          ghostState.data_[6] = 2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState.data_[5] = layer * ghostState.data_[5] - wVars.tke_;
            ghostState.data_[6] = layer * ghostState.data_[6] - wVars.sdr_;
          }
        }
      } else {  // low-Re wall treatment
        // don't need turbulent contribution b/c eddy viscosity is 0 at wall
        const auto t = this->Temperature(eqnState);
        const auto mu = trans->EffectiveViscosity(t);
        const auto kappa = trans->Conductivity(mu, t, thermo);
        // 2x wall distance as gradient length
        const auto tGhost =
            this->Temperature(eqnState) - qWall / kappa * 2.0 * wallDist;
        ghostState.data_[0] = eqnState->DensityTP(tGhost, ghostState.P());
        // numerical BCs for pressure, same as boundary state
      }
    } else {  // default is adiabatic -----------------------------------------
      if (bcData->IsWallLaw()) {
        wallLaw wl(bcData->VonKarmen(), bcData->WallConstant(), *this, wallDist,
                   inputVars.IsRANS());
        wVars = wl.AdiabaticBCs(normArea, velWall, eqnState, thermo, trans,
                                turb, isLower);

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState.data_[5] = 2.0 * wVars.tke_ - this->Tke();
          ghostState.data_[6] = 2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState.data_[5] = layer * ghostState.data_[5] - wVars.tke_;
            ghostState.data_[6] = layer * ghostState.data_[6] - wVars.sdr_;
          }
        }
      }
      // numerical BCs for pressure, density - same as boundary state
    }

    // turbulence bcs
    // for wall law, turbulence bcs are already calculated, unless low Re model
    // should be used
    if (inputVars.IsRANS() && (!bcData->IsWallLaw() || wVars.SwitchToLowRe())) {
      // tke at cell center is set to opposite of tke at boundary cell center
      // so that tke at face will be zero
      ghostState.data_[5] = -1.0 * this->Tke();

      const auto nuW =
          trans->Viscosity(this->Temperature(eqnState)) / this->Rho();
      const auto wWall = trans->NondimScaling() * trans->NondimScaling() *
                         60.0 * nuW / (wallDist * wallDist * turb->WallBeta());
      ghostState.data_[6] = 2.0 * wWall - this->Omega();

      if (layer > 1) {
        ghostState.data_[6] = layer * ghostState.data_[6] - wWall;
      }
    }

  // characteristic boundary condition
  // -------------------------------------------------------------------------
  // this is a characteristic based boundary condition which is appropriate for
  // subsonic and supersonic flow. It automatically switches between inflow
  // and outflow as needed. It also automatically determines the number of
  // characteristics that need to be specified at the boundary (subsonic vs
  // supersonic)
  } else if (bcType == "characteristic") {
    const auto & bcData = inputVars.BCData(tag);
    // freestream variables
    const auto freeVel = bcData->Velocity();
    const primVars freeState(bcData->Density(), freeVel, bcData->Pressure());

    // internal variables
    const auto velIntNorm = this->Velocity().DotProd(normArea);
    const auto SoSInt = this->SoS(thermo, eqnState);
    const auto machInt = fabs(velIntNorm) / SoSInt;

    if (machInt >= 1.0 && velIntNorm < 0.0) {  // supersonic inflow
      // -----------------------------------------------------
      // characteristics all go into the domain, so use freestream values for
      // both riemann invariants
      ghostState = freeState;

      // assign farfield conditions to turbulence variables
      if (inputVars.IsRANS()) {
        ghostState.ApplyFarfieldTurbBC(freeVel,
                                       bcData->TurbulenceIntensity(),
                                       bcData->EddyViscosityRatio(), trans,
                                       eqnState, turb);
      }

    } else if (machInt >= 1.0 && velIntNorm >= 0.0) {  // supersonic outflow
      // ----------------------------------------------
      // characteristics all leave the domain, so use interior values for both
      // riemann invariants
      ghostState = (*this);
    } else if (machInt < 1.0 && velIntNorm < 0.0) {  // subsonic inflow
      // ----------------------------------------------
      // characteristics go in both directions, use interior values for plus
      // characteristic and freestream values for minus characteristic
      const auto rhoSoSInt = this->Rho() * SoSInt;
      const auto velDiff = freeState.Velocity() - this->Velocity();

      // plus characteristic
      ghostState.data_[4] = 0.5 * (freeState.P() + this->P() -
                                   rhoSoSInt * normArea.DotProd(velDiff));
      const auto deltaPressure = freeState.P() - ghostState.P();

      // minus characteristic
      ghostState.data_[0] = freeState.Rho() - deltaPressure / (SoSInt * SoSInt);
      ghostState.data_[1] =
          freeState.U() - normArea.X() * deltaPressure / rhoSoSInt;
      ghostState.data_[2] =
          freeState.V() - normArea.Y() * deltaPressure / rhoSoSInt;
      ghostState.data_[3] =
          freeState.W() - normArea.Z() * deltaPressure / rhoSoSInt;

      // assign farfield conditions to turbulence variables
      if (inputVars.IsRANS()) {
        ghostState.ApplyFarfieldTurbBC(freeVel,
                                       bcData->TurbulenceIntensity(),
                                       bcData->EddyViscosityRatio(), trans,
                                       eqnState, turb);
      }

    } else if (machInt < 1.0 && velIntNorm >= 0.0) {  // subsonic outflow
      // ----------------------------------------------------------
      // characteristics go in both directions, use interior values for plus
      // characteristic and freestream values for minus characteristic
      const auto rhoSoSInt = this->Rho() * SoSInt;
      const auto deltaPressure = this->P() - freeState.P();

      // plus characteristic
      ghostState.data_[0] = this->Rho() - deltaPressure / (SoSInt * SoSInt);
      ghostState.data_[1] =
          this->U() + normArea.X() * deltaPressure / rhoSoSInt;
      ghostState.data_[2] =
          this->V() + normArea.Y() * deltaPressure / rhoSoSInt;
      ghostState.data_[3] =
          this->W() + normArea.Z() * deltaPressure / rhoSoSInt;
      ghostState.data_[4] = freeState.P();  // minus characteristic

      // numerical bcs for turbulence variables

    } else {
      cerr << "ERROR: flow condition for characteristic BC is not recognized!"
           << endl;
      cerr << "Interior state: " << (*this) << endl;
      cerr << "Ghost state: " << ghostState << endl;
      exit(EXIT_FAILURE);
    }

    // extrapolate from boundary to ghost cell
    ghostState = 2.0 * ghostState - (*this);

    if (layer > 1) {  // extrapolate to get ghost state at deeper layers
      ghostState = layer * ghostState - (*this);

      // assign farfield conditions to turbulence variables
      if (inputVars.IsRANS()) {
        ghostState.ApplyFarfieldTurbBC(freeVel,
                                       bcData->TurbulenceIntensity(),
                                       bcData->EddyViscosityRatio(), trans,
                                       eqnState, turb);
      }
    }

  // supersonic inflow boundary condition
  // -------------------------------------------------------------------------
  // this boundary condition enforces the entire state as the specified
  // freestream state
  } else if (bcType == "supersonicInflow") {
    const auto & bcData = inputVars.BCData(tag);
    // physical boundary conditions - fix everything
    const auto vel = bcData->Velocity();

    ghostState.data_[0] = bcData->Density();
    ghostState.data_[1] = vel.X();
    ghostState.data_[2] = vel.Y();
    ghostState.data_[3] = vel.Z();
    ghostState.data_[4] = bcData->Pressure();

    // assign farfield conditions to turbulence variables
    if (inputVars.IsRANS()) {
      ghostState.ApplyFarfieldTurbBC(vel, bcData->TurbulenceIntensity(),
                                     bcData->EddyViscosityRatio(), trans,
                                     eqnState, turb);
    }

    // supersonic outflow boundary condition
    // --------------------------------------------------------------------------
    // this boundary condition enforces the entire state as extrapolated from
    // the
    // interior (zeroth order extrapolation)
  } else if (bcType == "supersonicOutflow") {
    // do nothing and return boundary state -- numerical BCs for all
    if (layer > 1) {  // extrapolate to get ghost state at deeper layers
      ghostState = layer * ghostState - (*this);
    }

  // stagnation inlet boundary condition
  // --------------------------------------------------------------------------
  // this boundary condition is appropriate for subsonic flow. It is
  // particularly well suited for internal flows. Implementation from Blazek
  } else if (bcType == "stagnationInlet") {
    const auto & bcData = inputVars.BCData(tag);

    const auto t = this->Temperature(eqnState);
    const auto g = thermo->Gamma(t) - 1.0;
    // calculate outgoing riemann invarient
    const auto rNeg = this->Velocity().DotProd(normArea) -
        2.0 * this->SoS(thermo, eqnState) / g;

    // calculate SoS on boundary
    const auto cosTheta = -1.0 * this->Velocity().DotProd(normArea) /
        this->Velocity().Mag();
    const auto stagSoSsq = pow(this->SoS(thermo, eqnState), 2.0) +
        0.5 * g * this->Velocity().MagSq();

    const auto sosB = -1.0 * rNeg * g / (g * cosTheta * cosTheta + 2.0) *
        (1.0 + cosTheta * sqrt((g * cosTheta * cosTheta + 2.0) *
                               stagSoSsq / (g * rNeg * rNeg) - 0.5 * g));
    const auto tb = bcData->StagnationTemperature() * (sosB * sosB / stagSoSsq);
    const auto pb = bcData->StagnationPressure() *
                    pow(sosB * sosB / stagSoSsq, thermo->Gamma(t) / g);
    const auto vbMag = sqrt(2.0 / g * (bcData->StagnationTemperature() - tb));

    ghostState.data_[0] = eqnState->DensityTP(tb, pb);
    ghostState.data_[1] = vbMag * bcData->Direction().X();
    ghostState.data_[2] = vbMag * bcData->Direction().Y();
    ghostState.data_[3] = vbMag * bcData->Direction().Z();
    ghostState.data_[4] = pb;

    // assign farfield conditions to turbulence variables
    if (inputVars.IsRANS()) {
      ghostState.ApplyFarfieldTurbBC(ghostState.Velocity(),
                                     bcData->TurbulenceIntensity(),
                                     bcData->EddyViscosityRatio(), trans,
                                     eqnState, turb);
    }

    // extrapolate from boundary to ghost cell
    ghostState = 2.0 * ghostState - (*this);

    if (layer > 1) {  // extrapolate to get ghost state at deeper layers
      ghostState = layer * ghostState - (*this);

      // assign farfield conditions to turbulence variables
      if (inputVars.IsRANS()) {
        ghostState.ApplyFarfieldTurbBC(ghostState.Velocity(),
                                       bcData->TurbulenceIntensity(),
                                       bcData->EddyViscosityRatio(), trans,
                                       eqnState, turb);
      }
    }

  // pressure outlet boundary condition
  // -----------------------------------------------------------------------
  // this boundary condition is appropriate for subsonic flow. Implementation
  // from Blazek
  } else if (bcType == "pressureOutlet") {
    const auto & bcData = inputVars.BCData(tag);

    // nondimensional pressure from input file
    const auto pb = bcData->Pressure();

    const auto SoSInt = this->SoS(thermo, eqnState);
    const auto rhoSoSInt = this->Rho() * SoSInt;
    
    if (bcData->IsNonreflecting()) {
      // calculate LODI terms
      const auto deltaVel =
          (this->Velocity() - stateN.Velocity()).DotProd(normArea);
      constexpr auto sigma = 0.25;
      const auto rhoN = stateN.Rho();
      const auto sosN = stateN.SoS(thermo, eqnState);
      const auto rhoSoSN = rhoN * sosN;
      const auto length = bcData->LengthScale();
      const auto k = sigma * sosN * (1.0 - maxMach * maxMach) / length;

      // calculate transverse terms
      const auto beta = avgMach;
      const auto pGradT = pressGrad - pressGrad.DotProd(normArea) * normArea;
      const auto velT =
          stateN.Velocity() - stateN.Velocity().DotProd(normArea) * normArea;
      const auto velGradT = velGrad.RemoveComponent(normArea);
      const auto dVelN_dTrans = velGradT.LinearCombination(normArea);
      const auto dVelT_dTrans = velGradT.Sum() - dVelN_dTrans.SumElem();
      const auto gamma = thermo->Gamma(stateN.Temperature(eqnState));

      const auto trans = -0.5 * (velT.DotProd(pGradT - rhoSoSN * dVelN_dTrans) +
                                 gamma * stateN.P() * dVelT_dTrans);

      ghostState.data_[4] =
          (stateN.P() + rhoSoSN * deltaVel + dt * k * pb - dt * beta * trans) /
          (1.0 + dt * k);
    } else {
      ghostState.data_[4] = pb;
    }

    const auto deltaPressure = this->P() - ghostState.P();
    ghostState.data_[0] = this->Rho() - deltaPressure / (SoSInt * SoSInt);
    ghostState.data_[1] = this->U() + normArea.X() * deltaPressure / rhoSoSInt;
    ghostState.data_[2] = this->V() + normArea.Y() * deltaPressure / rhoSoSInt;
    ghostState.data_[3] = this->W() + normArea.Z() * deltaPressure / rhoSoSInt;

    // numerical bcs for turbulence variables

    // check for supersonic flow
    if (ghostState.Velocity().DotProd(normArea) /
            ghostState.SoS(thermo, eqnState) >=
        1.0) {
      ghostState = *this;
    }

    // extrapolate from boundary to ghost cell
    ghostState = 2.0 * ghostState - (*this);

    if (layer > 1) {  // extrapolate to get ghost state at deeper layers
      ghostState = layer * ghostState - (*this);
    }

  // connection boundary condition
  // --------------------------------------------------------------------------
  // this boundary condition is appropriate for point matched interfaces between
  // physical blocks or processor blocks
  } else if (bcType == "interblock" || "periodic") {
    // do nothing -- assign interior state to ghost state (already done)
    // for second layer of ghost cells interior state should be 2nd interior
    // cell

  } else {
    cerr << "ERROR: Error in primVars::GetGhostState ghost state for BC type "
         << bcType << " is not supported!" << endl;
    cerr << "surface is " << surf << " and layer is " << layer << endl;
    exit(EXIT_FAILURE);
  }

  return ghostState;
}

// member function to take in a genArray of updates to the conservative
// variables, and update the primative variables with it.
// this is used in the implicit solver
primVars primVars::UpdateWithConsVars(const unique_ptr<eos> &eqnState,
                                      const unique_ptr<thermodynamic> &thermo,
                                      const genArray &du,
                                      const unique_ptr<turbModel> &turb) const {
  // eqnState -- equation of state
  // du -- updates to conservative variables
  // turb -- turbulence model

  // convert primative to conservative and update
  const auto consUpdate = this->ConsVars(eqnState, thermo) + du;

  return primVars(consUpdate, false, eqnState, thermo, turb);
}

bool primVars::IsZero() const {
  auto nonzero = false;
  for (auto &var : data_) {
    if (var != 0.0) {
      nonzero = true;
      break;
    }
  }
  return !nonzero;
}

// member function to apply farfield turbulence boundary conditions
// using the method in STAR-CCM+ involving turbulence intensity and
// eddy viscosity ratio
void primVars::ApplyFarfieldTurbBC(const vector3d<double> &vel,
                                   const double &turbInten,
                                   const double &viscRatio,
                                   const unique_ptr<transport> &trans,
                                   const unique_ptr<eos> &eqnState,
                                   const unique_ptr<turbModel> &turb) {
  // vel -- reference velocity (nondimensionalized)
  // turbInten -- turbulence intensity at farfield
  // viscRatio -- eddy viscosity ratio at farfield
  // trans -- viscous transport model
  // eqnState -- equation of state
  // turb --  turbulence model

  data_[5] = 1.5 * pow(turbInten * vel.Mag(), 2.0);
  data_[6] = data_[0] * data_[5] /
      (viscRatio * trans->Viscosity(this->Temperature(eqnState)));
  this->LimitTurb(turb);
}

void primVars::LimitTurb(const unique_ptr<turbModel> &turb) {
  // Adjust turbulence variables to be above minimum if necessary
  data_[5] = max(data_[5], turb->TkeMin());
  data_[6] = max(data_[6], turb->OmegaMin());
}

/*Function to return the inviscid spectral radius for one direction (i, j, or k)
given a cell state, equation of state, and 2 face area vectors

L = 0.5 * (A1 + A2) * (|Vn| + SoS)

In the above equation L is the spectral radius in either the i, j, or k
direction. A1 and A2 are the two face areas in that direction. Vn is the
cell velocity normal to that direction. SoS is the speed of sound at the cell
 */
double primVars::InvCellSpectralRadius(const unitVec3dMag<double> &fAreaL,
                                       const unitVec3dMag<double> &fAreaR,
                                       const unique_ptr<thermodynamic> &thermo,
                                       const unique_ptr<eos> &eqnState) const {
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
  return (fabs(this->Velocity().DotProd(normAvg)) +
          this->SoS(thermo, eqnState)) *
         fMag;
}

double primVars::InvFaceSpectralRadius(const unitVec3dMag<double> &fArea,
                                       const unique_ptr<thermodynamic> &thermo,
                                       const unique_ptr<eos> &eqnState) const {
  // fArea -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state

  // return spectral radius
  return 0.5 * fArea.Mag() * (fabs(this->Velocity().DotProd(fArea.UnitVector()))
                              + this->SoS(thermo, eqnState));
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
double primVars::ViscCellSpectralRadius(
    const unitVec3dMag<double> &fAreaL, const unitVec3dMag<double> &fAreaR,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<eos> &eqnState,
    const unique_ptr<transport> &trans, const double &vol, const double &mu,
    const double &mut, const unique_ptr<turbModel> &turb) const {
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
  const auto t = this->Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * this->Rho()), thermo->Gamma(t) / this->Rho());
  // viscous term
  const auto viscTerm = trans->NondimScaling() *
      (mu / thermo->Prandtl(t) +  mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return maxTerm * viscTerm * fMag * fMag / vol;
}

double primVars::ViscFaceSpectralRadius(
    const unitVec3dMag<double> &fArea, const unique_ptr<thermodynamic> &thermo,
    const unique_ptr<eos> &eqnState, const unique_ptr<transport> &trans,
    const double &dist, const double &mu, const double &mut,
    const unique_ptr<turbModel> &turb) const {
  // fArea -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // dist -- distacne from cell center to cell center
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model

  const auto t = this->Temperature(eqnState);
  const auto maxTerm =
      max(4.0 / (3.0 * this->Rho()), thermo->Gamma(t) / this->Rho());
  // viscous term
  const auto viscTerm = trans->NondimScaling() *
      (mu / thermo->Prandtl(t) +  mut / turb->TurbPrandtlNumber());

  // return viscous spectral radius
  return fArea.Mag() / dist * maxTerm * viscTerm;
}

double primVars::CellSpectralRadius(
    const unitVec3dMag<double> &fAreaL, const unitVec3dMag<double> &fAreaR,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<eos> &eqnState,
    const unique_ptr<transport> &trans, const double &vol, const double &mu,
    const double &mut, const unique_ptr<turbModel> &turb,
    const bool &isViscous) const {
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

  auto specRad = this->InvCellSpectralRadius(fAreaL, fAreaR, thermo, eqnState);

  if (isViscous) {
    // factor 2 2 because viscous spectral radius is not halved (Blazek 6.53)
    specRad += 2.0 *
               this->ViscCellSpectralRadius(fAreaL, fAreaR, thermo, eqnState,
                                            trans, vol, mu, mut, turb);
  }
  return specRad;
}

double primVars::FaceSpectralRadius(const unitVec3dMag<double> &fArea,
                                    const unique_ptr<thermodynamic> &thermo,
                                    const unique_ptr<eos> &eqnState,
                                    const unique_ptr<transport> &trans,
                                    const double &dist, const double &mu,
                                    const double &mut,
                                    const unique_ptr<turbModel> &turb,
                                    const bool &isViscous) const {
  // fAreaL -- face area
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // trans -- viscous transport model
  // dist -- distance from cell center to cell center
  // mu -- laminar viscosity
  // mut -- turbulent viscosity
  // turb -- turbulence model
  // isViscous -- flag that is true if simulation is viscous

  auto specRad = this->InvFaceSpectralRadius(fArea, thermo, eqnState);

  if (isViscous) {
    specRad += this->ViscFaceSpectralRadius(fArea, thermo, eqnState, trans,
                                            dist, mu, mut, turb);
  }
  return specRad;
}


// function to calculate the Roe averaged state
primVars RoeAveragedState(const primVars &left, const primVars &right) {
  // compute Rho averaged quantities
  // density ratio
  const auto denRatio = sqrt(right.Rho() / left.Rho());
  // Roe averaged density
  const auto rhoR = left.Rho() * denRatio;
  // Roe averaged velocities - u, v, w
  const auto uR = (left.U() + denRatio * right.U()) / (1.0 + denRatio);
  const auto vR = (left.V() + denRatio * right.V()) / (1.0 + denRatio);
  const auto wR = (left.W() + denRatio * right.W()) / (1.0 + denRatio);

  // Roe averaged pressure
  const auto pR = (left.P() + denRatio * right.P()) / (1.0 + denRatio);

  // Roe averaged tke
  const auto kR = (left.Tke() + denRatio * right.Tke()) / (1.0 + denRatio);
  // Roe averaged specific dissipation (omega)
  const auto omR = (left.Omega() + denRatio * right.Omega()) / (1.0 + denRatio);

  return primVars(rhoR, uR, vR, wR, pR, kR, omR);
}

// return element by element squared values
primVars primVars::Squared() const {
  return (*this) * (*this);
}

// return element by element absolute value
primVars primVars::Abs() const {
  auto abs = *this;
  for (auto &val : abs.data_) {
    val = fabs(val);
  }
  return abs;
}

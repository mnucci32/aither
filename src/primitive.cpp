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
#include "primitive.hpp"
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

primitive::primitive(const conserved &cons, const unique_ptr<eos> &eqnState,
                     const unique_ptr<thermodynamic> &thermo,
                     const unique_ptr<turbModel> &turb) {
  // cons -- array of conserved variables
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // turb -- turbulence model

  *this = primitive(cons.Size(), cons.NumSpecies());

  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    (*this)[ii] = cons.RhoN(ii);
  }
  
  const auto rho = cons.Rho();
  (*this)[this->MomentumXIndex()] = cons.RhoU() / rho;
  (*this)[this->MomentumYIndex()] = cons.RhoV() / rho;
  (*this)[this->MomentumZIndex()] = cons.RhoW() / rho;
  
  const auto energy = cons.RhoE() / rho;
  (*this)[this->EnergyIndex()] =
      eqnState->PressFromEnergy(thermo, rho, energy, this->Velocity().Mag());

  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] = cons.RhoTurbN(ii) / rho;
  }

  // Adjust turbulence variables to be above minimum if necessary
  this->LimitTurb(turb);
}

// member function to initialize a state with nondimensional values
void primitive::NondimensionalInitialize(const unique_ptr<eos> &eqnState,
                                        const input &inp,
                                        const unique_ptr<transport> &trans,
                                        const int &parBlock,
                                        const unique_ptr<turbModel> &turb) {
  // get initial condition state for parent block
  auto ic = inp.ICStateForBlock(parBlock);
  auto mf = ic.MassFractions();

  for (auto pi = this->begin(), auto mi = std::begin(mf);
       ii < this->begin() + this->NumSpecies(); ++pi, ++mi) {
    *pi = *mi * ic.Density();
  }

  (*this)[this->MomentumXIndex()] = ic.Velocity().X();
  (*this)[this->MomentumYIndex()] = ic.Velocity().Y();
  (*this)[this->MomentumZIndex()] = ic.Velocity().Z();

  (*this)[this->EnergyIndex()] = ic.Pressure();
  
  if (this->HasTurbulenceData()) {
    // Initialize turbulence quantities based off of specified turublence
    // intensity and eddy viscosity ratio. This is the default for
    // STAR-CCM+
    this->ApplyFarfieldTurbBC(this->Velocity(), ic.TurbulenceIntensity(),
                              ic.EddyViscosityRatio(), trans, eqnState, turb);
  }
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const primitive &prim) {
  for (auto rr = 0; rr < prim.Size(); rr++) {
    os << prim[rr] << endl;
  }
  return os;
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
primitive primitive::GetGhostState(
    const string &bcType, const vector3d<double> &areaVec,
    const double &wallDist, const int &surf, const input &inputVars,
    const int &tag, const unique_ptr<eos> &eqnState,
    const unique_ptr<thermodynamic> &thermo, const unique_ptr<transport> &trans,
    const unique_ptr<turbModel> &turb, wallVars &wVars, const int &layer,
    const double &dt, const primitive &stateN, const vector3d<double> &pressGrad,
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

  // the instance of primitive being acted upon should be the interior cell
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

    ghostState[ghostState.MomentumXIndex()] = ghostVel.X();
    ghostState[ghostState.MomentumYIndex()] = ghostVel.Y();
    ghostState[ghostState.MomentumZIndex()] = ghostVel.Z();

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
    ghostState[ghostState.MomentumXIndex()] = ghostVel.X();
    ghostState[ghostState.MomentumYIndex()] = ghostVel.Y();
    ghostState[ghostState.MomentumZIndex()] = ghostVel.Z();

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
          const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
          for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
            ghostState[ii] = rho * this->MassFraction(ii);
          }
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
          const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
          for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
            ghostState[ii] = rho * this->MassFraction(ii);
          }
        }

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState[ghostState.TurbulenceIndex()] =
              2.0 * wVars.tke_ - this->Tke();
          ghostState[ghostState.TurbulenceIndex() + 1] =
              2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState[ghostState.TurbulenceIndex()] =
                layer * ghostState[ghostState.TurbulenceIndex()] - wVars.tke_;
            ghostState[ghostState.TurbulenceIndex() + 1] =
                layer * ghostState[ghostState.TurbulenceIndex() + 1] -
                wVars.sdr_;
          }
        }
      } else {  // low-Re wall treatment
        const auto tGhost = 2.0 * tWall - this->Temperature(eqnState);
        const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
        for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
          ghostState[ii] = rho * this->MassFraction(ii);
        }
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
          const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
          for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
            ghostState[ii] = rho * this->MassFraction(ii);
          }
        } else {
          // use wall law wall temperature to get ghost cell density
          const auto tGhost =
              2.0 * wVars.temperature_ - this->Temperature(eqnState);
          const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
          for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
            ghostState[ii] = rho * this->MassFraction(ii);
          }
        }

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState[ghostState.TurbulenceIndex()] =
              2.0 * wVars.tke_ - this->Tke();
          ghostState[ghostState.TurbulenceIndex() + 1] =
              2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState[ghostState.TurbulenceIndex()] =
                layer * ghostState[ghostState.TurbulenceIndex()] - wVars.tke_;
            ghostState[ghostState.TurbulenceIndex() + 1] =
                layer * ghostState[ghostState.TurbulenceIndex() + 1] -
                wVars.sdr_;
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
        const auto rho = eqnState->DensityTP(tGhost, ghostState.P());
        for (auto ii = 0; ii < ghostState.NumSpecies(); ++ii) {
          ghostState[ii] = rho * this->MassFraction(ii);
        }
        // numerical BCs for pressure, same as boundary state
      }
    } else {  // default is adiabatic -----------------------------------------
      if (bcData->IsWallLaw()) {
        wallLaw wl(bcData->VonKarmen(), bcData->WallConstant(), *this, wallDist,
                   inputVars.IsRANS());
        wVars = wl.AdiabaticBCs(normArea, velWall, eqnState, thermo, trans,
                                turb, isLower);

        if (inputVars.IsRANS() && !wVars.SwitchToLowRe()) {
          ghostState[ghostState.TurbulenceIndex()] =
              2.0 * wVars.tke_ - this->Tke();
          ghostState[ghostState.TurbulenceIndex() + 1] =
              2.0 * wVars.sdr_ - this->Omega();
          if (layer > 1) {
            ghostState[ghostState.TurbulenceIndex()] =
                layer * ghostState[ghostState.TurbulenceIndex()] - wVars.tke_;
            ghostState[ghostState.TurbulenceIndex() + 1] =
                layer * ghostState[ghostState.TurbulenceIndex() + 1] -
                wVars.sdr_;
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
      ghostState[ghostState.TurbulenceIndex()] = -1.0 * this->Tke();

      const auto nuW =
          trans->Viscosity(this->Temperature(eqnState)) / this->Rho();
      const auto wWall = trans->NondimScaling() * trans->NondimScaling() *
                         60.0 * nuW / (wallDist * wallDist * turb->WallBeta());
      ghostState[ghostState.TurbulenceIndex() + 1] =
          2.0 * wWall - this->Omega();

      if (layer > 1) {
        ghostState[ghostState.TurbulenceIndex() + 1] =
            layer * ghostState[ghostState.TurbulenceIndex() + 1] - wWall;
      }
    }

    // characteristic boundary condition
    // -------------------------------------------------------------------------
    // this is a characteristic based boundary condition which is appropriate
    // for subsonic and supersonic flow. It automatically switches between
    // inflow and outflow as needed. It also automatically determines the number
    // of characteristics that need to be specified at the boundary (subsonic vs
    // supersonic)
  } else if (bcType == "characteristic") {
    const auto & bcData = inputVars.BCData(tag);
    // freestream variables
    const auto freeVel = bcData->Velocity();
    const primitive freeState(bcData->Density(), freeVel, bcData->Pressure());

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
      ghostState[ghostState.EnergyIndex()] = 0.5 * (freeState.P() + this->P() -
                                   rhoSoSInt * normArea.DotProd(velDiff));
      const auto deltaPressure = freeState.P() - ghostState.P();

      // minus characteristic
      ghostState.data_[0] = freeState.Rho() - deltaPressure / (SoSInt * SoSInt);
      ghostState[ghostState.MomentumXIndex()] =
          freeState.U() - normArea.X() * deltaPressure / rhoSoSInt;
      ghostState[ghostState.MomentumYIndex()] =
          freeState.V() - normArea.Y() * deltaPressure / rhoSoSInt;
      ghostState[ghostState.MomentumZIndex()] =
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
      ghostState[ghostState.MomentumXIndex()] =
          this->U() + normArea.X() * deltaPressure / rhoSoSInt;
      ghostState[ghostState.MomentumYIndex()] =
          this->V() + normArea.Y() * deltaPressure / rhoSoSInt;
      ghostState[ghostState.MomentumZIndex()] =
          this->W() + normArea.Z() * deltaPressure / rhoSoSInt;
      ghostState[ghostState.EnergyIndex()] = freeState.P();  // minus characteristic

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

  // inlet boundary condition
  // -------------------------------------------------------------------------
  } else if (bcType == "inlet") {
    const auto & bcData = inputVars.BCData(tag);
    // freestream variables
    const auto freeVel = bcData->Velocity();
    const primitive freeState(bcData->Density(), freeVel, bcData->Pressure());

    // internal variables
    const auto velIntNorm = this->Velocity().DotProd(normArea);
    const auto SoSInt = this->SoS(thermo, eqnState);
    const auto machInt = fabs(velIntNorm) / SoSInt;

    if (machInt >= 1.0) {  // supersonic inflow
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
    } else {  // subsonic inflow
      // ----------------------------------------------
      // characteristics go in both directions, use interior values for plus
      // characteristic and freestream values for minus characteristic
      const auto rhoSoSInt = this->Rho() * SoSInt;
      const auto velDiff = freeState.Velocity() - this->Velocity();

      // plus characteristic
      ghostState[ghostState.EnergyIndex()] =
          0.5 *
          (freeState.P() + this->P() - rhoSoSInt * normArea.DotProd(velDiff));

      if (bcData->IsNonreflecting()) {
        // minus characteristic
        // calculate LODI terms
        constexpr auto sigma = 0.25;
        const auto rhoN = stateN.Rho();
        const auto sosN = stateN.SoS(thermo, eqnState);
        const auto rhoSoSN = rhoN * sosN;
        const auto deltaPressure = ghostState.P() - stateN.P();
        const auto length = bcData->LengthScale();
        const auto alphaR = sigma / (sosN * length);

        ghostState.data_[0] = (rhoN + dt * alphaR * freeState.Rho() +
                               deltaPressure / (sosN * sosN)) /
                              (1.0 + dt * alphaR);

        const auto alpha = sigma * sosN / length;
        const auto k = alpha * (1.0 - maxMach * maxMach);
        auto vel = (stateN.Velocity() + dt * k * freeState.Velocity() -
                    normArea * deltaPressure / rhoSoSN) /
                   (1.0 + dt * k);

        ghostState[ghostState.MomentumXIndex()] = vel.X();
        ghostState[ghostState.MomentumYIndex()] = vel.Y();
        ghostState[ghostState.MomentumZIndex()] = vel.Z();

      } else {
        const auto deltaPressure = freeState.P() - ghostState.P();

        // minus characteristic
        ghostState.data_[0] =
            freeState.Rho() - deltaPressure / (SoSInt * SoSInt);
        ghostState[ghostState.MomentumXIndex()] =
            freeState.U() - normArea.X() * deltaPressure / rhoSoSInt;
        ghostState[ghostState.MomentumYIndex()] =
            freeState.V() - normArea.Y() * deltaPressure / rhoSoSInt;
        ghostState[ghostState.MomentumZIndex()] =
            freeState.W() - normArea.Z() * deltaPressure / rhoSoSInt;
      }

      // assign farfield conditions to turbulence variables
      if (inputVars.IsRANS()) {
        ghostState.ApplyFarfieldTurbBC(freeVel,
                                       bcData->TurbulenceIntensity(),
                                       bcData->EddyViscosityRatio(), trans,
                                       eqnState, turb);
      }
    }

    // extrapolate from boundary to ghost cell
    ghostState = 2.0 * ghostState - (*this);

    if (layer > 1) {  // extrapolate to get ghost state at deeper layers
      ghostState = layer * ghostState - (*this);
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
    ghostState[ghostState.MomentumXIndex()] = vel.X();
    ghostState[ghostState.MomentumYIndex()] = vel.Y();
    ghostState[ghostState.MomentumZIndex()] = vel.Z();
    ghostState[ghostState.EnergyIndex()] = bcData->Pressure();

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
    ghostState[ghostState.MomentumXIndex()] = vbMag * bcData->Direction().X();
    ghostState[ghostState.MomentumYIndex()] = vbMag * bcData->Direction().Y();
    ghostState[ghostState.MomentumZIndex()] = vbMag * bcData->Direction().Z();
    ghostState[ghostState.EnergyIndex()] = pb;

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

      ghostState[ghostState.EnergyIndex()] =
          (stateN.P() + rhoSoSN * deltaVel + dt * k * pb - dt * beta * trans) /
          (1.0 + dt * k);
    } else {
      ghostState[ghostState.EnergyIndex()] = pb;
    }

    const auto deltaPressure = this->P() - ghostState.P();
    ghostState.data_[0] = this->Rho() - deltaPressure / (SoSInt * SoSInt);
    ghostState[ghostState.MomentumXIndex()] = this->U() + normArea.X() * deltaPressure / rhoSoSInt;
    ghostState[ghostState.MomentumYIndex()] = this->V() + normArea.Y() * deltaPressure / rhoSoSInt;
    ghostState[ghostState.MomentumZIndex()] = this->W() + normArea.Z() * deltaPressure / rhoSoSInt;

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
    cerr << "ERROR: Error in primitive::GetGhostState ghost state for BC type "
         << bcType << " is not supported!" << endl;
    cerr << "surface is " << surf << " and layer is " << layer << endl;
    exit(EXIT_FAILURE);
  }

  return ghostState;
}

// member function to take in a genArray of updates to the conservative
// variables, and update the primitive variables with it.
// this is used in the implicit solver
primitive primitive::UpdateWithConsVars(
    const unique_ptr<eos> &eqnState, const unique_ptr<thermodynamic> &thermo,
    const conserved &du, const unique_ptr<turbModel> &turb) const {
  // eqnState -- equation of state
  // du -- updates to conservative variables
  // turb -- turbulence model

  // convert primitive to conservative and update
  const auto consUpdate = this->ConsVars(eqnState, thermo) + du;
  return primitive(consUpdate, eqnState, thermo, turb);
}

// member function to apply farfield turbulence boundary conditions
// using the method in STAR-CCM+ involving turbulence intensity and
// eddy viscosity ratio
void primitive::ApplyFarfieldTurbBC(const vector3d<double> &vel,
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

  (*this)[this->TurbulenceIndex()] = 1.5 * pow(turbInten * vel.Mag(), 2.0);
  (*this)[this->TurbulenceIndex() + 1] = this->Rho() * this->Tke() /
      (viscRatio * trans->Viscosity(this->Temperature(eqnState)));
  this->LimitTurb(turb);
}

void primitive::LimitTurb(const unique_ptr<turbModel> &turb) {
  // Adjust turbulence variables to be above minimum if necessary
  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] =
        max((*this)[TurbulenceIndex() + ii], turb->TurbMinN(ii));
  }
}

/*Function to return the inviscid spectral radius for one direction (i, j, or k)
given a cell state, equation of state, and 2 face area vectors

L = 0.5 * (A1 + A2) * (|Vn| + SoS)

In the above equation L is the spectral radius in either the i, j, or k
direction. A1 and A2 are the two face areas in that direction. Vn is the
cell velocity normal to that direction. SoS is the speed of sound at the cell
 */
double primitive::InvCellSpectralRadius(const unitVec3dMag<double> &fAreaL,
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

double primitive::InvFaceSpectralRadius(const unitVec3dMag<double> &fArea,
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
double primitive::ViscCellSpectralRadius(
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

double primitive::ViscFaceSpectralRadius(
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

double primitive::CellSpectralRadius(
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

double primitive::FaceSpectralRadius(const unitVec3dMag<double> &fArea,
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
primitive RoeAveragedState(const primitive &left, const primitive &right) {
  // compute Rho averaged quantities
  primitive rhoState(left.Size(), left.NumSpecies());
  // density ratio
  const auto denRatio = sqrt(right.Rho() / left.Rho());
  // Roe averaged density
  for (auto ii = 0; ii < rhoState.NumSpecies(); ++ii) {
    rhoState[ii] = left.RhoN(ii) * denRatio;
  }
  // Roe averaged velocities - u, v, w
  rhoState[rhoState->MomentumXIndex()] =
      (left.U() + denRatio * right.U()) / (1.0 + denRatio);
  rhoState[rhoState->MomentumYIndex()] =
      (left.V() + denRatio * right.V()) / (1.0 + denRatio);
  rhoState[rhoState->MomentumZIndex()] =
      (left.W() + denRatio * right.W()) / (1.0 + denRatio);

  // Roe averaged pressure
  rhoState[rhoState->EnergyIndex()] =
      (left.P() + denRatio * right.P()) / (1.0 + denRatio);

  // Roe averaged turbulence variables
  for (auto ii = 0; ii < rhoState.NumTurbulence(); ++ii) {
    rhoState[rhoState->TurbulenceIndex() + ii] =
        (left.TurbN(ii) + denRatio * right.TurbN(ii)) / (1.0 + denRatio);
  }

  return rhoState;
}

// return element by element absolute value
primitive primitive::Abs() const {
  auto abs = *this;
  std::transform(std::begin(*this), std::end(*this), std::begin(abs),
                 [](const double &val) { return fabs(val); });
  return abs;
}

/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include "procBlock.hpp"
#include "plot3d.hpp"              // plot3d
#include "eos.hpp"                 // idealGas
#include "inviscidFlux.hpp"        // inviscidFlux
#include "viscousFlux.hpp"         // viscousFlux
#include "input.hpp"               // inputVars
#include "turbulence.hpp"
#include "gradients.hpp"
#include "slices.hpp"
#include "source.hpp"
#include "resid.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::to_string;
using std::max;
using std::min;

// constructors for procBlock class
procBlock::procBlock() {
  numGhosts_ = 0;
  parBlock_ = 0;
  rank_ = 0;
  globalPos_ = 0;
  localPos_ = 0;

  // size vectors holding cell center and cell face values appropriately
  state_ = multiArray3d<primVars>(1, 1, 1);

  center_ = multiArray3d<vector3d<double> >(1, 1, 1);
  fAreaI_ = multiArray3d<unitVec3dMag<double> >(1, 1, 1);
  fAreaJ_ = multiArray3d<unitVec3dMag<double> >(1, 1, 1);
  fAreaK_ = multiArray3d<unitVec3dMag<double> >(1, 1, 1);
  fCenterI_ = multiArray3d<vector3d<double> >(1, 1, 1);
  fCenterJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
  fCenterK_ = multiArray3d<vector3d<double> >(1, 1, 1);

  residual_ = multiArray3d<genArray>(1, 1, 1);

  vol_ = multiArray3d<double>(1, 1, 1);
  avgWaveSpeed_ = multiArray3d<double>(1, 1, 1);
  dt_ = multiArray3d<double>(1, 1, 1);
  wallDist_ = multiArray3d<double>(1, 1, 1);

  bc_ = boundaryConditions();
}

// constructor -- assign passed state to initialize state vector
procBlock::procBlock(const primVars &inputState, const plot3dBlock &blk,
                     const int &numBlk, const int &numG,
                     const boundaryConditions &bound, const int &pos,
                     const int &r, const int &lpos) {
  // inputState -- state_ to initialize block with (primative)
  // blk -- plot3d block of which this procBlock is a subset of
  // numBlk -- the block number of blk (the parent block)
  // numG -- number of ghost cells
  // bound -- boundary conditions for block
  // pos -- global position of block, an identifying number unique to this block
  // r -- processor rank_ that procBlock should be on
  // lpos -- local position of block on processor

  numGhosts_ = numG;
  parBlock_ = numBlk;

  rank_ = r;
  globalPos_ = pos;
  localPos_ = lpos;

  bc_ = bound;

  // pad stored variable vectors with ghost cells
  state_ = PadWithGhosts(multiArray3d<primVars>(blk.NumI() - 1, blk.NumJ() - 1,
                                                blk.NumK() - 1, inputState),
                         numGhosts_);

  vol_ = PadWithGhosts(blk.Volume(), numGhosts_);
  center_ = PadWithGhosts(blk.Centroid(), numGhosts_);
  fAreaI_ = PadWithGhosts(blk.FaceAreaI(), numGhosts_);
  fAreaJ_ = PadWithGhosts(blk.FaceAreaJ(), numGhosts_);
  fAreaK_ = PadWithGhosts(blk.FaceAreaK(), numGhosts_);
  fCenterI_ = PadWithGhosts(blk.FaceCenterI(), numGhosts_);
  fCenterJ_ = PadWithGhosts(blk.FaceCenterJ(), numGhosts_);
  fCenterK_ = PadWithGhosts(blk.FaceCenterK(), numGhosts_);

  avgWaveSpeed_ = multiArray3d<double>(blk.NumI() - 1, blk.NumJ() - 1,
                                       blk.NumK() - 1);
  dt_ = multiArray3d<double>(blk.NumI() - 1, blk.NumJ() - 1, blk.NumK() - 1);
  wallDist_ = multiArray3d<double>(blk.NumI() - 1, blk.NumJ() - 1,
                                   blk.NumK() - 1, DEFAULTWALLDIST);
  residual_ = multiArray3d<genArray>(blk.NumI() - 1, blk.NumJ() - 1,
                                     blk.NumK() - 1);
}

// constructor -- allocate space for procBlock
procBlock::procBlock(const int &ni, const int &nj, const int &nk,
                     const int &numG) {
  // ni -- i-dimension (cell)
  // nj -- j-dimension (cell)
  // nk -- k-dimension (cell)
  // numG -- number of ghost cell layers

  numGhosts_ = numG;
  parBlock_ = 0;

  rank_ = 0;
  globalPos_ = 0;
  localPos_ = 0;

  bc_ = boundaryConditions();

  // pad stored variable vectors with ghost cells
  state_ = multiArray3d<primVars>(ni + 2 * numG, nj + 2 * numG, nk + 2 * numG);
  center_ = multiArray3d<vector3d<double> >(ni + 2 * numG, nj + 2 * numG,
                                            nk + 2 * numG);
  fAreaI_ = multiArray3d<unitVec3dMag<double> >(ni + 2 * numG + 1,
                                                nj + 2 * numG, nk + 2 * numG);
  fAreaJ_ = multiArray3d<unitVec3dMag<double> >(ni + 2 * numG,
                                                nj + 2 * numG + 1,
                                                nk + 2 * numG);
  fAreaK_ = multiArray3d<unitVec3dMag<double> >(ni + 2 * numG,
                                                nj + 2 * numG,
                                                nk + 2 * numG + 1);
  fCenterI_ = multiArray3d<vector3d<double> >(ni + 2 * numG + 1, nj + 2 * numG,
                                              nk + 2 * numG);
  fCenterJ_ = multiArray3d<vector3d<double> >(ni + 2 * numG, nj + 2 * numG + 1,
                                              nk + 2 * numG);
  fCenterK_ = multiArray3d<vector3d<double> >(ni + 2 * numG, nj + 2 * numG,
                                              nk + 2 * numG + 1);
  residual_ = multiArray3d<genArray>(ni, nj, nk);
  vol_ = multiArray3d<double>(ni + 2 * numG, nj + 2 * numG, nk + 2 * numG);
  avgWaveSpeed_ = multiArray3d<double>(ni, nj, nk);
  dt_ = multiArray3d<double>(ni, nj, nk);
  wallDist_ = multiArray3d<double>(ni, nj, nk);
}

// member function to add a member of the inviscid flux class to the residual
void procBlock::AddToResidual(const inviscidFlux &flux, const int &ii,
                              const int &jj, const int &kk) {
  // flux -- inviscid flux to add to residual
  // ii -- i-location of residual to add to
  // jj -- j-location of residual to add to
  // kk -- k-location of residual to add to

  (*this).residual_(ii, jj, kk)[0] += flux.RhoVel();
  (*this).residual_(ii, jj, kk)[1] += flux.RhoVelU();
  (*this).residual_(ii, jj, kk)[2] += flux.RhoVelV();
  (*this).residual_(ii, jj, kk)[3] += flux.RhoVelW();
  (*this).residual_(ii, jj, kk)[4] += flux.RhoVelH();
  (*this).residual_(ii, jj, kk)[5] += flux.RhoVelK();
  (*this).residual_(ii, jj, kk)[6] += flux.RhoVelO();
}

// member function to add a member of the viscous flux class to the residual_
void procBlock::AddToResidual(const viscousFlux &flux, const int &ii,
                              const int &jj, const int &kk) {
  // flux -- inviscid flux to add to residual_
  // ii -- location of residual_ to add to
  // jj -- j-location of residual to add to
  // kk -- k-location of residual to add to

  (*this).residual_(ii, jj, kk)[1] += flux.MomX();
  (*this).residual_(ii, jj, kk)[2] += flux.MomY();
  (*this).residual_(ii, jj, kk)[3] += flux.MomZ();
  (*this).residual_(ii, jj, kk)[4] += flux.Engy();
  (*this).residual_(ii, jj, kk)[5] += flux.MomK();
  (*this).residual_(ii, jj, kk)[6] += flux.MomO();
}

// member function to add a member of the inviscid source class to the residual
void procBlock::AddToResidual(const source &src, const int &ii,
                              const int &jj, const int &kk) {
  // src -- source to add to residual
  // ii -- location of residual to add to
  // jj -- j-location of residual to add to
  // kk -- k-location of residual to add to

  (*this).residual_(ii, jj, kk)[0] += src.SrcMass();
  (*this).residual_(ii, jj, kk)[1] += src.SrcMomX();
  (*this).residual_(ii, jj, kk)[2] += src.SrcMomY();
  (*this).residual_(ii, jj, kk)[3] += src.SrcMomZ();
  (*this).residual_(ii, jj, kk)[4] += src.SrcEngy();
  (*this).residual_(ii, jj, kk)[5] += src.SrcTke();
  (*this).residual_(ii, jj, kk)[6] += src.SrcOmg();
}

//---------------------------------------------------------------------
// function declarations

/* Function to calculate the inviscid fluxes on the i-faces. All phyiscal
(non-ghost) i-faces are looped over. The left and right states are
calculated, and then the flux at the face is calculated. The flux at the
face contributes to the residual of the cells to the left and right of
the face. This contribution from the flux is added to the residuals and
the wave speed is accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Ui       -->   Ui+1   |
  |            |            |
  |____________|____________|
Ui-1/2       Ui+1/2       Ui+3/2

Using the above diagram, the flux is calculated at face Ui+1/2. Since the area
vector at the face always points from lower indices to higher indices it
points from Ui to Ui+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Ui, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Ui+1.

The spectral radius in the i-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Ui+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.
*/
void procBlock::CalcInvFluxI(const idealGas &eqnState, const input &inp) {
  // eqnState -- equation of state
  // inp -- all input variables


  // loop over all physical i-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaI_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaI_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.g <
               fAreaI_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // use constant reconstruction (first order)
        if (inp.OrderOfAccuracy() == "first") {
          primVars faceStateLower = state_(ii.g - 1, jj.g, kk).FaceReconConst();
          primVars faceStateUpper = state_(ii.g, jj.g, kk).FaceReconConst();
        } else {  // second order accuracy -- use MUSCL extrapolation
          // length of second upwind, first upwind, and downwind cells in
          // i-direction
          double upwind2L = fCenterI_(ii.g - 1, jj.g, kk.g).
              Distance(fCenterI_(ii.g - 2, jj.g, kk.g));
          double upwindL = fCenterI_(ii.g, jj.g, kk.g).
              Distance(fCenterI_(ii.g - 1, jj.g, kk.g));
          double downwindL = fCenterI_(ii.g, jj.g, kk.g).
              Distance(fCenterI_(ii.g + 1, jj.g, kk.g));

          primVars faceStateLower = state_(ii.g - 1, jj.g, kk.g).FaceReconMUSCL(
              state_(ii.g - 2, jj.g, kk.g), state_(ii.g, jj.g, kk.g),
              inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL);

          // length of second upwind, first upwind, and downwind cells in
          // i-direction
          double upwind2U = fCenterI_(ii.g + 1, jj.g, kk.g).
              Distance(fCenterI_(ii.g + 2, jj.g, kk.g));
          double upwindU = fCenterI_(ii.g, jj.g, kk.g).
              Distance(fCenterI_(ii.g + 1, jj.g, kk.g));
          double downwindU = fCenterI_(ii.g, jj.g, kk.g).
              Distance(fCenterI_(ii.g - 1, jj.g, kk.g));

          primVars faceStateUpper = state_(ii.g, jj.g, kk.g).FaceReconMUSCL(
              state_(ii.g + 1, jj.g, kk.g), state_(ii.g - 1, jj.g, kk.g),
              inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU);
        }
        // calculate Roe flux at face
        inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper,
                                        eqnState,
                                        (*this).FAreaUnitI(ii.g, jj.g, kk.g));

        // area vector points from left to right, so add to left cell, subtract
        // from right cell
        // at left boundary there is no left cell to add to
        if (ii.g > numGhosts_) {
          (*this).AddToResidual(tempFlux * (*this).FAreaMagI(ii.g, jj.g, kk.g),
                                ii.p - 1, jj.p, kk.p);
        }
        // at right boundary there is no right cell to add to
        if (ii.g < fAreaI_.NumI() - 2 * numGhosts_ - 1) {
          (*this).AddToResidual(-1.0 * tempFlux *
                                (*this).FAreaMagI(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);
          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_(ii.p, jj.p, kk.p) +=
              CellSpectralRadius((*this).FAreaI(ii.g, jj.g, kk.g),
                                 (*this).FAreaI(ii.g + 1, jj.g, kk.g),
                                 state_(ii.g, jj.g, kk.g), eqnState);
        }
      }
    }
  }
}

/* Function to calculate the inviscid fluxes on the j-faces. All physical
(non-ghost) j-faces are looped over. The left and right states are
calculated, and then the flux at the face is calculated. The flux at the
face contributes to the residual of the cells to the left and right of the
face. This contribution from the flux is added to the residuals and the wave
speed is accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Uj       -->   Uj+1   |
  |            |            |
  |____________|____________|
Uj-1/2       Uj+1/2       Uj+3/2

Using the above diagram, the flux is calculated at face Uj+1/2. Since the area
vector at the face always points from lower indices to higher indices it points
from Uj to Uj+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Uj, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Uj+1.

The spectral radius in the j-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Uj+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.
*/
void procBlock::CalcInvFluxJ(const idealGas &eqnState, const input &inp) {
  // eqnState -- equation of state
  // inp -- all input variables

  // loop over all physical j-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaJ_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaJ_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.g <
               fAreaJ_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // use constant reconstruction (first order)
        if (inp.OrderOfAccuracy() == "first") {
          primVars faceStateLower = state_(ii.g, jj.g - 1, kk).FaceReconConst();
          primVars faceStateUpper = state_(ii.g, jj.g, kk).FaceReconConst();
        } else {  // second order accuracy -- use MUSCL extrapolation
          // length of second upwind, first upwind, and downwind cells in
          // j-direction
          double upwind2L = fCenterJ_(ii.g, jj.g - 1, kk.g).
              Distance(fCenterJ_(ii.g, jj.g - 2, kk.g));
          double upwindL = fCenterJ_(ii.g, jj.g, kk.g).
              Distance(fCenterJ_(ii.g, jj.g - 1, kk.g));
          double downwindL = fCenterJ_(ii.g, jj.g, kk.g).
              Distance(fCenterJ_(ii.g, jj.g + 1, kk.g));

          primVars faceStateLower = state_(ii.g, jj.g - 1, kk.g).FaceReconMUSCL(
              state_(ii.g, jj.g - 2, kk.g), state_(ii.g, jj.g, kk.g),
              inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL);

          // length of second upwind, first upwind, and downwind cells in
          // j-direction
          double upwind2U = fCenterJ_(ii.g, jj.g + 1, kk.g).
              Distance(fCenterJ_(ii.g, jj.g + 2, kk.g));
          double upwindU = fCenterJ_(ii.g, jj.g, kk.g).
              Distance(fCenterJ_(ii.g, jj.g + 1, kk.g));
          double downwindU = fCenterJ_(ii.g, jj.g, kk.g).
              Distance(fCenterJ_(ii.g, jj.g - 1, kk.g));

          primVars faceStateUpper = state_(ii.g, jj.g, kk.g).FaceReconMUSCL(
              state_(ii.g, jj.g + 1, kk.g), state_(ii.g, jj.g - 1, kk.g),
              inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU);
        }

        // calculate Roe flux at face
        inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper,
                                        eqnState,
                                        (*this).FAreaUnitJ(ii.g, jj.g, kk.g));

        // area vector points from left to right, so add to left cell, subtract
        // from right cell
        // at left boundary no left cell to add to
        if (jj.g > numGhosts_) {
          (*this).AddToResidual(tempFlux * (*this).FAreaMagJ(ii.g, jj.g, kk.g),
                                ii.p, jj.p - 1, kk.p);
        }
        // at right boundary no right cell to add to
        if (jj.g < fAreaJ_.NumJ() - 2 * numGhosts_ - 1) {
          (*this).AddToResidual(-1.0 * tempFlux *
                                (*this).FAreaMagJ(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);

          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_(ii.p, jj.p, kk.p) +=
              CellSpectralRadius((*this).FAreaJ(ii.g, jj.g, kk.g),
                                 (*this).FAreaJ(ii.g, jj.g + 1, kk.g),
                                 state_(ii.g, jj.g, kk.g), eqnState);
        }
      }
    }
  }
}

/* Function to calculate the inviscid fluxes on the k-faces. All phyiscal
(non-ghost) k-faces are looped over. The left and right states are calculated,
and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This
contribution from the flux is added to the residuals and the wave speed is
accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Uk       -->   Uk+1   |
  |            |            |
  |____________|____________|
Uk-1/2       Uk+1/2       Uk+3/2

Using the above diagram, the flux is calculated at face Uk+1/2. Since the area
vector at the face always points from lower indices to higher indices it points
from Uk to Uk+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Uk, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Uk+1.

The spectral radius in the k-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Uk+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.
*/
void procBlock::CalcInvFluxK(const idealGas &eqnState, const input &inp) {
  // eqnState -- equation of state
  // inp -- all input variables

  // loop over all physical k-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaK_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaK_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} kk = {0, numGhosts_}; ii.g <
               fAreaK_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // use constant reconstruction (first order)
        if (inp.OrderOfAccuracy() == "first") {
          faceStateLower = state_(ii.g, jj.g, kk.g - 1).FaceReconConst();
          faceStateUpper = state_(ii.g, jj.g, kk.g).FaceReconConst();
        } else {  // second order accuracy -- use MUSCL extrapolation
          // length of second upwind, first upwind, and downwind cells in
          // j-direction
          double upwind2L = fCenterK_(ii.g, jj.g, kk.g - 1).
              Distance(fCenterK_(ii.g, jj.g, kk.g - 2));
          double upwindL = fCenterK_(ii.g, jj.g, kk.g).
              Distance(fCenterK_(ii.g, jj.g, kk.g - 1));
          double downwindL = fCenterK_(ii.g, jj.g, kk.g).
              Distance(fCenterK_(ii.g, jj.g, kk.g + 1));

          primVars faceStateLower = state_(ii.g, jj.g, kk.g - 1).FaceReconMUSCL(
              state_(ii.g, jj.g, kk.g - 2), state_(ii.g, jj.g, kk.g),
              inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL);

          // length of second upwind, first upwind, and downwind cells in
          // j-direction
          double upwind2U = fCenterK_(ii.g, jj.g, kk.g + 1).
              Distance(fCenterK_(ii.g, jj.g, kk.g + 2));
          double upwindU = fCenterK_(ii.g, jj.g, kk.g).
              Distance(fCenterK_(ii.g, jj.g, kk.g + 1));
          double downwindU = fCenterK_(ii.g, jj.g, kk.g).
              Distance(fCenterK_(ii.g, jj.g, kk.g - 1));

          primVars faceStateUpper = state_(ii.g, jj.g, kk.g).FaceReconMUSCL(
              state_(ii.g, jj.g, kk.g + 1), state_(ii.g, jj.g, kk.g - 1),
              inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU);
        }

        // calculate Roe flux at face
        inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper,
                                        eqnState,
                                        (*this).FAreaUnitK(ii.g, jj.g, kk.g));

        // area vector points from left to right, so add to left cell, subtract
        // from right cell
        // at left boundary no left cell to add to
        if (kk.g > numGhosts_) {
          (*this).AddToResidual(tempFlux *
                                (*this).FAreaMagK(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p - 1);
        }
        // at right boundary no right cell to add to
        if (kk.g < fAreaK_.NumK() - 2 * numGhosts_ - 1) {
          (*this).AddToResidual(-1.0 * tempFlux *
                                (*this).FAreaMagK(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);

          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_(ii.p, jj.p, kk.p) +=
              CellSpectralRadius((*this).FAreaK(ii.g, jj.g, kk.g),
                                 (*this).FAreaK(ii.g, jj.g, kk.g + 1),
                                 state_(ii.g, jj.g, kk.g), eqnState);
        }
      }
    }
  }
}

/* Member function to calculate the local time step. (i,j,k) are cell indices.
The following equation is used:

dt = CFL * V / (Lci + Lcj + Lck + C * (Lvi + Lvj + Lvk)) (Blazek 6.18)

In the above equation dt is the time step, CFL is the CFL number, V is the cell
volume, Lci, Lcj, Lck are the convective (inviscid) spectral radii in the i, j,
and k directions, C is a constant (typical value b/w 1 and 4), and Lvi, Lvj,
Lvk are the viscous spectral radii. This function is only used when the time
step isn't explicitly defined by the user.
*/
void procBlock::CalcCellDt(const int &ii, const int &jj, const int &kk,
                           const double &cfl) {
  // ii -- i index of cell
  // jj -- j index of cell
  // kk -- k index of cell
  // cfl -- cfl number

  // use nondimensional time
  dt_(ii, jj, kk) = cfl * (vol_(ii + numGhosts_, jj + numGhosts_,
                                kk + numGhosts_) /
                           avgWaveSpeed_(ii, jj, kk));
}

/* Member function to calculate the time step for all cells in the procBlock. If
the time step is user specified assign that time step (after
nondimensionalization) to dt variable. If time step is to be determined using CFL
number, call function to do so.
*/
void procBlock::CalcBlockTimeStep(const input &inputVars, const double &aRef) {
  // inputVars -- all input variables
  // aRef -- reference speed of sound (used for time non dimensionalization)

  // loop over all physical cells - no ghost cells for dt variable
  for (int kk = 0; kk < (*this).NumK(); kk++) {
    for (int jj = 0; jj < (*this).NumJ(); jj++) {
      for (int ii = 0; ii < (*this).NumI(); ii++) {
        // dt specified, use global time stepping
        if (inputVars.Dt() > 0.0) {
          // nondimensional time
          dt_(ii, jj, kk) = inputVars.Dt() * aRef / inputVars.LRef();

        // cfl specified, use local time stepping
        } else if (inputVars.CFL() > 0.0) {
          (*this).CalcCellDt(ii, jj, kk, inputVars.CFL());
        } else {
          cerr << "ERROR: Neither dt or cfl was specified!" << endl;
          exit(0);
        }
      }
    }
  }
}

/* Member function to update the procBlock to advance to a new time step. For
explicit methods it calls the appropriate explicit method to update. For
implicit methods it uses the correction du and calls the implicit updater.
*/
void procBlock::UpdateBlock(const input &inputVars, const int &impFlag,
                            const idealGas &eos, const double &aRef,
                            const multiArray3d<genArray> &du, genArray &l2,
                            resid &linf) {
  // inputVars -- all input variables
  // impFlag -- flag to determine if simulation is to be solved via explicit or
  // implicit time stepping
  // eos -- equation of state
  // aRef -- reference speed of sound (for nondimensionalization)
  // bb
  // du -- updates to conservative variables (only used in implicit solver)
  // l2 -- l-2 norm of residual
  // linf -- l-infinity norm of residual

  // if not runge-kutta 4 step method for time integration
  if (inputVars.TimeIntegration() != "rk4") {
    // loop over all physical cells
    for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
             (*this).NumK(); kk.g++, kk.p++) {
      for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
               (*this).NumJ(); jj.g++, jj.p++) {
        for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
                 (*this).NumI(); ii.g++, ii.p++) {
          // explicit euler time integration
          if (inputVars.TimeIntegration() == "explicitEuler") {
            (*this).ExplicitEulerTimeAdvance(eos, ii.g, jj.g, kk.g,
                                             ii.p, jj.p, kk.p);
          } else if (impFlag) {  // if implicit use update (du)
            (*this).ImplicitTimeAdvance(du(ii.p, jj.p, kk.p), eos,
                                        ii.g, jj.g, kk.g);
          }

          // accumulate l2 norm of residual
          l2 = l2 + residual_(ii.p, jj.p, kk.p) * residual_(ii.p, jj.p, kk.p);

          // if any residual is larger than previous residual, a new linf
          // residual is found
          for (int ll = 0; ll < NUMVARS; ll++) {
            if ((*this).Residual(ii.p, jj.p, kk.p, ll) > linf.Linf()) {
              linf.UpdateMax((*this).Residual(ii.p, jj.p, kk.p, ll),
                             parentBlock_, ii.p, jj.p, kk.p, ll + 1);
            }
          }
        }
      }
    }
  // using min storage rk4 method
  } else if (inputVars.TimeIntegration() == "rk4") {
    // save state and local time step at time n
    multiArray3d<primVars> stateN = state_;
    multiArray3d<double> dtN = dt_;

    // loop over rk stages
    for (int rr = 0; rr < 4; rr++) {
      // loop over all physical cells
      for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
               (*this).NumK(); kk.g++, kk.p++) {
        for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
                 (*this).NumJ(); jj.g++, jj.p++) {
          for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
                   (*this).NumI(); ii.g++, ii.p++) {
            // advance 1 RK stage
            (*this).RK4TimeAdvance(stateN(ii.g, jj.g, kk.g), eos,
                                   dtN(ii.p, jj.p, kk.p), ii.g, jj.g, kk.g,
                                   ii.p, jj.p, kk.p, rr);

            // at last stage
            // accumulate l2 norm of residual
            if (rr == 3) {
              l2 = l2 + residual_(ii.p, jj.p, kk.p) *
                  residual_(ii.p, jj.p, kk.p);

              for (int ll = 0; ll < NUMVARS; ll++) {
                if ((*this).Residual(ii.p, jj.p, kk.p, ll) > linf.Linf()) {
                  linf.UpdateMax((*this).Residual(ii.p, jj.p, kk.p, ll),
                                 parentBlock_, ii.p, jj.p, kk.p, ll + 1);
                }
              }
            }
          }
        }
      }
      // for multistage RK4 method, calculate fluxes and residuals again
      if (rr < 3) {  // no need to calculate fluxes after final RK interation
        // UPDATE NEEDED -- have to calculate grads, visc fluxes, source terms
        // again. Maybe BCs as well
        (*this).CalcInvFluxI(eos, inputVars);
        (*this).CalcInvFluxJ(eos, inputVars);
        (*this).CalcInvFluxK(eos, inputVars);
        (*this).CalcBlockTimeStep(inputVars, aRef);
      }
    }
  } else {
    cerr << "ERROR: Time integration scheme " << inputVars.TimeIntegration()
         << " is not recognized!" << endl;
  }
}

/* Member function to advance the state_ vector to time n+1 using explicit Euler
method. The following equation is used:

 Un+1 = Un - dt_/V * R

Un is the conserved variables at time n, Un+1 is the conserved variables at time
n+1, dt_ is the cell's time step, V is the cell's volume, and R is the cell's
residual.
 */
void procBlock::ExplicitEulerTimeAdvance(const idealGas &eqnState,
                                         const int &ig, const int &jg,
                                         const int &kg, const int &ip,
                                         const int &jp, const int &kp) {
  // eqnState -- equation of state
  // ig -- i-location of cell (including ghost cells)
  // jg -- j-location of cell (including ghost cells)
  // kg -- k-location of cell (including ghost cells)
  // ip -- i-location of cell (without ghost cells)
  // jp -- j-location of cell (without ghost cells)
  // kp -- k-location of cell (without ghost cells)

  // Get conserved variables for current state (time n)
  genArray consVars = state_(ig, jg, kg).ConsVars(eqnState);
  // calculate updated conserved variables
  consVars = consVars - dt_(ip, jp, kp) / vol_(ig, jg, kg) *
      residual_(ip, jp, kp);

  // calculate updated primative variables and update state
  state_[ig, jg, kg] = primVars(consVars, false, eqnState);
}

// member function to advance the state_ vector to time n+1 (for implicit
// methods)
void procBlock::ImplicitTimeAdvance(const genArray &du,
                                    const idealGas &eqnState, const int &ii,
                                    const int &jj, const int &kk) {
  // du -- update for a specific cell (to move from time n to n+1)
  // eqnState -- equation of state
  // ii -- i-location of cell (with ghosts)
  // jj -- j-location of cell (with ghosts)
  // kk -- k-location of cell (with ghosts)

  // calculate updated state (primative variables)
  state_(ii, jj, kk) = state_(ii, jj, kk).UpdateWithConsVars(eqnState, du);
}

/*member function to advance the state_ vector to time n+1 using 4th order
(minimum storage) Runge-Kutta method (2nd order accurate)

 Un+1 = Un - dt/V * alpha * R

Un is the conserved variables at time n, Un+1 is the conserved variables at time
n+1, dt_ is the cell's time step, V is the cell's volume, alpha is the runge-kutta
coefficient, and R is the cell's residual.
 */
void procBlock::RK4TimeAdvance(const primVars &currState,
                               const idealGas &eqnState, const int &ig,
                               const int &jg, const int &kg, const int &ip,
                               const int &jp, const int &kp, const int &rk) {
  // currState -- current state_ (including steps within RK4) (primative)
  // eqnState -- equation of state
  // ig -- i-location of cell (including ghost cells)
  // jg -- j-location of cell (including ghost cells)
  // kg -- k-location of cell (including ghost cells)
  // ip -- i-location of cell (without ghost cells)
  // jp -- j-location of cell (without ghost cells)
  // kp -- k-location of cell (without ghost cells)
  // rk -- runge-kutta step number

  // runge-kutta step coefficients (low storage 4 step)
  double alpha[4] = {0.25, 1.0 / 3.0, 0.5, 1.0};

  // update conserved variables
  genArray consVars = currState.ConsVars(eqnState) -
      dt_(ip, jp, kp) / vol_(ig, jg, kg) * alpha[rk] * residual_(ip, jp, kp);

  // calculate updated primative variables
  state_(ig, jg, kg) = primVars tempState(consVars, false, eqnState);
}

// member function to reset the residual and wave speed back to zero after an
// iteration. This is done because the residual and wave
// speed are accumulated over many function calls.
void procBlock::ResetResidWS() {
  // loop over all physical cells - no ghost cells in residual variable
  for (int kk = 0; kk < (*this).NumK(); kk++) {
    for (int jj = 0; jj < (*this).NumJ(); jj++) {
      for (int ii = 0; ii < (*this).NumI(); ii++) {
        // reset residual
        residual_(ii, jj, kk) = genArray(0.0);

        // reset wave speed
        avgWaveSpeed_(ii, jj, kk) = 0.0;
      }
    }
  }
}

/* Member function to add the cell volume divided by the cell time step to the
main diagonal of the time m minus time n term.

dU/dt = V/t * [ ((1 + zeta) * FD - zeta * BD) / ((1 + theta) * FD )] * Un = -Rn

The above equation shows the governing equations written in the Beam & Warming
format for time integration. U is the vector of conserved variables
where n represents the time step. Theta and zeta are Beam & Warming parameters,
t is the time step, V is the cell volume, and R is the residual.
FD and BD are the forward and backward difference operators respectively. These
opererators operate in the time domain. For example FD(U) =
Un+1 - Un and BD(U) = Un - Un-1. Solving the above equation for FD(Qn) we get
the following:

FD(Un) = (-t * Rn - t * theta * FD(Rn) + zeta * V * FD(Un-1)) / ((1 + zeta) * V)

FD(Rn) requires us to know the residual at time n+1, but this is unknown. To
bypass this difficulty we linearize the residual_ using a Taylor series
expansion about time n. Rn+1 = Rn + J*FD(Un) where J is the flux jacobian dR/dU.
Rearranging the above equation we get the following:

[J + (1+zeta)*V/(t*theta)] * FD(Un) = -Rn/theta + zeta*V/(t*theta) * FD(Un-1)

The above equation shows that the time m minus time n term (FD(Un)) requires a
(1+zeta)V/(t*theta) term multiplied by it. That is the purpose of this
function.
*/
multiArray3d<genArray> procBlock::AddVolTime(const multiArray3d<genArray> &m,
                                             const multiArray3d<genArray> &n,
                                             const double &theta,
                                             const double &zeta) const {
  // m -- solution for block at time m
  // n -- solution for block at time n
  // theta -- Beam & Warming coefficient theta for time integration
  // zeta -- Beam & Warming coefficient zeta for time integration

  // initialize a vector to hold the returned values
  multiArray3d<genArray> mMinusN(m.NumI(), m.NumJ(), m.NumK());

  // loop over all physical cells
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
           (*this).NumK(); kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
             (*this)NumJ(); jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
               (*this).NumI(); ii.g++, ii.p++) {
        double I = (vol_(ii.g, jj.g, kk.g) * (1.0 + zeta)) /
            (dt_(ii.p, jj.p, kk.p) * theta);
        mMinusN(ii.p, jj.p, kk.p) =
            I * (m(ii.p, jj.p, kk.p) - n(ii.p, jj.p, kk.p));
      }
    }
  }
  return mMinusN;
}

/* Member function to calculate the delta n-1 term for the implicit bdf2 solver.

dU/dt_ = V/t * [ ((1 + zeta) * FD - zeta * BD) / ((1 + theta) * FD )] * Un = -Rn

The above equation shows the governing equations written in the Beam & Warming
format for time integration. U is the vector of conserved variables
where n represents the time step. Theta and zeta are Beam & Warming parameters,
t is the time step, V is the cell volume, and R is the residual.
FD and BD are the forward and backward difference operators respectively. These
opererators operate in the time domain. For example FD(U) =
Un+1 - Un and BD(U) = Un - Un-1. Solving the above equation for FD(Qn) we get
the following:

FD(Un) = (-t * Rn - t * theta * FD(Rn) + zeta * V * FD(Un-1)) / ((1 + zeta) * V)

FD(Rn) requires us to know the residual at time n+1, but this is unknown. To
bypass this difficulty we linearize the residual using a Taylor series
expansion about time n. Rn+1 = Rn + J*FD(Un) where J is the flux jacobian dR/dU.
Rearranging the above equation we get the following:

[J + (1+zeta)*V/(t*theta)] * FD(Un) = -Rn/theta + zeta*V/(t*theta) * FD(Un-1)

The above equation shows that the time n minus time n-1 term (FD(Un-1)) requires
a zeta*V/(t*theta) term multiplied by it. That is the purpose of this
function.
*/
void procBlock::DeltaNMinusOne(multiArray3d<genArray> &solDeltaNm1,
                               const multiArray3d<genArray> &solTimeN,
                               const idealGas &eqnState, const double &theta,
                               const double &zeta) {
  // solDeltaNm1 -- The solution at time n minus the solution at time n-1. (Un -
  // Un-1) (output)
  // solTimeN -- The solution at time n
  // eqnState -- equation of state_
  // theta -- Beam & Warming coefficient theta for time integration
  // zeta -- Beam & Warming coefficient zeta for time integration

  // loop over physical cells
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
           (*this).NumK(); kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
             (*this)NumJ(); jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
               (*this).NumI(); ii.g++, ii.p++) {
        double coeff = (vol_(ii.g, jj.g, kk.g) * zeta) /
            (dt_(ii.p, jj.p, kk.p) * theta);
        solDeltaNm1(ii.p, jj.p, kk.p) = coeff *
            (state_(ii.g, jj.g, kk.g).ConsVars(eqnState) -
             solTimeN(ii.p, jj.p, kk.));
      }
    }
  }
}

// Member function to return a copy of the conserved variables. This is useful
// for "saving" the solution at a specific time.
// It is used in the implicit solver.
multiArray3d<genArray> procBlock::GetCopyConsVars(
    const idealGas &eqnState) const {
  // initialize array to proper size (no ghost cells)
  multiArray3d<genArray> consVars((*this).NumI(), (*this).NumJ(),
                                  (*this).NumK());

  // loop over physical cells
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
           consVars.NumK(); kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
             consVars.NumJ(); jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
               consVars.NumI(); ii.g++, ii.p++) {
        // convert state to conservative variables
        consVars(ii.p, jj.p, kk.p) =
            state_(ii.g, jj.g, kk.g).ConsVars(eqnState);
      }
    }
  }
  return consVars;
}

/* Member function to calculate update to solution implicitly using Lower-Upper
Symmetric Gauss Seidel (LUSGS) method.

Un+1 = Un - t/V * Rn+1

The above equation shows a simple first order implicit method to calculate the
solution at the next time step (n+1). The equation shows that this method
requires the residual (R) at time n+1 which is unknown. In the equation, t is
the time step, and V is the volume. Since the residual at n+1 in unknown, it
must be linearized about time n as shown below.

Rn+1 = Rn + dRn/dUn * FD(Un)

In the above equation FD is the forward difference operator in time (FD(Un) =
Un+1 - Un). The derivative of the residual term can be further simplified as
below.

Rn = (sum over all faces) Fni

In the above equation n refers to the time level and i refers to the face index.
The summation over all faces operator will now be abbreviated as (SF).
Substituting the second and third equations into the first and rearranging we
get the following.

[d(SF)Fni/dUnj + V/t] * FD(Un) = -Rn

In the above equation the index j refers to all of the cells in the stencil
going into the calculation of the flux at face i. The above equation
can be simplified to A*x=b. The matrix A is an MxM block matrix where M is the
number of cells in the block. Each block is an LxL block where L
is the number of equations being solved for. The sparsity of A depends on the
stencil used in flux calculation. In 3D for a first order simulation
the matrix A is block pentadiagonal. For a second order approximation it would
have 13 diagonals. This increases the storage requirements so in
practice a first order approximation is used. The order of accuracy is
determined by the residual calculation. The accuracy of the matrix A helps
with convergence. A poorer approximation will eventually get the correct answer
with enough iteration (defect correction). Fully implicit methods
calculate and store the flux jacobians needed to populate the matrix A. The
LUSGS method does not do this and instead calculates an approximate
flux jacobian "on-the-fly" so there is no need for storage. Because an
approximate flux jacobian is being used, there is no need for a highly
accurate linear solver. Therefore the Symmetric Gauss-Seidel (SGS) method is a
good candidate, and only one iteration is needed. This approximate
flux jacobian along with the SGS linear solver form the basics of the LUSGS
method (Jameson & Yoon). The LUSGS method begins by factoring the
matrix A as shown below.

A = (D + L) * D^-1 * (D + U)

In the above equation D is the diagonal of A, L is the lower triangular portion
of A, and U is the upper triangular portion of A. This allows the equation A*x=b
to be solved in one SGS sweep as shown below.

Forward sweep:  (D + L) * FD(Un*) = -Rn
Backward sweep: (D + U) * FD(Un) = D * FD(Un*)

Another key component of the LUSGS scheme is to sweep along hyperplanes.
Hyperplanes are planes of i+j+k=constant within a plot3d block. The diagram
below shows a 2D example of a block reordered to sweep along hyperplanes
           ____ ____ ____ ____ ____ ____ ____ ____
          | 20 | 26 | 32 | 37 | 41 | 44 | 46 | 47 |
          |____|____|____|____|____|____|____|____|
          | 14 | 19 | 25 | 31 | 36 | 40 | 43 | 45 |
          |____|____|____|____|____|____|____|____|
          | 9  | 13 | 18 | 24 | 30 | 35 | 39 | 42 |
   A=     |____|____|____|____|____|____|____|____|
          | 5  | 8  | 12 | 17 | 23 | 29 | 34 | 38 |
          |____|____|____|____|____|____|____|____|
          | 2  | 4  | 7  | 11 | 16 | 22 | 28 | 33 |
          |____|____|____|____|____|____|____|____|
          | 0  | 1  | 3  | 6  | 10 | 15 | 21 | 27 |
          |____|____|____|____|____|____|____|____|

This is advantageous because on the forward sweep the lower matrix L can be
calculated with data from time n+1 because all of the cells contributing
to L would already have been updated. The same is true with the upper matrix U
for the backward sweep. This removes the need for any storage of the
matrix. For a given location (say A12) the matrix L would be constructed of A8
and A7, and the matrix U would be constructed of A18 and A17. This
requires the product of the flux jacobian multiplied with the update (FD(Un)) to
be calculated at these locations. For the LUSGS method the flux
jacobians are approximated as follows:

A * S = 0.5 * (Ac * S + K * I)

A is the flux jacobian, S is the face area, Ac is the convective flux jacobian
(dF/dU), K is the spectral radius multiplied by a factor, and I is
the identity matrix. The addition of the spectrial radius improves diagonal
dominance which improves stability at the cost of convergence. When the
factor multiplied by K is 1, the method is SGS, when it is < 1, it is successive
overrelaxation. Reducing the factor improves convergence but hurts
stability. The product of the approximate flux jacobian with the update is
calculated as shown below.

A * S * FD(Unj) = 0.5 * (Ac * FD(Unj) * S + K * I * FD(Unj)) = 0.5 * ( dFi/dUnj
* FD(Unj) * S + K * I * FD(Unj)) = 0.5 * (dFi * S + K * I * FD(Unj))

The above equation shows that all that is needed to calculate the RHS of A*x=b
is the update of the convective flux, and the update to the convervative
variabes (FD(Unj)) which is known due to sweeping along hyperplanes.

For viscous simulations, the viscous contribution to the spectral radius K is
used, and everything else remains the same.
 */
double procBlock::LUSGS(const vector<vector3d<int> > &reorder,
                        multiArray3d<genArray> &x,
                        const multiArray3d<genArray> &solTimeMmN,
                        const multiArray3d<genArray> &solDeltaNm1,
                        const idealGas &eqnState, const input &inp,
                        const sutherland &suth, const turbModel *turb) const {
  // reorder -- order of cells to visit (this should be ordered in hyperplanes)
  // x -- correction - added to solution at time n to get to time n+1 (assumed
  // to be zero to start)
  // solTimeMmn -- solution at time m minus n
  // solDeltaNm1 -- solution at time n minus solution at time n-1
  // eqnState -- equation of state_
  // inp -- all input variables
  // suth -- method to get temperature varying viscosity (Sutherland's law)
  // turb -- turbulence model

  double thetaInv = 1.0 / inp.Theta();

  // initialize genArray to zero
  genArray initial(0.0);

  // initialize L and U matrices
  multiArray3d<genArray> U((*this).NumI(), (*this).NumJ(), (*this).NumK(),
                           initial);
  multiArray3d<genArray> L((*this).NumI(), (*this).NumJ(), (*this).NumK(),
                           initial);

  //--------------------------------------------------------------------
  // forward sweep over all physical cells
  for (int ii = 0; ii < (*this).NumCells(); ii++) {
    // indices for variables without ghost cells
    int ip = reorder[ii].X();
    int jp = reorder[ii].Y();
    int kp = reorder[ii].Z();
    // indices for variables with ghost cells
    int ig = reorder[ii].X() - numGhosts_;
    int jg = reorder[ii].Y() - numGhosts_;
    int kg = reorder[ii].Z() - numGhosts_;

    // if i lower diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip - 1, jp, kp, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaI(ig - 1, jg, kg),
                             (*this).FAreaI(ig, jg, kg),
                             state_(ig - 1, jg, kg).
                             UpdateWithConsVars(eqnState,
                                                x(ip - 1, jp, kp)), eqnState);

      // if viscous add viscous contribution to spectral radius
      if (inp.EquationSet() != "euler") {
        specRad +=
            ViscCellSpectralRadius((*this).FAreaI(ig - 1, jg, kg),
                                   (*this).FAreaI(ig, jg, kg),
                                   state_(ig - 1, jg, kg).
                                   UpdateWithConsVars(eqnState,
                                                      x(ip - 1, jp, kp)),
                                   eqnState, suth, vol_(ig - 1, jg, kg),
                                   turb->EddyViscNoLim(state_(ig - 1, jg, kg)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig - 1, jg, kg), eqnState, (*this).FAreaUnitI(ig, jg, kg),
          x(ip - 1, jp, kp));

      // update L matrix
      L(ip, jp, kp) = L(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagI(ig, jg, kg) * fluxChange + inp.MatrixRelaxation() *
           specRad * x(ip - 1, jp, kp));
    }
    // if j lower diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip, jp - 1, kp, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaJ(ig, jg - 1, kg),
                             (*this).FAreaJ(ig, jg, kg),
                             state_(ig, jg - 1, kg).
                             UpdateWithConsVars(eqnState,
                                                x(ip, jp - 1, kp)), eqnState);

      // if viscous add viscous contribution to spectral radius
      if (inp.EquationSet() != "euler") {
        specRad +=
            ViscCellSpectralRadius(
                (*this).FAreaJ(ig, jg - 1, kg), (*this).FAreaJ(ig, jg, kg),
                state_(ig, jg - 1, kg).
                UpdateWithConsVars(eqnState, x(ip, jp - 1, kp)),
                eqnState, suth, vol_(ig, jg - 1, kg),
                turb->EddyViscNoLim(state_(ig, jg - 1, kg)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig, jg - 1, kg), eqnState, (*this).FAreaUnitJ(ig, jg, kg),
          x(ip, jp - 1, kp));

      // update L matrix
      L(ip, jp, kp) = L(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagJ(ig, jg, kg) * fluxChange + inp.MatrixRelaxation() *
           specRad * x(ip, jp - 1, kp));
    }
    // if k lower diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip, jp, kp - 1, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaK(ig, jg, kg - 1),
                             (*this).FAreaK(ig, jg, kg),
                             state_(ig, jg, kg - 1).
                             UpdateWithConsVars(eqnState,
                                                x(ip, jp, kp - 1)), eqnState);

      // if viscous add viscous contribution to spectral radius
      if (inp.EquationSet() != "euler") {
        specRad +=
            ViscCellSpectralRadius(
                (*this).FAreaK(ig, jg, kg - 1), (*this).FAreaK(ig, jg, kg),
                state_(ig, jg, kg - 1).
                UpdateWithConsVars(eqnState, x(ip, jp, kp - 1)),
                eqnState, suth, vol_(ig, jg, kg - 1),
                turb->EddyViscNoLim(state_(ig, jg, kg - 1)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig, jg, kg - 1), eqnState, (*this).FAreaUnitK(ig, jg, kg),
          x(ip, jp, kp - 1));

      // update L matrix
      L(ip, jp, kp) = L(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagK(ig, jg, kg) * fluxChange + inp.MatrixRelaxation() *
           specRad * x(ip, jp, kp - 1));
    }

    // add dual time stepping contribution to main diagonal
    double diagTimeVol = (vol_(ig, jg, kg) * (1.0 + inp.Zeta())) /
                         (dt_(ip, jp, kp) * inp.Theta());
    if (inp.DualTimeCFL() > 0.0) {  // use dual time stepping
      double tau = avgWaveSpeed_(ip, jp, kp) /
                   inp.DualTimeCFL();  // equal to volume / tau
      diagTimeVol += tau;
    }

    double AiiInv = 1.0 / ((avgWaveSpeed_(ip, jp, kp) + diagTimeVol) *
                           inp.MatrixRelaxation());

    // calculate intermediate update
    // normal at lower boundaries needs to be reversed, so add instead
    // of subtract L
    x(ip, jp, kp) = AiiInv * (-1.0 * thetaInv * residual_(ip, jp, kp) -
                              solDeltaNm1(ip, jp, kp) - solTimeMmN(ip, jp, kp) +
                              L(ip, jp, kp));
  }  // end forward sweep

  //----------------------------------------------------------------------
  // backward sweep over all physical cells
  for (int ii = (*this).NumCells() - 1; ii >= 0; ii--) {
    // indices for variables without ghost cells
    int ip = reorder[ii].X();
    int jp = reorder[ii].Y();
    int kp = reorder[ii].Z();
    // indices for variables with ghost cells
    int ig = reorder[ii].X() - numGhosts_;
    int jg = reorder[ii].Y() - numGhosts_;
    int kg = reorder[ii].Z() - numGhosts_;

    // if i upper diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip + 1, jp, kp, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaI(ig + 2, jg, kg),
                             (*this).FAreaI(ig + 1, jg, kg),
                             state_(ig + 1, jg, kg).
                             UpdateWithConsVars(eqnState,
                                                x(ip + 1, jp, kp)), eqnState);

      if (inp.EquationSet() != "euler") {  // viscous
        specRad +=
            ViscCellSpectralRadius((*this).FAreaI(ig + 2, jg, kg),
                                   (*this).FAreaI(ig + 1, jg, kg),
                                   state_(ig + 1, jg, kg).
                                   UpdateWithConsVars(eqnState,
                                                      x(ip + 1, jp, kp)),
                                   eqnState, suth, vol_(ig + 1, jg, kg),
                                   turb->EddyViscNoLim(state_(ig + 1, jg, kg)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig + 1, jg, kg), eqnState, (*this).FAreaUnitI(ig + 1, jg, kg),
          x(ip + 1, jp, kp));

      // update U matrix
      U(ip, jp, kp) = U(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagI(ig + 1, jg, kg) * fluxChange -
           inp.MatrixRelaxation() * specRad * x(ip + 1, jp, kp));
    }
    // if j upper diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip, jp + 1, kp, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaJ(ig, jg + 2, kg),
                             (*this).FAreaJ(ig, jg + 1, kg),
                             state_(ig, jg + 1, kg).
                             UpdateWithConsVars(eqnState,
                                                x(ip, jp + 1, kp)), eqnState);

      if (inp.EquationSet() != "euler") {  // viscous
        specRad +=
            ViscCellSpectralRadius((*this).FAreaJ(ig, jg + 2, kg),
                                   (*this).FAreaJ(ig, jg + 1, kg),
                                   state_(ig, jg + 1, kg).
                                   UpdateWithConsVars(eqnState,
                                                      x(ip, jp + 1, kp)),
                                   eqnState, suth, vol_(ig, jg + 1, kg),
                                   turb->EddyViscNoLim(state_(ig, jg + 1, kg)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig, jg + 1, kg), eqnState, (*this).FAreaUnitJ(ig, jg + 1, kg),
          x(ip, jp + 1, kp));

      // update U matrix
      U(ip, jp, kp) = U(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagJ(ig, jg + 1, kg) * fluxChange -
           inp.MatrixRelaxation() * specRad * x(ip, jp + 1, kp));
    }
    // if k upper diagonal cell is in physical location there is a contribution
    // from it
    if ((*this).IsPhysical(ip, jp, kp + 1, false)) {
      // at given face location, call function to calculate spectral radius,
      // since values are constant throughout cell, cell center values are used
      double specRad =
          CellSpectralRadius((*this).FAreaK(ig, jg, kg + 2),
                             (*this).FAreaK(ig, jg, kg + 1),
                             state_(ig, jg, kg + 1).
                             UpdateWithConsVars(eqnState,
                                                x(ip, jp, kp + 1)), eqnState);

      if (inp.EquationSet() != "euler") {  // viscous
        specRad +=
            ViscCellSpectralRadius((*this).FAreaK(ig, jg, kg + 2),
                                   (*this).FAreaK(ig, jg, kg + 1),
                                   state_(ig, jg, kg + 1).
                                   UpdateWithConsVars(eqnState,
                                                      x(ip, jp, kp + 1)),
                                   eqnState, suth, vol_(ig, jg, kg + 1),
                                   turb->EddyViscNoLim(state_(ig, jg, kg + 1)));
      }

      // at given face location, call function to calculate convective flux
      // change
      genArray fluxChange = ConvectiveFluxUpdate(
          state_(ig, jg, kg + 1), eqnState, (*this).FAreaUnitK(ig, jg, kg + 1),
          x(ip, jp, kp + 1));

      // update U matrix
      U(ip, jp, kp) = U(ip, jp, kp) + 0.5 *
          ((*this).FAreaMagK(ig, jg, kg + 1) * fluxChange -
           inp.MatrixRelaxation() * specRad * x(ip, jp, kp + 1));
    }

    // add dual time stepping contribution to main diagonal
    double diagTimeVol = (vol_(ig, jg, kg) * (1.0 + inp.Zeta())) /
                         (dt_(ip, jp, kp) * inp.Theta());
    if (inp.DualTimeCFL() > 0.0) {  // use dual time stepping
      double tau = avgWaveSpeed_(ip, jp, kp) /
                   inp.DualTimeCFL();  // equal to volume / tau
      diagTimeVol += tau;
    }

    double AiiInv = 1.0 / ((avgWaveSpeed_(ip, jp, kp) + diagTimeVol) *
                           inp.MatrixRelaxation());

    // calculate update
    x(ip, jp, kp) = x(ip, jp, kp) - AiiInv * U(ip, jp, kp);
  }  // end backward sweep

  //-------------------------------------------------------------------
  // calculate residual
  // initialize LUSGS residual vector
  genArray l2Resid(0.0);

  // loop over physical cells, can use natural order for speed
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
           (*this).NumK(); kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
             (*this).NumJ(); jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
               (*this).NumI(); ii.g++, ii.p++) {
        // calclate dual time stepping contribution
        double diagTimeVol = (vol_(ii.g, jj.g, kk.g) * (1.0 + inp.Zeta())) /
            (dt_(ii.p, jj.p, kk.p) * inp.Theta());
        if (inp.DualTimeCFL() > 0.0) {  // use dual time stepping
          double tau = avgWaveSpeed_(ii.p, jj.p, kk.p) /
              inp.DualTimeCFL();  // equal to volume / tau
          diagTimeVol += tau;
        }

        double Aii = (avgWaveSpeed_(ii.p, jj.p, kk.p) + diagTimeVol) *
            inp.MatrixRelaxation();

        // normal at lower boundaries needs to be reversed, so add instead of
        // subtract L
        genArray resid = -1.0 * thetaInv * residual_(ii.p, jj.p, kk.p) +
            solDeltaNm1(ii.p, jj.p, kk.p) + solTimeMmN(ii.p, jj.p, kk.p) - Aii *
            x(ii.p, jj.p, kk.p) + L(ii.p, jj.p, kk.p) - U(ii.p, jj.p, kk.p);
        l2Resid = l2Resid + resid * resid;
      }
    }
  }
  return l2Resid.Sum();
}

/*Function to return the inviscid spectral radius for one direction (i, j, or k)
given a cell state, equation of state, and 2 face area vectors

L = 0.5 * (A1 + A2) * (|Vn| + SoS)

In the above equation L is the spectral radius in either the i, j, or k
direction. A1 and A2 are the two face areas in that direction. Vn is the
cell velocity normal to that direction. SoS is the speed of sound at the cell
 */
double CellSpectralRadius(const unitVec3dMag<double> &fAreaL,
                          const unitVec3dMag<double> &fAreaR,
                          const primVars &state, const idealGas &eqnState) {
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // state -- state at cell center (primative)
  // eqnState -- equation of state

  // normalize face areas
  vector3d<double> normAvg = (0.5 * (fAreaL.UnitVector() +
                                     fAreaR.UnitVector())).Normalize();
  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());  // average area magnitude

  // return spectral radius
  return (fabs(state.Velocity().DotProd(normAvg)) + state.SoS(eqnState)) *
         fMag;
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
double ViscCellSpectralRadius(const unitVec3dMag<double> &fAreaL,
                              const unitVec3dMag<double> &fAreaR,
                              const primVars &state, const idealGas &eqnState,
                              const sutherland &suth, const double &vol,
                              const double &eddyVisc) {
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // state -- state at cell center (primative)
  // eqnState -- equation of state
  // suth -- method to the temperature varying visosity and Prandtl number
  // (Sutherland's law)
  // vol -- cell volume

  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());  // average area magnitude
  double maxTerm =
      max(4.0 / (3.0 * state.Rho()), eqnState.Gamma() / state.Rho());
  double mu = (suth.Viscosity(state.Temperature(eqnState)) +
               eddyVisc) * suth.NondimScaling();  // viscosity at cell center
  double viscTerm = mu / eqnState.Prandtl();

  // return viscous spectral radius
  return maxTerm * viscTerm * fMag * fMag / vol;
}

// function to reconstruct cell variables to the face using central
// differences
template <class T>
T FaceReconCentral(const T &velU, const T &velD, const vector3d<double> &pU,
                   const vector3d<double> &pD, const vector3d<double> &pF) {
  // velU -- velocity at the cell center of the upwind cell
  // velD -- velocity at the cell center of the downwind cell
  // pU -- position of the cell center of the upwind cell
  // pD -- position of the cell center of the downwind cell
  // pF -- position of the face center of the face on which the reconstruction
  // is happening

  // distance from cell center to cell center
  double cen2cen = pU.Distance(pD);
  // distance from upwind cell center to cell face
  double up2face = pU.Distance(pF);

  // reconstruct with central difference
  return velD * (up2face / cen2cen) + velU * (1.0 - (up2face / cen2cen));
}

/* Function to pad a vector with a specified number of ghost cells
           ___ ___ ___ ___ ___ ___ ___ ___
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|

In the above diagram, the cells marked with an "X" represent physical cells. The
entire diagram represents the block (in 2D) padded with 2 layers of ghost cells.
The cells marked with "G" are regualar ghost cells. The cells marked with "E" are
ghost cells located along one of the 12 edges that form a plot3d block. In 3D
there are also "corner" cells located at the 8 corners that form the plot3d block.
These cells are not used though. There is a place in the vector for them to make
accessing the padded vector of cells the same as for a plot3d block without ghost
cells.
*/
template <class T>
multiArray3d<T> PadWithGhosts(const multiArray3d<T> &var,
                              const int &numGhosts) {
  // var -- vector of variables to pad (no ghost cells included)
  // numGhosts -- number of layers of ghost cells to pad var with

  // initialize added array
  multiArray3d<T> padBlk(var.NumI() + 2 * numGhosts, var.NumJ() + 2 * numGhosts,
                         var.NumK() + 2 * numGhosts);

  padBlk.Insert(numGhosts, varNumI() + numGhosts - 1, numGhosts, var.NumJ() +
                numGhosts - 1, numGhosts, var.NumK() + numGhosts - 1, var);
  return padBlk;
}

/* Function to calculate the velocity gradient at the cell center using the
Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a scalar is calculated using the
Green-Gauss method. U is a scalar. Aij is the area at face i (component j).
V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out
of the control volume.
 */
tensor<double> CalcVelGradGG(
    const vector3d<double> &vil, const vector3d<double> &viu,
    const vector3d<double> &vjl, const vector3d<double> &vju,
    const vector3d<double> &vkl, const vector3d<double> &vku,
    const vector3d<double> &ail, const vector3d<double> &aiu,
    const vector3d<double> &ajl, const vector3d<double> &aju,
    const vector3d<double> &akl, const vector3d<double> &aku,
    const double &vol) {
  // vil -- velocity vector at the i-lower face of the cell at which the
  // velocity gradient is being calculated
  // viu -- velocity vector at the i-upper face of the cell at which the
  // velocity gradient is being calculated
  // vjl -- velocity vector at the j-lower face of the cell at which the
  // velocity gradient is being calculated
  // vju -- velocity vector at the j-upper face of the cell at which the
  // velocity gradient is being calculated
  // vkl -- velocity vector at the k-lower face of the cell at which the
  // velocity gradient is being calculated
  // vku -- velocity vector at the k-upper face of the cell at which the
  // velocity gradient is being calculated

  // ail -- area vector at the lower i-face of the cell at which the velocity
  // gradient is being calculated
  // aiu -- area vector at the upper i-face of the cell at which the velocity
  // gradient is being calculated
  // ajl -- area vector at the lower j-face of the cell at which the velocity
  // gradient is being calculated
  // aju -- area vector at the upper j-face of the cell at which the velocity
  // gradient is being calculated
  // akl -- area vector at the lower k-face of the cell at which the velocity
  // gradient is being calculated
  // aku -- area vector at the upper k-face of the cell at which the velocity
  // gradient is being calculated

  // vol -- cell volume

  tensor<double> temp;
  double invVol = 1.0 / vol;

  // define velocity gradient tensor
  // convention is for area vector to point out of cell, so lower values are
  // negative, upper are positive
  temp.SetXX(invVol *
             (viu.X() * aiu.X() - vil.X() * ail.X() + vju.X() * aju.X() -
              vjl.X() * ajl.X() + vku.X() * aku.X() - vkl.X() * akl.X()));
  temp.SetXY(invVol *
             (viu.Y() * aiu.X() - vil.Y() * ail.X() + vju.Y() * aju.X() -
              vjl.Y() * ajl.X() + vku.Y() * aku.X() - vkl.Y() * akl.X()));
  temp.SetXZ(invVol *
             (viu.Z() * aiu.X() - vil.Z() * ail.X() + vju.Z() * aju.X() -
              vjl.Z() * ajl.X() + vku.Z() * aku.X() - vkl.Z() * akl.X()));

  temp.SetYX(invVol *
             (viu.X() * aiu.Y() - vil.X() * ail.Y() + vju.X() * aju.Y() -
              vjl.X() * ajl.Y() + vku.X() * aku.Y() - vkl.X() * akl.Y()));
  temp.SetYY(invVol *
             (viu.Y() * aiu.Y() - vil.Y() * ail.Y() + vju.Y() * aju.Y() -
              vjl.Y() * ajl.Y() + vku.Y() * aku.Y() - vkl.Y() * akl.Y()));
  temp.SetYZ(invVol *
             (viu.Z() * aiu.Y() - vil.Z() * ail.Y() + vju.Z() * aju.Y() -
              vjl.Z() * ajl.Y() + vku.Z() * aku.Y() - vkl.Z() * akl.Y()));

  temp.SetZX(invVol *
             (viu.X() * aiu.Z() - vil.X() * ail.Z() + vju.X() * aju.Z() -
              vjl.X() * ajl.Z() + vku.X() * aku.Z() - vkl.X() * akl.Z()));
  temp.SetZY(invVol *
             (viu.Y() * aiu.Z() - vil.Y() * ail.Z() + vju.Y() * aju.Z() -
              vjl.Y() * ajl.Z() + vku.Y() * aku.Z() - vkl.Y() * akl.Z()));
  temp.SetZZ(invVol *
             (viu.Z() * aiu.Z() - vil.Z() * ail.Z() + vju.Z() * aju.Z() -
              vjl.Z() * ajl.Z() + vku.Z() * aku.Z() - vkl.Z() * akl.Z()));

  return temp;
}

/* Function to calculate the gradient of a scalar at the cell center using the
Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a scalar is calculated using the
Green-Gauss method. U is a scalar. Aij is the area at face i (component j).
V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out
of the control volume.
 */
vector3d<double> CalcScalarGradGG(
    const double &til, const double &tiu, const double &tjl, const double &tju,
    const double &tkl, const double &tku, const vector3d<double> &ail,
    const vector3d<double> &aiu, const vector3d<double> &ajl,
    const vector3d<double> &aju, const vector3d<double> &akl,
    const vector3d<double> &aku, const double &vol) {
  // til -- scalar value at the lower face of the cell at which the scalar
  // gradient is being calculated
  // tiu -- scalar value at the upper face of the cell at which the scalar
  // gradient is being calculated
  // tjl -- scalar value at the lower face of the cell at which the scalar
  // gradient is being calculated
  // tju -- scalar value at the upper face of the cell at which the scalar
  // gradient is being calculated
  // tkl -- scalar value at the lower face of the cell at which the scalar
  // gradient is being calculated
  // tku -- scalar value at the upper face of the cell at which the scalar
  // gradient is being calculated

  // ail -- area vector at the lower face of the cell at which the scalar
  // gradient is being calculated
  // aiu -- area vector at the upper face of the cell at which the scalar
  // gradient is being calculated
  // ajl -- area vector at the lower face of the cell at which the scalar
  // gradient is being calculated
  // aju -- area vector at the upper face of the cell at which the scalar
  // gradient is being calculated
  // akl -- area vector at the lower face of the cell at which the scalar
  // gradient is being calculated
  // aku -- area vector at the upper face of the cell at which the scalar
  // gradient is being calculated

  // vol -- cell volume

  vector3d<double> temp;
  double invVol = 1.0 / vol;

  // define scalar gradient vector
  // convention is for area vector to point out of cell, so lower values are
  // negative, upper are positive
  temp.SetX(invVol * (tiu * aiu.X() - til * ail.X() + tju * aju.X() -
                      tjl * ajl.X() + tku * aku.X() - tkl * akl.X()));
  temp.SetY(invVol * (tiu * aiu.Y() - til * ail.Y() + tju * aju.Y() -
                      tjl * ajl.Y() + tku * aku.Y() - tkl * akl.Y()));
  temp.SetZ(invVol * (tiu * aiu.Z() - til * ail.Z() + tju * aju.Z() -
                      tjl * ajl.Z() + tku * aku.Z() - tkl * akl.Z()));

  return temp;
}

/* Function to calculate the viscous fluxes on the i-faces. All phyiscal
(non-ghost) i-faces are looped over. The left and right states are
calculated, and then the flux at the face is calculated. The flux at the
face contributes to the residual of the cells to the left and right of the
face. This contribution from the flux is added to the residuals and the wave
speed is accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Ui       -->   Ui+1   |
  |            |            |
  |____________|____________|
Ui-1/2       Ui+1/2       Ui+3/2

Using the above diagram, the flux is calculated at face Ui+1/2. Since the area
vector at the face always points from lower indices to higher indices it points
from Ui to Ui+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Ui, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Ui+1.

The spectral radius in the i-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Ui+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by
constructing an alternative control volume centered around that face as show below.
  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done.
In this example the gradient is being calculated at the face between cell Ui,j
and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the
alternative control volume. The grid cells themselves cannot be used as the
control volume for the gradients (averaging values at adjacent cells to get
gradients at the face) because this leads to odd/even decoupling. The face areas,
volumes, and states at the center_ of the faces are needed for the alternative
control volume. The left and right sides of the alternate control volume pass
through the center_ of cells Ui,j and Ui+1,j respectively. Therefore these values
are used for the face states. The left and right face areas are calculated as the
average of the face areas of the cells that they split. For example, the left face
area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j). The top and bottom sides
of the alternative control volume pass through 4 cells each. Therefore the value
of the state at the face center is the average of these four states. For example
the state at the top face is calculated as 0.25 * (Ui,j + Ui+1,j + Ui,j+1 +
Ui+1,j+1). The face areas of the top and bottom sides are calculated as the
average of the 2 face areas that each one passes through. For example, the top
face area is calculated as 0.5 * (Ai,j+1/2 + Ai+1,j+1/2). In three dimensions
each gradient calculation touches the values at 10 cells (6 shown and 4 more
in/out of the page). The stencil for the gradients of all faces in a cell
touches 15 cells. The gradient calculation with this stencil uses the "edge"
ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxI(const sutherland &suth, const idealGas &eqnState,
                              const input &inp, const gradients &grads,
                              const turbModel *turb) {
  // suth -- method to get viscosity as a function of temperature (Sutherland's
  // law)
  // eqnState -- equation of state
  // inp -- all input variables
  // grads -- class holding gradients at face for velocity, temperature, tke,
  // and omega
  // turb -- turbulence model

  // coefficient for viscous spectral radii
  double vCoeff = 1.0;

  // loop over all physical i-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaI_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaI_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.g <
               fAreaI_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // Get velocity at face
        vector3d<double> vel =
            FaceReconCentral(state_(ii.g - 1, jj.g, kk.g).Velocity(),
                             state_(ii.g, jj.g, kk.g).Velocity(),
                             center_(ii.g - 1, jj.g, kk.g),
                             center_(ii.g, jj.g, kk.g),
                             fCenterI_(ii.g, jj.g, kk.g));

        // Get state at face
        primVars state =
            FaceReconCentral(state_(ii.g - 1, jj.g, kk.g),
                             state_(ii.g, jj.g, kk.g),
                             center_(ii.g - 1, jj.g, kk.g),
                             center_(ii.g, jj.g, kk.g),
                             fCenterI_(ii.g, jj.g, kk.g));

        // Get viscosity at face
        double mu =
            FaceReconCentral(
                suth.EffectiveViscosity(state_(ii.g - 1, jj.g, kk.g).
                                        Temperature(eqnState)),
                suth.EffectiveViscosity(state_(ii.g, jj.g, kk.g).
                                        Temperature(eqnState)),
                center_(ii.g - 1, jj.g, kk.g), center_(ii.g, jj.g, kk.g),
                fCenterI_(ii.g, jj.g, kk.g));

        double eddyVisc =
            FaceReconCentral(
            turb->EddyVisc(state_(ii.g - 1, jj.g, kk.g),
                           grads.VelGradI(ii.p, jj.p, kk.p), suth),
            turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                           grads.VelGradI(ii.p, jj.p, kk.p), suth),
            center_(ii.g - 1, jj.g, kk.g), center_(ii.g, jj.g, kk.g),
            fCenterI_(ii.g, jj.g, kk.g));
        eddyVisc *= suth.NondimScaling();  // effective viscosity

        // calculate viscous flux
        vector3d<double> tkeGrad, omegaGrad;
        if (inp.IsTurbulent()) {
          tkeGrad = grads.TkeGradI(ii.p, jj.p, kk.p);
          omegaGrad = grads.OmegaGradI(ii.p, jj.p, kk.p);
        }
        viscousFlux tempViscFlux(grads.VelGradI(ii.p, jj.p, kk.p), vel, mu,
                                 eddyVisc, suth, eqnState,
                                 grads.TempGradI(ii.p, jj.p, kk.p),
                                 (*this).FAreaUnitI(ii.g, jj.g, kk.g), tkeGrad,
                                 omegaGrad, turb, state);

        // area vector points from left to right, so add to left cell, subtract
        // from right cell but viscous fluxes are subtracted from inviscid
        // fluxes, so sign is reversed
        // at left boundary there is no left cell to add to
        if (ii.g > numGhosts_) {
          (*this).AddToResidual(-1.0 * tempViscFlux *
                                (*this).FAreaMagI(ii.g, jj.g, kk.g),
                                ii.p - 1, jj.p, kk.p);
        }
        // at right boundary ther eis to right cell to add to
        if (ii.g < fAreaI_.NumI() - 2 *numGhosts_ - 1) {
          (*this).AddToResidual(tempViscFlux *
                                (*this).FAreaMagI(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);

          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_[ii.p, jj.p, kk.p] += vCoeff *
              ViscCellSpectralRadius(
                  (*this).FAreaI(ii.g, jj.g, kk.g),
                  (*this).FAreaI(ii.g + 1, jj.g, kk.g),
                  state_(ii.g, jj.g, kk.g), eqnState, suth,
                  vol_(ii.g, jj.g, kk.g),
                  turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                                 grads.VelGradI(ii.p, jj.p, kk.p), suth));
        }
      }
    }
  }
}

/* Function to calculate the viscous fluxes on the j-faces. All phyiscal
(non-ghost) j-faces are looped over. The left and right states are calculated,
and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This
contribution from the flux is added to the residuals and the wave speed is
accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Uj       -->   Uj+1   |
  |            |            |
  |____________|____________|
Uj-1/2       Uj+1/2       Uj+3/2

Using the above diagram, the flux is calculated at face Uj+1/2. Since the area
vector at the face always points from lower indices to higher indices it points
from Uj to Uj+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Uj, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Uj+1.

The spectral radius in the j-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Uj+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by
constructing an alternative control volume centered around that face as shown
below.

  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done.
In this example the gradient is being calculated at the face between cell Ui,j
and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the
alternative control volume. The grid cells themselves cannot be used as the
control volume for the gradients (averaging values at adjacent cells to get
gradients at the face) because this leads to odd/even decoupling. The face
areas, volumes, and states at the center_ of the faces are needed for the
alternative control volume. The left and right sides of the alternate control
volume pass through the center of cells Ui,j and Ui+1,j respectively. Therefore
these values are used for the face states. The left and right face areas are
calculated as the average of the face areas of the cells that they split. For
example, the left face area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j).
The top and bottom sides of the alternative control volume pass through 4 cells
each. Therefore the value of the state at the face center is the average of
these four states. For example the state at the top face is calculated as 0.25 *
(Ui,j + Ui+1,j + Ui,j+1, Ui+1,j+1). The face areas of the top and bottom sides
are calculated as the average of the 2 face areas that each one passes through.
For example, the top face area is calculated as 0.5 * (Ai,j+1/2 + Ai+1,j+1/2).
In three dimensions each gradient calculation touches the values at 10 cells
(6 shown and 4 more in/out of the page). The stencil for the gradients of all
faces in a cell touches 15 cells. The gradient calculation with this stencil uses
the "edge" ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxJ(const sutherland &suth, const idealGas &eqnState,
                              const input &inp, const gradients &grads,
                              const turbModel *turb) {
  // suth -- method to get viscosity as a function of temperature (Sutherland's
  // law)
  // eqnState -- equation of state_
  // inp -- all input variables
  // grads -- class holding gradients at face for velocity, temperature, tke,
  // and omega
  // turb -- turbulence model

  // coefficient for viscous spectral radii
  double vCoeff = 1.0;


  // loop over all physical j-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaJ_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaJ_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.g <
               fAreaJ_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // Get velocity at face
        vector3d<double> vel =
            FaceReconCentral(state_(ii.g, jj.g - 1, kk.g).Velocity(),
                             state_(ii.g, jj.g, kk.g).Velocity(),
                             center_(ii.g, jj.g - 1, kk.g),
                             center_(ii.g, jj.g, kk.g),
                             fCenterJ_(ii.g, jj.g, kk.g));

        // Get velocity at face
        primVars state =
            FaceReconCentral(state_(ii.g, jj.g - 1, kk.g),
                             state_(ii.g, jj.g, kk.g),
                             center_(ii.g, jj.g - 1, kk.g),
                             center_(ii.g, jj.g, kk.g),
                             fCenterJ_(ii.g, jj.g, kk.g));

        // Get viscosity at face
        double mu =
            FaceReconCentral(
                suth.EffectiveViscosity(state_(ii.g, jj.g - 1, kk.g).
                                        Temperature(eqnState)),
                suth.EffectiveViscosity(state_(ii.g, jj.g, kk.g).
                                        Temperature(eqnState)),
                center_(ii.g, jj.g - 1, kk.g), center_(ii.g, jj.g, kk.g),
                fCenterJ_(ii.g, jj.g, kk.g));

        double eddyVisc =
            FaceReconCentral(
                turb->EddyVisc(state_(ii.g, jj.g - 1, kk.g),
                               grads.VelGradJ(ii.p, jj.p, kk.p), suth),
                turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                               grads.VelGradJ(ii.p, jj.p, kk.p), suth),
                             center_(ii.g, jj.g - 1, kk.g),
                             center_(ii.g, jj.g, kk.g),
                             fCenterJ_(ii.g, jj.g, kk.g));
        // effective viscosity (due to nondimensionalization)
        eddyVisc *= suth.NondimScaling();

        // calculate viscous flux
        vector3d<double> tkeGrad, omegaGrad;
        if (inp.IsTurbulent()) {
          tkeGrad = grads.TkeGradJ(ii.p, jj.p, kk.p);
          omegaGrad = grads.OmegaGradJ(ii.p, jj.p, kk.p);
        }
        viscousFlux tempViscFlux(grads.VelGradJ(ii.p, jj.p, kk.p), vel, mu,
                                 eddyVisc, suth, eqnState,
                                 grads.TempGradJ(ii.p, jj.p, kk.p),
                                 (*this).FAreaUnitJ(ii.g, jj.g, kk.g), tkeGrad,
                                 omegaGrad, turb, state);

        // area vector points from left to right, so add to left cell, subtract
        // from right cell but viscous fluxes are subtracted from inviscid
        // fluxes, so sign is reversed
        // at left boundary there is no left cell to add to
        if (jj.g > numGhosts_) {
          (*this).AddToResidual(-1.0 * tempViscFlux *
                                (*this).FAreaMagJ(ii.g, jj.g, kk.g),
                                ii.p, jj.p - 1, kk.p);
        }
        // at right boundary there is no right cell to add to
        if (jj.g < fAreaJ_.NumJ() - 2 * numGhosts_ - 1) {
          (*this).AddToResidual(tempViscFlux *
                                (*this).FAreaMagJ(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);

          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_(ii.p, jj.p, kk.p) += vCoeef *
              ViscCellSpectralRadius(
                  (*this).FAreaJ(ii.g, jj.g, kk.g),
                  (*this).FAreaJ(ii.g, jj.g + 1, kk.g),
                  state_(ii.g, jj.g, kk.g), eqnState, suth,
                  vol_(ii.g, jj.g, kk.g),
                  turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                                 grads.VelGradJ(ii.p, jj.p, kk.p), suth));
        }
      }
    }
  }
}

/* Function to calculate the viscous fluxes on the k-faces. All phyiscal
(non-ghost) k-faces are looped over. The left and right states are calculated,
and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This
contribution from the flux is added to the residuals and the wave speed is
accumulated as well.
  ___________________________
  |            |            |
  |            |            |
  |   Uk       -->   Uk+1   |
  |            |            |
  |____________|____________|
Uk-1/2       Uk+1/2       Uk+3/2

Using the above diagram, the flux is calculated at face Uk+1/2. Since the area
vector at the face always points from lower indices to higher indices it points
from Uk to Uk+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the
residual at Uk, but the sign needs to be flipped when adding the contribution of
the flux to the residual at Uk+1.

The spectral radius in the k-direction is also calculated. Since this is done on
a cell basis instead of a face bases, it is only calculated for the upper cell
(Uk+1 in this case). The spectral radius is added to the average wave speed
variable and is eventually used in the time step calculation if the time step
isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by
constructing an alternative control volume centered around that face as shown
below.
  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done.
In this example the gradient is being calculated at the face between cell Ui,j
and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the
alternative control volume. The grid cells themselves cannot be used as the
control volume for the gradients (averaging values at adjacent cells to get
gradients at the face) because this leads to odd/even decoupling. The face
areas, volumes, and states at the center_ of the faces are needed for the
alternative control volume. The left and right sides of the alternate control
volume pass through the center_ of cells Ui,j and Ui+1,j respectively. Therefore
these values are used for the face states. The left and right face areas are
calculated as the average of the face areas of the cells that they split. For
example, the left face area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j).
The top and bottom sides of the alternative control volume pass through 4 cells
each. Therefore the value of the state_ at the face center_ is the average of
these four states. For example the state_ at the top face is calculated as 0.25 *
(Ui,j + Ui+1,j + Ui,j+1, Ui+1,j+1). The face areas of the top and bottom sides
are calculated as the average of the 2 face areas that each one passes through.
For example, the top face area is calculated as 0.5 * (Ai,j+1/2 + Ai+1,j+1/2).
In three dimensions each gradient calculation touches the values at 10 cells
(6 shown and 4 more in/out of the page). The stencil for the gradients of all
faces in a cell touches 15 cells. The gradient calculation with this stencil uses
the "edge" ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxK(const sutherland &suth, const idealGas &eqnState,
                              const input &inp, const gradients &grads,
                              const turbModel *turb) {
  // suth -- method to get viscosity as a function of temperature (Sutherland's
  // law)
  // eqnState -- equation of state_
  // inp -- all input variables
  // grads -- class holding gradients at face for velocity, temperature, tke,
  // and omega
  // turb -- turbulence model

  // coefficient for viscous spectral radii
  double vCoeff = 1.0;

  // loop over all physical k-faces
  // in struct p is for physical index, g is for index with ghosts
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.g <
           fAreaK_.NumK() - numGhosts_; kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.g <
             fAreaK_.NumJ() - numGhosts_; jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.g <
               fAreaK_.NumI() - numGhosts_; ii.g++, ii.p++) {
        // Get velocity at face
        vector3d<double> vel =
            FaceReconCentral(state_(ii.g, jj.g, kk.g - 1).Velocity(),
                             state_(ii.g, jj.g, kk.g).Velocity(),
                             center_(ii.g, jj.g, kk.g - 1),
                             center_(ii.g, jj.g, kk.g),
                             fCenterK_(ii.g, jj.g, kk.g));

        // Get velocity at face
        primVars state =
            FaceReconCentral(state_(ii.g, jj.g, kk.g - 1),
                             state_(ii.g, jj.g, kk.g),
                             center_(ii.g, jj.g, kk.g - 1),
                             center_(ii.g, jj.g, kk.g),
                             fCenterK_(ii.g, jj.g, kk.g));

        // Get viscosity at face
        double mu =
            FaceReconCentral(
                suth.EffectiveViscosity(state_(ii.g, jj.g, kk.g - 1).
                                        Temperature(eqnState)),
                suth.EffectiveViscosity(state_(ii.g, jj.g, kk.g).
                                        Temperature(eqnState)),
                center_(ii.g, jj.g, kk.g - 1), center_(ii.g, jj.g, kk.g),
                fCenterK_(ii.g, jj.g, kk.g));

        double eddyVisc =
            FaceReconCentral(
            turb->EddyVisc(state_(ii.g, jj.g, kk.g - 1),
                           grads.VelGradK(ii.p, jj.p, kk.p), suth),
            turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                           grads.VelGradK(ii.p, jj.p, kk.p), suth),
            center_(ii.g, jj.g, kk.g - 1), center_(ii.g, jj.g, kk.g),
            fCenterK_(ii.g, jj.g, kk.g));
        // effective viscosity (due to nondimensionalization)
        eddyVisc *= suth.NondimScaling();

        // calculate viscous flux
        vector3d<double> tkeGrad, omegaGrad;
        if (inp.IsTurbulent()) {
          tkeGrad = grads.TkeGradK(ii.p, jj.p, kk.p);
          omegaGrad = grads.OmegaGradK(ii.p, jj.p, kk.p);
        }
        viscousFlux tempViscFlux(grads.VelGradK(ii.p, jj.p, kk.p), vel, mu,
                                 eddyVisc, suth, eqnState,
                                 grads.TempGradK(ii.p, jj.p, kk.p),
                                 (*this).FAreaUnitK(ii.g, jj.g, kk.g), tkeGrad,
                                 omegaGrad, turb, state);

        // area vector points from left to right, so add to left cell, subtract
        // from right cell but viscous fluxes are subtracted from inviscid
        // fluxes, so sign is reversed
        // at left boundary there is no left cell to add to
        if (kk.g > numGhosts_) {
          (*this).AddToResidual(-1.0 * tempViscFlux *
                                (*this).FAreaMagK(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p - 1);
        }
        // at right boundary there is no right cell to add to
        if (kk.g < fAreaK_.NumK() - 2 * numGhosts_ - 1) {
          (*this).AddToResidual(tempViscFlux *
                                (*this).FAreaMagK(ii.g, jj.g, kk.g),
                                ii.p, jj.p, kk.p);

          // calculate component of wave speed. This is done on a cell by cell
          // basis, so only at the upper faces
          avgWaveSpeed_(ii.p, jj.p, kk.p) += vCoeff *
              ViscCellSpectralRadius(
                  (*this).FAreaK(ii.g, jj.g, kk.g),
                  (*this).FAreaK(ii.g, jj.g, kk.g + 1),
                  state_(ii.g, jj.g, kk.g), eqnState, suth,
                  vol_(ii.g, jj.g, kk.g),
                  turb->EddyVisc(state_(ii.g, jj.g, kk.g),
                                 grads.VelGradK(ii.p, jj.p, kk.p), suth));
        }
      }
    }
  }
}

/* Member function to assign geometric quantities such as volume, face area,
cell centroid, and face center_ to ghost cells. This assigns values for
regular ghost cells and "edge" ghost cells. "Corner" cells are left with no
value as they are not used.
           ____ ____ ____ ____ ____ ____ ____ ____
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|

In the above diagram where X represents the physical cells, cells marked G
(regular ghost cells) and E ("edge" ghost cells) are assigned geometric
values. G1 represents the first layer of ghost cells and G2 represents the
second layer.
*/
void procBlock::AssignGhostCellsGeom() {
  // inp -- all input variables

  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * numGhosts_;
  int jmaxG = jmax + 2 * numGhosts_;
  int kmaxG = kmax + 2 * numGhosts_;

  //------------------------------------------------------------------------
  // loop over physical I faces and assign values for regular ghost cells
  // ------------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(0, jj - numGhosts_,
                                              kk - numGhosts_, "il");
      string bcNameU = (*this).BC().GetBCName(imax, jj - numGhosts_,
                                              kk - numGhosts_, "iu");

      if (bcNameL == "interblock" && bcNameU == "interblock" && imax < 2) {
        cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot "
                "have interblock BCs surrounding a block that is 1 cell thick!"
             << endl;
        exit(0);
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at lower i-boundary
        int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first ghost cell at lower
        // i-boundary
        int lFaceG1_il = GetLowerFaceI(1, jj, kk, imaxG, jmaxG);
        int lFaceG1_jl = GetLowerFaceJ(1, jj, kk, imaxG, jmaxG);
        int lFaceG1_kl = GetLowerFaceK(1, jj, kk, imaxG, jmaxG);

        // location of second ghost cell at lower i-boundary
        int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first ghost cell at lower
        // i-boundary
        int lFaceG2_il = GetLowerFaceI(0, jj, kk, imaxG, jmaxG);
        int lFaceG2_jl = GetLowerFaceJ(0, jj, kk, imaxG, jmaxG);
        int lFaceG2_kl = GetLowerFaceK(0, jj, kk, imaxG, jmaxG);

        // location of first interior cell at lower i-boundary
        int cellLowIn1 = GetLoc1D(numGhosts_, jj, kk, imaxG, jmaxG);
        // location of upper i and lower j, k faces of first interior cell at
        // lower i-boundary
        int lFaceIn1_iu =
            GetUpperFaceI(numGhosts_, jj, kk, imaxG, jmaxG);
        int lFaceIn1_jl =
            GetLowerFaceJ(numGhosts_, jj, kk, imaxG, jmaxG);
        int lFaceIn1_kl =
            GetLowerFaceK(numGhosts_, jj, kk, imaxG, jmaxG);

        // location of second interior cell at lower i-boundary
        int cellLowIn2 =
            GetLoc1D(numGhosts_ + 1, jj, kk, imaxG, jmaxG);
        // location of upper i and lower j, k faces of second interior cell at
        // lower i-boundary
        int lFaceIn2_iu =
            GetUpperFaceI(numGhosts_ + 1, jj, kk, imaxG, jmaxG);
        int lFaceIn2_jl =
            GetLowerFaceJ(numGhosts_ + 1, jj, kk, imaxG, jmaxG);
        int lFaceIn2_kl =
            GetLowerFaceK(numGhosts_ + 1, jj, kk, imaxG, jmaxG);

        // location of lower i-boundary face
        int lFaceB = GetLowerFaceI(numGhosts_, jj, kk, imaxG, jmaxG);

        // Assign volume
        // -----------------------------------------------------------------
        // mirror volume values from adjacent cells across i-boundary
        // first layer of ghost cells
        (*this).vol_[cellLowG1] = vol_(cellLowIn1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellLowG2] = vol_(cellLowIn1);
        } else {
          (*this).vol_[cellLowG2] = vol_(cellLowIn2);
        }

        // Assign face areas
        // -----------------------------------------------------------------
        // mirror face area values from adjacent cells across i-boundary
        // first layer of ghost cells
        (*this).fAreaI_[lFaceG1_il] = (*this).FAreaI(lFaceIn1_iu);
        (*this).fAreaJ_[lFaceG1_jl] = (*this).FAreaJ(lFaceIn1_jl);
        (*this).fAreaK_[lFaceG1_kl] = (*this).FAreaK(lFaceIn1_kl);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn1_iu);
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn1_jl);
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn1_kl);
        } else {
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn2_iu);
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn2_jl);
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn2_kl);
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at lower i-boundary
          int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG);
          int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

          // location of upper j-face for interior cells at lower i-boundary
          int lFaceIn1_ju =
              GetUpperFaceJ(numGhosts_, jj, kk, imaxG, jmaxG);
          int lFaceIn2_ju =
              GetUpperFaceJ(numGhosts_ + 1, jj, kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaJ_[lFaceG1_ju] = (*this).FAreaJ(lFaceIn1_ju);

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaJ_[lFaceG2_ju] = (*this).FAreaJ(lFaceIn1_ju);
          } else {
            (*this).fAreaJ_[lFaceG2_ju] = (*this).FAreaJ(lFaceIn2_ju);
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at lower i-boundary
          int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG);
          int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at lower i-boundary
          int lFaceIn1_ku =
              GetUpperFaceK(numGhosts_, jj, kk, imaxG, jmaxG);
          int lFaceIn2_ku =
              GetUpperFaceK(numGhosts_ + 1, jj, kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaK_[lFaceG1_ku] = (*this).FAreaK(lFaceIn1_ku);

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaK_[lFaceG2_ku] = (*this).FAreaK(lFaceIn1_ku);
          } else {
            (*this).fAreaK_[lFaceG2_ku] = (*this).FAreaK(lFaceIn2_ku);
          }
        }

        // Assign cell centroid
        // ------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu);
        (*this).center_[cellLowG1] = center_(cellLowIn1) + dist2Move;

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu));
          (*this).center_[cellLowG2] = center_(cellLowG1) + dist2Move;
        } else {
          dist2Move = fCenterI_(lFaceB) - fCenterI_(lFaceIn2_iu);
          (*this).center_[cellLowG2] = center_(cellLowIn1) + dist2Move;
        }

        // Assign face centers
        // ------------------------------------------------------------------
        // face center is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu);
        (*this).fCenterI_[lFaceG1_il] = fCenterI_(lFaceB) + dist2Move;
        (*this).fCenterJ_[lFaceG1_jl] =
            fCenterJ_(lFaceIn1_jl) + dist2Move;
        (*this).fCenterK_[lFaceG1_kl] =
            fCenterK_(lFaceIn1_kl) + dist2Move;

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu));
          (*this).fCenterI_[lFaceG2_il] =
              fCenterI_(lFaceG1_il) + dist2Move;
          (*this).fCenterJ_[lFaceG2_jl] =
              fCenterJ_(lFaceG1_jl) + dist2Move;
          (*this).fCenterK_[lFaceG2_kl] =
              fCenterK_(lFaceG1_kl) + dist2Move;
        } else {
          dist2Move = fCenterI_(lFaceB) - fCenterI_(lFaceIn2_iu);
          (*this).fCenterI_[lFaceG2_il] = fCenterI_(lFaceB) + dist2Move;
          (*this).fCenterJ_[lFaceG2_jl] =
              fCenterJ_(lFaceIn1_jl) + dist2Move;
          (*this).fCenterK_[lFaceG2_kl] =
              fCenterK_(lFaceIn1_kl) + dist2Move;
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // face areas too

          // location of upper j-face for ghost cells at lower i-boundary
          int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG);
          int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

          // location of upper j-face for interior cell at lower i-boundary
          int lFaceIn1_ju =
              GetUpperFaceJ(numGhosts_, jj, kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu);
          (*this).fCenterJ_[lFaceG1_ju] =
              fCenterJ_(lFaceIn1_ju) + dist2Move;

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterI_(lFaceB) -
                               fCenterI_(lFaceIn1_iu));
            (*this).fCenterJ_[lFaceG2_ju] =
                fCenterJ_(lFaceG1_ju) + dist2Move;
          } else {
            dist2Move =
                fCenterI_(lFaceB) - fCenterI_(lFaceIn2_iu);
            (*this).fCenterJ_[lFaceG2_ju] =
                fCenterJ_(lFaceIn1_ju) + dist2Move;
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // face areas too

          // location of upper k-face for ghost cells at lower i-boundary
          int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG);
          int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at lower i-boundary
          int lFaceIn1_ku =
              GetUpperFaceK(numGhosts_, jj, kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterI_(lFaceB) - fCenterI_(lFaceIn1_iu);
          (*this).fCenterK_[lFaceG1_ku] =
              fCenterK_(lFaceIn1_ku) + dist2Move;

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterI_(lFaceB) -
                               fCenterI_(lFaceIn1_iu));
            (*this).fCenterK_[lFaceG2_ku] =
                fCenterK_(lFaceG1_ku) + dist2Move;
          } else {
            dist2Move =
                fCenterI_(lFaceB) - fCenterI_(lFaceIn2_iu);
            (*this).fCenterK_[lFaceG2_ku] =
                fCenterK_(lFaceIn1_ku) + dist2Move;
          }
        }
      }

      // upper surface ----------------------------------------------------
      if (bcNameU != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at upper i-boundary
        int cellUpG1 = GetLoc1D(imaxG - 2, jj, kk, imaxG, jmaxG);
        // location of upper i and lower j, k faces of first ghost cell at upper
        // i-boundary
        int uFaceG1_iu = GetUpperFaceI(imaxG - 2, jj, kk, imaxG, jmaxG);
        int uFaceG1_jl = GetLowerFaceJ(imaxG - 2, jj, kk, imaxG, jmaxG);
        int uFaceG1_kl = GetLowerFaceK(imaxG - 2, jj, kk, imaxG, jmaxG);

        // location of second ghost cell at upper i-boundary
        int cellUpG2 = GetLoc1D(imaxG - 1, jj, kk, imaxG, jmaxG);
        // location of upper i and lower j, k faces of second ghost cell at
        // upper i-boundary
        int uFaceG2_iu = GetUpperFaceI(imaxG - 1, jj, kk, imaxG, jmaxG);
        int uFaceG2_jl = GetLowerFaceJ(imaxG - 1, jj, kk, imaxG, jmaxG);
        int uFaceG2_kl = GetLowerFaceK(imaxG - 1, jj, kk, imaxG, jmaxG);

        // location of first interior cell at upper i-boundary
        int cellUpIn1 =
            GetLoc1D(imaxG - 1 - numGhosts_, jj, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first interior cell at upper
        // i-boundary
        int uFaceIn1_il = GetLowerFaceI(imaxG - 1 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);
        int uFaceIn1_jl = GetLowerFaceJ(imaxG - 1 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);
        int uFaceIn1_kl = GetLowerFaceK(imaxG - 1 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);

        // location of second interior cell at upper i-boundary
        int cellUpIn2 =
            GetLoc1D(imaxG - 2 - numGhosts_, jj, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of second interior cell at upper
        // i-boundary
        int uFaceIn2_il = GetLowerFaceI(imaxG - 2 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);
        int uFaceIn2_jl = GetLowerFaceJ(imaxG - 2 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);
        int uFaceIn2_kl = GetLowerFaceK(imaxG - 2 - numGhosts_, jj, kk,
                                        imaxG, jmaxG);

        // location of upper i-boundary face
        int uFaceB = GetUpperFaceI(imaxG - 1 - numGhosts_, jj, kk,
                                   imaxG, jmaxG);

        // Assign volume
        // --------------------------------------------------------------------
        // mirror volume values from adjacent cells across i-boundary
        // first layer of ghost cells
        (*this).vol_[cellUpG1] = vol_(cellUpIn1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellUpG2] = vol_(cellUpIn1);
        } else {
          (*this).vol_[cellUpG2] = vol_(cellUpIn2);
        }

        // Assign face areas
        // --------------------------------------------------------------------
        // mirror face area values from adjacent cells across i-boundary
        // first layer of ghost cells
        (*this).fAreaI_[uFaceG1_iu] = (*this).FAreaI(uFaceIn1_il);
        (*this).fAreaJ_[uFaceG1_jl] = (*this).FAreaJ(uFaceIn1_jl);
        (*this).fAreaK_[uFaceG1_kl] = (*this).FAreaK(uFaceIn1_kl);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn1_il);
          (*this).fAreaJ_[uFaceG2_jl] = (*this).FAreaJ(uFaceIn1_jl);
          (*this).fAreaK_[uFaceG2_kl] = (*this).FAreaK(uFaceIn1_kl);
        } else {
          (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn2_il);
          (*this).fAreaJ_[uFaceG2_jl] = (*this).FAreaJ(uFaceIn2_jl);
          (*this).fAreaK_[uFaceG2_kl] = (*this).FAreaK(uFaceIn2_kl);
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at upper i-boundary
          int uFaceG1_ju = GetUpperFaceJ(imaxG - 2, jj, kk, imaxG, jmaxG);
          int uFaceG2_ju = GetUpperFaceJ(imaxG - 1, jj, kk, imaxG, jmaxG);

          // location of upper j-face for interior cells at upper i-boundary
          int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);
          int uFaceIn2_ju = GetUpperFaceJ(imaxG - 2 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaJ_[uFaceG1_ju] = (*this).FAreaJ(uFaceIn1_ju);

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn1_ju);
          } else {
            (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn2_ju);
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at upper i-boundary
          int uFaceG1_ku = GetUpperFaceK(imaxG - 2, jj, kk, imaxG, jmaxG);
          int uFaceG2_ku = GetUpperFaceK(imaxG - 1, jj, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at upper i-boundary
          int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);
          int uFaceIn2_ku = GetUpperFaceK(imaxG - 2 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaK_[uFaceG1_ku] = (*this).FAreaK(uFaceIn1_ku);

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn1_ku);
          } else {
            (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn2_ku);
          }
        }

        // Assign cell centroid
        // ---------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il);
        (*this).center_[cellUpG1] = center_(cellUpIn1) + dist2Move;

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il));
          (*this).center_[cellUpG2] = center_(cellUpG1) + dist2Move;
        } else {
          dist2Move = fCenterI_(uFaceB) - fCenterI_(uFaceIn2_il);
          (*this).center_[cellUpG2] = center_(cellUpIn1) + dist2Move;
        }

        // Assign face centers
        // ---------------------------------------------------------------------
        // face center_ is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il);
        (*this).fCenterI_[uFaceG1_iu] = fCenterI_(uFaceB) + dist2Move;
        (*this).fCenterJ_[uFaceG1_jl] =
            fCenterJ_(uFaceIn1_jl) + dist2Move;
        (*this).fCenterK_[uFaceG1_kl] =
            fCenterK_(uFaceIn1_kl) + dist2Move;

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il));
          (*this).fCenterI_[uFaceG2_iu] =
              fCenterI_(uFaceG1_iu) + dist2Move;
          (*this).fCenterJ_[uFaceG2_jl] =
              fCenterJ_(uFaceG1_jl) + dist2Move;
          (*this).fCenterK_[uFaceG2_kl] =
              fCenterK_(uFaceG1_kl) + dist2Move;
        } else {
          dist2Move = fCenterI_(uFaceB) - fCenterI_(uFaceIn2_il);
          (*this).fCenterI_[uFaceG2_iu] = fCenterI_(uFaceB) + dist2Move;
          (*this).fCenterJ_[uFaceG2_jl] =
              fCenterJ_(uFaceIn1_jl) + dist2Move;
          (*this).fCenterK_[uFaceG2_kl] =
              fCenterK_(uFaceIn1_kl) + dist2Move;
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // face areas too

          // location of upper j-face for ghost cells at lower i-boundary
          int uFaceG1_ju = GetUpperFaceJ(imaxG - 2, jj, kk, imaxG, jmaxG);
          int uFaceG2_ju = GetUpperFaceJ(imaxG - 1, jj, kk, imaxG, jmaxG);

          // location of upper j-face for interior cell at upper i-boundary
          int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il);
          (*this).fCenterJ_[uFaceG1_ju] =
              fCenterJ_(uFaceIn1_ju) + dist2Move;

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterI_(uFaceB) -
                               fCenterI_(uFaceIn1_il));
            (*this).fCenterJ_[uFaceG2_ju] =
                fCenterJ_(uFaceG1_ju) + dist2Move;
          } else {
            dist2Move =
                fCenterI_(uFaceB) - fCenterI_(uFaceIn2_il);
            (*this).fCenterJ_[uFaceG2_ju] =
                fCenterJ_(uFaceIn1_ju) + dist2Move;
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // face areas too

          // location of upper j-face for ghost cells at upper i-boundary
          int uFaceG1_ku = GetUpperFaceK(imaxG - 2, jj, kk, imaxG, jmaxG);
          int uFaceG2_ku = GetUpperFaceK(imaxG - 1, jj, kk, imaxG, jmaxG);

          // location of upper j-face for interior cells at upper i-boundary
          int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - numGhosts_, jj,
                                          kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterI_(uFaceB) - fCenterI_(uFaceIn1_il);
          (*this).fCenterK_[uFaceG1_ku] =
              fCenterK_(uFaceIn1_ku) + dist2Move;

          // second layer of ghost cells
          if (imax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterI_(uFaceB) -
                               fCenterI_(uFaceIn1_il));
            (*this).fCenterK_[uFaceG2_ku] =
                fCenterK_(uFaceG1_ku) + dist2Move;
          } else {
            dist2Move =
                fCenterI_(uFaceB) - fCenterI_(uFaceIn2_il);
            (*this).fCenterK_[uFaceG2_ku] =
                fCenterK_(uFaceIn1_ku) + dist2Move;
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // loop over physical J faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(ii - numGhosts_, 0,
                                              kk - numGhosts_, "jl");
      string bcNameU = (*this).BC().GetBCName(ii - numGhosts_, jmax,
                                              kk - numGhosts_, "ju");

      if (bcNameL == "interblock" && bcNameU == "interblock" && jmax < 2) {
        cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot "
                "have interblock BCs surrounding a block that is 1 cell thick!"
             << endl;
        exit(0);
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at lower j-boundary
        int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first ghost cell at lower
        // j-boundary
        int lFaceG1_il = GetLowerFaceI(ii, 1, kk, imaxG, jmaxG);
        int lFaceG1_jl = GetLowerFaceJ(ii, 1, kk, imaxG, jmaxG);
        int lFaceG1_kl = GetLowerFaceK(ii, 1, kk, imaxG, jmaxG);

        // location of second ghost cell at lower i-boundary
        int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first ghost cell at lower
        // j-boundary
        int lFaceG2_il = GetLowerFaceI(ii, 0, kk, imaxG, jmaxG);
        int lFaceG2_jl = GetLowerFaceJ(ii, 0, kk, imaxG, jmaxG);
        int lFaceG2_kl = GetLowerFaceK(ii, 0, kk, imaxG, jmaxG);

        // location of first interior cell at lower i-boundary
        int cellLowIn1 = GetLoc1D(ii, numGhosts_, kk, imaxG, jmaxG);
        // location of upper j and lower i, k faces of first interior cell at
        // lower j-boundary
        int lFaceIn1_il =
            GetLowerFaceI(ii, numGhosts_, kk, imaxG, jmaxG);
        int lFaceIn1_ju =
            GetUpperFaceJ(ii, numGhosts_, kk, imaxG, jmaxG);
        int lFaceIn1_kl =
            GetLowerFaceK(ii, numGhosts_, kk, imaxG, jmaxG);

        // location of second interior cell at lower i-boundary
        int cellLowIn2 =
            GetLoc1D(ii, numGhosts_ + 1, kk, imaxG, jmaxG);
        // location of upper j and lower i, k faces of second interior cell at
        // lower j-boundary
        int lFaceIn2_il =
            GetLowerFaceI(ii, numGhosts_ + 1, kk, imaxG, jmaxG);
        int lFaceIn2_ju =
            GetUpperFaceJ(ii, numGhosts_ + 1, kk, imaxG, jmaxG);
        int lFaceIn2_kl =
            GetLowerFaceK(ii, numGhosts_ + 1, kk, imaxG, jmaxG);

        // location of lower j-boundary face
        int lFaceB = GetLowerFaceJ(ii, numGhosts_, kk, imaxG, jmaxG);

        // Assign volume
        // --------------------------------------------------------------------
        // mirror volume values from adjacent cells across j-boundary
        // first layer of ghost cells
        (*this).vol_[cellLowG1] = vol_(cellLowIn1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellLowG2] = vol_(cellLowIn1);
        } else {
          (*this).vol_[cellLowG2] = vol_(cellLowIn2);
        }

        // Assign face areas
        // --------------------------------------------------------------------
        // mirror face area values from adjacent cells across j-boundary
        // first layer of ghost cells
        (*this).fAreaJ_[lFaceG1_jl] = (*this).FAreaJ(lFaceIn1_ju);
        (*this).fAreaI_[lFaceG1_il] = (*this).FAreaI(lFaceIn1_il);
        (*this).fAreaK_[lFaceG1_kl] = (*this).FAreaK(lFaceIn1_kl);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn1_ju);
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn1_il);
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn1_kl);
        } else {
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn2_ju);
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn2_il);
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn2_kl);
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at lower j-boundary
          int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG);
          int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

          // location of upper i-face for interior cells at lower j-boundary
          int lFaceIn1_iu =
              GetUpperFaceI(ii, numGhosts_, kk, imaxG, jmaxG);
          int lFaceIn2_iu =
              GetUpperFaceI(ii, numGhosts_ + 1, kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaI_[lFaceG1_iu] = (*this).FAreaI(lFaceIn1_iu);

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaI_[lFaceG2_iu] = (*this).FAreaI(lFaceIn1_iu);
          } else {
            (*this).fAreaI_[lFaceG2_iu] = (*this).FAreaI(lFaceIn2_iu);
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at lower j-boundary
          int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG);
          int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at upper j-boundary
          int lFaceIn1_ku =
              GetUpperFaceK(ii, numGhosts_, kk, imaxG, jmaxG);
          int lFaceIn2_ku =
              GetUpperFaceK(ii, numGhosts_ + 1, kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaK_[lFaceG1_ku] = (*this).FAreaK(lFaceIn1_ku);

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaK_[lFaceG2_ku] = (*this).FAreaK(lFaceIn1_ku);
          } else {
            (*this).fAreaK_[lFaceG2_ku] = (*this).FAreaK(lFaceIn2_ku);
          }
        }

        // Assign cell centroid
        // --------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju);
        (*this).center_[cellLowG1] = center_(cellLowIn1) + dist2Move;

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju));
          (*this).center_[cellLowG2] = center_(cellLowG1) + dist2Move;
        } else {
          dist2Move = fCenterJ_(lFaceB) - fCenterJ_(lFaceIn2_ju);
          (*this).center_[cellLowG2] = center_(cellLowIn1) + dist2Move;
        }

        // Assign face centers
        // -------------------------------------------------------------------
        // face center_ is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju);
        (*this).fCenterJ_[lFaceG1_jl] = fCenterJ_(lFaceB) + dist2Move;
        (*this).fCenterI_[lFaceG1_il] =
            fCenterI_(lFaceIn1_il) + dist2Move;
        (*this).fCenterK_[lFaceG1_kl] =
            fCenterK_(lFaceIn1_kl) + dist2Move;

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju));
          (*this).fCenterJ_[lFaceG2_jl] =
              fCenterJ_(lFaceG1_jl) + dist2Move;
          (*this).fCenterI_[lFaceG2_il] =
              fCenterI_(lFaceG1_il) + dist2Move;
          (*this).fCenterK_[lFaceG2_kl] =
              fCenterK_(lFaceG1_kl) + dist2Move;
        } else {
          dist2Move = fCenterJ_(lFaceB) - fCenterJ_(lFaceIn2_ju);
          (*this).fCenterJ_[lFaceG2_jl] = fCenterJ_(lFaceB) + dist2Move;
          (*this).fCenterI_[lFaceG2_il] =
              fCenterI_(lFaceIn1_il) + dist2Move;
          (*this).fCenterK_[lFaceG2_kl] =
              fCenterK_(lFaceIn1_kl) + dist2Move;
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at lower j-boundary
          int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG);
          int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

          // location of upper i-face for interior cells at lower j-boundary
          int lFaceIn1_iu =
              GetUpperFaceI(ii, numGhosts_, kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju);
          (*this).fCenterI_[lFaceG1_iu] =
              fCenterI_(lFaceIn1_iu) + dist2Move;

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterJ_(lFaceB) -
                               fCenterJ_(lFaceIn1_ju));
            (*this).fCenterI_[lFaceG2_iu] =
                fCenterI_(lFaceG1_iu) + dist2Move;
          } else {
            dist2Move =
                fCenterJ_(lFaceB) - fCenterJ_(lFaceIn2_ju);
            (*this).fCenterI_[lFaceG2_iu] =
                fCenterI_(lFaceIn1_iu) + dist2Move;
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at lower j-boundary
          int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG);
          int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at lower j-boundary
          int lFaceIn1_ku =
              GetUpperFaceK(ii, numGhosts_, kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterJ_(lFaceB) - fCenterJ_(lFaceIn1_ju);
          (*this).fCenterK_[lFaceG1_ku] =
              fCenterK_(lFaceIn1_ku) + dist2Move;

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterJ_(lFaceB) -
                               fCenterJ_(lFaceIn1_ju));
            (*this).fCenterK_[lFaceG2_ku] =
                fCenterK_(lFaceG1_ku) + dist2Move;
          } else {
            dist2Move =
                fCenterJ_(lFaceB) - fCenterJ_(lFaceIn2_ju);
            (*this).fCenterK_[lFaceG2_ku] =
                fCenterK_(lFaceIn1_ku) + dist2Move;
          }
        }
      }

      // upper surface ----------------------------------------------------
      if (bcNameU != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at upper j-boundary
        int cellUpG1 = GetLoc1D(ii, jmaxG - 2, kk, imaxG, jmaxG);
        // location of upper j and lower i, k faces of first ghost cell at lower
        // j-boundary
        int uFaceG1_il = GetLowerFaceI(ii, jmaxG - 2, kk, imaxG, jmaxG);
        int uFaceG1_ju = GetUpperFaceJ(ii, jmaxG - 2, kk, imaxG, jmaxG);
        int uFaceG1_kl = GetLowerFaceK(ii, jmaxG - 2, kk, imaxG, jmaxG);

        // location of second ghost cell at upper j-boundary
        int cellUpG2 = GetLoc1D(ii, jmaxG - 1, kk, imaxG, jmaxG);
        // location of upper j and lower i, k faces of second ghost cell at
        // lower j-boundary
        int uFaceG2_il = GetLowerFaceI(ii, jmaxG - 1, kk, imaxG, jmaxG);
        int uFaceG2_ju = GetUpperFaceJ(ii, jmaxG - 1, kk, imaxG, jmaxG);
        int uFaceG2_kl = GetLowerFaceK(ii, jmaxG - 1, kk, imaxG, jmaxG);

        // location of first interior cell at upper j-boundary
        int cellUpIn1 =
            GetLoc1D(ii, jmaxG - 1 - numGhosts_, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of first interior cell at lower
        // j-boundary
        int uFaceIn1_il = GetLowerFaceI(ii, jmaxG - 1 - numGhosts_, kk,
                                        imaxG, jmaxG);
        int uFaceIn1_jl = GetLowerFaceJ(ii, jmaxG - 1 - numGhosts_, kk,
                                        imaxG, jmaxG);
        int uFaceIn1_kl = GetLowerFaceK(ii, jmaxG - 1 - numGhosts_, kk,
                                        imaxG, jmaxG);

        // location of second interior cell at upper j-boundary
        int cellUpIn2 =
            GetLoc1D(ii, jmaxG - 2 - numGhosts_, kk, imaxG, jmaxG);
        // location of lower i, j, k faces of second interior cell at lower
        // j-boundary
        int uFaceIn2_il = GetLowerFaceI(ii, jmaxG - 2 - numGhosts_, kk,
                                        imaxG, jmaxG);
        int uFaceIn2_jl = GetLowerFaceJ(ii, jmaxG - 2 - numGhosts_, kk,
                                        imaxG, jmaxG);
        int uFaceIn2_kl = GetLowerFaceK(ii, jmaxG - 2 - numGhosts_, kk,
                                        imaxG, jmaxG);

        // location of upper j-boundary face
        int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - numGhosts_, kk,
                                   imaxG, jmaxG);

        // Assign volume
        // --------------------------------------------------------------------
        // mirror volume values from adjacent cells across j-boundary
        // first layer of ghost cells
        (*this).vol_[cellUpG1] = vol_(cellUpIn1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellUpG2] = vol_(cellUpIn1);
        } else {
          (*this).vol_[cellUpG2] = vol_(cellUpIn2);
        }

        // Assign face areas
        // --------------------------------------------------------------------
        // mirror face area values from adjacent cells across j-boundary
        // first layer of ghost cells
        (*this).fAreaJ_[uFaceG1_ju] = (*this).FAreaJ(uFaceIn1_jl);
        (*this).fAreaI_[uFaceG1_il] = (*this).FAreaI(uFaceIn1_il);
        (*this).fAreaK_[uFaceG1_kl] = (*this).FAreaK(uFaceIn1_kl);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn1_jl);
          (*this).fAreaI_[uFaceG2_il] = (*this).FAreaI(uFaceIn1_il);
          (*this).fAreaK_[uFaceG2_kl] = (*this).FAreaK(uFaceIn1_kl);
        } else {
          (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn2_jl);
          (*this).fAreaI_[uFaceG2_il] = (*this).FAreaI(uFaceIn2_il);
          (*this).fAreaK_[uFaceG2_kl] = (*this).FAreaK(uFaceIn2_kl);
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at upper j-boundary
          int uFaceG1_iu = GetUpperFaceI(ii, jmaxG - 2, kk, imaxG, jmaxG);
          int uFaceG2_iu = GetUpperFaceI(ii, jmaxG - 1, kk, imaxG, jmaxG);

          // location of upper i-face for interior cells at upper j-boundary
          int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - numGhosts_,
                                          kk, imaxG, jmaxG);
          int uFaceIn2_iu = GetUpperFaceI(ii, jmaxG - 2 - numGhosts_,
                                          kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaI_[uFaceG1_iu] = (*this).FAreaI(uFaceIn1_iu);

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn1_iu);
          } else {
            (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn2_iu);
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at upper j-boundary
          int uFaceG1_ku = GetUpperFaceK(ii, jmaxG - 2, kk, imaxG, jmaxG);
          int uFaceG2_ku = GetUpperFaceK(ii, jmaxG - 1, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at upper j-boundary
          int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - numGhosts_,
                                          kk, imaxG, jmaxG);
          int uFaceIn2_ku = GetUpperFaceK(ii, jmaxG - 2 - numGhosts_,
                                          kk, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaK_[uFaceG1_ku] = (*this).FAreaK(uFaceIn1_ku);

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn1_ku);
          } else {
            (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn2_ku);
          }
        }

        // Assign cell centroid
        // -------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl);
        (*this).center_[cellUpG1] = center_(cellUpIn1) + dist2Move;

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl));
          (*this).center_[cellUpG2] = center_(cellUpG1) + dist2Move;
        } else {
          dist2Move = fCenterJ_(uFaceB) - fCenterJ_(uFaceIn2_jl);
          (*this).center_[cellUpG2] = center_(cellUpIn1) + dist2Move;
        }

        // Assign face centers
        // -------------------------------------------------------------------
        // face center_ is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl);
        (*this).fCenterJ_[uFaceG1_ju] = fCenterJ_(uFaceB) + dist2Move;
        (*this).fCenterI_[uFaceG1_il] =
            fCenterI_(uFaceIn1_il) + dist2Move;
        (*this).fCenterK_[uFaceG1_kl] =
            fCenterK_(uFaceIn1_kl) + dist2Move;

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl));
          (*this).fCenterJ_[uFaceG2_ju] =
              fCenterJ_(uFaceG1_ju) + dist2Move;
          (*this).fCenterI_[uFaceG2_il] =
              fCenterI_(uFaceG1_il) + dist2Move;
          (*this).fCenterK_[uFaceG2_kl] =
              fCenterK_(uFaceG1_kl) + dist2Move;
        } else {
          dist2Move = fCenterJ_(uFaceB) - fCenterJ_(uFaceIn2_jl);
          (*this).fCenterJ_[uFaceG2_ju] = fCenterJ_(uFaceB) + dist2Move;
          (*this).fCenterI_[uFaceG2_il] =
              fCenterI_(uFaceIn1_il) + dist2Move;
          (*this).fCenterK_[uFaceG2_kl] =
              fCenterK_(uFaceIn1_kl) + dist2Move;
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at upper j-boundary
          int uFaceG1_iu = GetUpperFaceI(ii, jmaxG - 2, kk, imaxG, jmaxG);
          int uFaceG2_iu = GetUpperFaceI(ii, jmaxG - 1, kk, imaxG, jmaxG);

          // location of upper i-face for interior cells at upper j-boundary
          int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - numGhosts_,
                                          kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl);
          (*this).fCenterI_[uFaceG1_iu] =
              fCenterI_(uFaceIn1_iu) + dist2Move;

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterJ_(uFaceB) -
                               fCenterJ_(uFaceIn1_jl));
            (*this).fCenterI_[uFaceG2_iu] =
                fCenterI_(uFaceG1_iu) + dist2Move;
          } else {
            dist2Move =
                fCenterJ_(uFaceB) - fCenterJ_(uFaceIn2_jl);
            (*this).fCenterI_[uFaceG2_iu] =
                fCenterI_(uFaceIn1_iu) + dist2Move;
          }
        }

        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // k-face areas too

          // location of upper k-face for ghost cells at upper j-boundary
          int uFaceG1_ku = GetUpperFaceK(ii, jmaxG - 2, kk, imaxG, jmaxG);
          int uFaceG2_ku = GetUpperFaceK(ii, jmaxG - 1, kk, imaxG, jmaxG);

          // location of upper k-face for interior cells at upper j-boundary
          int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - numGhosts_,
                                          kk, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterJ_(uFaceB) - fCenterJ_(uFaceIn1_jl);
          (*this).fCenterK_[uFaceG1_ku] =
              fCenterK_(uFaceIn1_ku) + dist2Move;

          // second layer of ghost cells
          if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterJ_(uFaceB) -
                               fCenterJ_(uFaceIn1_jl));
            (*this).fCenterK_[uFaceG2_ku] =
                fCenterK_(uFaceG1_ku) + dist2Move;
          } else {
            dist2Move =
                fCenterJ_(uFaceB) - fCenterJ_(uFaceIn2_jl);
            (*this).fCenterK_[uFaceG2_ku] =
                fCenterK_(uFaceIn1_ku) + dist2Move;
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // loop over physical K faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, 0, "kl");
      string bcNameU = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, kmax, "ku");

      if (bcNameL == "interblock" && bcNameU == "interblock" && kmax < 2) {
        cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot "
                "have interblock BCs surrounding a block that is 1 cell thick!"
             << endl;
        exit(0);
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at lower k-boundary
        int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
        // location of lower i, j, k faces of first ghost cell at lower
        // k-boundary
        int lFaceG1_il = GetLowerFaceI(ii, jj, 1, imaxG, jmaxG);
        int lFaceG1_jl = GetLowerFaceJ(ii, jj, 1, imaxG, jmaxG);
        int lFaceG1_kl = GetLowerFaceK(ii, jj, 1, imaxG, jmaxG);

        // location of second ghost cell at lower k-boundary
        int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);
        // location of lower i, j, k faces of second ghost cell at lower
        // k-boundary
        int lFaceG2_il = GetLowerFaceI(ii, jj, 0, imaxG, jmaxG);
        int lFaceG2_jl = GetLowerFaceJ(ii, jj, 0, imaxG, jmaxG);
        int lFaceG2_kl = GetLowerFaceK(ii, jj, 0, imaxG, jmaxG);

        // location of first interior cell at lower k-boundary
        int cellLowIn1 = GetLoc1D(ii, jj, numGhosts_, imaxG, jmaxG);
        // location of upper k and lower i, j faces of first interior cell at
        // lower k-boundary
        int lFaceIn1_il =
            GetLowerFaceI(ii, jj, numGhosts_, imaxG, jmaxG);
        int lFaceIn1_jl =
            GetLowerFaceJ(ii, jj, numGhosts_, imaxG, jmaxG);
        int lFaceIn1_ku =
            GetUpperFaceK(ii, jj, numGhosts_, imaxG, jmaxG);

        // location of second interior cell at lower k-boundary
        int cellLowIn2 =
            GetLoc1D(ii, jj, numGhosts_ + 1, imaxG, jmaxG);
        // location of upper k and lower i, j faces of second interior cell at
        // lower k-boundary
        int lFaceIn2_il =
            GetLowerFaceI(ii, jj, numGhosts_ + 1, imaxG, jmaxG);
        int lFaceIn2_jl =
            GetLowerFaceJ(ii, jj, numGhosts_ + 1, imaxG, jmaxG);
        int lFaceIn2_ku =
            GetUpperFaceK(ii, jj, numGhosts_ + 1, imaxG, jmaxG);

        // location of lower k-boundary face
        int lFaceB = GetLowerFaceK(ii, jj, numGhosts_, imaxG, jmaxG);

        // Assign volume
        // -------------------------------------------------------------------
        // mirror volume values from adjacent cells across k-boundary
        // first layer of ghost cells
        (*this).vol_[cellLowG1] = vol_(cellLowIn1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellLowG2] = vol_(cellLowIn1);
        } else {
          (*this).vol_[cellLowG2] = vol_(cellLowIn2);
        }

        // Assign face areas
        // -------------------------------------------------------------------
        // mirror face area values from adjacent cells across k-boundary
        // first layer of ghost cells
        (*this).fAreaK_[lFaceG1_kl] = (*this).FAreaK(lFaceIn1_ku);
        (*this).fAreaI_[lFaceG1_il] = (*this).FAreaI(lFaceIn1_il);
        (*this).fAreaJ_[lFaceG1_jl] = (*this).FAreaJ(lFaceIn1_jl);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn1_ku);
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn1_il);
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn1_jl);
        } else {
          (*this).fAreaK_[lFaceG2_kl] = (*this).FAreaK(lFaceIn2_ku);
          (*this).fAreaI_[lFaceG2_il] = (*this).FAreaI(lFaceIn2_il);
          (*this).fAreaJ_[lFaceG2_jl] = (*this).FAreaJ(lFaceIn2_jl);
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at lower k-boundary
          int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG);
          int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

          // location of upper i-face for interior cells at lower k-boundary
          int lFaceIn1_iu =
              GetUpperFaceI(ii, jj, numGhosts_, imaxG, jmaxG);
          int lFaceIn2_iu =
              GetUpperFaceI(ii, jj, numGhosts_ + 1, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaI_[lFaceG1_iu] = (*this).FAreaI(lFaceIn1_iu);

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaI_[lFaceG2_iu] = (*this).FAreaI(lFaceIn1_iu);
          } else {
            (*this).fAreaI_[lFaceG2_iu] = (*this).FAreaI(lFaceIn2_iu);
          }
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at lower k-boundary
          int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG);
          int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

          // location of upper j-face for interior cells at lower k-boundary
          int lFaceIn1_ju =
              GetUpperFaceJ(ii, jj, numGhosts_, imaxG, jmaxG);
          int lFaceIn2_ju =
              GetUpperFaceJ(ii, jj, numGhosts_ + 1, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaJ_[lFaceG1_ju] = (*this).FAreaJ(lFaceIn1_ju);

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaJ_[lFaceG2_ju] = (*this).FAreaJ(lFaceIn1_ju);
          } else {
            (*this).fAreaJ_[lFaceG2_ju] = (*this).FAreaJ(lFaceIn2_ju);
          }
        }

        // Assign cell centroid
        // --------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku);
        (*this).center_[cellLowG1] = center_(cellLowIn1) + dist2Move;

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku));
          (*this).center_[cellLowG2] = center_(cellLowG1) + dist2Move;
        } else {
          dist2Move = fCenterK_(lFaceB) - fCenterK_(lFaceIn2_ku);
          (*this).center_[cellLowG2] = center_(cellLowIn1) + dist2Move;
        }

        // Assign face centers
        // --------------------------------------------------------------------
        // face center_ is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku);
        (*this).fCenterK_[lFaceG1_kl] = fCenterK_(lFaceB) + dist2Move;
        (*this).fCenterI_[lFaceG1_il] =
            fCenterI_(lFaceIn1_il) + dist2Move;
        (*this).fCenterJ_[lFaceG1_jl] =
            fCenterJ_(lFaceIn1_jl) + dist2Move;

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku));
          (*this).fCenterK_[lFaceG2_kl] =
              fCenterK_(lFaceG1_kl) + dist2Move;
          (*this).fCenterI_[lFaceG2_il] =
              fCenterI_(lFaceG1_il) + dist2Move;
          (*this).fCenterJ_[lFaceG2_jl] =
              fCenterJ_(lFaceG1_jl) + dist2Move;
        } else {
          dist2Move = fCenterK_(lFaceB) - fCenterK_(lFaceIn2_ku);
          (*this).fCenterK_[lFaceG2_kl] = fCenterK_(lFaceB) + dist2Move;
          (*this).fCenterI_[lFaceG2_il] =
              fCenterI_(lFaceIn1_il) + dist2Move;
          (*this).fCenterJ_[lFaceG2_jl] =
              fCenterJ_(lFaceIn1_jl) + dist2Move;
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at lower k-boundary
          int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG);
          int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

          // location of upper i-face for interior cells at lower k-boundary
          int lFaceIn1_iu =
              GetUpperFaceI(ii, jj, numGhosts_, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku);
          (*this).fCenterI_[lFaceG1_iu] =
              fCenterI_(lFaceIn1_iu) + dist2Move;

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterK_(lFaceB) -
                               fCenterK_(lFaceIn1_ku));
            (*this).fCenterI_[lFaceG2_iu] =
                fCenterI_(lFaceG1_iu) + dist2Move;
          } else {
            dist2Move =
                fCenterK_(lFaceB) - fCenterK_(lFaceIn2_ku);
            (*this).fCenterI_[lFaceG2_iu] =
                fCenterI_(lFaceIn1_iu) + dist2Move;
          }
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at lower k-boundary
          int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG);
          int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

          // location of upper j-face for interior cells at lower k-boundary
          int lFaceIn1_ju =
              GetUpperFaceJ(ii, jj, numGhosts_, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterK_(lFaceB) - fCenterK_(lFaceIn1_ku);
          (*this).fCenterJ_[lFaceG1_ju] =
              fCenterJ_(lFaceIn1_ju) + dist2Move;

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterK_(lFaceB) -
                               fCenterK_(lFaceIn1_ku));
            (*this).fCenterJ_[lFaceG2_ju] =
                fCenterJ_(lFaceG1_ju) + dist2Move;
          } else {
            dist2Move =
                fCenterK_(lFaceB) - fCenterK_(lFaceIn2_ku);
            (*this).fCenterJ_[lFaceG2_ju] =
                fCenterJ_(lFaceIn1_ju) + dist2Move;
          }
        }
      }

      // upper surface
      // -----------------------------------------------------------------
      if (bcNameU != "interblock") {  // only supply geometry values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of first ghost cell at upper k-boundary
        int cellUpG1 = GetLoc1D(ii, jj, kmaxG - 2, imaxG, jmaxG);
        // location of upper k and lower i, j faces of first ghost cell at lower
        // k-boundary
        int uFaceG1_il = GetLowerFaceI(ii, jj, kmaxG - 2, imaxG, jmaxG);
        int uFaceG1_jl = GetLowerFaceJ(ii, jj, kmaxG - 2, imaxG, jmaxG);
        int uFaceG1_ku = GetUpperFaceK(ii, jj, kmaxG - 2, imaxG, jmaxG);

        // location of second ghost cell at upper k-boundary
        int cellUpG2 = GetLoc1D(ii, jj, kmaxG - 1, imaxG, jmaxG);
        // location of upper k and lower i, j faces of second ghost cell at
        // lower k-boundary
        int uFaceG2_il = GetLowerFaceI(ii, jj, kmaxG - 1, imaxG, jmaxG);
        int uFaceG2_jl = GetLowerFaceJ(ii, jj, kmaxG - 1, imaxG, jmaxG);
        int uFaceG2_ku = GetUpperFaceK(ii, jj, kmaxG - 1, imaxG, jmaxG);

        // location of first interior cell at upper k-boundary
        int cellUpIn1 =
            GetLoc1D(ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);
        // location of lower i, j, k faces of first interior cell at lower
        // k-boundary
        int uFaceIn1_il = GetLowerFaceI(ii, jj, kmaxG - 1 - numGhosts_,
                                        imaxG, jmaxG);
        int uFaceIn1_jl = GetLowerFaceJ(ii, jj, kmaxG - 1 - numGhosts_,
                                        imaxG, jmaxG);
        int uFaceIn1_kl = GetLowerFaceK(ii, jj, kmaxG - 1 - numGhosts_,
                                        imaxG, jmaxG);

        // location of second interior cell at upper k-boundary
        int cellUpIn2 =
            GetLoc1D(ii, jj, kmaxG - 2 - numGhosts_, imaxG, jmaxG);
        // location of lower i, j, k faces of second interior cell at lower
        // k-boundary
        int uFaceIn2_il = GetLowerFaceI(ii, jj, kmaxG - 2 - numGhosts_,
                                        imaxG, jmaxG);
        int uFaceIn2_jl = GetLowerFaceJ(ii, jj, kmaxG - 2 - numGhosts_,
                                        imaxG, jmaxG);
        int uFaceIn2_kl = GetLowerFaceK(ii, jj, kmaxG - 2 - numGhosts_,
                                        imaxG, jmaxG);

        // location of upper k-boundary face
        int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - numGhosts_,
                                   imaxG, jmaxG);

        // Assign volume
        // --------------------------------------------------------------------
        // mirror volume values from adjacent cells across k-boundary
        // first layer of ghost cells
        (*this).vol_[cellUpG1] = vol_(cellUpIn1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).vol_[cellUpG2] = vol_(cellUpIn1);
        } else {
          (*this).vol_[cellUpG2] = vol_(cellUpIn2);
        }

        // Assign face areas
        // --------------------------------------------------------------------
        // mirror face area values from adjacent cells across k-boundary
        // first layer of ghost cells
        (*this).fAreaK_[uFaceG1_ku] = (*this).FAreaK(uFaceIn1_kl);
        (*this).fAreaI_[uFaceG1_il] = (*this).FAreaI(uFaceIn1_il);
        (*this).fAreaJ_[uFaceG1_jl] = (*this).FAreaJ(uFaceIn1_jl);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn1_kl);
          (*this).fAreaI_[uFaceG2_il] = (*this).FAreaI(uFaceIn1_il);
          (*this).fAreaJ_[uFaceG2_jl] = (*this).FAreaJ(uFaceIn1_jl);
        } else {
          (*this).fAreaK_[uFaceG2_ku] = (*this).FAreaK(uFaceIn2_kl);
          (*this).fAreaI_[uFaceG2_il] = (*this).FAreaI(uFaceIn2_il);
          (*this).fAreaJ_[uFaceG2_jl] = (*this).FAreaJ(uFaceIn2_jl);
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at upper k-boundary
          int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG - 2, imaxG, jmaxG);
          int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG - 1, imaxG, jmaxG);

          // location of upper i-face for interior cells at upper k-boundary
          int uFaceIn1_iu = GetUpperFaceI(
              ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);
          int uFaceIn2_iu = GetUpperFaceI(
              ii, jj, kmaxG - 2 - numGhosts_, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaI_[uFaceG1_iu] = (*this).FAreaI(uFaceIn1_iu);

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn1_iu);
          } else {
            (*this).fAreaI_[uFaceG2_iu] = (*this).FAreaI(uFaceIn2_iu);
          }
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at upper k-boundary
          int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG - 2, imaxG, jmaxG);
          int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG - 1, imaxG, jmaxG);

          // location of upper j-face for interior cells at upper k-boundary
          int uFaceIn1_ju = GetUpperFaceJ(
              ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);
          int uFaceIn2_ju = GetUpperFaceJ(
              ii, jj, kmaxG - 2 - numGhosts_, imaxG, jmaxG);

          // mirror face area values from adjacent cells
          // first layer of ghost cells
          (*this).fAreaJ_[uFaceG1_ju] = (*this).FAreaJ(uFaceIn1_ju);

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn1_ju);
          } else {
            (*this).fAreaJ_[uFaceG2_ju] = (*this).FAreaJ(uFaceIn2_ju);
          }
        }

        // Assign cell centroid
        // --------------------------------------------------------------------
        // cell centroid is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        vector3d<double> dist2Move =
            fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl);
        (*this).center_[cellUpG1] = center_(cellUpIn1) + dist2Move;

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl));
          (*this).center_[cellUpG2] = center_(cellUpG1) + dist2Move;
        } else {
          dist2Move = fCenterK_(uFaceB) - fCenterK_(uFaceIn2_kl);
          (*this).center_[cellUpG2] = center_(cellUpIn1) + dist2Move;
        }

        // Assign face centers
        // --------------------------------------------------------------------
        // face center_ is moved interior cell width in the boundary normal
        // direction
        // first layer of ghost cells
        dist2Move = fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl);
        (*this).fCenterK_[uFaceG1_ku] = fCenterK_(uFaceB) + dist2Move;
        (*this).fCenterI_[uFaceG1_il] =
            fCenterI_(uFaceIn1_il) + dist2Move;
        (*this).fCenterJ_[uFaceG1_jl] =
            fCenterJ_(uFaceIn1_jl) + dist2Move;

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          dist2Move =
              2.0 * (fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl));
          (*this).fCenterK_[uFaceG2_ku] =
              fCenterK_(uFaceG1_ku) + dist2Move;
          (*this).fCenterI_[uFaceG2_il] =
              fCenterI_(uFaceG1_il) + dist2Move;
          (*this).fCenterJ_[uFaceG2_jl] =
              fCenterJ_(uFaceG1_jl) + dist2Move;
        } else {
          dist2Move = fCenterK_(uFaceB) - fCenterK_(uFaceIn2_kl);
          (*this).fCenterK_[uFaceG2_ku] = fCenterK_(uFaceB) + dist2Move;
          (*this).fCenterI_[uFaceG2_il] =
              fCenterI_(uFaceIn1_il) + dist2Move;
          (*this).fCenterJ_[uFaceG2_jl] =
              fCenterJ_(uFaceIn1_jl) + dist2Move;
        }

        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // i-face areas too

          // location of upper i-face for ghost cells at upper k-boundary
          int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG - 2, imaxG, jmaxG);
          int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG - 1, imaxG, jmaxG);

          // location of upper i-face for interior cells at upper k-boundary
          int uFaceIn1_iu = GetUpperFaceI(
              ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl);
          (*this).fCenterI_[uFaceG1_iu] =
              fCenterI_(uFaceIn1_iu) + dist2Move;

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterK_(uFaceB) -
                               fCenterK_(uFaceIn1_kl));
            (*this).fCenterI_[uFaceG2_iu] =
                fCenterI_(uFaceG1_iu) + dist2Move;
          } else {
            dist2Move =
                fCenterK_(uFaceB) - fCenterK_(uFaceIn2_kl);
            (*this).fCenterI_[uFaceG2_iu] =
                fCenterI_(uFaceIn1_iu) + dist2Move;
          }
        }

        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // j-face areas too

          // location of upper j-face for ghost cells at upper k-boundary
          int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG - 2, imaxG, jmaxG);
          int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG - 1, imaxG, jmaxG);

          // location of upper j-face for interior cells at upper k-boundary
          int uFaceIn1_ju = GetUpperFaceJ(
              ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);

          // face center_ is moved interior cell width in the boundary normal
          // direction
          // first layer of ghost cells
          dist2Move = fCenterK_(uFaceB) - fCenterK_(uFaceIn1_kl);
          (*this).fCenterJ_[uFaceG1_ju] =
              fCenterJ_(uFaceIn1_ju) + dist2Move;

          // second layer of ghost cells
          if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
            dist2Move = 2.0 * (fCenterK_(uFaceB) -
                               fCenterK_(uFaceIn1_kl));
            (*this).fCenterJ_[uFaceG2_ju] =
                fCenterJ_(uFaceG1_ju) + dist2Move;
          } else {
            dist2Move =
                fCenterK_(uFaceB) - fCenterK_(uFaceIn2_kl);
            (*this).fCenterJ_[uFaceG2_ju] =
                fCenterJ_(uFaceIn1_ju) + dist2Move;
          }
        }
      }
    }
  }

  // fill ghost cell edge lines with geometric values
  // (*this).AssignGhostCellsGeomEdge(inp);
}

/* Member function to assign geometric quantities such as volume, face area,
cell centroid, and face center to ghost cells located on the 12 block edges.
Assumes AssignGhostCellsGeom has already been run.

           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J

In the above diagram the cells marked X represent physical cells. Cells marked
G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each
corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells
(1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells
it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to
the values of the ghost cells they tough at level "e". The values at edge cell
2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it
touches.
*/
void procBlock::AssignGhostCellsGeomEdge() {
  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * numGhosts_;
  int jmaxG = jmax + 2 * numGhosts_;
  int kmaxG = kmax + 2 * numGhosts_;

  //-------------------------------------------------------------------------
  // loop over edges at lower and upper j sides of block - this will include 4
  // edges that run in the i-direction -------------------------------
  // edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this
  // loop
  for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int j1, k1, j2, k2, je, ke;

      int gfe_j1_k1_jl, gfe_j1_k1_kl;
      int gfe_j1_k2_jl, gfe_j1_k2_kl;
      int gfe_j2_k1_jl, gfe_j2_k1_kl;
      int gfe_j2_k2_jl, gfe_j2_k2_kl;

      int gf_j1_ke_jl, gf_j1_ke_kl, gf_j1_ke_ku;
      int gf_j2_ke_jl, gf_j2_ke_kl;
      int gf_je_k1_jl, gf_je_k1_ju;
      int gf_je_k2_jl, gf_je_k2_kl;

      string bc_J, bc_K, surfJ, surfK;

      if (cc == 0) {  // at jl/kl edge - ghost cells are in the lower direction
                      // of both j and k, so use GetLowerFace for both
        j2 = 0;
        j1 = 1;
        je = numGhosts_;

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, j# = j
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of j line of cells, on first layer
        // of k line of cells
        gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG);
        gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of j line of cells, on second
        // layer of k line of cells
        gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG);
        gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on first
        // layer of k line of cells
        gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG);
        gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on second
        // layer of k line of cells
        gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG);
        gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG);

        // ghost face, on first layer of j line of cells, on non-edge layer of k
        // line of cells
        gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);

        // ghost face, on second layer of j line of cells, on non-edge layer of
        // k line of cells
        gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first layer of k
        // line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on second layer of
        // k line of cells
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);
        gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG);

        surfJ = "jl";
        surfK = "kl";

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at jl/ku edge - ghost cells are in the lower
                             // direction of j and upper direction of k, so use
                             // GetLowerFace for J
        j2 = 0;
        j1 = 1;
        je = numGhosts_;

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, j# = j
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of j line of cells, on first layer
        // of k line of cells
        gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG);
        gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of j line of cells, on second
        // layer of k line of cells
        gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG);
        gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on first
        // layer of k line of cells
        gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG);
        gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on second
        // layer of k line of cells
        gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG);
        gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG);

        // ghost face, on first layer of j line of cells, on non-edge layer of k
        // line of cells
        gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);

        // ghost face, on second layer of j line of cells, on non-edge layer of
        // k line of cells
        gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first layer of k
        // line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on second layer of
        // k line of cells
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);
        gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG);

        surfJ = "jl";
        surfK = "ku";

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at ju/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of j so use
                             // GetLowerFace for k
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, j# = j
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of j line of cells, on first layer
        // of k line of cells
        gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG);
        gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of j line of cells, on second
        // layer of k line of cells
        gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG);
        gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on first
        // layer of k line of cells
        gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG);
        gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on second
        // layer of k line of cells
        gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG);
        gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG);

        // ghost face, on first layer of j line of cells, on non-edge layer of k
        // line of cells
        gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);

        // ghost face, on second layer of j line of cells, on non-edge layer of
        // k line of cells
        gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first layer of k
        // line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on second layer of
        // k line of cells
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);
        gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG);

        surfJ = "ju";
        surfK = "kl";

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_ + 1,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at ju/ku edge - ghost cells are in the upper
                             // direction of both j and k, use GetUpperFace for
                             // both
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, j# = j
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of j line of cells, on first layer
        // of k line of cells
        gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG);
        gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of j line of cells, on second
        // layer of k line of cells
        gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG);
        gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on first
        // layer of k line of cells
        gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG);
        gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of j line of cells, on second
        // layer of k line of cells
        gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG);
        gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG);

        // ghost face, on first layer of j line of cells, on non-edge layer of k
        // line of cells
        gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);

        // ghost face, on second layer of j line of cells, on non-edge layer of
        // k line of cells
        gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first layer of k
        // line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on second layer of
        // k line of cells
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);
        gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG);

        surfJ = "ju";
        surfK = "ku";

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_ + 1,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if bouth
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_J != "interblock" && bc_K != "interblock") {
        // cell indices and remaining face indices
        int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
        int gfe_j1_k1_il = GetLowerFaceI(
            ii, j1, k1, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // j line of cells, on first layer of k
                                        // line of cells

        int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of j line of cells, on second layer
                                          // of k line of cells
        int gfe_j1_k2_il = GetLowerFaceI(
            ii, j1, k2, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // j line of cells, on second layer of k
                                        // line of cells

        int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
        int gfe_j2_k1_il = GetLowerFaceI(
            ii, j2, k1, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of j line of cells, on first layer of
                                        // k line of cells

        int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of j line of cells, on second layer
                                          // of k line of cells
        int gfe_j2_k2_il = GetLowerFaceI(
            ii, j2, k2, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of j line of cells, on second layer
                                        // of k line of cells

        int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of j
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gf_j1_ke_il = GetLowerFaceI(
            ii, j1, ke, imaxG, jmaxG);  // ghost face, on first layer of j line
                                        // of cells, on non-edge layer of k line
                                        // of cells

        int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of j
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gf_j2_ke_il = GetLowerFaceI(
            ii, j2, ke, imaxG, jmaxG);  // ghost face, on second layer of j line
                                        // of cells, on non-edge layer of k line
                                        // of cells

        int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // j line of cells,
                                                            // on first layer of
                                                            // k line of cells
        int gf_je_k1_il = GetLowerFaceI(
            ii, je, k1, imaxG, jmaxG);  // ghost face, on non-edge layer of j
                                        // line of cells, on first layer of k
                                        // line of cells

        int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of j
                                         // line of cells, on second layer of k
                                         // line of cells
        int gf_je_k2_il = GetLowerFaceI(
            ii, je, k2, imaxG, jmaxG);  // ghost face, on non-edge layer of j
                                        // line of cells, on second layer of k
                                        // line of cells

        // volume
        // ----------------------------------------------------------------
        (*this).vol_[gce_j1_k1] =
            0.5 * (vol_(gc_j1_ke) + vol_(gc_je_k1));
        (*this).vol_[gce_j2_k1] = vol_(gc_j2_ke);
        (*this).vol_[gce_j1_k2] = vol_(gc_je_k2);
        (*this).vol_[gce_j2_k2] =
            0.5 * (vol_(gc_j2_ke) + vol_(gc_je_k2));

        // face areas
        // ----------------------------------------------------------------
        (*this).fAreaI_[gfe_j1_k1_il] =
            0.5 * ((*this).FAreaI(gf_je_k1_il) + (*this).FAreaI(gf_j1_ke_il));
        (*this).fAreaJ_[gfe_j1_k1_jl] = (*this).FAreaJ(gf_je_k1_jl);
        (*this).fAreaK_[gfe_j1_k1_kl] = (*this).FAreaK(gf_j1_ke_kl);

        (*this).fAreaI_[gfe_j1_k2_il] = (*this).FAreaI(gf_je_k2_il);
        (*this).fAreaJ_[gfe_j1_k2_jl] = (*this).FAreaJ(gf_je_k2_jl);
        (*this).fAreaK_[gfe_j1_k2_kl] = (*this).FAreaK(gf_j1_ke_kl);

        (*this).fAreaI_[gfe_j2_k1_il] = (*this).FAreaI(gf_j2_ke_il);
        (*this).fAreaJ_[gfe_j2_k1_jl] = (*this).FAreaJ(gf_je_k1_jl);
        (*this).fAreaK_[gfe_j2_k1_kl] = (*this).FAreaK(gf_j2_ke_kl);

        (*this).fAreaI_[gfe_j2_k2_il] =
            0.5 * ((*this).FAreaI(gf_j2_ke_il) + (*this).FAreaI(gf_je_k2_il));
        (*this).fAreaJ_[gfe_j2_k2_jl] = (*this).FAreaJ(gf_je_k2_jl);
        (*this).fAreaK_[gfe_j2_k2_kl] = (*this).FAreaK(gf_j2_ke_kl);

        // centroids
        // ----------------------------------------------------------------
        // edge centroid is moved distance of cell width normal to face dividing
        // regular and edge ghost cells
        vector3d<double> dist2MoveK =
            fCenterK_(gf_j1_ke_kl) - fCenterK_(gf_j1_ke_ku);
        vector3d<double> dist2MoveJ =
            fCenterJ_(gf_je_k1_jl) - fCenterJ_(gf_je_k1_ju);
        (*this).center_[gce_j1_k1] = center_(gc_j1_ke) + dist2MoveK;
        (*this).center_[gce_j2_k1] = center_(gc_j2_ke) + dist2MoveK;
        (*this).center_[gce_j1_k2] = center_(gc_je_k2) + dist2MoveJ;
        (*this).center_[gce_j2_k2] =
            center_(gc_je_k2) + 2.0 * dist2MoveJ;

        // face centers
        // ----------------------------------------------------------------
        // edge face centers are moved distance of cell width normal to face
        // dividing regular and edge ghost cells
        (*this).fCenterI_[gfe_j1_k1_il] =
            fCenterI_(gf_j1_ke_il) + dist2MoveK;
        (*this).fCenterJ_[gfe_j1_k1_jl] =
            fCenterJ_(gf_j1_ke_jl) + dist2MoveK;
        (*this).fCenterK_[gfe_j1_k1_kl] =
            fCenterK_(gf_j1_ke_kl) + dist2MoveK;

        (*this).fCenterI_[gfe_j2_k1_il] =
            fCenterI_(gf_j2_ke_il) + dist2MoveK;
        (*this).fCenterJ_[gfe_j2_k1_jl] =
            fCenterJ_(gf_j2_ke_jl) + dist2MoveK;
        (*this).fCenterK_[gfe_j2_k1_kl] =
            fCenterK_(gf_j2_ke_kl) + dist2MoveK;

        (*this).fCenterI_[gfe_j1_k2_il] =
            fCenterI_(gf_je_k2_il) + dist2MoveJ;
        (*this).fCenterJ_[gfe_j1_k2_jl] =
            fCenterJ_(gf_je_k2_jl) + dist2MoveJ;
        (*this).fCenterK_[gfe_j1_k2_kl] =
            fCenterK_(gf_je_k2_kl) + dist2MoveJ;

        (*this).fCenterI_[gfe_j2_k2_il] =
            fCenterI_(gf_je_k2_il) + 2.0 * dist2MoveJ;
        (*this).fCenterJ_[gfe_j2_k2_jl] =
            fCenterJ_(gf_je_k2_jl) + 2.0 * dist2MoveJ;
        (*this).fCenterK_[gfe_j2_k2_kl] =
            fCenterK_(gf_je_k2_kl) + 2.0 * dist2MoveJ;

        // this is only done at the end of the i loop
        if (ii == imax - 1 + numGhosts_) {  // at end of i-line of
                                                     // cells assign cell upper
                                                     // face areas too

          int gfe_j1_k1_2il = GetLowerFaceI(
              ii - 1, j1, k1, imaxG, jmaxG);  // ghost face on edge, on first
                                              // layer of j line of cells, on
                                              // first layer of k line of cells
          vector3d<double> dist2MoveI =
              fCenterI_(gfe_j1_k1_il) -
              (*this)
                  .FCenterI(gfe_j1_k1_2il);  // i-width of cell adjacent to edge

          int gfe_j1_k1_iu = GetUpperFaceI(
              ii, j1, k1, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
          int gfe_j1_k2_iu = GetUpperFaceI(
              ii, j1, k2, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of j line of cells, on second layer
                                          // of k line of cells
          int gfe_j2_k1_iu = GetUpperFaceI(
              ii, j2, k1, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
          int gfe_j2_k2_iu = GetUpperFaceI(
              ii, j2, k2, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of j line of cells, on second layer
                                          // of k line of cells

          // face areas
          (*this).fAreaI_[gfe_j1_k1_iu] = (*this).FAreaI(gfe_j1_k1_il);
          (*this).fAreaI_[gfe_j1_k2_iu] = (*this).FAreaI(gfe_j1_k2_il);
          (*this).fAreaI_[gfe_j2_k1_iu] = (*this).FAreaI(gfe_j2_k1_il);
          (*this).fAreaI_[gfe_j2_k2_iu] = (*this).FAreaI(gfe_j2_k2_il);

          // face centers
          (*this).fCenterI_[gfe_j1_k1_iu] =
              fCenterI_(gfe_j1_k1_il) + dist2MoveI;
          (*this).fCenterI_[gfe_j1_k2_iu] =
              fCenterI_(gfe_j1_k2_il) + dist2MoveI;
          (*this).fCenterI_[gfe_j2_k1_iu] =
              fCenterI_(gfe_j2_k1_il) + dist2MoveI;
          (*this).fCenterI_[gfe_j2_k2_iu] =
              fCenterI_(gfe_j2_k2_il) + dist2MoveI;
        }
      }
    }
  }

  //-------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the j-direction --------------------------------
  // edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this
  // loop
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, k1, i2, k2, ie, ke;

      int gfe_i1_k1_il, gfe_i1_k1_kl;
      int gfe_i1_k2_il, gfe_i1_k2_kl;
      int gfe_i2_k1_il, gfe_i2_k1_kl;
      int gfe_i2_k2_il, gfe_i2_k2_kl;

      int gf_i1_ke_il, gf_i1_ke_ku, gf_i1_ke_kl;
      int gf_i2_ke_il, gf_i2_ke_kl;
      int gf_ie_k1_il, gf_ie_k1_iu;
      int gf_ie_k2_il, gf_ie_k2_kl;

      string bc_I, bc_K, surfI, surfK;

      if (cc == 0) {  // at il/kl edge - ghost cells are in the lower direction
                      // of both i and k, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of k line of cells
        gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG);
        gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of k line of cells
        gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG);
        gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of k line of cells
        gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG);
        gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of k line of cells
        gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG);
        gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of k
        // line of cells
        gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // k line of cells
        gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of k
        // line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // k line of cells
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);
        gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);

        surfI = "il";
        surfK = "kl";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at il/ku edge - ghost cells are in the lower
                             // direction of i and upper direction of k, so use
                             // GetLowerFace for J
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of k line of cells
        gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG);
        gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of k line of cells
        gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG);
        gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of k line of cells
        gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG);
        gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of k line of cells
        gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG);
        gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of k
        // line of cells
        gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // k line of cells
        gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of k
        // line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // k line of cells
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);
        gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG);

        surfI = "il";
        surfK = "ku";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at iu/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of i so use
                             // GetLowerFace for k
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of k line of cells
        gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG);
        gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of k line of cells
        gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG);
        gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of k line of cells
        gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG);
        gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of k line of cells
        gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG);
        gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of k
        // line of cells
        gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // k line of cells
        gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of k
        // line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // k line of cells
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);
        gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);

        surfI = "iu";
        surfK = "kl";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at iu/ku edge - ghost cells are in the upper
                             // direction of both i and k, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, k# = k position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of k line of cells
        gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG);
        gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of k line of cells
        gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG);
        gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of k line of cells
        gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG);
        gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of k line of cells
        gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG);
        gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of k
        // line of cells
        gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // k line of cells
        gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of k
        // line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // k line of cells
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);
        gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG);

        surfI = "iu";
        surfK = "ku";

        // boundary conditioins at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if bouth
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_I != "interblock" && bc_K != "interblock") {
        // cell indices and remaining face indices
        int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gfe_i1_k1_jl = GetLowerFaceJ(
            i1, jj, k1, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // i line of cells, on first layer of k
                                        // line of cells

        int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of k line of cells
        int gfe_i1_k2_jl = GetLowerFaceJ(
            i1, jj, k2, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // i line of cells, on second layer of k
                                        // line of cells

        int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gfe_i2_k1_jl = GetLowerFaceJ(
            i2, jj, k1, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of i line of cells, on first layer of
                                        // k line of cells

        int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of k line of cells
        int gfe_i2_k2_jl = GetLowerFaceJ(
            i2, jj, k2, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of i line of cells, on second layer
                                        // of k line of cells

        int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gf_i1_ke_jl = GetLowerFaceJ(
            i1, jj, ke, imaxG, jmaxG);  // ghost face, on first layer of i line
                                        // of cells, on non-edge layer of k line
                                        // of cells

        int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gf_i2_ke_jl = GetLowerFaceJ(
            i2, jj, ke, imaxG, jmaxG);  // ghost face, on second layer of i line
                                        // of cells, on non-edge layer of k line
                                        // of cells

        int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // i line of cells,
                                                            // on first layer of
                                                            // k line of cells
        int gf_ie_k1_jl = GetLowerFaceJ(
            ie, jj, k1, imaxG, jmaxG);  // ghost face, on non-edge layer of i
                                        // line of cells, on first layer of k
                                        // line of cells

        int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of i
                                         // line of cells, on second layer of k
                                         // line of cells
        int gf_ie_k2_jl = GetLowerFaceJ(
            ie, jj, k2, imaxG, jmaxG);  // ghost face, on non-edge layer of i
                                        // line of cells, on second layer of k
                                        // line of cells

        // volume
        // ------------------------------------------------------------------
        (*this).vol_[gce_i1_k1] =
            0.5 * (vol_(gc_i1_ke) + vol_(gc_ie_k1));
        (*this).vol_[gce_i2_k1] = vol_(gc_i2_ke);
        (*this).vol_[gce_i1_k2] = vol_(gc_ie_k2);
        (*this).vol_[gce_i2_k2] =
            0.5 * (vol_(gc_i2_ke) + vol_(gc_ie_k2));

        // face areas
        // ------------------------------------------------------------------
        (*this).fAreaJ_[gfe_i1_k1_jl] =
            0.5 * ((*this).FAreaJ(gf_ie_k1_jl) + (*this).FAreaJ(gf_i1_ke_jl));
        (*this).fAreaI_[gfe_i1_k1_il] = (*this).FAreaI(gf_ie_k1_il);
        (*this).fAreaK_[gfe_i1_k1_kl] = (*this).FAreaK(gf_i1_ke_kl);

        (*this).fAreaJ_[gfe_i1_k2_jl] = (*this).FAreaJ(gf_ie_k2_jl);
        (*this).fAreaI_[gfe_i1_k2_il] = (*this).FAreaI(gf_ie_k2_il);
        (*this).fAreaK_[gfe_i1_k2_kl] = (*this).FAreaK(gf_i1_ke_kl);

        (*this).fAreaJ_[gfe_i2_k1_jl] = (*this).FAreaJ(gf_i2_ke_jl);
        (*this).fAreaI_[gfe_i2_k1_il] = (*this).FAreaI(gf_ie_k1_il);
        (*this).fAreaK_[gfe_i2_k1_kl] = (*this).FAreaK(gf_i2_ke_kl);

        (*this).fAreaJ_[gfe_i2_k2_jl] =
            0.5 * ((*this).FAreaJ(gf_i2_ke_jl) + (*this).FAreaJ(gf_ie_k2_jl));
        (*this).fAreaI_[gfe_i2_k2_il] = (*this).FAreaI(gf_ie_k2_il);
        (*this).fAreaK_[gfe_i2_k2_kl] = (*this).FAreaK(gf_i2_ke_kl);

        // centroids
        // ------------------------------------------------------------------
        // edge centroid is moved distance of cell width normal to face dividing
        // regular and edge ghost cells
        vector3d<double> dist2MoveK =
            fCenterK_(gf_i1_ke_kl) - fCenterK_(gf_i1_ke_ku);
        vector3d<double> dist2MoveI =
            fCenterI_(gf_ie_k1_il) - fCenterI_(gf_ie_k1_iu);
        (*this).center_[gce_i1_k1] = center_(gc_i1_ke) + dist2MoveK;
        (*this).center_[gce_i2_k1] = center_(gc_i2_ke) + dist2MoveK;
        (*this).center_[gce_i1_k2] = center_(gc_ie_k2) + dist2MoveI;
        (*this).center_[gce_i2_k2] =
            center_(gc_ie_k2) + 2.0 * dist2MoveI;

        // face centers
        // -----------------------------------------------------------------
        // edge face centers are moved distance of cell width normal to face
        // dividing regular and edge ghost cells
        (*this).fCenterI_[gfe_i1_k1_il] =
            fCenterI_(gf_i1_ke_il) + dist2MoveK;
        (*this).fCenterJ_[gfe_i1_k1_jl] =
            fCenterJ_(gf_i1_ke_jl) + dist2MoveK;
        (*this).fCenterK_[gfe_i1_k1_kl] =
            fCenterK_(gf_i1_ke_kl) + dist2MoveK;

        (*this).fCenterI_[gfe_i2_k1_il] =
            fCenterI_(gf_i2_ke_il) + dist2MoveK;
        (*this).fCenterJ_[gfe_i2_k1_jl] =
            fCenterJ_(gf_i2_ke_jl) + dist2MoveK;
        (*this).fCenterK_[gfe_i2_k1_kl] =
            fCenterK_(gf_i2_ke_kl) + dist2MoveK;

        (*this).fCenterI_[gfe_i1_k2_il] =
            fCenterI_(gf_ie_k2_il) + dist2MoveI;
        (*this).fCenterJ_[gfe_i1_k2_jl] =
            fCenterJ_(gf_ie_k2_jl) + dist2MoveI;
        (*this).fCenterK_[gfe_i1_k2_kl] =
            fCenterK_(gf_ie_k2_kl) + dist2MoveI;

        (*this).fCenterI_[gfe_i2_k2_il] =
            fCenterI_(gf_ie_k2_il) + 2.0 * dist2MoveI;
        (*this).fCenterJ_[gfe_i2_k2_jl] =
            fCenterJ_(gf_ie_k2_jl) + 2.0 * dist2MoveI;
        (*this).fCenterK_[gfe_i2_k2_kl] =
            fCenterK_(gf_ie_k2_kl) + 2.0 * dist2MoveI;

        // this is only done at the end of the j loop
        if (jj == jmax - 1 + numGhosts_) {  // at end of j-line of
                                                     // cells assign cell upper
                                                     // face areas too

          int gfe_i1_k1_2jl = GetLowerFaceJ(
              i1, jj - 1, k1, imaxG, jmaxG);  // ghost face on edge, on first
                                              // layer of j line of cells, on
                                              // first layer of k line of cells
          vector3d<double> dist2MoveJ =
              fCenterJ_(gfe_i1_k1_jl) -
              fCenterJ_(gfe_i1_k1_2jl);  // j-width of adjacent cell

          int gfe_i1_k1_ju = GetUpperFaceJ(
              i1, jj, k1, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
          int gfe_i1_k2_ju = GetUpperFaceJ(
              i1, jj, k2, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of k line of cells
          int gfe_i2_k1_ju = GetUpperFaceJ(
              i2, jj, k1, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
          int gfe_i2_k2_ju = GetUpperFaceJ(
              i2, jj, k2, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of k line of cells

          // face areas
          (*this).fAreaJ_[gfe_i1_k1_ju] = (*this).FAreaJ(gfe_i1_k1_jl);
          (*this).fAreaJ_[gfe_i1_k2_ju] = (*this).FAreaJ(gfe_i1_k2_jl);
          (*this).fAreaJ_[gfe_i2_k1_ju] = (*this).FAreaJ(gfe_i2_k1_jl);
          (*this).fAreaJ_[gfe_i2_k2_ju] = (*this).FAreaJ(gfe_i2_k2_jl);

          // face centers
          (*this).fCenterJ_[gfe_i1_k1_ju] =
              fCenterJ_(gfe_i1_k1_jl) + dist2MoveJ;
          (*this).fCenterJ_[gfe_i1_k2_ju] =
              fCenterJ_(gfe_i1_k2_jl) + dist2MoveJ;
          (*this).fCenterJ_[gfe_i2_k1_ju] =
              fCenterJ_(gfe_i2_k1_jl) + dist2MoveJ;
          (*this).fCenterJ_[gfe_i2_k2_ju] =
              fCenterJ_(gfe_i2_k2_jl) + dist2MoveJ;
        }
      }
    }
  }

  //-------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the k-direction --------------------------------
  // edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this
  // loop
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, j1, i2, j2, ie, je;

      int gfe_i1_j1_il, gfe_i1_j1_jl;
      int gfe_i1_j2_il, gfe_i1_j2_jl;
      int gfe_i2_j1_il, gfe_i2_j1_jl;
      int gfe_i2_j2_il, gfe_i2_j2_jl;

      int gf_i1_je_il, gf_i1_je_jl, gf_i1_je_ju;
      int gf_i2_je_il, gf_i2_je_jl;
      int gf_ie_j1_il, gf_ie_j1_iu;
      int gf_ie_j2_il, gf_ie_j2_jl;

      string bc_I, bc_J, surfI, surfJ;

      if (cc == 0) {  // at il/jl edge - ghost cells are in the lower direction
                      // of both i and j, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;

        j2 = 0;
        j1 = 1;
        je = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, j# = j position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of j line of cells
        gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG);
        gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of j line of cells
        gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG);
        gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of j line of cells
        gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG);
        gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of j line of cells
        gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG);
        gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of j
        // line of cells
        gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // j line of cells
        gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of j
        // line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // j line of cells
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);
        gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG);

        surfI = "il";
        surfJ = "jl";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 1) {  // at il/ju edge - ghost cells are in the lower
                             // direction of i and upper direction of j, so use
                             // GetLowerFace for I
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, j# = j position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of j line of cells
        gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG);
        gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of j line of cells
        gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG);
        gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of j line of cells
        gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG);
        gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of j line of cells
        gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG);
        gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of j
        // line of cells
        gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // j line of cells
        gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of j
        // line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // j line of cells
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);
        gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);

        surfI = "il";
        surfJ = "ju";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_ + 1,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 2) {  // at iu/jl edge - ghost cells are in the upper
                             // direction of i, and lower direction of j so use
                             // GetLowerFace for J
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;

        j2 = 0;
        j1 = 1;
        je = numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, j# = j position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of j line of cells
        gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG);
        gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of j line of cells
        gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG);
        gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of j line of cells
        gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG);
        gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of j line of cells
        gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG);
        gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of j
        // line of cells
        gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // j line of cells
        gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of j
        // line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // j line of cells
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);
        gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG);

        surfI = "iu";
        surfJ = "jl";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 3) {  // at iu/ju edge - ghost cells are in the upper
                             // direction of both i and j, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;

        // ghost edge face indices
        // naming convention - g = ghost, f = face, e = edge ghost cell, i# = i
        // position of location, j# = j position of location, jl = j lower face
        // ghost face on edge, on first layer of i line of cells, on first layer
        // of j line of cells
        gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG);
        gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on first layer of i line of cells, on second
        // layer of j line of cells
        gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG);
        gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on first
        // layer of j line of cells
        gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG);
        gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG);

        // ghost face on edge, on second layer of i line of cells, on second
        // layer of j line of cells
        gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG);
        gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG);

        // ghost face, on first layer of i line of cells, on non-edge layer of j
        // line of cells
        gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);

        // ghost face, on second layer of i line of cells, on non-edge layer of
        // j line of cells
        gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first layer of j
        // line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on second layer of
        // j line of cells
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);
        gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);

        surfI = "iu";
        surfJ = "ju";

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_ + 1,
                                      kk - numGhosts_, surfJ);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if bouth
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_I != "interblock" && bc_J != "interblock") {

        // ghost cell and remaining face indices
        int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of j line of cells
        int gfe_i1_j1_kl = GetLowerFaceK(
            i1, j1, kk, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // i line of cells, on first layer of j
                                        // line of cells

        int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of j line of cells
        int gfe_i1_j2_kl = GetLowerFaceK(
            i1, j2, kk, imaxG, jmaxG);  // ghost face on edge, on first layer of
                                        // i line of cells, on second layer of j
                                        // line of cells

        int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of j line of cells
        int gfe_i2_j1_kl = GetLowerFaceK(
            i2, j1, kk, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of i line of cells, on first layer of
                                        // j line of cells

        int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of j line of cells
        int gfe_i2_j2_kl = GetLowerFaceK(
            i2, j2, kk, imaxG, jmaxG);  // ghost face on edge, on second layer
                                        // of i line of cells, on second layer
                                        // of j line of cells

        int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // j line of cells
        int gf_i1_je_kl = GetLowerFaceK(
            i1, je, kk, imaxG, jmaxG);  // ghost face, on first layer of i line
                                        // of cells, on non-edge layer of j line
                                        // of cells

        int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // j line of cells
        int gf_i2_je_kl = GetLowerFaceK(
            i2, je, kk, imaxG, jmaxG);  // ghost face, on second layer of i line
                                        // of cells, on non-edge layer of j line
                                        // of cells

        int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // i line of cells,
                                                            // on first layer of
                                                            // j line of cells
        int gf_ie_j1_kl = GetLowerFaceK(
            ie, j1, kk, imaxG, jmaxG);  // ghost face, on non-edge layer of i
                                        // line of cells, on first layer of j
                                        // line of cells

        int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of i
                                         // line of cells, on second layer of j
                                         // line of cells
        int gf_ie_j2_kl = GetLowerFaceK(
            ie, j2, kk, imaxG, jmaxG);  // ghost face, on non-edge layer of i
                                        // line of cells, on second layer of j
                                        // line of cells

        // volume
        // --------------------------------------------------------------------
        (*this).vol_[gce_i1_j1] =
            0.5 * (vol_(gc_i1_je) + vol_(gc_ie_j1));
        (*this).vol_[gce_i2_j1] = vol_(gc_i2_je);
        (*this).vol_[gce_i1_j2] = vol_(gc_ie_j2);
        (*this).vol_[gce_i2_j2] =
            0.5 * (vol_(gc_i2_je) + vol_(gc_ie_j2));

        // face areas
        // --------------------------------------------------------------------
        (*this).fAreaK_[gfe_i1_j1_kl] =
            0.5 * ((*this).FAreaK(gf_ie_j1_kl) + (*this).FAreaK(gf_i1_je_kl));
        (*this).fAreaI_[gfe_i1_j1_il] = (*this).FAreaI(gf_ie_j1_il);
        (*this).fAreaJ_[gfe_i1_j1_jl] = (*this).FAreaJ(gf_i1_je_jl);

        (*this).fAreaK_[gfe_i1_j2_kl] = (*this).FAreaK(gf_ie_j2_kl);
        (*this).fAreaI_[gfe_i1_j2_il] = (*this).FAreaI(gf_ie_j2_il);
        (*this).fAreaJ_[gfe_i1_j2_jl] = (*this).FAreaJ(gf_i1_je_jl);

        (*this).fAreaK_[gfe_i2_j1_kl] = (*this).FAreaK(gf_i2_je_kl);
        (*this).fAreaI_[gfe_i2_j1_il] = (*this).FAreaI(gf_ie_j1_il);
        (*this).fAreaJ_[gfe_i2_j1_jl] = (*this).FAreaJ(gf_i2_je_jl);

        (*this).fAreaK_[gfe_i2_j2_kl] =
            0.5 * ((*this).FAreaK(gf_i2_je_kl) + (*this).FAreaK(gf_ie_j2_kl));
        (*this).fAreaI_[gfe_i2_j2_il] = (*this).FAreaI(gf_ie_j2_il);
        (*this).fAreaJ_[gfe_i2_j2_jl] = (*this).FAreaJ(gf_i2_je_jl);

        // centroids
        // --------------------------------------------------------------------
        // edge centroid is moved distance of cell width normal to face dividing
        // regular and edge ghost cells
        vector3d<double> dist2MoveJ =
            fCenterJ_(gf_i1_je_jl) - fCenterJ_(gf_i1_je_ju);
        vector3d<double> dist2MoveI =
            fCenterI_(gf_ie_j1_il) - fCenterI_(gf_ie_j1_iu);
        (*this).center_[gce_i1_j1] = center_(gc_i1_je) + dist2MoveJ;
        (*this).center_[gce_i2_j1] = center_(gc_i2_je) + dist2MoveJ;
        (*this).center_[gce_i1_j2] = center_(gc_ie_j2) + dist2MoveI;
        (*this).center_[gce_i2_j2] =
            center_(gc_ie_j2) + 2.0 * dist2MoveI;

        // face centers
        // --------------------------------------------------------------------
        // edge face centers are moved distance of cell width normal to face
        // dividing regular and edge ghost cells
        (*this).fCenterI_[gfe_i1_j1_il] =
            fCenterI_(gf_i1_je_il) + dist2MoveJ;
        (*this).fCenterJ_[gfe_i1_j1_jl] =
            fCenterJ_(gf_i1_je_jl) + dist2MoveJ;
        (*this).fCenterK_[gfe_i1_j1_kl] =
            fCenterK_(gf_i1_je_kl) + dist2MoveJ;

        (*this).fCenterI_[gfe_i2_j1_il] =
            fCenterI_(gf_i2_je_il) + dist2MoveJ;
        (*this).fCenterJ_[gfe_i2_j1_jl] =
            fCenterJ_(gf_i2_je_jl) + dist2MoveJ;
        (*this).fCenterK_[gfe_i2_j1_kl] =
            fCenterK_(gf_i2_je_kl) + dist2MoveJ;

        (*this).fCenterI_[gfe_i1_j2_il] =
            fCenterI_(gf_ie_j2_il) + dist2MoveI;
        (*this).fCenterJ_[gfe_i1_j2_jl] =
            fCenterJ_(gf_ie_j2_jl) + dist2MoveI;
        (*this).fCenterK_[gfe_i1_j2_kl] =
            fCenterK_(gf_ie_j2_kl) + dist2MoveI;

        (*this).fCenterI_[gfe_i2_j2_il] =
            fCenterI_(gf_ie_j2_il) + 2.0 * dist2MoveI;
        (*this).fCenterJ_[gfe_i2_j2_jl] =
            fCenterJ_(gf_ie_j2_jl) + 2.0 * dist2MoveI;
        (*this).fCenterK_[gfe_i2_j2_kl] =
            fCenterK_(gf_ie_j2_kl) + 2.0 * dist2MoveI;

        // this is only done at the end of the k loop
        if (kk == kmax - 1 + numGhosts_) {  // at end of k-line of
                                                     // cells assign cell upper
                                                     // face areas too

          int gfe_i1_j1_2kl = GetLowerFaceK(
              i1, j1, kk - 1, imaxG, jmaxG);  // ghost face on edge, on first
                                              // layer of j line of cells, on
                                              // first layer of k line of cells
          vector3d<double> dist2MoveK =
              fCenterK_(gfe_i1_j1_kl) -
              fCenterK_(gfe_i1_j1_2kl);  // k-width of adjacent cell

          int gfe_i1_j1_ku = GetUpperFaceK(
              i1, j1, kk, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of j line of cells
          int gfe_i1_j2_ku = GetUpperFaceK(
              i1, j2, kk, imaxG, jmaxG);  // ghost face on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of j line of cells
          int gfe_i2_j1_ku = GetUpperFaceK(
              i2, j1, kk, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of j line of cells
          int gfe_i2_j2_ku = GetUpperFaceK(
              i2, j2, kk, imaxG, jmaxG);  // ghost face on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of j line of cells

          // face areas
          (*this).fAreaK_[gfe_i1_j1_ku] = (*this).FAreaK(gfe_i1_j1_kl);
          (*this).fAreaK_[gfe_i1_j2_ku] = (*this).FAreaK(gfe_i1_j2_kl);
          (*this).fAreaK_[gfe_i2_j1_ku] = (*this).FAreaK(gfe_i2_j1_kl);
          (*this).fAreaK_[gfe_i2_j2_ku] = (*this).FAreaK(gfe_i2_j2_kl);

          // face centers
          (*this).fCenterK_[gfe_i1_j1_ku] =
              fCenterK_(gfe_i1_j1_kl) + dist2MoveK;
          (*this).fCenterK_[gfe_i1_j2_ku] =
              fCenterK_(gfe_i1_j2_kl) + dist2MoveK;
          (*this).fCenterK_[gfe_i2_j1_ku] =
              fCenterK_(gfe_i2_j1_kl) + dist2MoveK;
          (*this).fCenterK_[gfe_i2_j2_ku] =
              fCenterK_(gfe_i2_j2_kl) + dist2MoveK;
        }
      }
    }
  }
}

/* Member function to assign values for ghost cells for the inviscid flux
calculation. This function assigns values for regular ghost cells and "edge"
ghost cells. "Corner" cells are left with no value as they are not used.
           ____ ____ ____ ____ ____ ____ ____ ____
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|

In the above diagram where X represents the physical cells, cells marked G
(regular ghost cells) and E ("edge" ghost cells) are assigned geometric
values. G1 represents the first layer of ghost cells and G2 represents the
second layer.
*/
void procBlock::AssignInviscidGhostCells(const input &inp,
                                         const idealGas &eos,
                                         const sutherland &suth) {
  // inp -- all input variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity

  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * numGhosts_;
  int jmaxG = (*this).NumJ() + 2 * numGhosts_;
  int kmaxG = (*this).NumK() + 2 * numGhosts_;

  //------------------------------------------------------------------------
  // loop over physical I faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(0, jj - numGhosts_,
                                              kk - numGhosts_, "il");
      string bcNameU = (*this).BC().GetBCName(imax, jj - numGhosts_,
                                              kk - numGhosts_, "iu");

      // inviscid fluxes require different bc_ than viscous fluxes - treat all
      // walls as the same
      if (bcNameL == "viscousWall") {
        bcNameL = "slipWall";
      }
      if (bcNameU == "viscousWall") {
        bcNameU = "slipWall";
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at lower i-boundary
        int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
        int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

        // location of interior cells at lower i-boundary
        int cellLowIn1 = GetLoc1D(numGhosts_, jj, kk, imaxG, jmaxG);
        int cellLowIn2 =
            GetLoc1D(numGhosts_ + 1, jj, kk, imaxG, jmaxG);

        // location of lower i-boundary face
        int lFaceB = GetLowerFaceI(numGhosts_, jj, kk, imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitI(lFaceB), "il", inp, eos, suth, 1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          if (bcNameL == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
                bcNameL, (*this).FAreaUnitI(lFaceB), "il", inp, eos, suth, 1);
          } else {
            (*this).state_[cellLowG2] = state_(cellLowIn1).GetGhostState(
                bcNameL, (*this).FAreaUnitI(lFaceB), "il", inp, eos, suth, 2);
          }
        }
      }

      // upper surface
      // --------------------------------------------------------------------
      if (bcNameU != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at upper i-boundary
        int cellUpG1 = GetLoc1D(imaxG - 2, jj, kk, imaxG, jmaxG);
        int cellUpG2 = GetLoc1D(imaxG - 1, jj, kk, imaxG, jmaxG);

        // location of interior cells at upper i-boundary
        int cellUpIn1 =
            GetLoc1D(imaxG - 1 - numGhosts_, jj, kk, imaxG, jmaxG);
        int cellUpIn2 =
            GetLoc1D(imaxG - 2 - numGhosts_, jj, kk, imaxG, jmaxG);

        // location of upper i-boundary face
        int uFaceB = GetUpperFaceI(imaxG - 1 - numGhosts_, jj, kk,
                                   imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitI(uFaceB), "iu", inp, eos, suth, 1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          if (bcNameU == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
                bcNameU, (*this).FAreaUnitI(uFaceB), "iu", inp, eos, suth, 1);
          } else {
            (*this).state_[cellUpG2] = state_(cellUpIn1).GetGhostState(
                bcNameU, (*this).FAreaUnitI(uFaceB), "iu", inp, eos, suth, 2);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // loop over physical J faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(ii - numGhosts_, 0,
                                              kk - numGhosts_, "jl");
      string bcNameU = (*this).BC().GetBCName(ii - numGhosts_, jmax,
                                              kk - numGhosts_, "ju");

      // inviscid fluxes require different bc_ than viscous fluxes - treat all
      // walls as the same
      if (bcNameL == "viscousWall") {
        bcNameL = "slipWall";
      }
      if (bcNameU == "viscousWall") {
        bcNameU = "slipWall";
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at lower j-boundary
        int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
        int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

        // location of interior cells at lower j-boundary
        int cellLowIn1 = GetLoc1D(ii, numGhosts_, kk, imaxG, jmaxG);
        int cellLowIn2 =
            GetLoc1D(ii, numGhosts_ + 1, kk, imaxG, jmaxG);

        // location of lower j-boundary face
        int lFaceB = GetLowerFaceJ(ii, numGhosts_, kk, imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitJ(lFaceB), "jl", inp, eos, suth, 1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use once cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          if (bcNameL == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
                bcNameL, (*this).FAreaUnitJ(lFaceB), "jl", inp, eos, suth, 1);
          } else {
            (*this).state_[cellLowG2] = state_(cellLowIn1).GetGhostState(
                bcNameL, (*this).FAreaUnitJ(lFaceB), "jl", inp, eos, suth, 2);
          }
        }
      }

      // upper surface
      // ----------------------------------------------------------
      if (bcNameU != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at upper j-boundary
        int cellUpG1 = GetLoc1D(ii, jmaxG - 2, kk, imaxG, jmaxG);
        int cellUpG2 = GetLoc1D(ii, jmaxG - 1, kk, imaxG, jmaxG);

        // location of interior cells at upper j-boundary
        int cellUpIn1 =
            GetLoc1D(ii, jmaxG - 1 - numGhosts_, kk, imaxG, jmaxG);
        int cellUpIn2 =
            GetLoc1D(ii, jmaxG - 2 - numGhosts_, kk, imaxG, jmaxG);

        // location of upper j-boundary face
        int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - numGhosts_, kk,
                                   imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitJ(uFaceB), "ju", inp, eos, suth, 1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use once cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          if (bcNameU == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
                bcNameU, (*this).FAreaUnitJ(uFaceB), "ju", inp, eos, suth, 1);
          } else {
            (*this).state_[cellUpG2] = state_(cellUpIn1).GetGhostState(
                bcNameU, (*this).FAreaUnitJ(uFaceB), "ju", inp, eos, suth, 2);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // loop over physical K faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // name of boundary conditions at lower and upper boundaries
      string bcNameL = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, 0, "kl");
      string bcNameU = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, kmax, "ku");

      // inviscid fluxes require different bc_ than viscous fluxes - treat all
      // walls as the same
      if (bcNameL == "viscousWall") {
        bcNameL = "slipWall";
      }
      if (bcNameU == "viscousWall") {
        bcNameU = "slipWall";
      }

      // lower surface
      // ----------------------------------------------------------
      if (bcNameL != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at lower k-boundary
        int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
        int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

        // location of interior cells at lower k-boundary
        int cellLowIn1 = GetLoc1D(ii, jj, numGhosts_, imaxG, jmaxG);
        int cellLowIn2 =
            GetLoc1D(ii, jj, numGhosts_ + 1, imaxG, jmaxG);

        // location of lower k-boundary face
        int lFaceB = GetLowerFaceK(ii, jj, numGhosts_, imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitK(lFaceB), "kl", inp, eos, suth, 1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use once cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          if (bcNameL == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
                bcNameL, (*this).FAreaUnitK(lFaceB), "kl", inp, eos, suth, 1);
          } else {
            (*this).state_[cellLowG2] = state_(cellLowIn1).GetGhostState(
                bcNameL, (*this).FAreaUnitK(lFaceB), "kl", inp, eos, suth, 2);
          }
        }
      }

      // upper surface
      // ----------------------------------------------------------
      if (bcNameU != "interblock") {  // only supply state_ values for non
                                      // interblock BCs, for interblock, do
                                      // nothing

        // location of ghost cells at upper k-boundary
        int cellUpG1 = GetLoc1D(ii, jj, kmaxG - 2, imaxG, jmaxG);
        int cellUpG2 = GetLoc1D(ii, jj, kmaxG - 1, imaxG, jmaxG);

        // location of interior cells at upper k-boundary
        int cellUpIn1 =
            GetLoc1D(ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);
        int cellUpIn2 =
            GetLoc1D(ii, jj, kmaxG - 2 - numGhosts_, imaxG, jmaxG);

        // location of upper k-boundary face
        int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - numGhosts_,
                                   imaxG, jmaxG);

        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitK(uFaceB), "ku", inp, eos, suth, 1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use once cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          if (bcNameU == "slipWall") {  // if slipWall, reflect second interior
                                        // state_ over boundary face instead of
                                        // extrapolation
            (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
                bcNameU, (*this).FAreaUnitK(uFaceB), "ku", inp, eos, suth, 1);
          } else {
            (*this).state_[cellUpG2] = state_(cellUpIn1).GetGhostState(
                bcNameU, (*this).FAreaUnitK(uFaceB), "ku", inp, eos, suth, 2);
          }
        }
      }
    }
  }

  // assign values to edge ghost cells
  // (*this).AssignInviscidGhostCellsEdge(inp, eos, suth);
}

/* Member function to assign values to ghost cells located on the 12 block edges
for the inviscid flux calculation. Assumes AssignInviscidGhostCells has already
been run.
           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J

In the above diagram the cells marked X represent physical cells. Cells marked
G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each
corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells
(1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells
it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to
the values of the ghost cells they tough at level "e". The values at edge cell
2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it
touches. The exception to this rule occurs when either of the boundaries that meet
at the corner are wall boundaries (slipWall, viscousWall) and the other is not.
When this occurs the wall boundaries are "extended" into the ghost cells. This
implementation is described in Blazek.
*/
void procBlock::AssignInviscidGhostCellsEdge(const input &inp,
                                             const idealGas &eos,
                                             const sutherland &suth) {
  // inp -- all input variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity

  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * numGhosts_;
  int jmaxG = jmax + 2 * numGhosts_;
  int kmaxG = kmax + 2 * numGhosts_;

  //--------------------------------------------------------------------------
  // loop over edges at lower and upper j sides of block - this will include 4
  // edges that run in the i-direction --------------------------------
  // edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this
  // loop
  for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int j1, k1, j2, k2, je, ke, je2, ke2;

      int gf_j1_ke_kl;
      int gf_j2_ke_kl;
      int gf_je_k1_jl;
      int gf_je_k2_jl;

      string bc_J, bc_K;
      string surfJ, surfK;

      if (cc == 0) {  // at jl/kl edge - ghost cells are in the lower direction
                      // of both j and k, so use GetLowerFace for both
        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfJ = "jl";
        surfK = "kl";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at jl/ku edge - ghost cells are in the lower
                             // direction of j and upper direction of k, so use
                             // GetLowerFace for J
        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfJ = "jl";
        surfK = "ku";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at ju/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of j so use
                             // GetLowerFace for k
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfJ = "ju";
        surfK = "kl";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_ + 1,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at ju/ku edge - ghost cells are in the upper
                             // direction of both j and k, use GetUpperFace for
                             // both
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfJ = "ju";
        surfK = "ku";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_J = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_ + 1,
                                      ke - numGhosts_, surfJ);
        bc_K = (*this).BC().GetBCName(ii - numGhosts_,
                                      je - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if both
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_J != "interblock" && bc_K != "interblock") {
        int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
        int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of j line of cells, on second layer
                                          // of k line of cells
        int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of j line of cells, on first layer
                                          // of k line of cells
        int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of j line of cells, on second layer
                                          // of k line of cells

        int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of j
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of j
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // j line of cells,
                                                            // on first layer of
                                                            // k line of cells
        int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of j
                                         // line of cells, on second layer of k
                                         // line of cells

        int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG,
                                 jmaxG);  // ghost cell, on first layer of j
                                          // line of cells, on second non-edge
                                          // layer of k line of cells
        int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG,
                                 jmaxG);  // ghost cell, on second layer of j
                                          // line of cells, on second non-edge
                                          // layer of k line of cells
        int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of j line of cells, on first
                                          // layer of k line of cells
        int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of j line of cells, on second
                                          // layer of k line of cells

        // inviscid fluxes require different bc_ than viscous fluxes - treat all
        // walls as the same
        if (bc_J == "viscousWall") {
          bc_J = "slipWall";
        }
        if (bc_K == "viscousWall") {
          bc_K = "slipWall";
        }

        if (bc_J == "slipWall" &&
            !(bc_K == "slipWall")) {  // j surface is a wall, but k surface is
                                      // not - extend wall bc_
          (*this).state_[gce_j1_k1] = state_(gc_je_k1).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_je_k1_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_j1_k2] = state_(gc_je_k2).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_je_k2_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_j2_k1] = state_(gc_je2_k1).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_je_k1_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_j2_k2] = state_(gc_je2_k2).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_je_k2_jl), surfJ, inp, eos, suth, 1);
        } else if (!(bc_J == "slipWall") &&
                   bc_K == "slipWall") {  // k surface is a wall, but j surface
                                          // is not - extend wall bc_
          (*this).state_[gce_j1_k1] = state_(gc_j1_ke).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_j1_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_j2_k1] = state_(gc_j2_ke).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_j2_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_j1_k2] = state_(gc_j1_ke2).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_j1_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_j2_k2] = state_(gc_j2_ke2).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_j2_ke_kl), surfK, inp, eos, suth, 1);
        } else {  // both surfaces are walls or neither are walls - proceed as
                  // normal
          (*this).state_[gce_j1_k1] =
              0.5 * (state_(gc_j1_ke) + state_(gc_je_k1));
          (*this).state_[gce_j2_k1] = state_(gc_j2_ke);
          (*this).state_[gce_j1_k2] = state_(gc_je_k2);
          (*this).state_[gce_j2_k2] =
              0.5 * (state_(gc_j2_ke) + state_(gc_je_k2));
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the j-direction --------------------------------
  // edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this
  // loop
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, k1, i2, k2, ie, ke, ie2, ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_I, bc_K;
      string surfI, surfK;

      if (cc == 0) {  // at il/kl edge - ghost cells are in the lower direction
                      // of both i and k, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfI = "il";
        surfK = "kl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at il/ku edge - ghost cells are in the lower
                             // direction of i and upper direction of k, so use
                             // GetLowerFace for I
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ke;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfI = "il";
        surfK = "ku";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at iu/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of i so use
                             // GetLowerFace for k
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ke;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfI = "iu";
        surfK = "kl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at iu/ku edge - ghost cells are in the upper
                             // direction of both j and k, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfI = "iu";
        surfK = "ku";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditioins at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      jj - numGhosts_,
                                      ke - numGhosts_, surfI);
        bc_K = (*this).BC().GetBCName(ie - numGhosts_,
                                      jj - numGhosts_,
                                      ke - numGhosts_ + 1, surfK);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if bouth
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_I != "interblock" && bc_K != "interblock") {
        // location of ghost cells
        int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of k line of cells
        int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of k line of cells

        int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // k line of cells
        int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // i line of cells,
                                                            // on first layer of
                                                            // k line of cells
        int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of i
                                         // line of cells, on second layer of k
                                         // line of cells

        int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG,
                                 jmaxG);  // ghost cell, on first layer of i
                                          // line of cells, on second non-edge
                                          // layer of k line of cells
        int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG,
                                 jmaxG);  // ghost cell, on second layer of i
                                          // line of cells, on second non-edge
                                          // layer of k line of cells
        int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of i line of cells, on first
                                          // layer of k line of cells
        int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of i line of cells, on second
                                          // layer of k line of cells

        // inviscid fluxes require different bc_ than viscous fluxes - treat all
        // walls as the same
        if (bc_I == "viscousWall") {
          bc_I = "slipWall";
        }
        if (bc_K == "viscousWall") {
          bc_K = "slipWall";
        }

        if (bc_I == "slipWall" &&
            !(bc_K == "slipWall")) {  // i surface is a wall, but k surface is
                                      // not - extend wall bc_
          (*this).state_[gce_i1_k1] = state_(gc_ie_k1).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_k1_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i1_k2] = state_(gc_ie_k2).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_k2_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i2_k1] = state_(gc_ie2_k1).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_k1_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i2_k2] = state_(gc_ie2_k2).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_k2_il), surfI, inp, eos, suth, 1);
        } else if (!(bc_I == "slipWall") &&
                   bc_K == "slipWall") {  // k surface is a wall, but i surface
                                          // is not - extend wall bc_
          (*this).state_[gce_i1_k1] = state_(gc_i1_ke).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_i1_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_i2_k1] = state_(gc_i2_ke).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_i2_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_i1_k2] = state_(gc_i1_ke2).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_i1_ke_kl), surfK, inp, eos, suth, 1);
          (*this).state_[gce_i2_k2] = state_(gc_i2_ke2).GetGhostState(
              bc_K, (*this).FAreaUnitK(gf_i2_ke_kl), surfK, inp, eos, suth, 1);

        } else {  // both surfaces are walls or neither are walls - proceed as
                  // normal
          (*this).state_[gce_i1_k1] =
              0.5 * (state_(gc_i1_ke) + state_(gc_ie_k1));
          (*this).state_[gce_i2_k1] = state_(gc_i2_ke);
          (*this).state_[gce_i1_k2] = state_(gc_ie_k2);
          (*this).state_[gce_i2_k2] =
              0.5 * (state_(gc_i2_ke) + state_(gc_ie_k2));
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the k-direction --------------------------------
  // edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this
  // loop
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, j1, i2, j2, ie, je, ie2, je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_I, bc_J;
      string surfI, surfJ;

      if (cc == 0) {  // at il/jl edge - ghost cells are in the lower direction
                      // of both i and j, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        surfI = "il";
        surfJ = "jl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 1) {  // at il/ju edge - ghost cells are in the lower
                             // direction of i and upper direction of j, so use
                             // GetLowerFace for I
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }
        surfI = "il";
        surfJ = "ju";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_ + 1,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 2) {  // at iu/jl edge - ghost cells are in the lower
                             // direction of j, and upper direction of i so use
                             // GetLowerFace for j
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        surfI = "iu";
        surfJ = "jl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfJ);

      } else if (cc == 3) {  // at iu/ju edge - ghost cells are in the upper
                             // direction of both i and j, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        surfI = "iu";
        surfJ = "ju";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_I = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                      je - numGhosts_,
                                      kk - numGhosts_, surfI);
        bc_J = (*this).BC().GetBCName(ie - numGhosts_,
                                      je - numGhosts_ + 1,
                                      kk - numGhosts_, surfJ);
      }

      // if both corner boundaries aren't interblocks, get edge cells; if bouth
      // boundaries are interblocks, edge cells come from ghost cell exchange
      if (bc_I != "interblock" && bc_J != "interblock") {
        // location of ghost cells
        int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on first layer
                                          // of i line of cells, on second layer
                                          // of k line of cells
        int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on first layer
                                          // of k line of cells
        int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG,
                                 jmaxG);  // ghost cell on edge, on second layer
                                          // of i line of cells, on second layer
                                          // of k line of cells

        int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // first layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // j line of cells
        int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // second layer of i
                                                            // line of cells, on
                                                            // non-edge layer of
                                                            // j line of cells
        int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);  // ghost cell, on
                                                            // non-edge layer of
                                                            // i line of cells,
                                                            // on first layer of
                                                            // j line of cells
        int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG,
                                jmaxG);  // ghost cell, on non-edge layer of i
                                         // line of cells, on second layer of j
                                         // line of cells

        int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG,
                                 jmaxG);  // ghost cell, on first layer of i
                                          // line of cells, on second non-edge
                                          // layer of j line of cells
        int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG,
                                 jmaxG);  // ghost cell, on second layer of i
                                          // line of cells, on second non-edge
                                          // layer of j line of cells
        int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of i line of cells, on first
                                          // layer of j line of cells
        int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG,
                                 jmaxG);  // ghost cell, on second non-edge
                                          // layer of i line of cells, on second
                                          // layer of j line of cells

        // inviscid fluxes require different bc_ than viscous fluxes - treat all
        // walls as the same
        if (bc_I == "viscousWall") {
          bc_I = "slipWall";
        }
        if (bc_J == "viscousWall") {
          bc_J = "slipWall";
        }

        if (bc_I == "slipWall" &&
            !(bc_J == "slipWall")) {  // i surface is a wall, but k surface is
                                      // not - extend wall bc_
          (*this).state_[gce_i1_j1] = state_(gc_ie_j1).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_j1_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i1_j2] = state_(gc_ie_j2).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_j2_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i2_j1] = state_(gc_ie2_j1).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_j1_il), surfI, inp, eos, suth, 1);
          (*this).state_[gce_i2_j2] = state_(gc_ie2_j2).GetGhostState(
              bc_I, (*this).FAreaUnitI(gf_ie_j2_il), surfI, inp, eos, suth, 1);
        } else if (!(bc_I == "slipWall") &&
                   bc_J == "slipWall") {  // k surface is a wall, but i surface
                                          // is not - extend wall bc_
          (*this).state_[gce_i1_j1] = state_(gc_i1_je).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_i1_je_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_i2_j1] = state_(gc_i2_je).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_i2_je_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_i1_j2] = state_(gc_i1_je2).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_i1_je_jl), surfJ, inp, eos, suth, 1);
          (*this).state_[gce_i2_j2] = state_(gc_i2_je2).GetGhostState(
              bc_J, (*this).FAreaUnitJ(gf_i2_je_jl), surfJ, inp, eos, suth, 1);
        } else {  // both surfaces are walls or neither are walls - proceed as
                  // normal
          (*this).state_[gce_i1_j1] =
              0.5 * (state_(gc_i1_je) + state_(gc_ie_j1));
          (*this).state_[gce_i2_j1] = state_(gc_i2_je);
          (*this).state_[gce_i1_j2] = state_(gc_ie_j2);
          (*this).state_[gce_i2_j2] =
              0.5 * (state_(gc_i2_je) + state_(gc_ie_j2));
        }
      }
    }
  }
}

/* Member function to assign ghost cells for the viscous flow calculation. This
 function assumes AssignInviscidGhostCells has been run first as
 it only overwrites the ghost cells associated with the viscousWall boundary
 condition. It overwrites both regular and edge ghost cells.
*/
void procBlock::AssignViscousGhostCells(const input &inp, const idealGas &eos,
                                        const sutherland &suth) {
  // inp -- all input variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity

  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * numGhosts_;
  int jmaxG = (*this).NumJ() + 2 * numGhosts_;
  int kmaxG = (*this).NumK() + 2 * numGhosts_;

  //------------------------------------------------------------------------
  // loop over physical I faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
      // location of ghost cells at lower i-boundary
      int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

      // location of interior cells at lower i-boundary
      int cellLowIn1 = GetLoc1D(numGhosts_, jj, kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(numGhosts_ + 1, jj, kk, imaxG, jmaxG);

      // location of lower i-boundary face
      int lFaceB = GetLowerFaceI(numGhosts_, jj, kk, imaxG, jmaxG);

      // location of ghost cells at upper i-boundary
      int cellUpG1 = GetLoc1D(imaxG - 2, jj, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(imaxG - 1, jj, kk, imaxG, jmaxG);

      // location of interior cells at upper i-boundary
      int cellUpIn1 =
          GetLoc1D(imaxG - 1 - numGhosts_, jj, kk, imaxG, jmaxG);
      int cellUpIn2 =
          GetLoc1D(imaxG - 2 - numGhosts_, jj, kk, imaxG, jmaxG);

      // location of upper i-boundary face
      int uFaceB =
          GetUpperFaceI(imaxG - 1 - numGhosts_, jj, kk, imaxG, jmaxG);

      // boundary condition at lower boundary
      string bcNameL = (*this).BC().GetBCName(0, jj - numGhosts_,
                                              kk - numGhosts_, "il");

      // if viscous, overwrite regular ghost cell
      if (bcNameL == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitI(lFaceB), "il", inp, eos, suth, 1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
              bcNameL, (*this).FAreaUnitI(lFaceB), "il", inp, eos, suth, 2);
        }
      }

      // boundary condition at upper boundary
      string bcNameU = (*this).BC().GetBCName(imax, jj - numGhosts_,
                                              kk - numGhosts_, "iu");

      // if viscous, overwrite regular ghost cell
      if (bcNameU == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitI(uFaceB), "iu", inp, eos, suth, 1);

        // second layer of ghost cells
        if (imax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
              bcNameU, (*this).FAreaUnitI(uFaceB), "iu", inp, eos, suth, 2);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // loop over physical J faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // location of ghost cells at lower j-boundary
      int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

      // location of interior cells at lower j-boundary
      int cellLowIn1 = GetLoc1D(ii, numGhosts_, kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, numGhosts_ + 1, kk, imaxG, jmaxG);

      // location of lower j-boundary face
      int lFaceB = GetLowerFaceJ(ii, numGhosts_, kk, imaxG, jmaxG);

      // location of ghost cells at upper j-boundary
      int cellUpG1 = GetLoc1D(ii, jmaxG - 2, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jmaxG - 1, kk, imaxG, jmaxG);

      // location of interior cells at upper j-boundary
      int cellUpIn1 =
          GetLoc1D(ii, jmaxG - 1 - numGhosts_, kk, imaxG, jmaxG);
      int cellUpIn2 =
          GetLoc1D(ii, jmaxG - 2 - numGhosts_, kk, imaxG, jmaxG);

      // location of upper j-boundary face
      int uFaceB =
          GetUpperFaceJ(ii, jmaxG - 1 - numGhosts_, kk, imaxG, jmaxG);

      // boundary condition at lower boundary
      string bcNameL = (*this).BC().GetBCName(ii - numGhosts_, 0,
                                              kk - numGhosts_, "jl");

      // if viscous, overwrite regular ghost cell
      if (bcNameL == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitJ(lFaceB), "jl", inp, eos, suth, 1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
              bcNameL, (*this).FAreaUnitJ(lFaceB), "jl", inp, eos, suth, 2);
        }
      }

      // boundary condition at upper boundary
      string bcNameU = (*this).BC().GetBCName(ii - numGhosts_, jmax,
                                              kk - numGhosts_, "ju");

      // if viscous, overwrite regular ghost cell
      if (bcNameU == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitJ(uFaceB), "ju", inp, eos, suth, 1);

        // second layer of ghost cells
        if (jmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
              bcNameU, (*this).FAreaUnitJ(uFaceB), "ju", inp, eos, suth, 2);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // loop over physical K faces and assign values for regular ghost cells
  // -----------------------------------------------------------------------
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
      // location of ghost cells at lower k-boundary
      int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

      // location of interior cells at lower k-boundary
      int cellLowIn1 = GetLoc1D(ii, jj, numGhosts_, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, jj, numGhosts_ + 1, imaxG, jmaxG);

      // location of lower k-boundary face
      int lFaceB = GetLowerFaceK(ii, jj, numGhosts_, imaxG, jmaxG);

      // location of interior cells at lower k-boundary
      int cellUpG1 = GetLoc1D(ii, jj, kmaxG - 2, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jj, kmaxG - 1, imaxG, jmaxG);

      // location of interior cells at upper k-boundary
      int cellUpIn1 =
          GetLoc1D(ii, jj, kmaxG - 1 - numGhosts_, imaxG, jmaxG);
      int cellUpIn2 =
          GetLoc1D(ii, jj, kmaxG - 2 - numGhosts_, imaxG, jmaxG);

      // location of upper k-boundary face
      int uFaceB =
          GetUpperFaceK(ii, jj, kmax - 1 - numGhosts_, imaxG, jmaxG);

      // name of boundary condition at lower boundary
      string bcNameL = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, 0, "kl");

      // if viscous, overwrite regular ghost cells
      if (bcNameL == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellLowG1] = state_(cellLowIn1).GetGhostState(
            bcNameL, (*this).FAreaUnitK(lFaceB), "kl", inp, eos, suth, 1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellLowG2] = state_(cellLowG1);
        } else {
          (*this).state_[cellLowG2] = state_(cellLowIn2).GetGhostState(
              bcNameL, (*this).FAreaUnitK(lFaceB), "kl", inp, eos, suth, 2);
        }
      }

      // boundary condition at upper boundary
      string bcNameU = (*this).BC().GetBCName(
          ii - numGhosts_, jj - numGhosts_, kmax, "ku");

      // if viscous, overwrite regular ghost cells
      if (bcNameU == "viscousWall") {
        // first layer of ghost cells
        (*this).state_[cellUpG1] = state_(cellUpIn1).GetGhostState(
            bcNameU, (*this).FAreaUnitK(uFaceB), "ku", inp, eos, suth, 1);

        // second layer of ghost cells
        if (kmax < 2) {  // one cell thick - use one cell for both ghost cells
          (*this).state_[cellUpG2] = state_(cellUpG1);
        } else {
          (*this).state_[cellUpG2] = state_(cellUpIn2).GetGhostState(
              bcNameU, (*this).FAreaUnitK(uFaceB), "ku", inp, eos, suth, 2);
        }
      }
    }
  }

  // Assign edge ghost cells
  (*this).AssignViscousGhostCellsEdge(inp, eos, suth);
}

/* Member function to assign values to ghost cells located on the 12 block edges
for the viscous flux calculation. Assumes AssignViscousGhostCells has already
been run. Only overwrites edge ghost cells if one of the boundaries at the corner
is a viscousWall boundary condition.
           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J

In the above diagram the cells marked X represent physical cells. Cells marked
G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each
corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells
(1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells
it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to
the values of the ghost cells they tough at level "e". The values at edge cell
2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it
touches. The exception to this rule occurs when either of the boundaries that
meet at the corner are viscousWall boundaries and the other is not. When this
occurs the viscousWall boundaries are "extended" into the ghost cells. This
implementation is described in Blazek.
*/
void procBlock::AssignViscousGhostCellsEdge(const input &inp,
                                            const idealGas &eos,
                                            const sutherland &suth) {
  // inp -- all input variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity

  // max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  // max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * numGhosts_;
  int jmaxG = jmax + 2 * numGhosts_;
  int kmaxG = kmax + 2 * numGhosts_;

  //--------------------------------------------------------------------------
  // loop over edges at lower and upper j sides of block - this will include 4
  // edges that run in the i-direction --------------------------------
  // edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this
  // loop
  for (int ii = numGhosts_; ii < imax + numGhosts_; ii++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int j1, k1, j2, k2, je, ke, je2, ke2;

      int gf_j1_ke_kl;
      int gf_j2_ke_kl;
      int gf_je_k1_jl;
      int gf_je_k2_jl;

      string bc_jl, bc_kl;
      string surfJ, surfK;

      if (cc == 0) {  // at jl/kl edge - ghost cells are in the lower direction
                      // of both j and k, so use GetLowerFace for both
        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfJ = "jl";
        surfK = "kl";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_jl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_, surfJ);
        bc_kl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at jl/ku edge - ghost cells are in the lower
                             // direction of j and upper direction of k, so use
                             // GetLowerFace for J
        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfJ = "jl";
        surfK = "ku";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_jl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_, surfJ);
        bc_kl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at ju/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of j so use
                             // GetLowerFace for k
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfJ = "ju";
        surfK = "kl";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_jl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_ + 1,
                                       ke - numGhosts_, surfJ);
        bc_kl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at ju/ku edge - ghost cells are in the upper
                             // direction of both j and k, use GetUpperFace for
                             // both
        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfJ = "ju";
        surfK = "ku";

        // ghost face, on first/second layer of j line of cells, on non-edge
        // layer of k line of cells
        gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);
        gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of j line of cells, on first/second
        // layer of k line of cells
        gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);
        gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_jl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_ + 1,
                                       ke - numGhosts_, surfJ);
        bc_kl = (*this).BC().GetBCName(ii - numGhosts_,
                                       je - numGhosts_,
                                       ke - numGhosts_ + 1, surfK);
      }

      // location of ghost cells
      int gce_j1_k1 =
          GetLoc1D(ii, j1, k1, imaxG, jmaxG);  // ghost cell on edge, on first
                                               // layer of j line of cells, on
                                               // first layer of k line of cells
      int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG,
                               jmaxG);  // ghost cell on edge, on first layer of
                                        // j line of cells, on second layer of k
                                        // line of cells
      int gce_j2_k1 =
          GetLoc1D(ii, j2, k1, imaxG, jmaxG);  // ghost cell on edge, on second
                                               // layer of j line of cells, on
                                               // first layer of k line of cells
      int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG,
                               jmaxG);  // ghost cell on edge, on second layer
                                        // of j line of cells, on second layer
                                        // of k line of cells

      int gc_j1_ke =
          GetLoc1D(ii, j1, ke, imaxG, jmaxG);  // ghost cell, on first layer of
                                               // j line of cells, on non-edge
                                               // layer of k line of cells
      int gc_j2_ke =
          GetLoc1D(ii, j2, ke, imaxG, jmaxG);  // ghost cell, on second layer of
                                               // j line of cells, on non-edge
                                               // layer of k line of cells
      int gc_je_k1 =
          GetLoc1D(ii, je, k1, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of j line of cells, on first
                                               // layer of k line of cells
      int gc_je_k2 =
          GetLoc1D(ii, je, k2, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of j line of cells, on second
                                               // layer of k line of cells

      int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG,
                               jmaxG);  // ghost cell, on first layer of j line
                                        // of cells, on second non-edge layer of
                                        // k line of cells
      int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG,
                               jmaxG);  // ghost cell, on second layer of j line
                                        // of cells, on second non-edge layer of
                                        // k line of cells
      int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of j line of cells, on first layer of
                                        // k line of cells
      int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of j line of cells, on second layer
                                        // of k line of cells

      if (bc_jl == "viscousWall" &&
          !(bc_kl == "viscousWall")) {  // j surface is a viscous wall, but k
                                        // surface is not - extend wall bc_
        (*this).state_[gce_j1_k1] = state_(gc_je_k1).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_je_k1_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_j1_k2] = state_(gc_je_k2).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_je_k2_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_j2_k1] = state_(gc_je2_k1).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_je_k1_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_j2_k2] = state_(gc_je2_k2).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_je_k2_jl), surfJ, inp, eos, suth, 1);
      } else if (!(bc_jl == "viscousWall") &&
                 bc_kl == "viscousWall") {  // k surface is a viscous wall, but
                                            // j surface is not - extend wall
                                            // bc_
        (*this).state_[gce_j1_k1] = state_(gc_j1_ke).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_j1_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_j2_k1] = state_(gc_j2_ke).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_j2_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_j1_k2] = state_(gc_j1_ke2).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_j1_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_j2_k2] = state_(gc_j2_ke2).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_j2_ke_kl), surfK, inp, eos, suth, 1);
      } else if (bc_jl == "viscousWall" &&
                 bc_kl == "viscousWall") {  // both surfaces are viscous walls -
                                            // proceed as normal
        (*this).state_[gce_j1_k1] =
            0.5 * (state_(gc_j1_ke) + state_(gc_je_k1));
        (*this).state_[gce_j2_k1] = state_(gc_j2_ke);
        (*this).state_[gce_j1_k2] = state_(gc_je_k2);
        (*this).state_[gce_j2_k2] =
            0.5 * (state_(gc_j2_ke) + state_(gc_je_k2));
      }
      // if no boundary is a viscous wall, do nothing
    }
  }

  //-------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the j-direction --------------------------------
  // edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this
  // loop
  for (int jj = numGhosts_; jj < jmax + numGhosts_; jj++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, k1, i2, k2, ie, ke, ie2, ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_il, bc_kl;
      string surfI, surfK;

      if (cc == 0) {  // at il/kl edge - ghost cells are in the lower direction
                      // of both i and k, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfI = "il";
        surfK = "kl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfI);
        bc_kl = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfK);

      } else if (cc == 1) {  // at il/ku edge - ghost cells are in the lower
                             // direction of i and upper direction of k, so use
                             // GetLowerFace for I
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfI = "il";
        surfK = "ku";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfI);
        bc_kl = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_ + 1, surfK);

      } else if (cc == 2) {  // at iu/kl edge - ghost cells are in the lower
                             // direction of k, and upper direction of i so use
                             // GetLowerFace for k
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        k2 = 0;
        k1 = 1;
        ke = numGhosts_;
        if (kmax > 1) {
          ke2 = numGhosts_ + 1;
        } else {
          ke2 = ke;
        }

        surfI = "iu";
        surfK = "kl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfI);
        bc_kl = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfK);

      } else if (cc == 3) {  // at iu/ku edge - ghost cells are in the upper
                             // direction of both i and k, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        k2 = kmaxG - 1;
        k1 = kmaxG - 2;
        ke = kmax - 1 + numGhosts_;
        if (kmax > 1) {
          ke2 = kmax - 1 + numGhosts_ - 1;
        } else {
          ke2 = ke;
        }

        surfI = "iu";
        surfK = "ku";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of k line of cells
        gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);
        gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of k line of cells
        gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);
        gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                       jj - numGhosts_,
                                       ke - numGhosts_, surfI);
        bc_kl = (*this).BC().GetBCName(ie - numGhosts_,
                                       jj - numGhosts_,
                                       ke - numGhosts_ + 1, surfK);
      }

      // location of ghost cells
      int gce_i1_k1 =
          GetLoc1D(i1, jj, k1, imaxG, jmaxG);  // ghost cell on edge, on first
                                               // layer of i line of cells, on
                                               // first layer of k line of cells
      int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG,
                               jmaxG);  // ghost cell on edge, on first layer of
                                        // i line of cells, on second layer of k
                                        // line of cells
      int gce_i2_k1 =
          GetLoc1D(i2, jj, k1, imaxG, jmaxG);  // ghost cell on edge, on second
                                               // layer of i line of cells, on
                                               // first layer of k line of cells
      int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG,
                               jmaxG);  // ghost cell on edge, on second layer
                                        // of i line of cells, on second layer
                                        // of k line of cells

      int gc_i1_ke =
          GetLoc1D(i1, jj, ke, imaxG, jmaxG);  // ghost cell, on first layer of
                                               // i line of cells, on non-edge
                                               // layer of k line of cells
      int gc_i2_ke =
          GetLoc1D(i2, jj, ke, imaxG, jmaxG);  // ghost cell, on second layer of
                                               // i line of cells, on non-edge
                                               // layer of k line of cells
      int gc_ie_k1 =
          GetLoc1D(ie, jj, k1, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of i line of cells, on first
                                               // layer of k line of cells
      int gc_ie_k2 =
          GetLoc1D(ie, jj, k2, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of i line of cells, on second
                                               // layer of k line of cells

      int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG,
                               jmaxG);  // ghost cell, on first layer of i line
                                        // of cells, on second non-edge layer of
                                        // k line of cells
      int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG,
                               jmaxG);  // ghost cell, on second layer of i line
                                        // of cells, on second non-edge layer of
                                        // k line of cells
      int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of i line of cells, on first layer of
                                        // k line of cells
      int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of i line of cells, on second layer
                                        // of k line of cells

      if (bc_il == "viscousWall" &&
          !(bc_kl == "viscousWall")) {  // i surface is a viscous wall, but k
                                        // surface is not - extend wall bc_
        (*this).state_[gce_i1_k1] = state_(gc_ie_k1).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_k1_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i1_k2] = state_(gc_ie_k2).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_k2_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i2_k1] = state_(gc_ie2_k1).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_k1_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i2_k2] = state_(gc_ie2_k2).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_k2_il), surfI, inp, eos, suth, 1);
      } else if (!(bc_il == "viscousWall") &&
                 bc_kl == "viscousWall") {  // k surface is a viscous wall, but
                                            // i surface is not - extend wall
                                            // bc_
        (*this).state_[gce_i1_k1] = state_(gc_i1_ke).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_i1_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_i2_k1] = state_(gc_i2_ke).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_i2_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_i1_k2] = state_(gc_i1_ke2).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_i1_ke_kl), surfK, inp, eos, suth, 1);
        (*this).state_[gce_i2_k2] = state_(gc_i2_ke2).GetGhostState(
            bc_kl, (*this).FAreaUnitK(gf_i2_ke_kl), surfK, inp, eos, suth, 1);
      } else if (bc_il == "viscousWall" &&
                 bc_kl == "viscousWall") {  // both surfaces are viscous walls -
                                            // proceed as normal
        (*this).state_[gce_i1_k1] =
            0.5 * (state_(gc_i1_ke) + state_(gc_ie_k1));
        (*this).state_[gce_i2_k1] = state_(gc_i2_ke);
        (*this).state_[gce_i1_k2] = state_(gc_ie_k2);
        (*this).state_[gce_i2_k2] =
            0.5 * (state_(gc_i2_ke) + state_(gc_ie_k2));
      }
      // if neither surface is a wall then do nothing
    }
  }

  //--------------------------------------------------------------------------
  // loop over edges at lower and upper i sides of block - this will include 4
  // edges that run in the k-direction --------------------------------
  // edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this
  // loop
  for (int kk = numGhosts_; kk < kmax + numGhosts_; kk++) {
    for (int cc = 0; cc < 4; cc++) {  // loop over 4 edges
      int i1, j1, i2, j2, ie, je, ie2, je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_il, bc_jl;
      string surfI, surfJ;

      if (cc == 0) {  // at il/jl edge - ghost cells are in the lower direction
                      // of both i and j, so use GetLowerFace for both
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        surfI = "il";
        surfJ = "jl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfI);
        bc_jl = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfJ);

      } else if (cc == 1) {  // at il/ju edge - ghost cells are in the lower
                             // direction of i and upper direction of j, so use
                             // GetLowerFace for I
        i2 = 0;
        i1 = 1;
        ie = numGhosts_;
        if (imax > 1) {
          ie2 = numGhosts_ + 1;
        } else {
          ie2 = ie;
        }

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        surfI = "il";
        surfJ = "ju";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfI);
        bc_jl = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_ + 1,
                                       kk - numGhosts_, surfJ);

      } else if (cc == 2) {  // at iu/jl edge - ghost cells are in the lower
                             // direction of j, and upper direction of i so use
                             // GetLowerFace for J
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        j2 = 0;
        j1 = 1;
        je = numGhosts_;
        if (jmax > 1) {
          je2 = numGhosts_ + 1;
        } else {
          je2 = je;
        }

        surfI = "iu";
        surfJ = "jl";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfI);
        bc_jl = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfJ);

      } else if (cc == 3) {  // at iu/ju edge - ghost cells are in the upper
                             // direction of both i and j, use GetUpperFace for
                             // both
        i2 = imaxG - 1;
        i1 = imaxG - 2;
        ie = imax - 1 + numGhosts_;
        if (imax > 1) {
          ie2 = imax - 1 + numGhosts_ - 1;
        } else {
          ie2 = ie;
        }

        j2 = jmaxG - 1;
        j1 = jmaxG - 2;
        je = jmax - 1 + numGhosts_;
        if (jmax > 1) {
          je2 = jmax - 1 + numGhosts_ - 1;
        } else {
          je2 = je;
        }

        surfI = "iu";
        surfJ = "ju";

        // ghost face, on first/second layer of i line of cells, on non-edge
        // layer of j line of cells
        gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);
        gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);

        // ghost face, on non-edge layer of i line of cells, on first/second
        // layer of j line of cells
        gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);
        gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);

        // boundary conditions at corner
        bc_il = (*this).BC().GetBCName(ie - numGhosts_ + 1,
                                       je - numGhosts_,
                                       kk - numGhosts_, surfI);
        bc_jl = (*this).BC().GetBCName(ie - numGhosts_,
                                       je - numGhosts_ + 1,
                                       kk - numGhosts_, surfJ);
      }

      // location of ghost cells
      int gce_i1_j1 =
          GetLoc1D(i1, j1, kk, imaxG, jmaxG);  // ghost cell on edge, on first
                                               // layer of i line of cells, on
                                               // first layer of k line of cells
      int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG,
                               jmaxG);  // ghost cell on edge, on first layer of
                                        // i line of cells, on second layer of k
                                        // line of cells
      int gce_i2_j1 =
          GetLoc1D(i2, j1, kk, imaxG, jmaxG);  // ghost cell on edge, on second
                                               // layer of i line of cells, on
                                               // first layer of k line of cells
      int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG,
                               jmaxG);  // ghost cell on edge, on second layer
                                        // of i line of cells, on second layer
                                        // of k line of cells

      int gc_i1_je =
          GetLoc1D(i1, je, kk, imaxG, jmaxG);  // ghost cell, on first layer of
                                               // i line of cells, on non-edge
                                               // layer of j line of cells
      int gc_i2_je =
          GetLoc1D(i2, je, kk, imaxG, jmaxG);  // ghost cell, on second layer of
                                               // i line of cells, on non-edge
                                               // layer of j line of cells
      int gc_ie_j1 =
          GetLoc1D(ie, j1, kk, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of i line of cells, on first
                                               // layer of j line of cells
      int gc_ie_j2 =
          GetLoc1D(ie, j2, kk, imaxG, jmaxG);  // ghost cell, on non-edge layer
                                               // of i line of cells, on second
                                               // layer of j line of cells

      int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG,
                               jmaxG);  // ghost cell, on first layer of i line
                                        // of cells, on second non-edge layer of
                                        // j line of cells
      int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG,
                               jmaxG);  // ghost cell, on second layer of i line
                                        // of cells, on second non-edge layer of
                                        // j line of cells
      int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of i line of cells, on first layer of
                                        // j line of cells
      int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG,
                               jmaxG);  // ghost cell, on second non-edge layer
                                        // of i line of cells, on second layer
                                        // of j line of cells

      if (bc_il == "viscousWall" &&
          !(bc_jl == "viscousWall")) {  // i surface is a viscous wall, but k
                                        // surface is not - extend wall bc_
        (*this).state_[gce_i1_j1] = state_(gc_ie_j1).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_j1_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i1_j2] = state_(gc_ie_j2).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_j2_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i2_j1] = state_(gc_ie2_j1).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_j1_il), surfI, inp, eos, suth, 1);
        (*this).state_[gce_i2_j2] = state_(gc_ie2_j2).GetGhostState(
            bc_il, (*this).FAreaUnitI(gf_ie_j2_il), surfI, inp, eos, suth, 1);
      } else if (!(bc_il == "viscousWall") &&
                 bc_jl == "viscousWall") {  // j surface is a viscous wall, but
                                            // i surface is not - extend wall
                                            // bc_
        (*this).state_[gce_i1_j1] = state_(gc_i1_je).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_i1_je_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_i2_j1] = state_(gc_i2_je).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_i2_je_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_i1_j2] = state_(gc_i1_je2).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_i1_je_jl), surfJ, inp, eos, suth, 1);
        (*this).state_[gce_i2_j2] = state_(gc_i2_je2).GetGhostState(
            bc_jl, (*this).FAreaUnitJ(gf_i2_je_jl), surfJ, inp, eos, suth, 1);
      } else if (bc_il == "viscousWall" &&
                 bc_jl == "viscousWall") {  // both surfaces are viscous walls -
                                            // proceed as normal
        (*this).state_[gce_i1_j1] =
            0.5 * (state_(gc_i1_je) + state_(gc_ie_j1));
        (*this).state_[gce_i2_j1] = state_(gc_i2_je);
        (*this).state_[gce_i1_j2] = state_(gc_ie_j2);
        (*this).state_[gce_i2_j2] =
            0.5 * (state_(gc_i2_je) + state_(gc_ie_j2));
      }
      // if neither surface is a viscous wall then do nothing
    }
  }
}

/* Member function to determine where in padded plot3dBlock an index is located.
It takes in an i, j, k cell location and returns a boolean indicating
if the given i, j, k location corresponds to a physical cell location.
 */
bool procBlock::IsPhysical(const int &ii, const int &jj, const int &kk,
                           const bool &includeGhost) const {
  // ii -- i index of location to test
  // jj -- j index of location to test
  // kk -- k index of location to test
  // includeGhost -- flag to determine if inputs include ghost cells or not

  bool isPhysical = true;

  int offset = includeGhost ? numGhosts_ : 0;

  // if any of (i, j, & k) are outside of the limits of physical cells, location
  // is non-physical
  if ((ii < offset ||
       ii > (*this).NumI() - 1 + offset) ||
      (jj < offset ||
       jj > (*this).NumJ() - 1 + offset) ||
      (kk < offset ||
       kk > (*this).NumK() - 1 + offset)) {
    isPhysical = false;
  }

  return isPhysical;
}

/* Member function to determine where in padded plot3dBlock an index is located.
It takes in an i, j, k cell location and returns a boolean indicating
if the given i, j, k location corresponds to a corner location. Corner locations
are not used at all. This function is NOT USED in the code but is useful
for debugging purposes.
 */
bool procBlock::AtCorner(const int &ii, const int &jj, const int &kk,
                         const bool &includeGhost) const {
  // ii -- i index of location to test
  // jj -- j index of location to test
  // kk -- k index of location to test
  // includeGhost -- flag to determine if inputs include ghost cells or not

  bool atCorner = false;

  int offset = includeGhost ? numGhosts_ : 0;

  // if all (i, j, & k) are outside of the limits of physical cells, location is
  // a corner location
  if ((ii < offset ||
       ii > (*this).NumI() - 1 + offset) &&
      (jj < offset ||
       jj > (*this).NumJ() - 1 + offset) &&
      (kk < offset ||
       kk > (*this).NumK() - 1 + offset)) {
    atCorner = true;
  }

  return atCorner;
}

/* Member function to determine where in padded plot3dBlock an index is located.
It takes in an i, j, k cell location and returns a boolean indicating
if the given i, j, k location corresponds to a edge location. Edge locations are
used in the gradient calculations.
 */
bool procBlock::AtEdge(const int &ii, const int &jj, const int &kk,
                       const bool &includeGhost, string &dir) const {
  // ii -- i index of location to test
  // jj -- j index of location to test
  // kk -- k index of location to test
  // includeGhost -- flag to determine if inputs include ghost cells or not
  // dir -- direction that edge runs in

  bool atEdge = false;

  int offset = includeGhost ? numGhosts_ : 0;

  if ((ii >= offset &&
       ii < (*this).NumI() + offset) &&  // at i-edge - i in
                                         // physical cell range,
                                         // j/k at first level of
                                         // ghost cells
      (jj == offset - 1 ||
       jj == (*this).NumJ() + offset) &&
      (kk == offset - 1 ||
       kk == (*this).NumK() + offset)) {
    atEdge = true;
    dir = "i";
  } else if ((ii == offset - 1 ||
              ii == (*this).NumI() +
                        offset) &&  // at j-edge - j in physical
                                    // cell range, i/k at first
                                    // level of ghost cells
             (jj >= offset &&
              jj < (*this).NumJ() + offset) &&
             (kk == offset - 1 ||
              kk == (*this).NumK() + offset)) {
    atEdge = true;
    dir = "j";
  } else if ((ii == offset - 1 ||
              ii == (*this).NumI() +
                        offset) &&  // at k-edge - k in physical
                                    // cell range, i/j at first
                                    // level of ghost cells
             (jj == offset - 1 ||
              jj == (*this).NumJ() + offset) &&
             (kk >= offset &&
              kk < (*this).NumK() + offset)) {
    atEdge = true;
    dir = "k";
  }

  return atEdge;
}

/* Function to swap ghost cell geometry between two blocks at an interblock
boundary. Slices are removed from the physical cells (extending into ghost cells
at the edges) of one block and inserted into the ghost cells of its partner
block. The reverse is also true. The slices are taken in the coordinate system
orientation of their parent block.

   Interior Cells    Ghost Cells               Ghost Cells   Interior Cells
   ________ ______|________ _________       _______________|_______ _________
Ui-3/2   Ui-1/2   |    Uj+1/2    Uj+3/2  Ui-3/2    Ui-1/2  |    Uj+1/2    Uj+3/2
  |        |      |        |         |     |        |      |       |         |
  | Ui-1   |  Ui  |  Uj    |  Uj+1   |     |  Ui-1  |   Ui |  Uj   |  Uj+1   |
  |        |      |        |         |     |        |      |       |         |
  |________|______|________|_________|     |________|______|_______|_________|
                  |                                        |

The above diagram shows the resulting values after the ghost cell swap. The
logic ensures that the ghost cells at the interblock boundary exactly match
their partner block as if there were no separation in the grid.

Only 3 faces at each ghost cell need to be swapped (i.e. the lower face for Ui
is the upper face for Ui-1). At the end of a line (i-line, j-line or k-line),
both the upper and lower faces need to be swapped.
*/
void SwapSlice(interblock &inter, procBlock &blk1, procBlock &blk2,
               const bool &geom) {
  // inter -- interblock boundary information
  // blk1 -- first block involved in interblock boundary
  // blk2 -- second block involved in interblock boundary
  // geom -- boolean to determine whether to swap geometry or states

  // Get indices for slice coming from first block to swap
  int is1, ie1, js1, je1, ks1, ke1;

  // if at upper boundary no need to adjust for ghost cells as constant surface
  // is already at the interior cells when acounting for ghost cells
  // if at the lower boundary adjust the constant surface by the number of ghost
  // cells to get to the first interior cell
  int upLowFac = (inter.BoundaryFirst() % 2 == 0) ? 0 : blk1.NumGhosts();

  if (inter.BoundaryFirst() == 1 ||
      inter.BoundaryFirst() == 2) {  // direction 3 is i
    // extend min/maxes to cover ghost cells
    is1 = inter.ConstSurfaceFirst() + upLowFac;
    ie1 = is1 + blk1.NumGhosts() - 1;

    // direction 1 is j
    js1 = inter.Dir1StartFirst();
    je1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();

    // direction 2 is k
    ks1 = inter.Dir2StartFirst();
    ke1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  } else if (inter.BoundaryFirst() == 3 ||
             inter.BoundaryFirst() == 4) {  // direction 3 is j
    // extend min/maxes to cover ghost cells
    js1 = inter.ConstSurfaceFirst() + upLowFac;
    je1 = js1 + blk1.NumGhosts() - 1;

    // direction 1 is k
    ks1 = inter.Dir1StartFirst();
    ke1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();

    // direction 2 is i
    is1 = inter.Dir2StartFirst();
    ie1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  } else if (inter.BoundaryFirst() == 5 ||
             inter.BoundaryFirst() == 6) {  // direction 3 is k
    // extend min/maxes to cover ghost cells
    ks1 = inter.ConstSurfaceFirst() + upLowFac;
    ke1 = ks1 + blk1.NumGhosts() - 1;

    // direction 1 is i
    is1 = inter.Dir1StartFirst();
    ie1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();

    // direction 2 is j
    js1 = inter.Dir2StartFirst();
    je1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  } else {
    cerr << "ERROR: Error in procBlock::SwapSlice(). Surface boundary "
         << inter.BoundaryFirst() << " is not recognized!" << endl;
    exit(0);
  }

  // Get indices for slice coming from second block to swap
  int is2, ie2, js2, je2, ks2, ke2;

  // if at upper boundary no need to adjust for ghost cells as constant surface
  // is already at the interior cells when acounting for ghost cells
  // if at the lower boundary adjust the constant surface by the number of ghost
  // cells to get to the first interior cell
  upLowFac = (inter.BoundarySecond() % 2 == 0) ? 0 : blk2.NumGhosts();

  if (inter.BoundarySecond() == 1 ||
      inter.BoundarySecond() == 2) {  // direction 3 is i
    // extend min/maxes to cover ghost cells
    is2 = inter.ConstSurfaceSecond() + upLowFac;
    ie2 = is2 + blk2.NumGhosts() - 1;

    // direction 1 is j
    js2 = inter.Dir1StartSecond();
    je2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();

    // direction 2 is k
    ks2 = inter.Dir2StartSecond();
    ke2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  } else if (inter.BoundarySecond() == 3 ||
             inter.BoundarySecond() == 4) {  // direction 3 is j
    // extend min/maxes to cover ghost cells
    js2 = inter.ConstSurfaceSecond() + upLowFac;
    je2 = js2 + blk2.NumGhosts() - 1;

    // direction 1 is k
    ks2 = inter.Dir1StartSecond();
    ke2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();

    // direction 2 is i
    is2 = inter.Dir2StartSecond();
    ie2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  } else if (inter.BoundarySecond() == 5 ||
             inter.BoundarySecond() == 6) {  // direction 3 is k
    // extend min/maxes to cover ghost cells
    ks2 = inter.ConstSurfaceSecond() + upLowFac;
    ke2 = ks2 + blk2.NumGhosts() - 1;

    // direction 1 is i
    is2 = inter.Dir1StartSecond();
    ie2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();

    // direction 2 is j
    js2 = inter.Dir2StartSecond();
    je2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  } else {
    cerr << "ERROR: Error in procBlock::SwapSlice(). Surface boundary "
         << inter.BoundarySecond() << " is not recognized!" << endl;
    exit(0);
  }

  geomSlice geom1, geom2;
  stateSlice state1, state2;
  if (geom) {  // get geomSlices to swap
    geom1 = geomSlice(blk1, is1, ie1, js1, je1, ks1, ke1);
    geom2 = geomSlice(blk2, is2, ie2, js2, je2, ks2, ke2);
  } else {  // get stateSlices to swap
    state1 = stateSlice(blk1, is1, ie1, js1, je1, ks1, ke1);
    state2 = stateSlice(blk2, is2, ie2, js2, je2, ks2, ke2);
  }

  // change interblocks to work with slice and ghosts
  interblock inter1 = inter;
  interblock inter2 = inter;
  inter1.AdjustForSlice(false, blk1.NumGhosts());
  inter2.AdjustForSlice(true, blk2.NumGhosts());

  // put slices in proper blocks
  if (geom) {  // put geomSlices in procBlock
    // return vector determining if any of the 4 edges of the interblock need to
    // be updated for a "t" intersection
    vector<bool> adjEdge1 =
        blk1.PutGeomSlice(geom2, inter2, blk2.NumGhosts(), blk2.NumGhosts());
    vector<bool> adjEdge2 =
        blk2.PutGeomSlice(geom1, inter1, blk1.NumGhosts(), blk1.NumGhosts());

    // if an interblock border needs to be updated, update
    for (unsigned int ii = 0; ii < adjEdge1.size(); ii++) {
      if (adjEdge1[ii]) {
        inter.UpdateBorderFirst(ii);
      }
      if (adjEdge2[ii]) {
        inter.UpdateBorderSecond(ii);
      }
    }
  } else {  // put stateSlices in procBlock
    blk1.PutStateSlice(state2, inter2, blk2.NumGhosts(), blk2.NumGhosts());
    blk2.PutStateSlice(state1, inter1, blk1.NumGhosts(), blk1.NumGhosts());
  }
}

/* Function to swap slice using MPI. This is similar to the SwapSlice member
function, but is called when the neighboring procBlocks are on different
processors.
*/
void procBlock::SwapSliceMPI(const interblock &inter, const int &rank,
                             const MPI_Datatype &MPI_cellData) {
  // inter -- interblock boundary information
  // rank -- processor rank
  // MPI_cellData -- MPI datatype for passing primVars, genArray

  // Get indices for slice coming from block to swap
  int is, ie, js, je, ks, ke;

  if (rank == inter.RankFirst()) {  // local block is first in interblock
    // if at upper boundary no need to adjust for ghost cells as constant
    // surface is already at the interior cells when acounting for ghost cells
    // if at the lower boundary adjust the constant surface by the number of
    // ghost cells to get to the first interior cell
    int upLowFac = (inter.BoundaryFirst() % 2 == 0) ? 0 : numGhosts_;

    if (inter.BoundaryFirst() == 1 ||
        inter.BoundaryFirst() == 2) {  // direction 3 is i
      // extend min/maxes to cover ghost cells
      is = inter.ConstSurfaceFirst() + upLowFac;
      ie = is + numGhosts_ - 1;

      // direction 1 is j
      js = inter.Dir1StartFirst();
      je = inter.Dir1EndFirst() - 1 + 2 * numGhosts_;

      // direction 2 is k
      ks = inter.Dir2StartFirst();
      ke = inter.Dir2EndFirst() - 1 + 2 * numGhosts_;
    } else if (inter.BoundaryFirst() == 3 ||
               inter.BoundaryFirst() == 4) {  // direction 3 is j
      // extend min/maxes to cover ghost cells
      js = inter.ConstSurfaceFirst() + upLowFac;
      je = js + numGhosts_ - 1;

      // direction 1 is k
      ks = inter.Dir1StartFirst();
      ke = inter.Dir1EndFirst() - 1 + 2 * numGhosts_;

      // direction 2 is i
      is = inter.Dir2StartFirst();
      ie = inter.Dir2EndFirst() - 1 + 2 * numGhosts_;
    } else if (inter.BoundaryFirst() == 5 ||
               inter.BoundaryFirst() == 6) {  // direction 3 is k
      // extend min/maxes to cover ghost cells
      ks = inter.ConstSurfaceFirst() + upLowFac;
      ke = ks + numGhosts_ - 1;

      // direction 1 is i
      is = inter.Dir1StartFirst();
      ie = inter.Dir1EndFirst() - 1 + 2 * numGhosts_;

      // direction 2 is j
      js = inter.Dir2StartFirst();
      je = inter.Dir2EndFirst() - 1 + 2 * numGhosts_;
    } else {
      cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Surface boundary "
           << inter.BoundaryFirst() << " is not recognized!" << endl;
      exit(0);
    }
  // local block is second in interblock
  } else if (rank == inter.RankSecond()) {
    // if at upper boundary no need to adjust for ghost cells as constant
    // surface is already at the interior cells when acounting for ghost cells
    // if at the lower boundary adjust the constant surface by the number of
    // ghost cells to get to the first interior cell
    int upLowFac = (inter.BoundarySecond() % 2 == 0) ? 0 : numGhosts_;

    if (inter.BoundarySecond() == 1 ||
        inter.BoundarySecond() == 2) {  // direction 3 is i
      // extend min/maxes to cover ghost cells
      is = inter.ConstSurfaceSecond() + upLowFac;
      ie = is + numGhosts_ - 1;

      // direction 1 is j
      js = inter.Dir1StartSecond();
      je = inter.Dir1EndSecond() - 1 + 2 * numGhosts_;

      // direction 2 is k
      ks = inter.Dir2StartSecond();
      ke = inter.Dir2EndSecond() - 1 + 2 * numGhosts_;
    } else if (inter.BoundarySecond() == 3 ||
               inter.BoundarySecond() == 4) {  // direction 3 is j
      // extend min/maxes to cover ghost cells
      js = inter.ConstSurfaceSecond() + upLowFac;
      je = js + numGhosts_ - 1;

      // direction 1 is k
      ks = inter.Dir1StartSecond();
      ke = inter.Dir1EndSecond() - 1 + 2 * numGhosts_;

      // direction 2 is i
      is = inter.Dir2StartSecond();
      ie = inter.Dir2EndSecond() - 1 + 2 * numGhosts_;
    } else if (inter.BoundarySecond() == 5 ||
               inter.BoundarySecond() == 6) {  // direction 3 is k
      // extend min/maxes to cover ghost cells
      ks = inter.ConstSurfaceSecond() + upLowFac;
      ke = ks + numGhosts_ - 1;

      // direction 1 is i
      is = inter.Dir1StartSecond();
      ie = inter.Dir1EndSecond() - 1 + 2 * numGhosts_;

      // direction 2 is j
      js = inter.Dir2StartSecond();
      je = inter.Dir2EndSecond() - 1 + 2 * numGhosts_;
    } else {
      cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Surface boundary "
           << inter.BoundarySecond() << " is not recognized!" << endl;
      exit(0);
    }
  } else {
    cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Processor rank does "
            "not match either of interblock ranks!" << endl;
    exit(0);
  }

  // get local state slice to swap
  stateSlice state((*this), is, ie, js, je, ks, ke);

  // swap state slices with partner block
  state.PackSwapUnpackMPI(inter, MPI_cellData, rank);

  // change interblocks to work with slice and ghosts
  interblock interAdj = inter;

  // block to insert into is first in interblock
  if (rank == inter.RankSecond()) {
    interAdj.AdjustForSlice(false, numGhosts_);
  } else {  // block to insert into is second in interblock, so pass swapped
            // version
    interAdj.AdjustForSlice(true, numGhosts_);
  }

  // insert stateSlice into procBlock
  (*this).PutStateSlice(state, interAdj, numGhosts_,
                        numGhosts_);
}

/* Function to return a vector of location indicies for ghost cells at an
interblock boundary. The vector is formatted as shown below:

  vector = [i j k]

The vector will contain 3 entries corresponding to the i, j, and k locations of
either the first or second pair in the interblock, depending on what is
specified in the pairID variable. The indices returned will correspond to cell
locations and will take into account the orientation of the patches that
comprise the interblock with relation to each other.
*/
vector3d<int> GetSwapLoc(const int &l1, const int &l2, const int &l3,
                         const interblock &inter, const bool &pairID) {
  // l1 -- index of direction 1 within slice to insert
  // l2 -- index of direction 2 within slice to insert
  // l3 -- index of direction 3 within slice to insert
  // inter -- interblock boundary condition
  // pairID -- returning index for first or second block in interblock match

  // preallocate vector to return
  vector3d<int> loc;

  if (pairID) {  // working on first in pair ------------------------------
    // first patch in pair is calculated using orientation 1
    if (inter.Direction3First() == "i") {  // i-patch
      // get direction 1 length
      loc[1] = inter.Dir1StartFirst() + l1;  // direction 1 is j
      loc[2] = inter.Dir2StartFirst() + l2;  // direction 2 is k
      loc[0] = inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else if (inter.Direction3First() == "j") {  // j-patch
      // get direction 1 length
      loc[2] = inter.Dir1StartFirst() + l1;  // direction 1 is k
      loc[0] = inter.Dir2StartFirst() + l2;  // direction 2 is i
      loc[1] = inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else if (inter.Direction3First() == "k") {  // k-patch
      // get direction 1 length
      loc[0] = inter.Dir1StartFirst() + l1;  // direction 1 is i
      loc[1] = inter.Dir2StartFirst() + l2;  // direction 2 is j
      loc[2] = inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else {
      cerr << "ERROR: Error in procBlock:GetSwapLoc(). Boundary direction "
           << inter.Direction3First() << " is not recognized!" << endl;
      exit(0);
    }
  //--------------------------------------------------------------------------
  // need to use orientation for second in pair
  } else {  // working on second in pair ---------------------------------------
    if (inter.Direction3Second() == "i") {  // i-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is j (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[2] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is k (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[1] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is j -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[1] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 1 is k -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[2] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      loc[0] = inter.ConstSurfaceSecond() + l3;  // add l3 to get to ghost cells

    //-------------------------------------------------------------------------
    } else if (inter.Direction3Second() == "j") {  // j-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is k (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[0] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is i (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[2] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is k -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[2] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 2 is i -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[0] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      loc[1] = inter.ConstSurfaceSecond() + l3;  // add l3 to get to ghost cells

    //------------------------------------------------------------------------
    } else if (inter.Direction3Second() == "k") {  // k-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is i (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[1] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is j (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[0] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is i -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[0] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 2 is j -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[1] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      loc[2] = inter.ConstSurfaceSecond() + l3;  // add l3 to get to ghost cells

    //--------------------------------------------------------------------------
    } else {
      cerr << "ERROR: Error in procBlock.cpp:GetSwapLoc(). Boundary surface of "
           << inter.Direction3Second() << " is not recognized!" << endl;
      exit(0);
    }
  }

  return loc;
}

/* Function to populate ghost cells with proper cell states for inviscid flow
calculation. This function operates on the entire grid and uses interblock
boundaries to pass the correct data between grid blocks.
*/
void GetBoundaryConditions(vector<procBlock> &states, const input &inp,
                           const idealGas &eos, const sutherland &suth,
                           vector<interblock> &connections, const int &rank,
                           const MPI_Datatype &MPI_cellData) {
  // states -- vector of all procBlocks in the solution domain
  // inp -- all input variables
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // connections -- vector of interblock connections

  // loop over all blocks and assign inviscid ghost cells
  for (unsigned int ii = 0; ii < states.size(); ii++) {
    states[ii].AssignInviscidGhostCells(inp, eos, suth);
  }

  // loop over connections and swap ghost cells where needed
  for (unsigned int ii = 0; ii < connections.size(); ii++) {
    if (connections[ii].RankFirst() == rank &&
        connections[ii].RankSecond() == rank) {  // both sides of interblock
                                                  // are on this processor, swap
                                                  // w/o mpi
      SwapSlice(connections[ii], states[connections[ii].LocalBlockFirst()],
                states[connections[ii].LocalBlockSecond()], false);
    } else if (connections[ii].RankFirst() ==
               rank) {  // rank matches rank of first side of interblock,
                         // swap over mpi
      states[connections[ii].LocalBlockFirst()]
          .SwapSliceMPI(connections[ii], rank, MPI_cellData);

    } else if (connections[ii].RankSecond() ==
               rank) {  // rank matches rank of second side of interblock,
                         // swap over mpi
      states[connections[ii].LocalBlockSecond()]
          .SwapSliceMPI(connections[ii], rank, MPI_cellData);
    }
    // if rank doesn't match either side of interblock, then do nothing and
    // move on to the next interblock
  }

  // loop over all blocks and get ghost cell edge data
  for (unsigned int ii = 0; ii < states.size(); ii++) {
    states[ii].AssignInviscidGhostCellsEdge(inp, eos, suth);
  }
}

/* Member function to overwrite a section of a procBlock's geometry with a
geomSlice. The function uses the orientation supplied in the interblock to
orient the geomSlice relative to the procBlock. It assumes that the procBlock
is listed first, and the geomSlice second in the interblock data structure.
It returns a vector of 4 bools that are returned true if one of the 4 edges
of the interblock need to be updated because they border an interblock,
creating a possible "t" intersection.
     __________     _________
    |       |  |   |  |      |
    |       |  |   |  |  2   |
    |       |  |   |__|______|
    |       |  |   |X_|______|
    |  0    |--|    _________
    |       |  |   |X_|______|
    |       |  |   |  |      |
    |       |  |   |  |  1   |
    |_______|__|   |__|______|

The block configuration shown above shows a "t"intersection. If blocks 0/1 are
swapped first, block 1 gets all of its ghost cells including the edge cells
(marked X) from block 0, but block 0 knows that the interblock swapping with
block 1 borders another interblock (the one swapping with block 2), so it does
not take the edge ghost cells from block 1. This is because the ghost cells that
the edge of block 1 would supply will be supplied by the swap with block 2. If
next, blocks 1/2 swap, block 1 now has garbage in its edge ghost cells. This is
because it is the first swap for block 2, so its ghost cells are empty (filled
with garbage). If next blocks 0/2 swap, block 0 does not take the ghost cells
from block 2 for the same reason it did not take them from block 1. This leaves
blocks 0 and 2 with the correct ghost cells, but block 1 has garbage in its edge
ghost cells. This problem occurs if the 1/2 swap is done between the 0/1 and
0/2 swaps!

To solve this problem, in this situation the 1/2 swap is done without edge ghost
cells. This fixes the problem because the correct ghost cells block 1 got from
block 0 are not overwritten with the garbage ghost cells from block 2. The way
this is enforced is the following. Initially all ghost cell volumes are
initialized to 0. Therefore, if a cell from a geomSlice with a volume of 0 is
trying to be inserted into the procBlock, this cell has not been given its proper
boundary conditions and should be ignored. Subsequently, the interblock should be
updated so that in the future this cell is not inserted.
*/
vector<bool> procBlock::PutGeomSlice(const geomSlice &slice, interblock &inter,
                                     const int &d3, const int &numG) {
  // slice -- geomSlice to insert int procBlock
  // inter -- interblock data structure describing the patches and their
  // orientation
  // d3 -- distance of direction normal to patch to insert
  // numG -- number of ghost cells

  // check that number of cells to insert matches
  int blkCell = (inter.Dir1EndFirst() - inter.Dir1StartFirst()) *
                (inter.Dir2EndFirst() - inter.Dir2StartFirst()) * d3;
  if (blkCell != slice.NumCells()) {
    cerr << "ERROR: Error in procBlock::PutGeomSlice(). Number of cells being "
            "inserted does not match designated space to insert to." << endl;
    cerr << "Direction 1, 2, 3 of procBlock: "
         << inter.Dir1EndFirst() - inter.Dir1StartFirst() << ", "
         << inter.Dir2EndFirst() - inter.Dir2StartFirst() << ", " << d3 << endl;
    cerr << "Direction I, J, K of geomSlice: " << slice.NumI() << ", "
         << slice.NumJ() << ", " << slice.NumK() << endl;
    exit(0);
  }

  // get procBlock maxes
  int imaxB = (*this).NumI() + 2.0 * numGhosts_;
  int jmaxB = (*this).NumJ() + 2.0 * numGhosts_;

  // get slice maxes
  int imaxS = slice.NumI();
  int jmaxS = slice.NumJ();

  // adjust insertion indices if patch borders another interblock on the same
  // surface of the block
  int adjS1 = (inter.Dir1StartInterBorderFirst()) ? numG : 0;
  int adjE1 = (inter.Dir1EndInterBorderFirst()) ? numG : 0;
  int adjS2 = (inter.Dir2StartInterBorderFirst()) ? numG : 0;
  int adjE2 = (inter.Dir2EndInterBorderFirst()) ? numG : 0;
  vector<bool> adjEdge(4, false);  // initialize all return values to false

  // determine if area direction needs to be reversed
  double aFac3 =
      ((inter.BoundaryFirst() + inter.BoundarySecond()) % 2 == 0) ? -1.0 : 1.0;
  double aFac1 = (inter.Orientation() == 3 || inter.Orientation() == 4 ||
                  inter.Orientation() == 7 || inter.Orientation() == 8)
                     ? -1.0
                     : 1.0;
  double aFac2 = (inter.Orientation() == 5 || inter.Orientation() == 6 ||
                  inter.Orientation() == 7 || inter.Orientation() == 8)
                     ? -1.0
                     : 1.0;

  // loop over cells to insert
  for (int l3 = 0; l3 < d3; l3++) {
    for (int l2 = adjS2;
         l2 < (inter.Dir2EndFirst() - inter.Dir2StartFirst() - adjE2); l2++) {
      for (int l1 = adjS1;
           l1 < (inter.Dir1EndFirst() - inter.Dir1StartFirst() - adjE1); l1++) {
        // get block and slice indices
        vector3d<int> indB = GetSwapLoc(l1, l2, l3, inter, true);
        vector3d<int> indS = GetSwapLoc(l1, l2, l3, inter, false);

        // get cell locations
        int locB = GetLoc1D(indB[0], indB[1], indB[2], imaxB, jmaxB);
        int locS = GetLoc1D(indS[0], indS[1], indS[2], imaxS, jmaxS);

        // get i-face locations
        int IlowB = GetLowerFaceI(indB[0], indB[1], indB[2], imaxB, jmaxB);
        int IupB = GetUpperFaceI(indB[0], indB[1], indB[2], imaxB, jmaxB);

        int IlowS = GetLowerFaceI(indS[0], indS[1], indS[2], imaxS, jmaxS);
        int IupS = GetUpperFaceI(indS[0], indS[1], indS[2], imaxS, jmaxS);

        // get j-face locations
        int JlowB = GetLowerFaceJ(indB[0], indB[1], indB[2], imaxB, jmaxB);
        int JupB = GetUpperFaceJ(indB[0], indB[1], indB[2], imaxB, jmaxB);

        int JlowS = GetLowerFaceJ(indS[0], indS[1], indS[2], imaxS, jmaxS);
        int JupS = GetUpperFaceJ(indS[0], indS[1], indS[2], imaxS, jmaxS);

        // get k-face locations
        int KlowB = GetLowerFaceK(indB[0], indB[1], indB[2], imaxB, jmaxB);
        int KupB = GetUpperFaceK(indB[0], indB[1], indB[2], imaxB, jmaxB);

        int KlowS = GetLowerFaceK(indS[0], indS[1], indS[2], imaxS, jmaxS);
        int KupS = GetUpperFaceK(indS[0], indS[1], indS[2], imaxS, jmaxS);

        // don't overwrite with garbage from partner block that hasn't recieved
        // its ghost value yet (needed at "t" intersection)
        if (slice.Vol(locS) == 0.0) {
          // find out if on edge, if so save edge
          string edgeDir;

          if ((*this).AtEdge(indB[0], indB[1], indB[2], true,
                             edgeDir)) {  // at a block edge -- possible need to
                                          // adjust interblock

            int dir1, dir2;
            if (inter.Direction1First() == "i") {
              dir1 = 0;  // direction 1 is i
              dir2 = 1;  // direction 2 is j
            } else if (inter.Direction1First() == "j") {
              dir1 = 1;  // direction 1 is j
              dir2 = 2;  // direction 2 is k
            } else {
              dir1 = 2;  // direction 1 is k
              dir2 = 0;  // direction 2 is i
            }

            // find out edge direction
            if (edgeDir == inter.Direction1First()) {  // edge direction matches
                                                       // interblock direction 1
              if (indB[dir2] <
                  inter.Dir2StartFirst() +
                      numGhosts_) {  // adjust edge on lower dir2 side
                adjEdge[2] = true;
              } else {  // adjust edge on upper dir2 side
                adjEdge[3] = true;
              }
            } else if (edgeDir ==
                       inter.Direction2First()) {  // edge direction matches
                                                   // interblock direction 2
              if (indB[dir1] <
                  inter.Dir1StartFirst() +
                      numGhosts_) {  // adjust edge on lower dir1 side
                adjEdge[0] = true;
              } else {  // adjust edge on upper dir1 side
                adjEdge[1] = true;
              }
            } else {
              cerr << "ERROR: Error in procBlock::PutStateSlice(). Ghost cell "
                      "edge direction does not match interblock direction 1 or "
                      "2." << endl;
              exit(0);
            }
          }

        // volume is not 0, ok to overwrite variables
        } else {
          // swap cell data
          (*this).vol_[locB] = slice.Vol(locS);
          (*this).center_[locB] = slice.Center(locS);

          //-------------------------------------------------------------------
          // swap face data
          if (inter.Direction3First() == "i" &&
              inter.Direction3Second() ==
                  "i") {  // both patches i, i to i, j to j, k to k
            // swap face data for direction 3
            (*this).fCenterI_[IlowB] = slice.FCenterI(IlowS);
            (*this).fAreaI_[IlowB] = aFac3 * slice.FAreaI(IlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterI_[IupB] = slice.FCenterI(IupS);
              (*this).fAreaI_[IupB] = aFac3 * slice.FAreaI(IupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaJ(JlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterJ(JupS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterJ(JupS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaJ(JupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterJ(JlowS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaJ(JlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterK(KlowS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaK(KlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterK(KupS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterK(KupS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaK(KupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterK(KlowS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaK(KlowS);
              }
            }

          //----------------------------------------------------------------
          } else if (inter.Direction3First() == "j" &&
                   inter.Direction3Second() ==
                       "j") {  // both patches j, j to j, k to k, i to i
            // swap face data for direction 3
            (*this).fCenterJ_[JlowB] = slice.FCenterJ(JlowS);
            (*this).fAreaJ_[JlowB] = aFac3 * slice.FAreaJ(JlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterJ_[JupB] = slice.FCenterJ(JupS);
              (*this).fAreaJ_[JupB] = aFac3 * slice.FAreaJ(JupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterK(KlowS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaK(KlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterK(KupS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterK(KupS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaK(KupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterK(KlowS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaK(KlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterI(IlowS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaI(IlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterI(IupS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterI(IupS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaI(IupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterI(IlowS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaI(IlowS);
              }
            }

          //-------------------------------------------------------------
          } else if (inter.Direction3First() == "k" &&
                   inter.Direction3Second() ==
                       "k") {  // both patches k, k to k, i to i, j to j
            // swap face data for direction 3
            (*this).fCenterK_[KlowB] = slice.FCenterK(KlowS);
            (*this).fAreaK_[KlowB] = aFac3 * slice.FAreaK(KlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterK_[KupB] = slice.FCenterK(KupS);
              (*this).fAreaK_[KupB] = aFac3 * slice.FAreaK(KupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterI(IlowS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaI(IlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterI(IupS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterI(IupS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaI(IupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterI(IlowS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaI(IlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaJ(JlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterJ(JupS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterJ(JupS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaJ(JupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterJ(JlowS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaJ(JlowS);
              }
            }

          //---------------------------------------------------------------
          } else if (inter.Direction3First() == "i" &&
                   inter.Direction3Second() ==
                       "j") {  // patches are i/j  - i to j, j to k, k to i
            // swap face data for direction 3
            (*this).fCenterI_[IlowB] = slice.FCenterJ(JlowS);
            (*this).fAreaI_[IlowB] = aFac3 * slice.FAreaJ(JlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterI_[IupB] = slice.FCenterJ(JupS);
              (*this).fAreaI_[IupB] = aFac3 * slice.FAreaJ(JupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterK(KlowS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaK(KlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterK(KupS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterK(KupS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaK(KupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterK(KlowS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaK(KlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterI(IlowS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaI(IlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterI(IupS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterI(IupS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaI(IupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterI(IlowS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaI(IlowS);
              }
            }

          //------------------------------------------------------------------
          } else if (inter.Direction3First() == "i" &&
                   inter.Direction3Second() ==
                       "k") {  // patches are i/k  - i to k, j to i, k to j
            // swap face data for direction 3
            (*this).fCenterI_[IlowB] = slice.FCenterK(KlowS);
            (*this).fAreaI_[IlowB] = aFac3 * slice.FAreaK(KlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterI_[IupB] = slice.FCenterK(KupS);
              (*this).fAreaI_[IupB] = aFac3 * slice.FAreaK(KupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterI(IlowS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaI(IlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterI(IupS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterI(IupS);
              (*this).fAreaJ_[JlowB] = aFac1 * slice.FAreaI(IupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterJ_[JupB] = slice.FCenterI(IlowS);
                (*this).fAreaJ_[JupB] = aFac1 * slice.FAreaI(IlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaJ(JlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterJ(JupS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterJ(JupS);
              (*this).fAreaK_[KlowB] = aFac2 * slice.FAreaJ(JupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterK_[KupB] = slice.FCenterJ(JlowS);
                (*this).fAreaK_[KupB] = aFac2 * slice.FAreaJ(JlowS);
              }
            }

          //--------------------------------------------------------------
          } else if (inter.Direction3First() == "j" &&
                   inter.Direction3Second() ==
                       "i") {  // patches are j/i, j to i, k to j, i to k
            // swap face data for direction 3
            (*this).fCenterJ_[JlowB] = slice.FCenterI(IlowS);
            (*this).fAreaJ_[JlowB] = aFac3 * slice.FAreaI(IlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterJ_[JupB] = slice.FCenterI(IupS);
              (*this).fAreaJ_[JupB] = aFac3 * slice.FAreaI(IupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaJ(JlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterJ(JupS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterJ(JupS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaJ(JupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterJ(JlowS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaJ(JlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterK(KlowS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaK(KlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterK(KupS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterK(KupS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaK(KupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterK(KlowS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaK(KlowS);
              }
            }


          //-----------------------------------------------------------------
          } else if (inter.Direction3First() > "j" &&
                   inter.Direction3Second() ==
                       "k") {  // patches are j/k, j to k, k to i, i to j
            // swap face data for direction 3
            (*this).fCenterJ_[JlowB] = slice.FCenterK(KlowS);
            (*this).fAreaJ_[JlowB] = aFac3 * slice.FAreaK(KlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterJ_[JupB] = slice.FCenterJ(JupS);
              (*this).fAreaJ_[JupB] = aFac3 * slice.FAreaJ(JupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterK_[KlowB] = slice.FCenterI(IlowS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaI(IlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterI(IupS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterK_[KlowB] = slice.FCenterI(IupS);
              (*this).fAreaK_[KlowB] = aFac1 * slice.FAreaI(IupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterK_[KupB] = slice.FCenterI(IlowS);
                (*this).fAreaK_[KupB] = aFac1 * slice.FAreaI(IlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaJ(JlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterJ(JupS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterJ(JupS);
              (*this).fAreaI_[IlowB] = aFac2 * slice.FAreaJ(JupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterI_[IupB] = slice.FCenterJ(JlowS);
                (*this).fAreaI_[IupB] = aFac2 * slice.FAreaJ(JlowS);
              }
            }

          //------------------------------------------------------------------
          } else if (inter.Direction3First() == "k" &&
                   inter.Direction3Second() ==
                       "i") {  // patches are k/i, k to i, i to j, j to k
            // swap face data for direction 3
            (*this).fCenterK_[KlowB] = slice.FCenterI(IlowS);
            (*this).fAreaK_[KlowB] = aFac3 * slice.FAreaI(IlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterK_[KupB] = slice.FCenterI(IupS);
              (*this).fAreaK_[KupB] = aFac3 * slice.FAreaI(IupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterJ(JlowS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaJ(JlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterJ(JupS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaJ(JupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterJ(JupS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaJ(JupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterJ(JlowS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaJ(JlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterK(KlowS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaK(KlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterK(KupS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterK(KupS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaK(KupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterK(KlowS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaK(KlowS);
              }
            }

          //-------------------------------------------------------------------
          } else if (inter.Direction3First() == "k" &&
                   inter.Direction3Second() ==
                       "j") {  // patches are k/j, k to j, i to k, j to i
            // swap face data for direction 3
            (*this).fCenterK_[KlowB] = slice.FCenterJ(JlowS);
            (*this).fAreaK_[KlowB] = aFac3 * slice.FAreaJ(JlowS);

            if (l3 == (d3 - 1)) {  // at end of direction 3 line
              (*this).fCenterK_[KupB] = slice.FCenterJ(JupS);
              (*this).fAreaK_[KupB] = aFac3 * slice.FAreaJ(JupS);
            }

            // swap face data for direction 1
            if (aFac1 == 1.0) {
              (*this).fCenterI_[IlowB] = slice.FCenterK(KlowS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaK(KlowS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterK(KupS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaK(KupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterI_[IlowB] = slice.FCenterK(KupS);
              (*this).fAreaI_[IlowB] = aFac1 * slice.FAreaK(KupS);

              if (l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() -
                         1)) {  // at end of direction 1 line
                (*this).fCenterI_[IupB] = slice.FCenterK(KlowS);
                (*this).fAreaI_[IupB] = aFac1 * slice.FAreaK(KlowS);
              }
            }

            // swap face data for direction 2
            if (aFac2 == 1.0) {
              (*this).fCenterJ_[JlowB] = slice.FCenterI(IlowS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaI(IlowS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterI(IupS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaI(IupS);
              }
            } else {  // if direction is reversed, upper/lower faces need to be
                      // swapped
              (*this).fCenterJ_[JlowB] = slice.FCenterI(IupS);
              (*this).fAreaJ_[JlowB] = aFac2 * slice.FAreaI(IupS);

              if (l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() -
                         1)) {  // at end of direction 2 line
                (*this).fCenterJ_[JupB] = slice.FCenterI(IlowS);
                (*this).fAreaJ_[JupB] = aFac2 * slice.FAreaI(IlowS);
              }
            }

          //-------------------------------------------------------------------
          } else {
            cerr << "ERROR: Error in procBlock::PutGeomSlice(). Unable to swap "
                    "face quantities because behavior for interface with "
                    "boundary pair " << inter.BoundaryFirst() << ", "
                 << inter.BoundarySecond() << " is not defined." << endl;
            exit(0);
          }
        }
      }
    }
  }

  return adjEdge;
}

/* Member function to overwrite a section of a procBlock's states with a
stateSlice. The function uses the orientation supplied in the interblock to
orient the stateSlice relative to the procBlock. It assumes that the procBlock
is listed first, and the stateSlice second in the interblock data structure.
*/
void procBlock::PutStateSlice(const stateSlice &slice, const interblock &inter,
                              const int &d3, const int &numG) {
  // slice -- geomSlice to insert int procBlock
  // inter -- interblock data structure defining the patches and their
  // orientation
  // d3 -- distance of direction normal to patch to insert
  // numG -- number of ghost cells

  // check that number of cells to insert matches
  int blkCell = (inter.Dir1EndFirst() - inter.Dir1StartFirst()) *
                (inter.Dir2EndFirst() - inter.Dir2StartFirst()) * d3;
  if (blkCell != slice.NumCells()) {
    cerr << "ERROR: Error in procBlock::PutStateSlice(). Number of cells being "
            "inserted does not match designated space to insert to." << endl;
    cerr << "Direction 1, 2, 3 of procBlock: "
         << inter.Dir1EndFirst() - inter.Dir1StartFirst() << ", "
         << inter.Dir2EndFirst() - inter.Dir2StartFirst() << ", " << d3 << endl;
    cerr << "Direction I, J, K of geomSlice: " << slice.NumI() << ", "
         << slice.NumJ() << ", " << slice.NumK() << endl;
    exit(0);
  }

  // get block maxes
  int imaxB = (*this).NumI() + 2.0 * numGhosts_;
  int jmaxB = (*this).NumJ() + 2.0 * numGhosts_;

  // get slice maxes
  int imaxS = slice.NumI();
  int jmaxS = slice.NumJ();

  // adjust insertion indices if patch borders another interblock on the same
  // surface of the block
  int adjS1 = (inter.Dir1StartInterBorderFirst()) ? numG : 0;
  int adjE1 = (inter.Dir1EndInterBorderFirst()) ? numG : 0;
  int adjS2 = (inter.Dir2StartInterBorderFirst()) ? numG : 0;
  int adjE2 = (inter.Dir2EndInterBorderFirst()) ? numG : 0;

  // loop over cells to insert
  for (int l3 = 0; l3 < d3; l3++) {
    for (int l2 = adjS2;
         l2 < (inter.Dir2EndFirst() - inter.Dir2StartFirst() - adjE2); l2++) {
      for (int l1 = adjS1;
           l1 < (inter.Dir1EndFirst() - inter.Dir1StartFirst() - adjE1); l1++) {
        // get block and slice indices
        vector3d<int> indB = GetSwapLoc(l1, l2, l3, inter, true);
        vector3d<int> indS = GetSwapLoc(l1, l2, l3, inter, false);

        // get cell locations
        int locB = GetLoc1D(indB[0], indB[1], indB[2], imaxB, jmaxB);
        int locS = GetLoc1D(indS[0], indS[1], indS[2], imaxS, jmaxS);

        // swap cell data
        (*this).state_[locB] = slice.State(locS);
      }
    }
  }
}

/*Member function to pack and send procBlock geometry data to appropriate
 * processor. */
void procBlock::PackSendGeomMPI(const MPI_Datatype &MPI_cellData,
                                const MPI_Datatype &MPI_vec3d,
                                const MPI_Datatype &MPI_vec3dMag) const {
  // MPI_cellData -- MPI data type for cell data
  // MPI_vec3d -- MPI data type for a vector3d
  // MPI_vec3dMag -- MPI data type for a unitVect3dMag

  // determine size of buffer to send
  int sendBufSize = 0;
  int tempSize = 0;
  MPI_Pack_size(10, MPI_INT, MPI_COMM_WORLD,
                &tempSize);  // add size for ints in class procBlock
  sendBufSize += tempSize;
  MPI_Pack_size((*this).state_.size(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for states
  sendBufSize += tempSize;
  MPI_Pack_size((*this).center_.size(), MPI_vec3d, MPI_COMM_WORLD,
                &tempSize);  // add size for cell centers
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fAreaI_.size(), MPI_vec3dMag, MPI_COMM_WORLD,
                &tempSize);  // add size for face area I
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fAreaJ_.size(), MPI_vec3dMag, MPI_COMM_WORLD,
                &tempSize);  // add size for face area J
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fAreaK_.size(), MPI_vec3dMag, MPI_COMM_WORLD,
                &tempSize);  // add size for face area K
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fCenterI_.size(), MPI_vec3d, MPI_COMM_WORLD,
                &tempSize);  // add size for face center_ I
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fCenterJ_.size(), MPI_vec3d, MPI_COMM_WORLD,
                &tempSize);  // add size for face center_ J
  sendBufSize += tempSize;
  MPI_Pack_size((*this).fCenterK_.size(), MPI_vec3d, MPI_COMM_WORLD,
                &tempSize);  // add size for face center_ K
  sendBufSize += tempSize;
  MPI_Pack_size((*this).vol_.size(), MPI_DOUBLE, MPI_COMM_WORLD,
                &tempSize);  // add size for volumes
  sendBufSize += tempSize;
  MPI_Pack_size(3, MPI_INT, MPI_COMM_WORLD,
                &tempSize);  // add size for number of surfaces
  sendBufSize += tempSize;
  // 8x because iMin, iMax, jMin, jMax, kMin, kMax, tags, string sizes
  MPI_Pack_size((*this).bc_.NumSurfaces() * 8, MPI_INT, MPI_COMM_WORLD,
                &tempSize);  // add size for BCs
  sendBufSize += tempSize;

  int stringSize = 0;
  for (int jj = 0; jj < (*this).bc_.NumSurfaces(); jj++) {
    MPI_Pack_size(
        (*this).bc_.GetBCTypes(jj).size() + 1, MPI_CHAR, MPI_COMM_WORLD,
        &tempSize);  // add size for bc_ types (+1 for c_str end character)
    stringSize += tempSize;
  }
  sendBufSize += stringSize;

  char *sendBuffer = new char[sendBufSize];  // allocate buffer to pack data
                                             // into

  // pack data to send into buffer
  int position = 0;
  // int and vector data
  MPI_Pack(&(*this).numCells_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numVars_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numI_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numJ_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numK_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numGhosts_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlock_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).rank_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).localPos_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).globalPos_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).state_[0], (*this).state_.size(), MPI_cellData, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).center_[0], (*this).center_.size(), MPI_vec3d, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fAreaI_[0], (*this).fAreaI_.size(), MPI_vec3dMag,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fAreaJ_[0], (*this).fAreaJ_.size(), MPI_vec3dMag,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fAreaK_[0], (*this).fAreaK_.size(), MPI_vec3dMag,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fCenterI_[0], (*this).fCenterI_.size(), MPI_vec3d,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fCenterJ_[0], (*this).fCenterJ_.size(), MPI_vec3d,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).fCenterK_[0], (*this).fCenterK_.size(), MPI_vec3d,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).vol_[0], (*this).vol_.size(), MPI_DOUBLE, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);

  // pack boundary condition data
  (*this).bc_.PackBC(sendBuffer, sendBufSize, position);

  // send buffer to appropriate processor
  MPI_Send(sendBuffer, sendBufSize, MPI_PACKED, (*this).rank_, 2,
           MPI_COMM_WORLD);

  delete[] sendBuffer;  // deallocate buffer
}

void procBlock::RecvUnpackGeomMPI(const MPI_Datatype &MPI_cellData,
                                  const MPI_Datatype &MPI_vec3d,
                                  const MPI_Datatype &MPI_vec3dMag) {
  // MPI_cellData -- MPI data type for cell data
  // MPI_vec3d -- MPI data type for a vector3d
  // MPI_vec3dMag -- MPI data type for a unitVect3dMag

  MPI_Status status;  // allocate MPI_Status structure

  // probe message to get correct data size
  int recvBufSize = 0;
  MPI_Probe(ROOTP, 2, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_CHAR, &recvBufSize);  // use MPI_CHAR because
                                                   // sending buffer was
                                                   // allocated with chars

  char *recvBuffer = new char[recvBufSize];  // allocate buffer of correct size

  // receive message from ROOT
  MPI_Recv(recvBuffer, recvBufSize, MPI_PACKED, ROOTP, 2, MPI_COMM_WORLD,
           &status);

  // unpack procBlock INTs
  int position = 0;
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numCells_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numVars_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numI_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numJ_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numK_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numGhosts_, 1,
             MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).parBlock_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).rank_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).localPos_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).globalPos_, 1,
             MPI_INT, MPI_COMM_WORLD);

  // clean and resize the vectors in the class to
  (*this).CleanResizeVecs();

  // unpack vector data into allocated vectors
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).state_[0],
             (*this).state_.size(), MPI_cellData,
             MPI_COMM_WORLD);  // unpack states
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).center_[0],
             (*this).center_.size(), MPI_vec3d,
             MPI_COMM_WORLD);  // unpack cell centers
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fAreaI_[0],
             (*this).fAreaI_.size(), MPI_vec3dMag,
             MPI_COMM_WORLD);  // unpack face area I
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fAreaJ_[0],
             (*this).fAreaJ_.size(), MPI_vec3dMag,
             MPI_COMM_WORLD);  // unpack face area J
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fAreaK_[0],
             (*this).fAreaK_.size(), MPI_vec3dMag,
             MPI_COMM_WORLD);  // unpack face area K
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fCenterI_[0],
             (*this).fCenterI_.size(), MPI_vec3d,
             MPI_COMM_WORLD);  // unpack face center_ I
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fCenterJ_[0],
             (*this).fCenterJ_.size(), MPI_vec3d,
             MPI_COMM_WORLD);  // unpack face center_ J
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).fCenterK_[0],
             (*this).fCenterK_.size(), MPI_vec3d,
             MPI_COMM_WORLD);  // unpack face center_ K
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).vol_[0],
             (*this).vol_.size(), MPI_DOUBLE,
             MPI_COMM_WORLD);  // unpack volumes

  // unpack boundary conditions
  (*this).bc_.UnpackBC(recvBuffer, recvBufSize, position);

  delete[] recvBuffer;  // deallocate receiving buffer
}

/*Member function to zero and resize the vectors in a procBlock to their
 * appropriate size given the i, j, and k dimensions.*/
void procBlock::CleanResizeVecs() {
  int numCellsGhosts = ((*this).numI_ + 2 * (*this).numGhosts_) *
                       ((*this).numJ_ + 2 * (*this).numGhosts_) *
                       ((*this).numK_ + 2 * (*this).numGhosts_);
  int numFaceI = ((*this).numI_ + 2 * (*this).numGhosts_ + 1) *
                 ((*this).numJ_ + 2 * (*this).numGhosts_) *
                 ((*this).numK_ + 2 * (*this).numGhosts_);
  int numFaceJ = ((*this).numI_ + 2 * (*this).numGhosts_) *
                 ((*this).numJ_ + 2 * (*this).numGhosts_ + 1) *
                 ((*this).numK_ + 2 * (*this).numGhosts_);
  int numFaceK = ((*this).numI_ + 2 * (*this).numGhosts_) *
                 ((*this).numJ_ + 2 * (*this).numGhosts_) *
                 ((*this).numK_ + 2 * (*this).numGhosts_ + 1);

  (*this).state_ = vector<primVars>(numCellsGhosts);

  (*this).center_ = vector<vector3d<double> >(numCellsGhosts);

  (*this).fCenterI_ = vector<vector3d<double> >(numFaceI);
  (*this).fAreaI_ = vector<unitVec3dMag<double> >(numFaceI);

  (*this).fCenterJ_ = vector<vector3d<double> >(numFaceJ);
  (*this).fAreaJ_ = vector<unitVec3dMag<double> >(numFaceJ);

  (*this).fCenterK_ = vector<vector3d<double> >(numFaceK);
  (*this).fAreaK_ = vector<unitVec3dMag<double> >(numFaceK);

  (*this).residual_ = vector<genArray>((*this).numCells_);

  (*this).vol_ = vector<double>(numCellsGhosts);

  (*this).wallDist_ = vector<double>((*this).numCells_, DEFAULTWALLDIST);

  avgWaveSpeed_ = vector<double>((*this).numCells_);
  (*this).dt_ = vector<double>((*this).numCells_);
}

/*Member function to receive and unpack procBlock state data. This is used to
 * gather the solution on the ROOT processor to write out the solution. */
void procBlock::RecvUnpackSolMPI(const MPI_Datatype &MPI_cellData) {
  // MPI_cellData -- MPI data type for cell data

  MPI_Status status;  // allocate MPI_Status structure

  // probe message to get correct data size
  int recvBufSize = 0;
  MPI_Probe((*this).rank_, (*this).globalPos_, MPI_COMM_WORLD,
            &status);  // global position used as tag because each block has a
                       // unique one
  MPI_Get_count(&status, MPI_CHAR, &recvBufSize);  // use MPI_CHAR because
                                                   // sending buffer was
                                                   // allocated with chars

  char *recvBuffer = new char[recvBufSize];  // allocate buffer of correct size

  // receive message from non-ROOT
  MPI_Recv(recvBuffer, recvBufSize, MPI_PACKED, (*this).rank_,
           (*this).globalPos_, MPI_COMM_WORLD, &status);

  // unpack vector data into allocated vectors
  int position = 0;
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).state_[0],
             (*this).state_.size(), MPI_cellData,
             MPI_COMM_WORLD);  // unpack states
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).residual_[0],
             (*this).residual_.size(), MPI_cellData,
             MPI_COMM_WORLD);  // unpack residuals
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).dt_[0],
             (*this).dt_.size(), MPI_DOUBLE,
             MPI_COMM_WORLD);  // unpack time steps
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).wallDist_[0],
             (*this).wallDist_.size(), MPI_DOUBLE,
             MPI_COMM_WORLD);  // unpack wall distance
  MPI_Unpack(recvBuffer, recvBufSize, &position, &avgWaveSpeed_[0],
             avgWaveSpeed_.size(), MPI_DOUBLE,
             MPI_COMM_WORLD);  // unpack average wave speeds

  delete[] recvBuffer;  // deallocate receiving buffer
}

/*Member function to pack and send procBlock state data to the ROOT proecessor.
 * This is used to gather the solution on the ROOT processor to write out the
 * solution. */
void procBlock::PackSendSolMPI(const MPI_Datatype &MPI_cellData) const {
  // MPI_cellData -- MPI data type for cell data

  // determine size of buffer to send
  int sendBufSize = 0;
  int tempSize = 0;
  MPI_Pack_size((*this).state_.size(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for states
  sendBufSize += tempSize;
  MPI_Pack_size((*this).residual_.size(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for residuals
  sendBufSize += tempSize;
  MPI_Pack_size((*this).dt_.size(), MPI_DOUBLE, MPI_COMM_WORLD,
                &tempSize);  // add size for time steps
  sendBufSize += tempSize;
  MPI_Pack_size((*this).wallDist_.size(), MPI_DOUBLE, MPI_COMM_WORLD,
                &tempSize);  // add size for wall distance
  sendBufSize += tempSize;
  MPI_Pack_size(avgWaveSpeed_.size(), MPI_DOUBLE, MPI_COMM_WORLD,
                &tempSize);  // add size for average wave speed
  sendBufSize += tempSize;

  char *sendBuffer = new char[sendBufSize];  // allocate buffer to pack data
                                             // into

  // pack data to send into buffer
  int position = 0;
  MPI_Pack(&(*this).state_[0], (*this).state_.size(), MPI_cellData, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).residual_[0], (*this).residual_.size(), MPI_cellData,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).dt_[0], (*this).dt_.size(), MPI_DOUBLE, sendBuffer,
           sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).wallDist_[0], (*this).wallDist_.size(), MPI_DOUBLE,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&avgWaveSpeed_[0], avgWaveSpeed_.size(), MPI_DOUBLE,
           sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);

  // send buffer to appropriate processor
  MPI_Send(sendBuffer, sendBufSize, MPI_PACKED, ROOTP, (*this).GlobalPos(),
           MPI_COMM_WORLD);

  delete[] sendBuffer;  // deallocate buffer
}

/* Member function to split a procBlock along a plane defined by a direction and
an index. The calling instance will retain the lower portion of the split,
and the returned instance will retain the upper portion of the split.
*/
procBlock procBlock::Split(const string &dir, const int &ind, const int &num,
                           vector<boundarySurface> &alteredSurf) {
  // dir -- plane to split along, either i, j, or k
  // ind -- index (face) to split at (w/o counting ghost cells)
  // num -- new block number
  // alteredSurf -- vector of surfaces whose partners will need to be altered
  // after this split

  int iMax = (*this).NumI() + 2 * numGhosts_;
  int jMax = (*this).NumJ() + 2 * numGhosts_;
  int kMax = (*this).NumK() + 2 * numGhosts_;

  boundaryConditions bound1 = (*this).BC();
  boundaryConditions bound2 =
      bound1.Split(dir, ind, (*this).ParentBlock(), num, alteredSurf);

  if (dir == "i") {  // split along i-plane
    int numI2 = (*this).NumI() - ind;
    int numI1 = (*this).NumI() - numI2;

    procBlock blk1(numI1, (*this).NumJ(), (*this).NumK(), numGhosts_);
    procBlock blk2(numI2, (*this).NumJ(), (*this).NumK(), numGhosts_);

    blk1.parBlock_ = (*this).ParentBlock();
    blk2.parBlock_ = (*this).ParentBlock();

    int iMax1 = numI1 + 2 * numGhosts_;
    int iMax2 = numI2 + 2 * numGhosts_;

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG = GetLoc1D(
              ii - numGhosts_, jj - numGhosts_,
              kk - numGhosts_, (*this).NumI(), (*this).NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //--------------------------------------------------------------
          // this portion of parent block overlaps with upper split
          if (ii >= ind) {
            int loc2 = GetLoc1D(ii - ind, jj, kk, iMax2, jMax);
            int loc2NG = GetLoc1D(
                ii - ind - numGhosts_, jj - numGhosts_,
                kk - numGhosts_, numI2, (*this).NumJ());

            int fLowI2 = GetLowerFaceI(ii - ind, jj, kk, iMax2, jMax);
            int fLowJ2 = GetLowerFaceJ(ii - ind, jj, kk, iMax2, jMax);
            int fLowK2 = GetLowerFaceK(ii - ind, jj, kk, iMax2, jMax);

            // assign cell variables
            blk2.state_[loc2] = (*this).state_[loc];
            blk2.vol_[loc2] = (*this).vol_[loc];
            blk2.center_[loc2] = (*this).center_[loc];

            if (ii >= (ind + numGhosts_) &&
                ii < (iMax - numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cells
              blk2.avgWaveSpeed_[loc2NG] = avgWaveSpeed_[locNG];
              blk2.dt_[loc2NG] = (*this).dt_[locNG];
              blk2.wallDist_[loc2NG] = (*this).wallDist_[locNG];
              blk2.residual_[loc2NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk2.fAreaI_[fLowI2] = (*this).fAreaI_[fLowI];
            blk2.fAreaJ_[fLowJ2] = (*this).fAreaJ_[fLowJ];
            blk2.fAreaK_[fLowK2] = (*this).fAreaK_[fLowK];

            blk2.fCenterI_[fLowI2] = (*this).fCenterI_[fLowI];
            blk2.fCenterJ_[fLowJ2] = (*this).fCenterJ_[fLowJ];
            blk2.fCenterK_[fLowK2] = (*this).fCenterK_[fLowK];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI2 = GetUpperFaceI(ii - ind, jj, kk, iMax2, jMax);

              blk2.fAreaI_[fUpI2] = (*this).fAreaI_[fUpI];
              blk2.fCenterI_[fUpI2] = (*this).fCenterI_[fUpI];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ2 = GetUpperFaceJ(ii - ind, jj, kk, iMax2, jMax);

              blk2.fAreaJ_[fUpJ2] = (*this).fAreaJ_[fUpJ];
              blk2.fCenterJ_[fUpJ2] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK2 = GetUpperFaceK(ii - ind, jj, kk, iMax2, jMax);

              blk2.fAreaK_[fUpK2] = (*this).fAreaK_[fUpK];
              blk2.fCenterK_[fUpK2] = (*this).fCenterK_[fUpK];
            }
          }

          //------------------------------------------------------------------
          if (ii < ind + 2 * numGhosts_) {  // this portion of parent
                                                     // block overlaps with
                                                     // lower split

            int loc1 = GetLoc1D(ii, jj, kk, iMax1, jMax);
            int loc1NG =
                GetLoc1D(ii - numGhosts_, jj - numGhosts_,
                         kk - numGhosts_, numI1, (*this).NumJ());

            int fLowI1 = GetLowerFaceI(ii, jj, kk, iMax1, jMax);
            int fLowJ1 = GetLowerFaceJ(ii, jj, kk, iMax1, jMax);
            int fLowK1 = GetLowerFaceK(ii, jj, kk, iMax1, jMax);

            // assign cell variables
            blk1.state_[loc1] = (*this).state_[loc];
            blk1.vol_[loc1] = (*this).vol_[loc];
            blk1.center_[loc1] = (*this).center_[loc];

            if (ii >= numGhosts_ && ii < (ind + numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cell
              blk1.avgWaveSpeed_[loc1NG] = avgWaveSpeed_[locNG];
              blk1.dt_[loc1NG] = (*this).dt_[locNG];
              blk1.wallDist_[loc1NG] = (*this).wallDist_[locNG];
              blk1.residual_[loc1NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk1.fAreaI_[fLowI1] = (*this).fAreaI_[fLowI];
            blk1.fAreaJ_[fLowJ1] = (*this).fAreaJ_[fLowJ];
            blk1.fAreaK_[fLowK1] = (*this).fAreaK_[fLowK];

            blk1.fCenterI_[fLowI1] = (*this).fCenterI_[fLowI];
            blk1.fCenterJ_[fLowJ1] = (*this).fCenterJ_[fLowJ];
            blk1.fCenterK_[fLowK1] = (*this).fCenterK_[fLowK];

            if (ii == ind + numGhosts_ -
                          1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI1 = GetUpperFaceI(ii, jj, kk, iMax1, jMax);

              blk1.fAreaI_[fUpI1] = (*this).fAreaI_[fUpI];
              blk1.fCenterI_[fUpI1] = (*this).fCenterI_[fUpI];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ1 = GetUpperFaceJ(ii, jj, kk, iMax1, jMax);

              blk1.fAreaJ_[fUpJ1] = (*this).fAreaJ_[fUpJ];
              blk1.fCenterJ_[fUpJ1] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK1 = GetUpperFaceK(ii, jj, kk, iMax1, jMax);

              blk1.fAreaK_[fUpK1] = (*this).fAreaK_[fUpK];
              blk1.fCenterK_[fUpK1] = (*this).fCenterK_[fUpK];
            }
          }
        }
      }
    }

    blk1.bc_ = bound1;
    (*this) = blk1;
    blk2.bc_ = bound2;
    return blk2;

  } else if (dir == "j") {  // split along j-plane
    int numJ2 = (*this).NumJ() - ind;
    int numJ1 = (*this).NumJ() - numJ2;

    procBlock blk1((*this).NumI(), numJ1, (*this).NumK(), numGhosts_);
    procBlock blk2((*this).NumI(), numJ2, (*this).NumK(), numGhosts_);

    blk1.parBlock_ = (*this).ParentBlock();
    blk2.parBlock_ = (*this).ParentBlock();

    int jMax1 = numJ1 + 2 * numGhosts_;
    int jMax2 = numJ2 + 2 * numGhosts_;

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG = GetLoc1D(
              ii - numGhosts_, jj - numGhosts_,
              kk - numGhosts_, (*this).NumI(), (*this).NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //-------------------------------------------------------------------
          // this portion of parent block overlaps with upper split
          if (jj >= ind) {
            int loc2 = GetLoc1D(ii, jj - ind, kk, iMax, jMax2);
            int loc2NG = GetLoc1D(
                ii - numGhosts_, jj - ind - numGhosts_,
                kk - numGhosts_, (*this).NumI(), numJ2);

            int fLowI2 = GetLowerFaceI(ii, jj - ind, kk, iMax, jMax2);
            int fLowJ2 = GetLowerFaceJ(ii, jj - ind, kk, iMax, jMax2);
            int fLowK2 = GetLowerFaceK(ii, jj - ind, kk, iMax, jMax2);

            // assign cell variables
            blk2.state_[loc2] = (*this).state_[loc];
            blk2.vol_[loc2] = (*this).vol_[loc];
            blk2.center_[loc2] = (*this).center_[loc];

            if (jj >= (ind + numGhosts_) &&
                jj < (jMax - numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cells
              blk2.avgWaveSpeed_[loc2NG] = avgWaveSpeed_[locNG];
              blk2.dt_[loc2NG] = (*this).dt_[locNG];
              blk2.wallDist_[loc2NG] = (*this).wallDist_[locNG];
              blk2.residual_[loc2NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk2.fAreaI_[fLowI2] = (*this).fAreaI_[fLowI];
            blk2.fAreaJ_[fLowJ2] = (*this).fAreaJ_[fLowJ];
            blk2.fAreaK_[fLowK2] = (*this).fAreaK_[fLowK];

            blk2.fCenterI_[fLowI2] = (*this).fCenterI_[fLowI];
            blk2.fCenterJ_[fLowJ2] = (*this).fCenterJ_[fLowJ];
            blk2.fCenterK_[fLowK2] = (*this).fCenterK_[fLowK];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI2 = GetUpperFaceI(ii, jj - ind, kk, iMax, jMax2);

              blk2.fAreaI_[fUpI2] = (*this).fAreaI_[fUpI];
              blk2.fCenterI_[fUpI2] = (*this).fCenterI_[fUpI];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ2 = GetUpperFaceJ(ii, jj - ind, kk, iMax, jMax2);

              blk2.fAreaJ_[fUpJ2] = (*this).fAreaJ_[fUpJ];
              blk2.fCenterJ_[fUpJ2] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK2 = GetUpperFaceK(ii, jj - ind, kk, iMax, jMax2);

              blk2.fAreaK_[fUpK2] = (*this).fAreaK_[fUpK];
              blk2.fCenterK_[fUpK2] = (*this).fCenterK_[fUpK];
            }
          }

          //------------------------------------------------------------------
          if (jj < ind + 2 * numGhosts_) {  // this portion of parent
                                                     // block overlaps with
                                                     // lower split

            int loc1 = GetLoc1D(ii, jj, kk, iMax, jMax1);
            int loc1NG =
                GetLoc1D(ii - numGhosts_, jj - numGhosts_,
                         kk - numGhosts_, (*this).NumI(), numJ1);

            int fLowI1 = GetLowerFaceI(ii, jj, kk, iMax, jMax1);
            int fLowJ1 = GetLowerFaceJ(ii, jj, kk, iMax, jMax1);
            int fLowK1 = GetLowerFaceK(ii, jj, kk, iMax, jMax1);

            // assign cell variables
            blk1.state_[loc1] = (*this).state_[loc];
            blk1.vol_[loc1] = (*this).vol_[loc];
            blk1.center_[loc1] = (*this).center_[loc];

            if (jj >= numGhosts_ && jj < (ind + numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cell
              blk1.avgWaveSpeed_[loc1NG] = avgWaveSpeed_[locNG];
              blk1.dt_[loc1NG] = (*this).dt_[locNG];
              blk1.wallDist_[loc1NG] = (*this).wallDist_[locNG];
              blk1.residual_[loc1NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk1.fAreaI_[fLowI1] = (*this).fAreaI_[fLowI];
            blk1.fAreaJ_[fLowJ1] = (*this).fAreaJ_[fLowJ];
            blk1.fAreaK_[fLowK1] = (*this).fAreaK_[fLowK];

            blk1.fCenterI_[fLowI1] = (*this).fCenterI_[fLowI];
            blk1.fCenterJ_[fLowJ1] = (*this).fCenterJ_[fLowJ];
            blk1.fCenterK_[fLowK1] = (*this).fCenterK_[fLowK];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI1 = GetUpperFaceI(ii, jj, kk, iMax, jMax1);

              blk1.fAreaI_[fUpI1] = (*this).fAreaI_[fUpI];
              blk1.fCenterI_[fUpI1] = (*this).fCenterI_[fUpI];
            }

            if (jj == ind + numGhosts_ -
                          1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ1 = GetUpperFaceJ(ii, jj, kk, iMax, jMax1);

              blk1.fAreaJ_[fUpJ1] = (*this).fAreaJ_[fUpJ];
              blk1.fCenterJ_[fUpJ1] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK1 = GetUpperFaceK(ii, jj, kk, iMax, jMax1);

              blk1.fAreaK_[fUpK1] = (*this).fAreaK_[fUpK];
              blk1.fCenterK_[fUpK1] = (*this).fCenterK_[fUpK];
            }
          }
        }
      }
    }

    blk1.bc_ = bound1;
    (*this) = blk1;
    blk2.bc_ = bound2;
    return blk2;

  } else if (dir == "k") {  // split along k-plane
    int numK2 = (*this).NumK() - ind;
    int numK1 = (*this).NumK() - numK2;

    procBlock blk1((*this).NumI(), (*this).NumJ(), numK1, numGhosts_);
    procBlock blk2((*this).NumI(), (*this).NumJ(), numK2, numGhosts_);

    blk1.parBlock_ = (*this).ParentBlock();
    blk2.parBlock_ = (*this).ParentBlock();

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG = GetLoc1D(
              ii - numGhosts_, jj - numGhosts_,
              kk - numGhosts_, (*this).NumI(), (*this).NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //---------------------------------------------------------------
          // this portion of parent block overlaps with upper split
          if (kk >= ind) {
            int loc2 = GetLoc1D(ii, jj, kk - ind, iMax, jMax);
            int loc2NG = GetLoc1D(
                ii - numGhosts_, jj - numGhosts_,
                kk - ind - numGhosts_, (*this).NumI(), (*this).NumJ());

            int fLowI2 = GetLowerFaceI(ii, jj, kk - ind, iMax, jMax);
            int fLowJ2 = GetLowerFaceJ(ii, jj, kk - ind, iMax, jMax);
            int fLowK2 = GetLowerFaceK(ii, jj, kk - ind, iMax, jMax);

            // assign cell variables
            blk2.state_[loc2] = (*this).state_[loc];
            blk2.vol_[loc2] = (*this).vol_[loc];
            blk2.center_[loc2] = (*this).center_[loc];

            if (kk >= (ind + numGhosts_) &&
                kk < (kMax - numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_)) {  // physical cells
              blk2.avgWaveSpeed_[loc2NG] = avgWaveSpeed_[locNG];
              blk2.dt_[loc2NG] = (*this).dt_[locNG];
              blk2.wallDist_[loc2NG] = (*this).wallDist_[locNG];
              blk2.residual_[loc2NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk2.fAreaI_[fLowI2] = (*this).fAreaI_[fLowI];
            blk2.fAreaJ_[fLowJ2] = (*this).fAreaJ_[fLowJ];
            blk2.fAreaK_[fLowK2] = (*this).fAreaK_[fLowK];

            blk2.fCenterI_[fLowI2] = (*this).fCenterI_[fLowI];
            blk2.fCenterJ_[fLowJ2] = (*this).fCenterJ_[fLowJ];
            blk2.fCenterK_[fLowK2] = (*this).fCenterK_[fLowK];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI2 = GetUpperFaceI(ii, jj, kk - ind, iMax, jMax);

              blk2.fAreaI_[fUpI2] = (*this).fAreaI_[fUpI];
              blk2.fCenterI_[fUpI2] = (*this).fCenterI_[fUpI];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ2 = GetUpperFaceJ(ii, jj, kk - ind, iMax, jMax);

              blk2.fAreaJ_[fUpJ2] = (*this).fAreaJ_[fUpJ];
              blk2.fCenterJ_[fUpJ2] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK2 = GetUpperFaceK(ii, jj, kk - ind, iMax, jMax);

              blk2.fAreaK_[fUpK2] = (*this).fAreaK_[fUpK];
              blk2.fCenterK_[fUpK2] = (*this).fCenterK_[fUpK];
            }
          }

          //-------------------------------------------------------------------
          if (kk < ind + 2 * numGhosts_) {  // this portion of parent
                                                     // block overlaps with
                                                     // lower split

            int loc1 = GetLoc1D(ii, jj, kk, iMax, jMax);
            int loc1NG = GetLoc1D(
                ii - numGhosts_, jj - numGhosts_,
                kk - numGhosts_, (*this).NumI(), (*this).NumJ());

            int fLowI1 = GetLowerFaceI(ii, jj, kk, iMax, jMax);
            int fLowJ1 = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
            int fLowK1 = GetLowerFaceK(ii, jj, kk, iMax, jMax);

            // assign cell variables
            blk1.state_[loc1] = (*this).state_[loc];
            blk1.vol_[loc1] = (*this).vol_[loc];
            blk1.center_[loc1] = (*this).center_[loc];

            if (kk >= numGhosts_ && kk < (ind + numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_)) {  // physical cell
              blk1.avgWaveSpeed_[loc1NG] = avgWaveSpeed_[locNG];
              blk1.dt_[loc1NG] = (*this).dt_[locNG];
              blk1.wallDist_[loc1NG] = (*this).wallDist_[locNG];
              blk1.residual_[loc1NG] = (*this).residual_[locNG];
            }

            // assign face variables
            blk1.fAreaI_[fLowI1] = (*this).fAreaI_[fLowI];
            blk1.fAreaJ_[fLowJ1] = (*this).fAreaJ_[fLowJ];
            blk1.fAreaK_[fLowK1] = (*this).fAreaK_[fLowK];

            blk1.fCenterI_[fLowI1] = (*this).fCenterI_[fLowI];
            blk1.fCenterJ_[fLowJ1] = (*this).fCenterJ_[fLowJ];
            blk1.fCenterK_[fLowK1] = (*this).fCenterK_[fLowK];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpI1 = GetUpperFaceI(ii, jj, kk, iMax, jMax);

              blk1.fAreaI_[fUpI1] = (*this).fAreaI_[fUpI];
              blk1.fCenterI_[fUpI1] = (*this).fCenterI_[fUpI];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJ1 = GetUpperFaceJ(ii, jj, kk, iMax, jMax);

              blk1.fAreaJ_[fUpJ1] = (*this).fAreaJ_[fUpJ];
              blk1.fCenterJ_[fUpJ1] = (*this).fCenterJ_[fUpJ];
            }

            if (kk == ind + numGhosts_ -
                          1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpK1 = GetUpperFaceK(ii, jj, kk, iMax, jMax);

              blk1.fAreaK_[fUpK1] = (*this).fAreaK_[fUpK];
              blk1.fCenterK_[fUpK1] = (*this).fCenterK_[fUpK];
            }
          }
        }
      }
    }

    blk1.bc_ = bound1;
    (*this) = blk1;
    blk2.bc_ = bound2;
    return blk2;

  } else {
    cerr << "ERROR: Error in procBlock::Split(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }
}

/* Member function to join a procBlock along a plane defined by a direction.
 * Assumes that calling instance is lower, and input instance is upper.*/
void procBlock::Join(const procBlock &blk, const string &dir,
                     vector<boundarySurface> &alteredSurf) {
  // blk -- block to join with
  // dir -- plane to split along, either i, j, or k
  // alteredSurf -- vector of surfaces whose partners will need to be altered
  // after this split

  if (dir == "i") {
    int newNumI = (*this).NumI() + blk.NumI();
    int newNumJ = (*this).NumJ();
    int newNumK = (*this).NumK();

    int iMax = newNumI + 2 * numGhosts_;
    int jMax = newNumJ + 2 * numGhosts_;
    int kMax = newNumK + 2 * numGhosts_;

    procBlock newBlk(newNumI, newNumJ, newNumK, numGhosts_);

    newBlk.bc_ = (*this).BC();
    // newBlk.bc_.Join(blk.BC(), dir, alteredSurf);

    int iMaxU = blk.NumI() + 2 * blk.NumGhosts();
    int iMaxL = (*this).NumI() + 2 * blk.NumGhosts();

    int ind = (*this).NumI();

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG =
              GetLoc1D(ii - newBlk.NumGhosts(), jj - newBlk.NumGhosts(),
                       kk - newBlk.NumGhosts(), newBlk.NumI(), newBlk.NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //------------------------------------------------------------------
          if (ii >= ind + newBlk.NumGhosts()) {  // this portion of parent block
                                                 // overlaps with upper split

            int locU = GetLoc1D(ii - ind, jj, kk, iMaxU, jMax);
            int locNGU =
                GetLoc1D(ii - ind - blk.NumGhosts(), jj - blk.NumGhosts(),
                         kk - blk.NumGhosts(), blk.NumI(), blk.NumJ());

            int fLowIU = GetLowerFaceI(ii - ind, jj, kk, iMaxU, jMax);
            int fLowJU = GetLowerFaceJ(ii - ind, jj, kk, iMaxU, jMax);
            int fLowKU = GetLowerFaceK(ii - ind, jj, kk, iMaxU, jMax);

            // assign cell variables
            newBlk.state_[loc] = blk.state_[locU];
            newBlk.vol_[loc] = blk.vol_[locU];
            newBlk.center_[loc] = blk.center_[locU];

            if (ii >= (ind + blk.NumGhosts()) &&
                ii < (iMax - blk.NumGhosts()) && jj >= blk.NumGhosts() &&
                jj < (jMax - blk.NumGhosts()) && kk >= blk.NumGhosts() &&
                kk < (kMax - blk.NumGhosts())) {  // physical cells

              newBlk.avgWaveSpeed_[locNG] = blk.avgWaveSpeed_[locNGU];
              newBlk.dt_[locNG] = blk.dt_[locNGU];
              newBlk.wallDist_[locNG] = blk.wallDist_[locNGU];
              newBlk.residual_[locNG] = blk.residual_[locNGU];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = blk.fAreaI_[fLowIU];
            newBlk.fAreaJ_[fLowJ] = blk.fAreaJ_[fLowJU];
            newBlk.fAreaK_[fLowK] = blk.fAreaK_[fLowKU];

            newBlk.fCenterI_[fLowI] = blk.fCenterI_[fLowIU];
            newBlk.fCenterJ_[fLowJ] = blk.fCenterJ_[fLowJU];
            newBlk.fCenterK_[fLowK] = blk.fCenterK_[fLowKU];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIU = GetUpperFaceI(ii - ind, jj, kk, iMaxU, jMax);

              newBlk.fAreaI_[fUpI] = blk.fAreaI_[fUpIU];
              newBlk.fCenterI_[fUpI] = blk.fCenterI_[fUpIU];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJU = GetUpperFaceJ(ii - ind, jj, kk, iMaxU, jMax);

              newBlk.fAreaJ_[fUpJ] = blk.fAreaJ_[fUpJU];
              newBlk.fCenterJ_[fUpJ] = blk.fCenterJ_[fUpJU];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKU = GetUpperFaceK(ii - ind, jj, kk, iMaxU, jMax);

              newBlk.fAreaK_[fUpK] = blk.fAreaK_[fUpKU];
              newBlk.fCenterK_[fUpK] = blk.fCenterK_[fUpKU];
            }

          //--------------------------------------------------------------
          } else {  // this portion of parent block overlaps with lower split
            int locL = GetLoc1D(ii, jj, kk, iMaxL, jMax);
            int locNGL = GetLoc1D(
                ii - numGhosts_, jj - numGhosts_,
                kk - numGhosts_, (*this).NumI(), (*this).NumJ());

            int fLowIL = GetLowerFaceI(ii, jj, kk, iMaxL, jMax);
            int fLowJL = GetLowerFaceJ(ii, jj, kk, iMaxL, jMax);
            int fLowKL = GetLowerFaceK(ii, jj, kk, iMaxL, jMax);

            // assign cell variables
            newBlk.state_[loc] = (*this).state_[locL];
            newBlk.vol_[loc] = (*this).vol_[locL];
            newBlk.center_[loc] = (*this).center_[locL];

            if (ii >= numGhosts_ && ii < (ind + numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cell

              newBlk.avgWaveSpeed_[locNG] = avgWaveSpeed_[locNGL];
              newBlk.dt_[locNG] = (*this).dt_[locNGL];
              newBlk.wallDist_[locNG] = (*this).wallDist_[locNGL];
              newBlk.residual_[locNG] = (*this).residual_[locNGL];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = (*this).fAreaI_[fLowIL];
            newBlk.fAreaJ_[fLowJ] = (*this).fAreaJ_[fLowJL];
            newBlk.fAreaK_[fLowK] = (*this).fAreaK_[fLowKL];

            newBlk.fCenterI_[fLowI] = (*this).fCenterI_[fLowIL];
            newBlk.fCenterJ_[fLowJ] = (*this).fCenterJ_[fLowJL];
            newBlk.fCenterK_[fLowK] = (*this).fCenterK_[fLowKL];

            if (ii == ind + numGhosts_ -
                          1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIL = GetUpperFaceI(ii, jj, kk, iMaxL, jMax);

              newBlk.fAreaI_[fUpI] = (*this).fAreaI_[fUpIL];
              newBlk.fCenterI_[fUpI] = (*this).fCenterI_[fUpIL];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJL = GetUpperFaceJ(ii, jj, kk, iMaxL, jMax);

              newBlk.fAreaJ_[fUpJ] = (*this).fAreaJ_[fUpJL];
              newBlk.fCenterJ_[fUpJ] = (*this).fCenterJ_[fUpJL];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKL = GetUpperFaceK(ii, jj, kk, iMaxL, jMax);

              newBlk.fAreaK_[fUpK] = (*this).fAreaK_[fUpKL];
              newBlk.fCenterK_[fUpK] = (*this).fCenterK_[fUpKL];
            }
          }
        }
      }
    }
    (*this) = newBlk;

  //------------------------------------------------------------------
  //------------------------------------------------------------------
  //------------------------------------------------------------------
  } else if (dir == "j") {
    int newNumI = (*this).NumI();
    int newNumJ = (*this).NumJ() + blk.NumJ();
    int newNumK = (*this).NumK();

    int iMax = newNumI + 2 * numGhosts_;
    int jMax = newNumJ + 2 * numGhosts_;
    int kMax = newNumK + 2 * numGhosts_;

    procBlock newBlk(newNumI, newNumJ, newNumK, numGhosts_);

    newBlk.bc_ = (*this).BC();
    // newBlk.bc_.Join(blk.BC(), dir, alteredSurf);

    int jMaxU = blk.NumJ() + 2 * blk.NumGhosts();
    int jMaxL = (*this).NumJ() + 2 * blk.NumGhosts();

    int ind = (*this).NumJ();

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG =
              GetLoc1D(ii - newBlk.NumGhosts(), jj - newBlk.NumGhosts(),
                       kk - newBlk.NumGhosts(), newBlk.NumI(), newBlk.NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //-------------------------------------------------------------------
          if (jj >= ind + newBlk.NumGhosts()) {  // this portion of parent block
                                                 // overlaps with upper split

            int locU = GetLoc1D(ii, jj - ind, kk, iMax, jMaxU);
            int locNGU =
                GetLoc1D(ii - blk.NumGhosts(), jj - ind - blk.NumGhosts(),
                         kk - blk.NumGhosts(), blk.NumI(), blk.NumJ());

            int fLowIU = GetLowerFaceI(ii, jj - ind, kk, iMax, jMaxU);
            int fLowJU = GetLowerFaceJ(ii, jj - ind, kk, iMax, jMaxU);
            int fLowKU = GetLowerFaceK(ii, jj - ind, kk, iMax, jMaxU);

            // assign cell variables
            newBlk.state_[loc] = blk.state_[locU];
            newBlk.vol_[loc] = blk.vol_[locU];
            newBlk.center_[loc] = blk.center_[locU];

            if (jj >= (ind + blk.NumGhosts()) &&
                jj < (jMax - blk.NumGhosts()) && ii >= blk.NumGhosts() &&
                ii < (iMax - blk.NumGhosts()) && kk >= blk.NumGhosts() &&
                kk < (kMax - blk.NumGhosts())) {  // physical cells

              newBlk.avgWaveSpeed_[locNG] = blk.avgWaveSpeed_[locNGU];
              newBlk.dt_[locNG] = blk.dt_[locNGU];
              newBlk.wallDist_[locNG] = blk.wallDist_[locNGU];
              newBlk.residual_[locNG] = blk.residual_[locNGU];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = blk.fAreaI_[fLowIU];
            newBlk.fAreaJ_[fLowJ] = blk.fAreaJ_[fLowJU];
            newBlk.fAreaK_[fLowK] = blk.fAreaK_[fLowKU];

            newBlk.fCenterI_[fLowI] = blk.fCenterI_[fLowIU];
            newBlk.fCenterJ_[fLowJ] = blk.fCenterJ_[fLowJU];
            newBlk.fCenterK_[fLowK] = blk.fCenterK_[fLowKU];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIU = GetUpperFaceI(ii, jj - ind, kk, iMax, jMaxU);

              newBlk.fAreaI_[fUpI] = blk.fAreaI_[fUpIU];
              newBlk.fCenterI_[fUpI] = blk.fCenterI_[fUpIU];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJU = GetUpperFaceJ(ii, jj - ind, kk, iMax, jMaxU);

              newBlk.fAreaJ_[fUpJ] = blk.fAreaJ_[fUpJU];
              newBlk.fCenterJ_[fUpJ] = blk.fCenterJ_[fUpJU];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKU = GetUpperFaceK(ii, jj - ind, kk, iMax, jMaxU);

              newBlk.fAreaK_[fUpK] = blk.fAreaK_[fUpKU];
              newBlk.fCenterK_[fUpK] = blk.fCenterK_[fUpKU];
            }

          //-----------------------------------------------------------------
          } else {  // this portion of parent block overlaps with lower split
            int locL = GetLoc1D(ii, jj, kk, iMax, jMaxL);
            int locNGL = GetLoc1D(
                ii - numGhosts_, jj - numGhosts_,
                kk - numGhosts_, (*this).NumI(), (*this).NumJ());

            int fLowIL = GetLowerFaceI(ii, jj, kk, iMax, jMaxL);
            int fLowJL = GetLowerFaceJ(ii, jj, kk, iMax, jMaxL);
            int fLowKL = GetLowerFaceK(ii, jj, kk, iMax, jMaxL);

            // assign cell variables
            newBlk.state_[loc] = (*this).state_[locL];
            newBlk.vol_[loc] = (*this).vol_[locL];
            newBlk.center_[loc] = (*this).center_[locL];

            if (jj >= numGhosts_ && jj < (ind + numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                kk >= numGhosts_ &&
                kk < (kMax - numGhosts_)) {  // physical cell

              newBlk.avgWaveSpeed_[locNG] = avgWaveSpeed_[locNGL];
              newBlk.dt_[locNG] = (*this).dt_[locNGL];
              newBlk.wallDist_[locNG] = (*this).wallDist_[locNGL];
              newBlk.residual_[locNG] = (*this).residual_[locNGL];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = (*this).fAreaI_[fLowIL];
            newBlk.fAreaJ_[fLowJ] = (*this).fAreaJ_[fLowJL];
            newBlk.fAreaK_[fLowK] = (*this).fAreaK_[fLowKL];

            newBlk.fCenterI_[fLowI] = (*this).fCenterI_[fLowIL];
            newBlk.fCenterJ_[fLowJ] = (*this).fCenterJ_[fLowJL];
            newBlk.fCenterK_[fLowK] = (*this).fCenterK_[fLowKL];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIL = GetUpperFaceI(ii, jj, kk, iMax, jMaxL);

              newBlk.fAreaI_[fUpI] = (*this).fAreaI_[fUpIL];
              newBlk.fCenterI_[fUpI] = (*this).fCenterI_[fUpIL];
            }

            if (jj == ind + numGhosts_ -
                          1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJL = GetUpperFaceJ(ii, jj, kk, iMax, jMaxL);

              newBlk.fAreaJ_[fUpJ] = (*this).fAreaJ_[fUpJL];
              newBlk.fCenterJ_[fUpJ] = (*this).fCenterJ_[fUpJL];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKL = GetUpperFaceK(ii, jj, kk, iMax, jMaxL);

              newBlk.fAreaK_[fUpK] = (*this).fAreaK_[fUpKL];
              newBlk.fCenterK_[fUpK] = (*this).fCenterK_[fUpKL];
            }
          }
        }
      }
    }
    (*this) = newBlk;

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  } else if (dir == "k") {
    int newNumI = (*this).NumI();
    int newNumJ = (*this).NumJ();
    int newNumK = (*this).NumK() + blk.NumK();

    int iMax = newNumI + 2 * numGhosts_;
    int jMax = newNumJ + 2 * numGhosts_;
    int kMax = newNumK + 2 * numGhosts_;

    procBlock newBlk(newNumI, newNumJ, newNumK, numGhosts_);

    newBlk.bc_ = (*this).BC();
    // newBlk.bc_.Join(blk.BC(), dir, alteredSurf);

    int ind = (*this).NumK();

    // loop over cell locations of of block
    for (int kk = 0; kk < kMax; kk++) {
      for (int jj = 0; jj < jMax; jj++) {
        for (int ii = 0; ii < iMax; ii++) {
          int loc = GetLoc1D(ii, jj, kk, iMax, jMax);
          int locNG =
              GetLoc1D(ii - newBlk.NumGhosts(), jj - newBlk.NumGhosts(),
                       kk - newBlk.NumGhosts(), newBlk.NumI(), newBlk.NumJ());

          int fLowI = GetLowerFaceI(ii, jj, kk, iMax, jMax);
          int fLowJ = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
          int fLowK = GetLowerFaceK(ii, jj, kk, iMax, jMax);

          //-------------------------------------------------------------------
          if (kk >= ind + newBlk.NumGhosts()) {  // this portion of parent block
                                                 // overlaps with upper split

            int locU = GetLoc1D(ii, jj, kk - ind, iMax, jMax);
            int locNGU =
                GetLoc1D(ii - blk.NumGhosts(), jj - blk.NumGhosts(),
                         kk - ind - blk.NumGhosts(), blk.NumI(), blk.NumJ());

            int fLowIU = GetLowerFaceI(ii, jj, kk - ind, iMax, jMax);
            int fLowJU = GetLowerFaceJ(ii, jj, kk - ind, iMax, jMax);
            int fLowKU = GetLowerFaceK(ii, jj, kk - ind, iMax, jMax);

            // assign cell variables
            newBlk.state_[loc] = blk.state_[locU];
            newBlk.vol_[loc] = blk.vol_[locU];
            newBlk.center_[loc] = blk.center_[locU];

            if (kk >= (ind + blk.NumGhosts()) &&
                kk < (kMax - blk.NumGhosts()) && ii >= blk.NumGhosts() &&
                ii < (iMax - blk.NumGhosts()) && jj >= blk.NumGhosts() &&
                jj < (jMax - blk.NumGhosts())) {  // physical cells

              newBlk.avgWaveSpeed_[locNG] = blk.avgWaveSpeed_[locNGU];
              newBlk.dt_[locNG] = blk.dt_[locNGU];
              newBlk.wallDist_[locNG] = blk.wallDist_[locNGU];
              newBlk.residual_[locNG] = blk.residual_[locNGU];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = blk.fAreaI_[fLowIU];
            newBlk.fAreaJ_[fLowJ] = blk.fAreaJ_[fLowJU];
            newBlk.fAreaK_[fLowK] = blk.fAreaK_[fLowKU];

            newBlk.fCenterI_[fLowI] = blk.fCenterI_[fLowIU];
            newBlk.fCenterJ_[fLowJ] = blk.fCenterJ_[fLowJU];
            newBlk.fCenterK_[fLowK] = blk.fCenterK_[fLowKU];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIU = GetUpperFaceI(ii, jj, kk - ind, iMax, jMax);

              newBlk.fAreaI_[fUpI] = blk.fAreaI_[fUpIU];
              newBlk.fCenterI_[fUpI] = blk.fCenterI_[fUpIU];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJU = GetUpperFaceJ(ii, jj, kk - ind, iMax, jMax);

              newBlk.fAreaJ_[fUpJ] = blk.fAreaJ_[fUpJU];
              newBlk.fCenterJ_[fUpJ] = blk.fCenterJ_[fUpJU];
            }

            if (kk == kMax - 1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKU = GetUpperFaceK(ii, jj, kk - ind, iMax, jMax);

              newBlk.fAreaK_[fUpK] = blk.fAreaK_[fUpKU];
              newBlk.fCenterK_[fUpK] = blk.fCenterK_[fUpKU];
            }

          //-----------------------------------------------------------------
          } else {  // this portion of parent block overlaps with lower split
            int locL = GetLoc1D(ii, jj, kk, iMax, jMax);
            int locNGL = GetLoc1D(
                ii - numGhosts_, jj - numGhosts_,
                kk - numGhosts_, (*this).NumI(), (*this).NumJ());

            int fLowIL = GetLowerFaceI(ii, jj, kk, iMax, jMax);
            int fLowJL = GetLowerFaceJ(ii, jj, kk, iMax, jMax);
            int fLowKL = GetLowerFaceK(ii, jj, kk, iMax, jMax);

            // assign cell variables
            newBlk.state_[loc] = (*this).state_[locL];
            newBlk.vol_[loc] = (*this).vol_[locL];
            newBlk.center_[loc] = (*this).center_[locL];

            if (kk >= numGhosts_ && kk < (ind + numGhosts_) &&
                ii >= numGhosts_ &&
                ii < (iMax - numGhosts_) &&
                jj >= numGhosts_ &&
                jj < (jMax - numGhosts_)) {  // physical cell

              newBlk.avgWaveSpeed_[locNG] = avgWaveSpeed_[locNGL];
              newBlk.dt_[locNG] = (*this).dt_[locNGL];
              newBlk.wallDist_[locNG] = (*this).wallDist_[locNGL];
              newBlk.residual_[locNG] = (*this).residual_[locNGL];
            }

            // assign face variables
            newBlk.fAreaI_[fLowI] = (*this).fAreaI_[fLowIL];
            newBlk.fAreaJ_[fLowJ] = (*this).fAreaJ_[fLowJL];
            newBlk.fAreaK_[fLowK] = (*this).fAreaK_[fLowKL];

            newBlk.fCenterI_[fLowI] = (*this).fCenterI_[fLowIL];
            newBlk.fCenterJ_[fLowJ] = (*this).fCenterJ_[fLowJL];
            newBlk.fCenterK_[fLowK] = (*this).fCenterK_[fLowKL];

            if (ii == iMax - 1) {  // at end of i-line assign upper face values
              int fUpI = GetUpperFaceI(ii, jj, kk, iMax, jMax);
              int fUpIL = GetUpperFaceI(ii, jj, kk, iMax, jMax);

              newBlk.fAreaI_[fUpI] = (*this).fAreaI_[fUpIL];
              newBlk.fCenterI_[fUpI] = (*this).fCenterI_[fUpIL];
            }

            if (jj == jMax - 1) {  // at end of j-line assign upper face values
              int fUpJ = GetUpperFaceJ(ii, jj, kk, iMax, jMax);
              int fUpJL = GetUpperFaceJ(ii, jj, kk, iMax, jMax);

              newBlk.fAreaJ_[fUpJ] = (*this).fAreaJ_[fUpJL];
              newBlk.fCenterJ_[fUpJ] = (*this).fCenterJ_[fUpJL];
            }

            if (kk == ind + numGhosts_ -
                          1) {  // at end of k-line assign upper face values
              int fUpK = GetUpperFaceK(ii, jj, kk, iMax, jMax);
              int fUpKL = GetUpperFaceK(ii, jj, kk, iMax, jMax);

              newBlk.fAreaK_[fUpK] = (*this).fAreaK_[fUpKL];
              newBlk.fCenterK_[fUpK] = (*this).fCenterK_[fUpKL];
            }
          }
        }
      }
    }
    (*this) = newBlk;
  } else {
    cerr << "ERROR: Error in procBlock::Join(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }
}

void procBlock::CalcGradsI(const int &ii, const int &jj, const int &kk,
                           const idealGas &eqnState, const bool &turbFlag,
                           tensor<double> &velGrad, vector3d<double> &tGrad,
                           vector3d<double> &tkeGrad,
                           vector3d<double> &omegaGrad) const {
  int imax = (*this).NumI() + 2 * numGhosts_ + 1;
  int jmax = (*this).NumJ() + 2 * numGhosts_;

  // location of current face (with ghost cells included)
  int loc = GetLoc1D(ii, jj, kk, imax, jmax);

  // location of faces in the upper and lower i-direction (with ghost cells
  // included)
  int fUpi = GetNeighborUpI(ii, jj, kk, imax, jmax);
  int fLowi = GetNeighborLowI(ii, jj, kk, imax, jmax);

  // location of j-faces in the upper and lower direction belonging to the cells
  // in the upper and lower i-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpjUpi = GetUpperFaceJ(ii, jj, kk, imax - 1, jmax);
  int fUpjLowi = GetUpperFaceJ(ii - 1, jj, kk, imax - 1, jmax);
  int fLowjUpi = GetLowerFaceJ(ii, jj, kk, imax - 1, jmax);
  int fLowjLowi = GetLowerFaceJ(ii - 1, jj, kk, imax - 1, jmax);

  // location of k-faces in the upper and lower direction belonging to the cells
  // in the upper and lower i-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpkUpi = GetUpperFaceK(ii, jj, kk, imax - 1, jmax);
  int fUpkLowi = GetUpperFaceK(ii - 1, jj, kk, imax - 1, jmax);
  int fLowkUpi = GetLowerFaceK(ii, jj, kk, imax - 1, jmax);
  int fLowkLowi = GetLowerFaceK(ii - 1, jj, kk, imax - 1, jmax);

  // location of cells in the upper and lower i-direction with respect to
  // baseline face (with ghost cells included)
  int iLow = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
  int iUp = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

  // calculate areas of faces in alternate control volume
  vector3d<double> aiu = 0.5 * ((*this).FAreaI(loc).Vector() +
                                (*this).FAreaI(fUpi).Vector());
  vector3d<double> ail = 0.5 * ((*this).FAreaI(loc).Vector() +
                                (*this).FAreaI(fLowi).Vector());

  vector3d<double> aju =
      0.5 * ((*this).FAreaJ(fUpjUpi).Vector() +
             (*this).FAreaJ(fUpjLowi).Vector());
  vector3d<double> ajl =
      0.5 * ((*this).FAreaJ(fLowjUpi).Vector() +
             (*this).FAreaJ(fLowjLowi).Vector());

  vector3d<double> aku =
      0.5 * ((*this).FAreaK(fUpkUpi).Vector() +
             (*this).FAreaK(fUpkLowi).Vector());
  vector3d<double> akl =
      0.5 * ((*this).FAreaK(fLowkUpi).Vector() +
             (*this).FAreaK(fLowkLowi).Vector());

  // calculate volume of alternate control volume
  double vol_ = 0.5 * (vol_(iLow) + vol_(iUp));

  // location of cells in the upper and lower j-direction and the upper
  // and lower i-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int jUpiUp = GetNeighborUpJ(ii, jj, kk, imax - 1, jmax);
  int jUpiLow = GetNeighborUpJ(ii - 1, jj, kk, imax - 1, jmax);
  int jLowiUp = GetNeighborLowJ(ii, jj, kk, imax - 1, jmax);
  int jLowiLow = GetNeighborLowJ(ii - 1, jj, kk, imax - 1, jmax);

  // location of cells in the upper and lower k-direction and the upper
  // and lower i-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int kUpiUp = GetNeighborUpK(ii, jj, kk, imax - 1, jmax);
  int kUpiLow = GetNeighborUpK(ii - 1, jj, kk, imax - 1, jmax);
  int kLowiUp = GetNeighborLowK(ii, jj, kk, imax - 1, jmax);
  int kLowiLow = GetNeighborLowK(ii - 1, jj, kk, imax - 1, jmax);

  // calculate average velocity on j and k faces of alternate control volume
  vector3d<double> vju =
      0.25 *
      (state_(iLow).Velocity() + state_(iUp).Velocity() +
       state_(jUpiUp).Velocity() + state_(jUpiLow).Velocity());
  vector3d<double> vjl =
      0.25 *
      (state_(iLow).Velocity() + state_(iUp).Velocity() +
       state_(jLowiUp).Velocity() + state_(jLowiLow).Velocity());

  vector3d<double> vku =
      0.25 *
      (state_(iLow).Velocity() + state_(iUp).Velocity() +
       state_(kUpiUp).Velocity() + state_(kUpiLow).Velocity());
  vector3d<double> vkl =
      0.25 *
      (state_(iLow).Velocity() + state_(iUp).Velocity() +
       state_(kLowiUp).Velocity() + state_(kLowiLow).Velocity());

  // Get velocity gradient at face
  velGrad = CalcVelGradGG(state_(iLow).Velocity(),
                          state_(iUp).Velocity(), vjl, vju, vkl, vku,
                          ail, aiu, ajl, aju, akl, aku, vol_);

  // calculate average temperature on j and k faces of alternate control volume
  double tju = 0.25 * (state_(iLow).Temperature(eqnState) +
                       state_(iUp).Temperature(eqnState) +
                       state_(jUpiUp).Temperature(eqnState) +
                       state_(jUpiLow).Temperature(eqnState));
  double tjl = 0.25 * (state_(iLow).Temperature(eqnState) +
                       state_(iUp).Temperature(eqnState) +
                       state_(jLowiUp).Temperature(eqnState) +
                       state_(jLowiLow).Temperature(eqnState));

  double tku = 0.25 * (state_(iLow).Temperature(eqnState) +
                       state_(iUp).Temperature(eqnState) +
                       state_(kUpiUp).Temperature(eqnState) +
                       state_(kUpiLow).Temperature(eqnState));
  double tkl = 0.25 * (state_(iLow).Temperature(eqnState) +
                       state_(iUp).Temperature(eqnState) +
                       state_(kLowiUp).Temperature(eqnState) +
                       state_(kLowiLow).Temperature(eqnState));

  // Get temperature gradient at face
  tGrad = CalcScalarGradGG(state_(iLow).Temperature(eqnState),
                           state_(iUp).Temperature(eqnState), tjl, tju,
                           tkl, tku, ail, aiu, ajl, aju, akl, aku, vol_);

  if (turbFlag) {
    // calculate average tke on j and k faces of alternate control volume
    double tkeju =
        0.25 * (state_(iLow).Tke() + state_(iUp).Tke() +
                state_(jUpiUp).Tke() + state_(jUpiLow).Tke());
    double tkejl =
        0.25 * (state_(iLow).Tke() + state_(iUp).Tke() +
                state_(jLowiUp).Tke() + state_(jLowiLow).Tke());

    double tkeku =
        0.25 * (state_(iLow).Tke() + state_(iUp).Tke() +
                state_(kUpiUp).Tke() + state_(kUpiLow).Tke());
    double tkekl =
        0.25 * (state_(iLow).Tke() + state_(iUp).Tke() +
                state_(kLowiUp).Tke() + state_(kLowiLow).Tke());

    // Get tke gradient at face
    tkeGrad = CalcScalarGradGG(state_(iLow).Tke(),
                               state_(iUp).Tke(), tkejl, tkeju, tkekl,
                               tkeku, ail, aiu, ajl, aju, akl, aku, vol_);

    // calculate average Omega on j and k faces of alternate control volume
    double omgju =
        0.25 * (state_(iLow).Omega() + state_(iUp).Omega() +
                state_(jUpiUp).Omega() + state_(jUpiLow).Omega());
    double omgjl =
        0.25 *
        (state_(iLow).Omega() + state_(iUp).Omega() +
         state_(jLowiUp).Omega() + state_(jLowiLow).Omega());

    double omgku =
        0.25 * (state_(iLow).Omega() + state_(iUp).Omega() +
                state_(kUpiUp).Omega() + state_(kUpiLow).Omega());
    double omgkl =
        0.25 *
        (state_(iLow).Omega() + state_(iUp).Omega() +
         state_(kLowiUp).Omega() + state_(kLowiLow).Omega());

    // Get omega gradient at face
    omegaGrad = CalcScalarGradGG(
        state_(iLow).Omega(), state_(iUp).Omega(), omgjl, omgju,
        omgkl, omgku, ail, aiu, ajl, aju, akl, aku, vol_);
  }
}

void procBlock::CalcGradsJ(const int &ii, const int &jj, const int &kk,
                           const idealGas &eqnState, const bool &turbFlag,
                           tensor<double> &velGrad, vector3d<double> &tGrad,
                           vector3d<double> &tkeGrad,
                           vector3d<double> &omegaGrad) const {
  int imax = (*this).NumI() + 2 * numGhosts_;
  int jmax = (*this).NumJ() + 2 * numGhosts_ + 1;

  // location of current face (with ghost cells included)
  int loc = GetLoc1D(ii, jj, kk, imax, jmax);

  // location of faces in the upper and lower j-direction (with ghost cells
  // included)
  int fUpj = GetNeighborUpJ(ii, jj, kk, imax, jmax);
  int fLowj = GetNeighborLowJ(ii, jj, kk, imax, jmax);

  // location of i-faces in the upper and lower direction belonging to the cells
  // in the upper and lower j-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpiUpj = GetUpperFaceI(ii, jj, kk, imax, jmax - 1);
  int fUpiLowj = GetUpperFaceI(ii, jj - 1, kk, imax, jmax - 1);
  int fLowiUpj = GetLowerFaceI(ii, jj, kk, imax, jmax - 1);
  int fLowiLowj = GetLowerFaceI(ii, jj - 1, kk, imax, jmax - 1);

  // location of k-faces in the upper and lower direction belonging to the cells
  // in the upper and lower j-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpkUpj = GetUpperFaceK(ii, jj, kk, imax, jmax - 1);
  int fUpkLowj = GetUpperFaceK(ii, jj - 1, kk, imax, jmax - 1);
  int fLowkUpj = GetLowerFaceK(ii, jj, kk, imax, jmax - 1);
  int fLowkLowj = GetLowerFaceK(ii, jj - 1, kk, imax, jmax - 1);

  // location of cells in the upper and lower j-direction with respect to
  // baseline face (with ghost cells included)
  int jLow = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
  int jUp = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

  // calculate areas of faces in alternate control volume
  vector3d<double> aju = 0.5 * ((*this).FAreaJ(loc).Vector() +
                                (*this).FAreaJ(fUpj).Vector());
  vector3d<double> ajl = 0.5 * ((*this).FAreaJ(loc).Vector() +
                                    (*this).FAreaJ(fLowj).Vector());

  vector3d<double> aiu =
      0.5 * ((*this).FAreaI(fUpiUpj).Vector() +
             (*this).FAreaI(fUpiLowj).Vector());
  vector3d<double> ail =
      0.5 * ((*this).FAreaI(fLowiUpj).Vector() +
             (*this).FAreaI(fLowiLowj).Vector());

  vector3d<double> aku =
      0.5 * ((*this).FAreaK(fUpkUpj).Vector() +
             (*this).FAreaK(fUpkLowj).Vector());
  vector3d<double> akl =
      0.5 * ((*this).FAreaK(fLowkUpj).Vector() +
             (*this).FAreaK(fLowkLowj).Vector());

  // calculate volume of alternate control volume
  double vol_ = 0.5 * (vol_(jLow) + vol_(jUp));

  // location of cells in the upper and lower i-direction and the upper
  // and lower j-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int iUpjUp = GetNeighborUpI(ii, jj, kk, imax, jmax - 1);
  int iUpjLow = GetNeighborUpI(ii, jj - 1, kk, imax, jmax - 1);
  int iLowjUp = GetNeighborLowI(ii, jj, kk, imax, jmax - 1);
  int iLowjLow = GetNeighborLowI(ii, jj - 1, kk, imax, jmax - 1);

  // location of cells in the upper and lower k-direction and the upper
  // and lower j-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int kUpjUp = GetNeighborUpK(ii, jj, kk, imax, jmax - 1);
  int kUpjLow = GetNeighborUpK(ii, jj - 1, kk, imax, jmax - 1);
  int kLowjUp = GetNeighborLowK(ii, jj, kk, imax, jmax - 1);
  int kLowjLow = GetNeighborLowK(ii, jj - 1, kk, imax, jmax - 1);

  // calculate average velocity on i and k faces of alternate control volume
  vector3d<double> viu =
      0.25 *
      (state_(jLow).Velocity() + state_(jUp).Velocity() +
       state_(iUpjUp).Velocity() + state_(iUpjLow).Velocity());
  vector3d<double> vil =
      0.25 *
      (state_(jLow).Velocity() + state_(jUp).Velocity() +
       state_(iLowjUp).Velocity() + state_(iLowjLow).Velocity());

  vector3d<double> vku =
      0.25 *
      (state_(jLow).Velocity() + state_(jUp).Velocity() +
       state_(kUpjUp).Velocity() + state_(kUpjLow).Velocity());
  vector3d<double> vkl =
      0.25 *
      (state_(jLow).Velocity() + state_(jUp).Velocity() +
       state_(kLowjUp).Velocity() + state_(kLowjLow).Velocity());

  // Get velocity gradient at face
  velGrad = CalcVelGradGG(vil, viu, state_(jLow).Velocity(),
                          state_(jUp).Velocity(), vkl, vku, ail, aiu,
                          ajl, aju, akl, aku, vol_);

  // calculate average temperature on i and k faces of alternate control volume
  double tiu = 0.25 * (state_(jLow).Temperature(eqnState) +
                       state_(jUp).Temperature(eqnState) +
                       state_(iUpjUp).Temperature(eqnState) +
                       state_(iUpjLow).Temperature(eqnState));
  double til = 0.25 * (state_(jLow).Temperature(eqnState) +
                       state_(jUp).Temperature(eqnState) +
                       state_(iLowjUp).Temperature(eqnState) +
                       state_(iLowjLow).Temperature(eqnState));

  double tku = 0.25 * (state_(jLow).Temperature(eqnState) +
                       state_(jUp).Temperature(eqnState) +
                       state_(kUpjUp).Temperature(eqnState) +
                       state_(kUpjLow).Temperature(eqnState));
  double tkl = 0.25 * (state_(jLow).Temperature(eqnState) +
                       state_(jUp).Temperature(eqnState) +
                       state_(kLowjUp).Temperature(eqnState) +
                       state_(kLowjLow).Temperature(eqnState));

  // Get temperature gradient at face
  tGrad = CalcScalarGradGG(til, tiu, state_(jLow).Temperature(eqnState),
                           state_(jUp).Temperature(eqnState), tkl, tku,
                           ail, aiu, ajl, aju, akl, aku, vol_);

  if (turbFlag) {
    // calculate average tke on i and k faces of alternate control volume
    double tkeiu =
        0.25 * (state_(jLow).Tke() + state_(jUp).Tke() +
                state_(iUpjUp).Tke() + state_(iUpjLow).Tke());
    double tkeil =
        0.25 * (state_(jLow).Tke() + state_(jUp).Tke() +
                state_(iLowjUp).Tke() + state_(iLowjLow).Tke());

    double tkeku =
        0.25 * (state_(jLow).Tke() + state_(jUp).Tke() +
                state_(kUpjUp).Tke() + state_(kUpjLow).Tke());
    double tkekl =
        0.25 * (state_(jLow).Tke() + state_(jUp).Tke() +
                state_(kLowjUp).Tke() + state_(kLowjLow).Tke());

    // Get temperature gradient at face
    tkeGrad = CalcScalarGradGG(tkeil, tkeiu, state_(jLow).Tke(),
                               state_(jUp).Tke(), tkekl, tkeku, ail, aiu,
                               ajl, aju, akl, aku, vol_);

    // calculate average omega on i and k faces of alternate control volume
    double omgiu =
        0.25 * (state_(jLow).Omega() + state_(jUp).Omega() +
                state_(iUpjUp).Omega() + state_(iUpjLow).Omega());
    double omgil =
        0.25 *
        (state_(jLow).Omega() + state_(jUp).Omega() +
         state_(iLowjUp).Omega() + state_(iLowjLow).Omega());

    double omgku =
        0.25 * (state_(jLow).Omega() + state_(jUp).Omega() +
                state_(kUpjUp).Omega() + state_(kUpjLow).Omega());
    double omgkl =
        0.25 *
        (state_(jLow).Omega() + state_(jUp).Omega() +
         state_(kLowjUp).Omega() + state_(kLowjLow).Omega());

    // Get temperature gradient at face
    omegaGrad = CalcScalarGradGG(omgil, omgiu, state_(jLow).Omega(),
                                 state_(jUp).Omega(), omgkl, omgku, ail,
                                 aiu, ajl, aju, akl, aku, vol_);
  }
}

void procBlock::CalcGradsK(const int &ii, const int &jj, const int &kk,
                           const idealGas &eqnState, const bool &turbFlag,
                           tensor<double> &velGrad, vector3d<double> &tGrad,
                           vector3d<double> &tkeGrad,
                           vector3d<double> &omegaGrad) const {
  int imax = (*this).NumI() + 2 * numGhosts_;
  int jmax = (*this).NumJ() + 2 * numGhosts_;

  // location of current face (with ghost cells included)
  int loc = GetLoc1D(ii, jj, kk, imax, jmax);

  // location of faces in the upper and lower k-direction (with ghost cells
  // included)
  int fUpk = GetNeighborUpK(ii, jj, kk, imax, jmax);
  int fLowk = GetNeighborLowK(ii, jj, kk, imax, jmax);

  // location of i-faces in the upper and lower direction belonging to the cells
  // in the upper and lower k-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpiUpk = GetUpperFaceI(ii, jj, kk, imax, jmax);
  int fUpiLowk = GetUpperFaceI(ii, jj, kk - 1, imax, jmax);
  int fLowiUpk = GetLowerFaceI(ii, jj, kk, imax, jmax);
  int fLowiLowk = GetLowerFaceI(ii, jj, kk - 1, imax, jmax);

  // location of j-faces in the upper and lower direction belonging to the cells
  // in the upper and lower k-direction of the current face (with ghost cells
  // included) - these are used in the gradient calculation to construct the
  // alternate control volume
  int fUpjUpk = GetUpperFaceJ(ii, jj, kk, imax, jmax);
  int fUpjLowk = GetUpperFaceJ(ii, jj, kk - 1, imax, jmax);
  int fLowjUpk = GetLowerFaceJ(ii, jj, kk, imax, jmax);
  int fLowjLowk = GetLowerFaceJ(ii, jj, kk - 1, imax, jmax);

  // location of cells in the upper and lower k-direction with respect to
  // baseline face (with ghost cells included)
  int kLow = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
  int kUp = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

  // calculate areas of faces in alternate control volume
  vector3d<double> aku = 0.5 * ((*this).FAreaK(loc).Vector() +
                                (*this).FAreaK(fUpk).Vector());
  vector3d<double> akl = 0.5 * ((*this).FAreaK(loc).Vector() +
                                (*this).FAreaK(fLowk).Vector());

  vector3d<double> aiu =
      0.5 * ((*this).FAreaI(fUpiUpk).Vector() +
             (*this).FAreaI(fUpiLowk).Vector());
  vector3d<double> ail =
      0.5 * ((*this).FAreaI(fLowiUpk).Vector() +
             (*this).FAreaI(fLowiLowk).Vector());

  vector3d<double> aju =
      0.5 * ((*this).FAreaJ(fUpjUpk).Vector() +
             (*this).FAreaJ(fUpjLowk).Vector());
  vector3d<double> ajl =
      0.5 * ((*this).FAreaJ(fLowjUpk).Vector()
             + (*this).FAreaJ(fLowjLowk).Vector());

  // calculate volume of alternate control volume
  double vol_ = 0.5 * (vol_(kLow) + vol_(kUp));

  // location of cells in the upper and lower k-direction and the upper
  // and lower i-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int iUpkUp = GetNeighborUpI(ii, jj, kk, imax, jmax);
  int iUpkLow = GetNeighborUpI(ii, jj, kk - 1, imax, jmax);
  int iLowkUp = GetNeighborLowI(ii, jj, kk, imax, jmax);
  int iLowkLow = GetNeighborLowI(ii, jj, kk - 1, imax, jmax);

  // location of cells in the upper and lower k-direction and the upper
  // and lower j-direction of the current face (with ghost cells included) -
  // these are used in the gradient calculation to construct the alternate
  // control volume
  int jUpkUp = GetNeighborUpK(ii, jj, kk, imax, jmax);
  int jUpkLow = GetNeighborUpK(ii, jj, kk - 1, imax, jmax);
  int jLowkUp = GetNeighborLowK(ii, jj, kk, imax, jmax);
  int jLowkLow = GetNeighborLowK(ii, jj, kk - 1, imax, jmax);

  // calculate average velocity on i and j faces of alternate control volume
  vector3d<double> viu =
      0.25 *
      (state_(kLow).Velocity() + state_(kUp).Velocity() +
       state_(iUpkUp).Velocity() + state_(iUpkLow).Velocity());
  vector3d<double> vil =
      0.25 *
      (state_(kLow).Velocity() + state_(kUp).Velocity() +
       state_(iLowkUp).Velocity() + state_(iLowkLow).Velocity());

  vector3d<double> vju =
      0.25 *
      (state_(kLow).Velocity() + state_(kUp).Velocity() +
       state_(jUpkUp).Velocity() + state_(jUpkLow).Velocity());
  vector3d<double> vjl =
      0.25 *
      (state_(kLow).Velocity() + state_(kUp).Velocity() +
       state_(jLowkUp).Velocity() + state_(jLowkLow).Velocity());

  // Get velocity gradient at face
  velGrad = CalcVelGradGG(vil, viu, vjl, vju, state_(kLow).Velocity(),
                          state_(kUp).Velocity(), ail, aiu, ajl, aju,
                          akl, aku, vol_);

  // calculate average temperature on i and j faces of alternate control volume
  double tiu = 0.25 * (state_(kLow).Temperature(eqnState) +
                       state_(kUp).Temperature(eqnState) +
                       state_(iUpkUp).Temperature(eqnState) +
                       state_(iUpkLow).Temperature(eqnState));
  double til = 0.25 * (state_(kLow).Temperature(eqnState) +
                       state_(kUp).Temperature(eqnState) +
                       state_(iLowkUp).Temperature(eqnState) +
                       state_(iLowkLow).Temperature(eqnState));

  double tju = 0.25 * (state_(kLow).Temperature(eqnState) +
                       state_(kUp).Temperature(eqnState) +
                       state_(jUpkUp).Temperature(eqnState) +
                       state_(jUpkLow).Temperature(eqnState));
  double tjl = 0.25 * (state_(kLow).Temperature(eqnState) +
                       state_(kUp).Temperature(eqnState) +
                       state_(jLowkUp).Temperature(eqnState) +
                       state_(jLowkLow).Temperature(eqnState));

  // Get temperature gradient at face
  tGrad = CalcScalarGradGG(til, tiu, tjl, tju,
                           state_(kLow).Temperature(eqnState),
                           state_(kUp).Temperature(eqnState), ail, aiu,
                           ajl, aju, akl, aku, vol_);

  if (turbFlag) {
    // calculate average tke on i and j faces of alternate control volume
    double tkeiu =
        0.25 * (state_(kLow).Tke() + state_(kUp).Tke() +
                state_(iUpkUp).Tke() + state_(iUpkLow).Tke());
    double tkeil =
        0.25 * (state_(kLow).Tke() + state_(kUp).Tke() +
                state_(iLowkUp).Tke() + state_(iLowkLow).Tke());

    double tkeju =
        0.25 * (state_(kLow).Tke() + state_(kUp).Tke() +
                state_(jUpkUp).Tke() + state_(jUpkLow).Tke());
    double tkejl =
        0.25 * (state_(kLow).Tke() + state_(kUp).Tke() +
                state_(jLowkUp).Tke() + state_(jLowkLow).Tke());

    // Get temperature gradient at face
    tkeGrad = CalcScalarGradGG(
        tkeil, tkeiu, tkejl, tkeju, state_(kLow).Tke(),
        state_(kUp).Tke(), ail, aiu, ajl, aju, akl, aku, vol_);

    // calculate average omega on i and j faces of alternate control volume
    double omgiu =
        0.25 * (state_(kLow).Omega() + state_(kUp).Omega() +
                state_(iUpkUp).Omega() + state_(iUpkLow).Omega());
    double omgil =
        0.25 *
        (state_(kLow).Omega() + state_(kUp).Omega() +
         state_(iLowkUp).Omega() + state_(iLowkLow).Omega());

    double omgju =
        0.25 * (state_(kLow).Omega() + state_(kUp).Omega() +
                state_(jUpkUp).Omega() + state_(jUpkLow).Omega());
    double omgjl =
        0.25 *
        (state_(kLow).Omega() + state_(kUp).Omega() +
         state_(jLowkUp).Omega() + state_(jLowkLow).Omega());

    // Get temperature gradient at face
    omegaGrad = CalcScalarGradGG(
        omgil, omgiu, omgjl, omgju, state_(kLow).Omega(),
        state_(kUp).Omega(), ail, aiu, ajl, aju, akl, aku, vol_);
  }
}

// Member function to calculate the source terms and add them to the residual
void procBlock::CalcSrcTerms(const gradients &grads, const sutherland &suth,
                             const turbModel *turb) {
  // grads -- gradients (vel, temp, tke, omega)
  // suth -- sutherland's law for viscosity
  // turb -- turbulence model

  // loop over all physical cells - no ghost cells needed for source terms
  for (struct {int p; int g;} kk = {0, numGhosts_}; kk.p <
           (*this).NumK(); kk.g++, kk.p++) {
    for (struct {int p; int g;} jj = {0, numGhosts_}; jj.p <
             (*this).NumJ(); jj.g++, jj.p++) {
      for (struct {int p; int g;} ii = {0, numGhosts_}; ii.p <
               (*this).NumI(); ii.g++, ii.p++) {
        // calculate turbulent source terms
        source src;
        src.CalcTurbSrc(turb, state_(loc), grads, suth,
                        ii.p, jj.p, kk.p);

        // add source terms to residual
        // multiply by -1 because residual is initially on opposite
        // side of equation
        (*this).AddToResidual(src * (-1.0 * vol_(ii.g, jj.g, kk.g)),
                              ii.p, jj.p, kk.p);
      }
    }
  }
}

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

#include <iostream>
#include "slices.hpp"
#include "procBlock.hpp"
#include "range.hpp"  // range

using std::cout;
using std::endl;
using std::cerr;

// constructors for geomSlice class
// constructor -- initialize state vector with dummy variables
geomSlice::geomSlice(const int &li, const int &lj, const int &lk,
                     const int &pblk) {
  // li -- size of direction i (cell)
  // lj -- size of direction j (cell)
  // lk -- size of direction k (cell)
  // pblk -- parent block that slice is coming from

  parBlock_ = pblk;

  center_ = {li, lj, lk, 0};
  fAreaI_ = {li + 1, lj, lk, 0};
  fAreaJ_ = {li, lj + 1, lk, 0};
  fAreaK_ = {li, lj, lk + 1, 0};
  fCenterI_ = {li + 1, lj, lk, 0};
  fCenterJ_ = {li, lj + 1, lk, 0};
  fCenterK_ = {li, lj, lk + 1, 0};

  vol_ = {li, lj, lk, 0};
}

/* constructor to get a slice (portion) of the geometry of a procBlock. The
geom slice remains in the orientation of its parent block unless
the revI, revJ, or revK flags are activated to reverse either the i, j, or k
directions respectively. This function essentially cuts out a portion
of a procBlock and returns a geomSlice with that portion of the block geometry.
*/
geomSlice::geomSlice(const procBlock &blk, const range &ir, const range &jr,
                     const range &kr, const bool revI, const bool revJ,
                     const bool revK) {
  // blk -- procBlock to extract slice from
  // ir -- range for i-direction
  // jr -- range for j-direction
  // kr -- range for k-direction
  // revI -- flag to reverse i direction of indices (default false)
  // revJ -- flag to reverse j direction of indices (default false)
  // revK -- flag to reverse k direction of indices (default false)

  // allocate dimensions
  const auto numI = ir.Size();
  const auto numJ = jr.Size();
  const auto numK = kr.Size();
  parBlock_ = blk.ParentBlock();

  // allocate size for vectors
  center_ = {numI, numJ, numK, 0};
  fAreaI_ = {numI + 1, numJ, numK, 0};
  fAreaJ_ = {numI, numJ + 1, numK, 0};
  fAreaK_ = {numI, numJ, numK + 1, 0};
  fCenterI_ = {numI + 1, numJ, numK, 0};
  fCenterJ_ = {numI, numJ + 1, numK, 0};
  fCenterK_ = {numI, numJ, numK + 1, 0};
  vol_ = {numI, numJ, numK, 0};

  // loop over all cells in slice and populate
  for (auto kk = vol_.StartK(), ks = kr.Start(); kk < vol_.EndK(); kk++) {
    for (auto jj = vol_.StartJ(), js = jr.Start(); jj < vol_.EndJ(); jj++) {
      for (auto ii = vol_.StartI(), is = ir.Start(); ii < vol_.EndI(); ii++) {
        // determine if direction needs to be reversed
        auto k = revK ? numK - 1 - kk : kk;
        auto kFac = revK ? -1.0 : 1.0;

        auto j = revJ ? numJ - 1 - jj : jj;
        auto jFac = revJ ? -1.0 : 1.0;

        auto i = revI ? numI - 1 - ii : ii;
        auto iFac = revI ? -1.0 : 1.0;

        // assign cell variables
        vol_(ii, jj, kk) = blk.Vol(is + i, js + j, ks + k);
        center_(ii, jj, kk) = blk.Center(is + i, js + j, ks + k);

        // assign i-face variables
        // if direction is reversed, upper/lower faces need to be swapped
        if (revI) {
          fAreaI_(ii, jj, kk) = iFac * blk.FAreaI(is + i + 1, js + j, ks + k);
          fCenterI_(ii, jj, kk) = blk.FCenterI(is + i + 1, js + j, ks + k);

          // at end of i-line assign upper face values too
          if (ii == numI - 1) {
            fAreaI_(ii + 1, jj, kk) = iFac * blk.FAreaI(is + i, js + j, ks + k);
            fCenterI_(ii + 1, jj, kk) = blk.FCenterI(is + i, js + j, ks + k);
          }
        } else {
          fAreaI_(ii, jj, kk) = iFac * blk.FAreaI(is + i, js + j, ks + k);
          fCenterI_(ii, jj, kk) = blk.FCenterI(is + i, js + j, ks + k);

          // at end of i-line assign upper face values too
          if (ii == numI - 1) {
            fAreaI_(ii + 1, jj, kk) = iFac * blk.FAreaI(is + i + 1, js + j,
                                                        ks + k);
            fCenterI_(ii + 1, jj, kk) = blk.FCenterI(is + i + 1, js + j,
                                                     ks + k);
          }
        }

        // assign j-face variables
        // if direction is reversed, upper/lower faces need to be swapped
        if (revJ) {
          fAreaJ_(ii, jj, kk) = jFac * blk.FAreaJ(is + i, js + j + 1, ks + k);
          fCenterJ_(ii, jj, kk) = blk.FCenterJ(is + i, js + j + 1, ks + k);

          // at end of j-line assign upper face values too
          if (jj == numJ - 1) {
            fAreaJ_(ii, jj + 1, kk) = jFac * blk.FAreaJ(is + i, js + j, ks + k);
            fCenterJ_(ii, jj + 1, kk) = blk.FCenterJ(is + i, js + j, ks + k);
          }
        } else {
          fAreaJ_(ii, jj, kk) = jFac * blk.FAreaJ(is + i, js + j, ks + k);
          fCenterJ_(ii, jj, kk) = blk.FCenterJ(is + i, js + j, ks + k);

          // at end of j-line assign upper face values too
          if (jj == numJ - 1) {
            fAreaJ_(ii, jj + 1, kk) = jFac * blk.FAreaJ(is + i, js + j + 1,
                                                        ks + k);
            fCenterJ_(ii, jj + 1, kk) = blk.FCenterJ(is + i, js + j + 1,
                                                     ks + k);
          }
        }

        // assign k-face variables
        // if direction is reversed, upper/lower faces need to be swapped
        if (revK) {
          fAreaK_(ii, jj, kk) = kFac * blk.FAreaK(is + i, js + j, ks + k + 1);
          fCenterK_(ii, jj, kk) = blk.FCenterK(is + i, js + j, ks + k + 1);

          // at end of k-line assign upper face values too
          if (kk == numK - 1) {
            fAreaK_(ii, jj, kk + 1) = kFac * blk.FAreaK(is + i, js + j, ks + k);
            fCenterK_(ii, jj, kk + 1) = blk.FCenterK(is + i, js + j, ks + k);
          }
        } else {
          fAreaK_(ii, jj, kk) = kFac * blk.FAreaK(is + i, js + j, ks + k);
          fCenterK_(ii, jj, kk) = blk.FCenterK(is + i, js + j, ks + k);

          // at end of k-line assign upper face values too
          if (kk == numK - 1) {
            fAreaK_(ii, jj, kk + 1) = kFac * blk.FAreaK(is + i, js + j,
                                                        ks + k + 1);
            fCenterK_(ii, jj, kk + 1) = blk.FCenterK(is + i, js + j,
                                                     ks + k + 1);
          }
        }
      }
    }
  }
}

void geomSlice::SameSizeResize(const int &numI, const int &numJ,
                               const int &numK, const int &ng) {
  center_.SameSizeResize(numI, numJ, numK);
  fAreaI_.SameSizeResize(numI + 1, numJ, numK);
  fAreaJ_.SameSizeResize(numI, numJ + 1, numK);
  fAreaK_.SameSizeResize(numI, numJ, numK + 1);
  fCenterI_.SameSizeResize(numI + 1, numJ, numK);
  fCenterJ_.SameSizeResize(numI, numJ + 1, numK);
  fCenterK_.SameSizeResize(numI, numJ, numK + 1);
  vol_.SameSizeResize(numI, numJ, numK);
}

void geomSlice::PackSwapUnpackMPI(const connection &inter,
                                  const MPI_Datatype &MPI_vec3d,
                                  const MPI_Datatype &MPI_vec3dMag,
                                  const int &rank, const int &tag) {
  // swap with mpi_send_recv_replace
  // pack data into buffer, but first get size
  auto bufSize = 0;
  auto tempSize = 0;
  // add size for centers
  MPI_Pack_size(center_.Size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  // add size for face areas
  MPI_Pack_size(fAreaI_.Size(), MPI_vec3dMag, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  MPI_Pack_size(fAreaJ_.Size(), MPI_vec3dMag, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  MPI_Pack_size(fAreaK_.Size(), MPI_vec3dMag, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  // add size for face centers
  MPI_Pack_size(fCenterI_.Size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  MPI_Pack_size(fCenterJ_.Size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  MPI_Pack_size(fCenterK_.Size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  // add size for volume
  MPI_Pack_size(vol_.Size(), MPI_DOUBLE, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;
  // add size for 5 ints for multiArray3d dims, ghost layers, and parent block
  MPI_Pack_size(5, MPI_INT, MPI_COMM_WORLD, &tempSize);
  bufSize += tempSize;

  // allocate buffer to pack data into
  // use unique_ptr to manage memory; use underlying pointer for MPI calls
  auto buffer = std::make_unique<char[]>(bufSize);
  auto *rawBuffer = buffer.get();

  // pack data into buffer
  auto numI = this->NumI();
  auto numJ = this->NumJ();
  auto numK = this->NumK();
  auto numGhosts = this->GhostLayers();
  auto position = 0;
  // pack ints
  MPI_Pack(&numI, 1, MPI_INT, rawBuffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numJ, 1, MPI_INT, rawBuffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numK, 1, MPI_INT, rawBuffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&numGhosts, 1, MPI_INT, rawBuffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&parBlock_, 1, MPI_INT, rawBuffer, bufSize, &position, MPI_COMM_WORLD);
  // pack center
  MPI_Pack(&(*std::begin(center_)), center_.Size(), MPI_vec3d, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  // pack face areas
  MPI_Pack(&(*std::begin(fAreaI_)), fAreaI_.Size(), MPI_vec3dMag, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(fAreaJ_)), fAreaJ_.Size(), MPI_vec3dMag, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(fAreaK_)), fAreaK_.Size(), MPI_vec3dMag, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  // pack face centers
  MPI_Pack(&(*std::begin(fCenterI_)), fCenterI_.Size(), MPI_vec3d, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(fCenterJ_)), fCenterJ_.Size(), MPI_vec3d, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*std::begin(fCenterK_)), fCenterK_.Size(), MPI_vec3d, rawBuffer,
           bufSize, &position, MPI_COMM_WORLD);
  // pack volume
  MPI_Pack(&(*std::begin(vol_)), vol_.Size(), MPI_DOUBLE, rawBuffer, bufSize,
           &position, MPI_COMM_WORLD);

  MPI_Status status;
  if (rank == inter.RankFirst()) {  // send/recv with second entry in connection
    MPI_Sendrecv_replace(rawBuffer, bufSize, MPI_PACKED, inter.RankSecond(),
                         tag, inter.RankSecond(), tag, MPI_COMM_WORLD, &status);
  } else {  // send/recv with first entry in connection
    MPI_Sendrecv_replace(rawBuffer, bufSize, MPI_PACKED, inter.RankFirst(), tag,
                         inter.RankFirst(), tag, MPI_COMM_WORLD, &status);
  }

  // put slice back into geomSlice
  position = 0;
  MPI_Unpack(rawBuffer, bufSize, &position, &numI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &numJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &numK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &numGhosts, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &parBlock_, 1, MPI_INT,
             MPI_COMM_WORLD);

  // resize slice
  this->SameSizeResize(numI, numJ, numK, numGhosts);
  // unpack center
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(center_)),
             center_.Size(), MPI_vec3d, MPI_COMM_WORLD);
  // unpack face areas
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fAreaI_)),
             fAreaI_.Size(), MPI_vec3dMag, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fAreaJ_)),
             fAreaJ_.Size(), MPI_vec3dMag, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fAreaK_)),
             fAreaK_.Size(), MPI_vec3dMag, MPI_COMM_WORLD);
  // unpack face centers
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fCenterI_)),
             fCenterI_.Size(), MPI_vec3d, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fCenterJ_)),
             fCenterJ_.Size(), MPI_vec3d, MPI_COMM_WORLD);
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(fCenterK_)),
             fCenterK_.Size(), MPI_vec3d, MPI_COMM_WORLD);
  // unpack volume
  MPI_Unpack(rawBuffer, bufSize, &position, &(*std::begin(vol_)), vol_.Size(),
             MPI_DOUBLE, MPI_COMM_WORLD);
}

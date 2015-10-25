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
#include "slices.hpp"
#include "procBlock.hpp"
#include "boundaryConditions.hpp"  // interblock

using std::cout;
using std::endl;
using std::cerr;

// constructors for geomSlice class
geomSlice::geomSlice() {
  parBlock_ = 0;

  center_ = multiArray3d<vector3d<double>>(1, 1, 1);
  fAreaI_ = multiArray3d<unitVec3dMag<double>>(1, 1, 1);
  fAreaJ_ = fAreaI_;
  fAreaK_ = fAreaI_;
  fCenterI_ = multiArray3d<vector3d<double>>(1, 1, 1);
  fCenterJ_ = fCenterI_;
  fCenterK_ = fCenterI_;

  vol_ = multiArray3d<double>(1, 1, 1);
}

// constructor -- initialize state_ vector with dummy variables
geomSlice::geomSlice(const int &li, const int &lj, const int &lk,
                     const int &pblk) {
  // li -- size of direction i (cell)
  // lj -- size of direction j (cell)
  // lk -- size of direction k (cell)
  // pblk -- parent block that slice is coming from

  parBlock_ = pblk;

  center_ = multiArray3d<vector3d<double>>(li, lj, lk);
  fAreaI_ = multiArray3d<unitVec3dMag<double>> (li + 1, lj, lk);
  fAreaJ_ = multiArray3d<unitVec3dMag<double>> (li, lj + 1, lk);
  fAreaK_ = multiArray3d<unitVec3dMag<double>> (li, lj, lk + 1);
  fCenterI_ = multiArray3d<vector3d<double>> (li + 1, lj, lk);
  fCenterJ_ = multiArray3d<vector3d<double>> (li, lj + 1, lk);
  fCenterK_ = multiArray3d<vector3d<double>> (li, lj, lk + 1);

  vol_ = multiArray3d<double>(li, lj, lk);
}

/* constructor to get a slice (portion) of the geometry of a procBlock. The
geom slice remains in the orientation of its parent block unless
the revI, revJ, or revK flags are activated to reverse either the i, j, or k
directions respectively. This function essentially cuts out a portion
of a procBlock and returns a geomSlice with that portion of the block geometry.
*/
geomSlice::geomSlice(const procBlock &blk, const int &is, const int &ie,
                    const int &js, const int &je, const int &ks, const int &ke,
                    const bool revI, const bool revJ, const bool revK) {
  // blk -- procBlock to extract slice from
  // is -- starting i-cell index for slice
  // ie -- ending i-cell index for slice
  // js -- starting j-cell index for slice
  // je -- ending j-cell index for slice
  // ks -- starting k-cell index for slice
  // ke -- ending k-cell index for slice
  // revI -- flag to reverse i direction of indices (default false)
  // revJ -- flag to reverse j direction of indices (default false)
  // revK -- flag to reverse k direction of indices (default false)

  // allocate dimensions
  int numI = ie - is + 1;
  int numJ = je - js + 1;
  int numK = ke - ks + 1;
  parBlock_ = blk.ParentBlock();

  // allocate size for vectors
  center_ = multiArray3d<vector3d<double>>(numI, numJ, numK);
  fAreaI_ = multiArray3d<unitVec3dMag<double>>(numI + 1, numJ, numK);
  fAreaJ_ = multiArray3d<unitVec3dMag<double>>(numI, numJ + 1, numK);
  fAreaK_ = multiArray3d<unitVec3dMag<double>>(numI, numJ, numK + 1);
  fCenterI_ = multiArray3d<vector3d<double>>(numI + 1, numJ, numK);
  fCenterJ_ = multiArray3d<vector3d<double>>(numI, numJ + 1, numK);
  fCenterK_ = multiArray3d<vector3d<double>>(numI, numJ, numK + 1);
  vol_ = multiArray3d<double>(numI, numJ, numK);

  // loop over all cells in slice and populate
  for (int kk = 0; kk < vol_.NumK(); kk++) {
    for (int jj = 0; jj < vol_.NumJ(); jj++) {
      for (int ii = 0; ii < vol_.NumI(); ii++) {
        // determine if direction needs to be reversed
        int k = revK ? numK - 1 - kk : kk;
        double kFac = revK ? -1.0 : 1.0;

        int j = revJ ? numJ - 1 - jj : jj;
        double jFac = revJ ? -1.0 : 1.0;

        int i = revI ? numI - 1 - ii : ii;
        double iFac = revI ? -1.0 : 1.0;

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

// constructors for stateSlice class
stateSlice::stateSlice() {
  parBlock_ = 0;
  state_ = multiArray3d<primVars>(1, 1, 1);
}
// constructor -- initialize state vector with dummy variables
stateSlice::stateSlice(const int &li, const int &lj, const int &lk,
                       const int &pblk) {
  // li -- size of direction i
  // lj -- size of direction j
  // lk -- size of direction k
  // pblk -- parent block that slice is coming from

  parBlock_ = pblk;
  state_ = multiArray3d<primVars>(li, lj, lk);
}

/*Member function to pack a stateslice into a buffer, swap it with its
 * interblock partner, and then unpack it into a stateslice.*/
void stateSlice::PackSwapUnpackMPI(const interblock &inter,
                                   const MPI_Datatype &MPI_cellData,
                                   const int &rank) {
  // inter -- interblock boundary for the swap
  // MPI_cellData -- MPI datatype to pass cell state_ data
  // rank -- processor rank

  // swap with mpi_send_recv_replace
  // pack data into buffer, but first get size
  int bufSize = 0;
  int tempSize = 0;
  MPI_Pack_size(this->NumCells(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for states
  bufSize += tempSize;
  // add size for ints in class stateSlice, and 3 ints for multiArray3d dims
  MPI_Pack_size(4, MPI_INT, MPI_COMM_WORLD,
                &tempSize);
  bufSize += tempSize;

  char *buffer = new char[bufSize];  // allocate buffer to pack data into

  // pack data into buffer
  int numI = this->NumI();
  int numJ = this->NumJ();
  int numK = this->NumK();
  int position = 0;
  MPI_Pack(&numI, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numJ, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numK, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&state_(0, 0, 0), this->NumCells(), MPI_cellData, buffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&parBlock_, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);

  MPI_Status status;
  if (rank == inter.RankFirst()) {  // send/recv with second entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankSecond(), 1,
                         inter.RankSecond(), 1, MPI_COMM_WORLD, &status);
  } else {  // send/recv with first entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankFirst(), 1,
                         inter.RankFirst(), 1, MPI_COMM_WORLD, &status);
  }

  // put slice back into stateSlice
  position = 0;
  MPI_Unpack(buffer, bufSize, &position, &numI, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numJ, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &numK, 1, MPI_INT,
             MPI_COMM_WORLD);
  // resize slice
  state_.SameSizeResize(numI, numJ, numK);

  MPI_Unpack(buffer, bufSize, &position, &state_(0, 0, 0),
             this->NumCells(), MPI_cellData, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &parBlock_, 1, MPI_INT,
             MPI_COMM_WORLD);

  delete[] buffer;
}


/* constructor to get a slice (portion) of the states of a procBlock. The
state slice remains in the orientation of its parent block unless
the revI, revJ, or revK flags are activated to reverse either the i, j, or k
directions respectively. This function essentially cuts out a portion
of a procBlock and returns a stateSlice with that portion of the block states.
*/
stateSlice::stateSlice(const procBlock &blk, const int &is, const int &ie,
                       const int &js, const int &je, const int &ks,
                       const int &ke, const bool revI, const bool revJ,
                       const bool revK) {
  // blk -- procBlock to get slice from
  // is -- starting i-cell index for slice
  // ie -- ending i-cell index for slice
  // js -- starting j-cell index for slice
  // je -- ending j-cell index for slice
  // ks -- starting k-cell index for slice
  // ke -- ending k-cell index for slice
  // revI -- flag to reverse i direction of indices (default false)
  // revJ -- flag to reverse j direction of indices (default false)
  // revK -- flag to reverse k direction of indices (default false)

  // allocate dimensions
  int numI = ie - is + 1;
  int numJ = je - js + 1;
  int numK = ke - ks + 1;
  parBlock_ = blk.ParentBlock();

  // allocate size for vectors
  state_ = multiArray3d<primVars>(numI, numJ, numK);

  // loop over all cells in slice and populate
  for (int kk = 0; kk < state_.NumK(); kk++) {
    for (int jj = 0; jj < state_.NumJ(); jj++) {
      for (int ii = 0; ii < state_.NumI(); ii++) {
        // determine if direction needs to be reversed
        int k = revK ? numK - 1 - kk : kk;
        int j = revJ ? numJ - 1 - jj : jj;
        int i = revI ? numI - 1 - ii : ii;

        // assign cell variables
        state_(ii, jj, kk) = blk.State(is + i, js + j, ks + k);
      }
    }
  }
}

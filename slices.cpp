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
#include <vector>
#include <string>
#include "slices.hpp"
#include "procBlock.hpp"
#include "boundaryConditions.hpp"  // interblock
#include "plot3d.hpp"              // location functions

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

// constructors for geomSlice class
geomSlice::geomSlice() {
  numCells_ = 1;
  numI_ = 1;
  numJ_ = 1;
  numK_ = 1;
  parBlock_ = 0;

  int numFaces = (numI_ + 1) * (numJ_) * (numK_);

  center_ = vector<vector3d<double> >(numCells_);
  fAreaI_ = vector<unitVec3dMag<double> >(numFaces);
  fAreaJ_ = fAreaI_;
  fAreaK_ = fAreaI_;
  fCenterI_ = vector<vector3d<double> >(numFaces);
  fCenterJ_ = fCenterI_;
  fCenterK_ = fCenterI_;

  vol_ = vector<double>(numCells_);
}

// constructor -- initialize state_ vector with dummy variables
geomSlice::geomSlice(const int &li, const int &lj, const int &lk,
                     const int &pblk) {
  // li -- size of direction i (cell)
  // lj -- size of direction j (cell)
  // lk -- size of direction k (cell)
  // pblk -- parent block that slice is coming from

  numI_ = li;
  numJ_ = lj;
  numK_ = lk;

  numCells_ = li * lj * lk;

  parBlock_ = pblk;

  int numIFaces = (numI_ + 1) * (numJ_) * (numK_);
  int numJFaces = (numI_) * (numJ_ + 1) * (numK_);
  int numKFaces = (numI_) * (numJ_) * (numK_ + 1);

  center_ = vector<vector3d<double> >(numCells_);
  fAreaI_ = vector<unitVec3dMag<double> > (numIFaces);
  fAreaJ_ = vector<unitVec3dMag<double> > (numJFaces);
  fAreaK_ = vector<unitVec3dMag<double> > (numKFaces);
  fCenterI_ = vector<vector3d<double> > (numIFaces);
  fCenterJ_ = vector<vector3d<double> > (numJFaces);
  fCenterK_ = vector<vector3d<double> > (numKFaces);

  vol_ = vector<double>(numCells_);
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
  numI_ = ie - is + 1;
  numJ_ = je - js + 1;
  numK_ = ke - ks + 1;
  numCells_ = numI_ * numJ_ * numK_;
  parBlock_ = blk.ParentBlock();

  // allocate size for vectors
  center_.resize(numCells_);
  fAreaI_.resize((numI_ + 1) * numJ_ * numK_);
  fAreaJ_.resize(numI_ * (numJ_ + 1) * numK_);
  fAreaK_.resize(numI_ * numJ_ * (numK_ + 1));
  fCenterI_.resize((numI_ + 1) * numJ_ * numK_);
  fCenterJ_.resize(numI_ * (numJ_ + 1) * numK_);
  fCenterK_.resize(numI_ * numJ_ * (numK_ + 1));
  vol_.resize(numCells_);

  // get parent block maxes
  int imaxPar = blk.NumI() + 2.0 * blk.NumGhosts();
  int jmaxPar = blk.NumJ() + 2.0 * blk.NumGhosts();

  // loop over all cells in slice and populate
  for (int kk = 0; kk < numK_; kk++) {
    for (int jj = 0; jj < numJ_; jj++) {
      for (int ii = 0; ii < numI_; ii++) {
        // determine if direction needs to be reversed
        int k = revK ? numK_ - 1 - kk : kk;
        double kFac = revK ? -1.0 : 1.0;

        int j = revJ ? numJ_ - 1 - jj : jj;
        double jFac = revJ ? -1.0 : 1.0;

        int i = revI ? numI_ - 1 - ii : ii;
        double iFac = revI ? -1.0 : 1.0;

        // cell locations
        int locPar = GetLoc1D(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int loc = GetLoc1D(ii, jj, kk, numI_, numJ_);

        // lower i-face locations
        int lowIPar = GetLowerFaceI(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int lowI = GetLowerFaceI(ii, jj, kk, numI_, numJ_);

        // upper i-face locations
        int upIPar = GetUpperFaceI(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int upI = GetUpperFaceI(ii, jj, kk, numI_, numJ_);

        // lower j-face locations
        int lowJPar = GetLowerFaceJ(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int lowJ = GetLowerFaceJ(ii, jj, kk, numI_, numJ_);

        // upper j-face locations
        int upJPar = GetUpperFaceJ(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int upJ = GetUpperFaceJ(ii, jj, kk, numI_, numJ_);

        // lower k-face locations
        int lowKPar = GetLowerFaceK(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int lowK = GetLowerFaceK(ii, jj, kk, numI_, numJ_);

        // upper k-face locations
        int upKPar = GetUpperFaceK(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int upK = GetUpperFaceK(ii, jj, kk, numI_, numJ_);

        // assign cell variables
        vol_[loc] = blk.Vol(locPar);
        center_[loc] = blk.Center(locPar);

        // assign i-face variables
        if (revI) {  // if direction is reversed, upper/lower faces need to be
                     // swapped
          fAreaI_[lowI] = iFac * blk.FAreaI(upIPar);
          fCenterI_[lowI] = blk.FCenterI(upIPar);

          if (ii ==
              numI_ - 1) {  // at end of i-line assign upper face values too
            fAreaI_[upI] = iFac * blk.FAreaI(lowIPar);
            fCenterI_[upI] = blk.FCenterI(lowIPar);
          }
        } else {
          fAreaI_[lowI] = iFac * blk.FAreaI(lowIPar);
          fCenterI_[lowI] = blk.FCenterI(lowIPar);

          if (ii ==
              numI_ - 1) {  // at end of i-line assign upper face values too
            fAreaI_[upI] = iFac * blk.FAreaI(upIPar);
            fCenterI_[upI] = blk.FCenterI(upIPar);
          }
        }

        // assign j-face variables
        if (revJ) {  // if direction is reversed, upper/lower faces need to be
                     // swapped
          fAreaJ_[lowJ] = jFac * blk.FAreaJ(upJPar);
          fCenterJ_[lowJ] = blk.FCenterJ(upJPar);

          if (jj ==
              numJ_ - 1) {  // at end of j-line assign upper face values too
            fAreaJ_[upJ] = jFac * blk.FAreaJ(lowJPar);
            fCenterJ_[upJ] = blk.FCenterJ(lowJPar);
          }
        } else {
          fAreaJ_[lowJ] = jFac * blk.FAreaJ(lowJPar);
          fCenterJ_[lowJ] = blk.FCenterJ(lowJPar);

          if (jj ==
              numJ_ - 1) {  // at end of j-line assign upper face values too
            fAreaJ_[upJ] = jFac * blk.FAreaJ(upJPar);
            fCenterJ_[upJ] = blk.FCenterJ(upJPar);
          }
        }

        // assign k-face variables
        if (revK) {  // if direction is reversed, upper/lower faces need to be
                     // swapped
          fAreaK_[lowK] = kFac * blk.FAreaK(upKPar);
          fCenterK_[lowK] = blk.FCenterK(upKPar);

          if (kk ==
              numK_ - 1) {  // at end of k-line assign upper face values too
            fAreaK_[upK] = kFac * blk.FAreaK(lowKPar);
            fCenterK_[upK] = blk.FCenterK(lowKPar);
          }
        } else {
          fAreaK_[lowK] = kFac * blk.FAreaK(lowKPar);
          fCenterK_[lowK] = blk.FCenterK(lowKPar);

          if (kk ==
              numK_ - 1) {  // at end of k-line assign upper face values too
            fAreaK_[upK] = kFac * blk.FAreaK(upKPar);
            fCenterK_[upK] = blk.FCenterK(upKPar);
          }
        }
      }
    }
  }
}

// constructors for stateSlice class
stateSlice::stateSlice() {
  numCells_ = 1;
  numI_ = 1;
  numJ_ = 1;
  numK_ = 1;
  parBlock_ = 0;

  state_ = vector<primVars>(numCells_);
}
// constructor -- initialize state vector with dummy variables
stateSlice::stateSlice(const int &li, const int &lj, const int &lk,
                       const int &pblk) {
  // li -- size of direction i
  // lj -- size of direction j
  // lk -- size of direction k
  // pblk -- parent block that slice is coming from

  numI_ = li;
  numJ_ = lj;
  numK_ = lk;

  numCells_ = li * lj * lk;

  parBlock_ = pblk;

  state_ = vector<primVars>(numCells_);
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
  MPI_Pack_size((*this).NumCells(), MPI_cellData, MPI_COMM_WORLD,
                &tempSize);  // add size for states
  bufSize += tempSize;
  MPI_Pack_size(5, MPI_INT, MPI_COMM_WORLD,
                &tempSize);  // add size for ints in class stateSlice
  bufSize += tempSize;

  char *buffer = new char[bufSize];  // allocate buffer to pack data into

  // pack data into buffer
  int position = 0;
  MPI_Pack(&(*this).state_[0], (*this).NumCells(), MPI_cellData, buffer,
           bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numCells_, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numI_, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numJ_, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).numK_, 1, MPI_INT, buffer, bufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlock_, 1, MPI_INT, buffer, bufSize, &position,
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
  MPI_Unpack(buffer, bufSize, &position, &(*this).state_[0], (*this).NumCells(),
             MPI_cellData, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numCells_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numI_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numJ_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numK_, 1, MPI_INT,
             MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlock_, 1, MPI_INT,
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
  numI_ = ie - is + 1;
  numJ_ = je - js + 1;
  numK_ = ke - ks + 1;
  numCells_ = numI_ * numJ_ * numK_;
  parBlock_ = blk.ParentBlock();

  // allocate size for vectors
  state_.resize(numCells_);

  // get parent block maxes
  int imaxPar = blk.NumI() + 2.0 * blk.NumGhosts();
  int jmaxPar = blk.NumJ() + 2.0 * blk.NumGhosts();

  // loop over all cells in slice and populate
  for (int kk = 0; kk < numK_; kk++) {
    for (int jj = 0; jj < numJ_; jj++) {
      for (int ii = 0; ii < numI_; ii++) {
        // determine if direction needs to be reversed
        int k = revK ? numK_ - 1 - kk : kk;
        int j = revJ ? numJ_ - 1 - jj : jj;
        int i = revI ? numI_ - 1 - ii : ii;

        // cell locations
        int locPar = GetLoc1D(is + i, js + j, ks + k, imaxPar, jmaxPar);
        int loc = GetLoc1D(ii, jj, kk, numI_, numJ_);

        // assign cell variables
        state_[loc] = blk.State(locPar);
      }
    }
  }
}

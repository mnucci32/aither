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

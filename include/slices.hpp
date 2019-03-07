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

#ifndef SLICESHEADERDEF  // only if the macro SLICESHEADERDEF is not
                            // defined execute these lines of code
#define SLICESHEADERDEF  // define the macro

/* This header contains the gradients class.

   The gradients class stores the gradient terms for the inviscid fluxes
   viscous fluxes, and source terms. */

#include <iostream>
#include "mpi.h"
#include "primitive.hpp"
#include "vector3d.hpp"
#include "multiArray3d.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::cerr;

// forward class declaration
class procBlock;
class connection;
class range;

class geomSlice {
  multiArray3d<vector3d<double>> center_;  // coordinates of cell center
  multiArray3d<unitVec3dMag<double>> fAreaI_;  // face area vector for i-faces
  multiArray3d<unitVec3dMag<double>> fAreaJ_;  // face area vector for j-faces
  multiArray3d<unitVec3dMag<double>> fAreaK_;  // face area vector for k-faces
  multiArray3d<vector3d<double>> fCenterI_;  // coordinates of i-face centers
  multiArray3d<vector3d<double>> fCenterJ_;  // coordinates of j-face centers
  multiArray3d<vector3d<double>> fCenterK_;  // coordinates of k-face centers
  multiArray3d<double> vol_;  // cell volume

  int parBlock_;  // parent block number

 public:
  // constructors
  geomSlice(const int &, const int &, const int &, const int &);
  geomSlice() : geomSlice(1, 1, 1, 0) {}
  geomSlice(const procBlock &, const range &, const range &, const range &,
            const bool = false, const bool = false, const bool = false);

  // move constructor and assignment operator
  geomSlice(geomSlice&&) noexcept = default;
  geomSlice& operator=(geomSlice&&) noexcept = default;

  // copy constructor and assignment operator
  geomSlice(const geomSlice&) = default;
  geomSlice& operator=(const geomSlice&) = default;

  // member functions
  int NumCells() const { return vol_.Size(); }
  int NumI() const { return vol_.NumI(); }
  int NumJ() const { return vol_.NumJ(); }
  int NumK() const { return vol_.NumK(); }
  int GhostLayers() const { return vol_.GhostLayers(); }
  int ParentBlock() const { return parBlock_; }

  double Vol(const int &ii, const int &jj, const int &kk) const {
    return vol_(ii, jj, kk);
  }
  vector3d<double> Center(const int &ii, const int &jj, const int &kk) const {
    return center_(ii, jj, kk);
  }
  unitVec3dMag<double> FAreaI(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaI_(ii, jj, kk);
  }
  unitVec3dMag<double> FAreaJ(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaJ_(ii, jj, kk);
  }
  unitVec3dMag<double> FAreaK(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaK_(ii, jj, kk);
  }
  vector3d<double> FCenterI(const int &ii, const int &jj, const int &kk) const {
    return fCenterI_(ii, jj, kk);
  }
  vector3d<double> FCenterJ(const int &ii, const int &jj, const int &kk) const {
    return fCenterJ_(ii, jj, kk);
  }
  vector3d<double> FCenterK(const int &ii, const int &jj, const int &kk) const {
    return fCenterK_(ii, jj, kk);
  }

  // destructor
  ~geomSlice() noexcept {}
};


// function definitions

#endif

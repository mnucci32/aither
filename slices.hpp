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

#ifndef SLICESHEADERDEF  // only if the macro SLICESHEADERDEF is not
                            // defined execute these lines of code
#define SLICESHEADERDEF  // define the macro

/* This header contains the gradients class.

   The gradients class stores the gradient terms for the inviscid fluxes
   viscous fluxes, and source terms. */

#include <vector>  // vector
#include <string>  // string
#include <iostream>
#include "primVars.hpp"
#include "vector3d.hpp"
#include "macros.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class procBlock;

class geomSlice {
  vector<vector3d<double> > center_;  // coordinates of cell center_
  vector<vector3d<double> > fAreaI_;  // face area vector for i-faces
  vector<vector3d<double> > fAreaJ_;  // face area vector for j-faces
  vector<vector3d<double> > fAreaK_;  // face area vector for k-faces
  vector<vector3d<double> > fCenterI_;  // coordinates of i-face centers
  vector<vector3d<double> > fCenterJ_;  // coordinates of j-face centers
  vector<vector3d<double> > fCenterK_;  // coordinates of k-face centers

  vector<double> vol_;  // cell volume

  int numCells_;  // number of cells in block
  int numI_;  // i-dimension of block (cells)
  int numJ_;  // j-dimension of block (cells)
  int numK_;  // k-dimension of block (cells)
  int parBlock_;  // parent block number

 public:
  // constructors
  geomSlice();
  geomSlice(const int &, const int &, const int &, const int &);
  geomSlice(const procBlock &, const int &, const int &, const int &,
            const int &, const int &, const int &, const bool = false,
            const bool = false, const bool = false);

  // member functions
  int NumCells() const { return numCells_; }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }
  int ParentBlock() const { return parBlock_; }

  double Vol(const int &ind) const { return vol_[ind]; }
  vector3d<double> Center(const int &ind) const { return center_[ind]; }
  vector3d<double> FAreaI(const int &ind) const { return fAreaI_[ind]; }
  vector3d<double> FAreaJ(const int &ind) const { return fAreaJ_[ind]; }
  vector3d<double> FAreaK(const int &ind) const { return fAreaK_[ind]; }
  vector3d<double> FCenterI(const int &ind) const { return fCenterI_[ind]; }
  vector3d<double> FCenterJ(const int &ind) const { return fCenterJ_[ind]; }
  vector3d<double> FCenterK(const int &ind) const { return fCenterK_[ind]; }

  // destructor
  ~geomSlice() {}
};

class stateSlice {
  vector<primVars> state_;  // cell states

  int numCells_;  // number of cells in block
  int numI_;  // i-dimension of block (cells)
  int numJ_;  // j-dimension of block (cells)
  int numK_;  // k-dimension of block (cells)
  int parBlock_;  // parent block number

 public:
  // constructors
  stateSlice();
  stateSlice(const int &, const int &, const int &, const int &);

  stateSlice(const procBlock &, const int &, const int &, const int &,
             const int &, const int &, const int &, const bool = false,
             const bool = false, const bool = false);

  // member functions
  int NumCells() const { return numCells_; }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }
  int ParentBlock() const { return parBlock_; }

  primVars State(const int &ind) const { return state_[ind]; }

  void PackSwapUnpackMPI(const interblock &, const MPI_Datatype &, const int &);

  // destructor
  ~stateSlice() {}
};

// function definitions

#endif

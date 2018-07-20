/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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

#ifndef GRIDLEVEL_HEADER_DEF
#define GRIDLEVEL_HEADER_DEF

#include <vector>
#include <string>
#include "procBlock.hpp"
#include "boundaryConditions.hpp"
#include "vector3d.hpp"
#include "matMultiArray3d.hpp"
#include "mpi.h"

using std::string;
using std::vector;

// forward class declaration
class plot3dBlock;
class input;
class decomposition;
class input;
class physics;
class residual;
class kdtree;

class gridLevel {
  vector<procBlock> blocks_;
  vector<connection> connections_;
  vector<matMultiArray3d> diagonal_;

 public:
  // Constructor
  gridLevel(const vector<plot3dBlock>& mesh,
            const vector<boundaryConditions>& bcs, const decomposition& decomp,
            const physics& phys, const vector<vector3d<int>>& origGridSizes,
            const string& restartFile, input& inp, residual& first);
  gridLevel(const int& numBlocks) : blocks_(numBlocks) {}
  gridLevel() : gridLevel(1) {}

  // move constructor and assignment operator
  gridLevel(gridLevel&&) noexcept = default;
  gridLevel& operator=(gridLevel&&) noexcept = default;

  // copy constructor and assignment operator
  gridLevel(const gridLevel&) = default;
  gridLevel& operator=(const gridLevel&) = default;

  // Member functions
  int NumBlocks() const { return blocks_.size(); }
  const vector<procBlock>& Blocks() const { return blocks_; }
  const procBlock& Block(const int &ii) const { return blocks_[ii]; }
  procBlock& Block(const int &ii) { return blocks_[ii]; }

  int NumConnections() const { return connections_.size(); }
  const vector<connection>& Connections() const { return connections_; }
  const connection& Connection(const int& ii) const { return connections_[ii]; }
  connection& Connection(const int& ii) { return connections_[ii]; }

  gridLevel SendGridLevel(const int& rank, const int& numProcBlock,
                          const MPI_Datatype& MPI_vec3d,
                          const MPI_Datatype& MPI_vec3dMag,
                          const MPI_Datatype& MPI_connection,
                          const input& inp) const;
  void GetGridLevel(const gridLevel& local, const int& rank,
                    const MPI_Datatype& MPI_uncoupledScalar,
                    const MPI_Datatype& MPI_vec3d,
                    const MPI_Datatype& MPI_tensorDouble, const input& inp);
  void CalcWallDistance(const kdtree& tree);
  void AssignSolToTimeN(const physics& phys);
  void AssignSolToTimeNm1();

  void CalcTimeStep(const input& inp);
  void ExplicitUpdate(const input& inp, const physics& phys, const int& mm,
                      residual& residL2, resid& residLinf);
  double ImplicitUpdate(const input& inp, const physics& phys, const int& mm,
                        residual& residL2, resid& residLinf, const int& rank);
  void ResizeMatrix(const input& inp, const int& numProcBlock);

  void GetBoundaryConditions(const input& inp, const physics& phys,
                             const int& rank);
  void SwapWallDist(const int& rank, const int& numGhosts);
  void SwapTurbVars(const int& rank, const int& numGhosts);
  void SwapEddyViscAndGradients(const int& rank,
                                const MPI_Datatype& MPI_tensorDouble,
                                const MPI_Datatype& MPI_vec3d,
                                const int& numGhosts);
  void CalcResidual(const physics& phys, const input& inp, const int& rank,
                    const MPI_Datatype& MPI_tensorDouble,
                    const MPI_Datatype& MPI_vec3d);
  void AuxillaryAndWidths(const physics& phys);

  // Destructor
  ~gridLevel() noexcept {}
};

#endif

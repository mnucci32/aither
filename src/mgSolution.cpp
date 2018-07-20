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

#include <iostream>     // cout
#include <cstdlib>      // exit()
#include <vector>
#include "mgSolution.hpp"
#include "gridLevel.hpp"
#include "parallel.hpp"
#include "input.hpp"
#include "boundaryConditions.hpp"
#include "plot3d.hpp"
#include "physicsModels.hpp"
#include "output.hpp"
#include "resid.hpp"
#include "vector3d.hpp"
#include "macros.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

// constructor
mgSolution::mgSolution(const int &numLevels) {
  MSG_ASSERT(numLevels > 0, "must at least have 1 multigrid level");
  solution_.reserve(numLevels);
}


// member functions
void mgSolution::ConstructFinestLevel(
    const vector<plot3dBlock>& mesh, const vector<boundaryConditions>& bcs,
    const decomposition& decomp, const physics& phys,
    const vector<vector3d<int>>& origGridSizes, const string& restartFile,
    input& inp, residual& first) {
  MSG_ASSERT(solution_.size() == 0U,
             "should only be called once to initialize");
  // inputs correspond to finest mesh
  solution_.emplace_back(mesh, bcs, decomp, phys, origGridSizes, restartFile,
                         inp, first);
}

mgSolution mgSolution::SendFinestGridLevel(const int& rank,
                                           const int& numProcBlock,
                                           const MPI_Datatype& MPI_vec3d,
                                           const MPI_Datatype& MPI_vec3dMag,
                                           const MPI_Datatype& MPI_connection,
                                           const input& inp) const {
  mgSolution local(inp.MultiGridLevels());
  local.solution_.emplace_back(this->Finest().SendGridLevel(
      rank, numProcBlock, MPI_vec3d, MPI_vec3dMag, MPI_connection, inp));
  return local;
}

void mgSolution::GetFinestGridLevel(const mgSolution& local, const int& rank,
                                    const MPI_Datatype& MPI_uncoupledScalar,
                                    const MPI_Datatype& MPI_vec3d,
                                    const MPI_Datatype& MPI_tensorDouble,
                                    const input& inp) {
  solution_[0].GetGridLevel(local.Finest(), rank, MPI_uncoupledScalar,
                            MPI_vec3d, MPI_tensorDouble, inp);
}

void mgSolution::ConstructMultigrids() {}

void mgSolution::AuxillaryAndWidths(const physics& phys) {
  for (auto& sol : solution_) {
    sol.AuxillaryAndWidths(phys);
  }
}

void mgSolution::CalcWallDistance(const kdtree &tree) {
  for (auto &sol : solution_) {
    sol.CalcWallDistance(tree);
  }
}

void mgSolution::SwapWallDist(const int& rank, const int& numGhosts) {
  for (auto &sol : solution_) {
    sol.SwapWallDist(rank, numGhosts);
  }
}

void mgSolution::ResizeMatrix(const input& inp, const int& numProcBlock) {
  for (auto &sol : solution_) {
    sol.ResizeMatrix(inp, numProcBlock);
  }
}

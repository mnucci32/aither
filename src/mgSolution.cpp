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
mgSolution::mgSolution(const int &numLevels, const int &cycle) {
  solution_.reserve(numLevels);
  mgCycleIndex_ = cycle;
}

// member functions
void mgSolution::ConstructFinestLevel(
    const vector<plot3dBlock>& mesh, const vector<boundaryConditions>& bcs,
    const decomposition& decomp, const physics& phys,
   const string& restartFile,
    input& inp, residual& first) {
  MSG_ASSERT(solution_.size() == 0U,
             "should only be called once to initialize");
  // inputs correspond to finest mesh
  // get original grid sizes (before decomposition)
  vector<vector3d<int>> gridSizes;
  gridSizes.reserve(mesh.size());
  for (const auto& msh : mesh) {
    gridSizes.emplace_back(msh.NumCellsI(), msh.NumCellsJ(), msh.NumCellsK());
  }

  solution_.emplace_back(mesh, bcs, decomp, phys, gridSizes, restartFile, inp,
                         first);
}

mgSolution mgSolution::SendFinestGridLevel(const int& rank,
                                           const int& numProcBlock,
                                           const MPI_Datatype& MPI_vec3d,
                                           const MPI_Datatype& MPI_vec3dMag,
                                           const MPI_Datatype& MPI_connection,
                                           const input& inp) const {
  mgSolution local(inp);
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

void mgSolution::ConstructMultigrids(const decomposition& decomp,
                                     const input& inp, const physics& phys) {
  const auto numLevels = solution_.capacity();
  while (solution_.size() < numLevels) {
    solution_.push_back(solution_.back().Coarsen(decomp, inp, phys));
  }
}

void mgSolution::AuxillaryAndWidths(const physics& phys) {
  for (auto& sol : solution_) {
    sol.AuxillaryAndWidths(phys);
  }
}

void mgSolution::StoreOldSolution(const input& inp, const physics& phys,
                                  const int& iter) {
  // Store time-n solution, for time integration methods that require it
  for (auto& sol : solution_) {
    if (inp.NeedToStoreTimeN()) {
      sol.AssignSolToTimeN(phys);
      if (!inp.IsRestart() && inp.IsMultilevelInTime() && iter == 0) {
        sol.AssignSolToTimeNm1();
      }
    }
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

// multigrid restriction - fine grid to coarse grid operator
void mgSolution::Restriction(
    const int& fi, const int& mm,
    const vector<blkMultiArray3d<varArray>>& matrixResid, const input& inp,
    const physics& phys, const int& rank, const MPI_Datatype& MPI_tensorDouble,
    const MPI_Datatype& MPI_vec3d) {
  MSG_ASSERT(fi >= 0 && fi < static_cast<int>(solution_.size() - 1),
             "index for restriction out of range");
  solution_[fi].Restriction(solution_[fi + 1], mm, matrixResid, inp, phys, rank,
                            MPI_tensorDouble, MPI_vec3d);
}

vector<blkMultiArray3d<varArray>> mgSolution::Relax(const int& ll,
                                                    const int& sweeps,
                                                    const physics& phys,
                                                    const input& inp,
                                                    const int& rank) {
  return solution_[ll].Relax(phys, inp, rank, sweeps);
}

void mgSolution::SubtractFromUpdate(
    const int& ll, const vector<blkMultiArray3d<varArray>>& coarseDu) {
  solution_[ll].SubtractFromUpdate(coarseDu);
}

// multigrid prolongation - coarse grid to fine grid operator
void mgSolution::Prolongation(const int &ci) {
  MSG_ASSERT(ci > 0 && ci < static_cast<int>(solution_.size()),
             "index for prolongation out of range");
  solution_[ci].Prolongation(solution_[ci - 1]);
}

double mgSolution::CycleAtLevel(const int& fl, const int& mm,
                                const physics& phys, const input& inp,
                                const int& rank,
                                const MPI_Datatype& MPI_tensorDouble,
                                const MPI_Datatype& MPI_vec3d) {
  // fl -- index for fine grid level
  MSG_ASSERT(fl >= 0 && fl < static_cast<int>(solution_.size()),
             "index for multigrid cycle out of range");
  // initialize matrix residual
  vector<blkMultiArray3d<varArray>> matrixResid;

  if (fl == this->NumGridLevels() - 1) {  // recursive base case
    matrixResid = this->Relax(fl, inp.MatrixSweeps(), phys, inp, rank);
  } else {
    // pre-relaxation sweeps
    matrixResid = this->Relax(fl, ceil(inp.MatrixSweeps()/2), phys, inp, rank);

    // coarse grid correction
    // restrict solution, matrix residual, and implicit update to coarse grid
    auto cl = fl + 1;
    this->Restriction(fl, mm, matrixResid, inp, phys, rank, MPI_tensorDouble,
                      MPI_vec3d);
    const auto coarseDu = solution_[cl].Update();

    // recursive call to next coarse level
    for (auto ii = 0; ii < mgCycleIndex_; ++ii) {
      this->CycleAtLevel(cl, mm, phys, inp, rank, MPI_tensorDouble, MPI_vec3d);
    }
    this->SubtractFromUpdate(cl, coarseDu);

    // interpolate coarse level correction and add to solution
    this->Prolongation(cl);

    // post-relaxation sweeps
    matrixResid = this->Relax(fl, ceil(inp.MatrixSweeps()/2), phys, inp, rank);
  }

  // calculate l2 norm of matrix residual
  auto l2Resid = 0.0;
  auto totalSize = 0;
  for (auto& mr : matrixResid) {
    mr *= mr;
    l2Resid += std::accumulate(std::begin(mr), std::end(mr), 0.0);
    totalSize += mr.Size();
  }
  return l2Resid / totalSize;
}

double mgSolution::ImplicitUpdate(const input& inp,
                                  const physics& phys, const int& mm,
                                  const int& rank,
                                  const MPI_Datatype& MPI_tensorDouble,
                                  const MPI_Datatype& MPI_vec3d,
                                  residual& residL2, resid& residLinf) {
  // inp -- input variables
  // phys -- physics models
  // mm -- nonlinear iteration
  // residL2 -- L2 residual
  // residLinf -- L infinity residual

  // initialize matrix error
  auto matrixError = 0.0;

  // add volume and time term and calculate inverse of main diagonal
  const auto fl = this->FinestIndex();
  solution_[fl].InvertDiagonal(inp);

  // initialize matrix update
  solution_[fl].InitializeMatrixUpdate(inp, phys);

  // Solve Ax=b with supported solver and multigrid
  matrixError =
      this->CycleAtLevel(fl, mm, phys, inp, rank, MPI_tensorDouble, MPI_vec3d);

  // Update blocks
  solution_[fl].UpdateBlocks(inp, phys, mm, residL2, residLinf);

  // Reset main diagonal on all levels
  for (auto level = fl; level < this->NumGridLevels(); ++level) {
    solution_[level].ResetDiagonal();
  }

  return matrixError;
}

double mgSolution::Iterate(const input& inp, const physics& phys,
                           const MPI_Datatype& MPI_tensorDouble,
                           const MPI_Datatype& MPI_vec3d, const int& mm,
                           const int& rank, residual& residL2,
                           resid& residLinf) {
  const auto fl = this->FinestIndex();
  // Get boundary conditions for all blocks
  solution_[fl].GetBoundaryConditions(inp, phys, rank);

  // Calculate residual (RHS)
  solution_[fl].CalcResidual(phys, inp, rank, MPI_tensorDouble, MPI_vec3d);

  // Calculate time step
  solution_[fl].CalcTimeStep(inp);

  auto matrixResid = 0.0;
  if (inp.IsImplicit()) {
    matrixResid = this->ImplicitUpdate(inp, phys, mm, rank, MPI_tensorDouble,
                                       MPI_vec3d, residL2, residLinf);
  } else {  // explicit time integration
    solution_[fl].ExplicitUpdate(inp, phys, mm, residL2, residLinf);
  }
  return matrixResid;
}

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

// multigrid restriction - fine grid to coarse grid operator
vector<blkMultiArray3d<varArray>> mgSolution::Restriction(
    const int& fi, const vector<blkMultiArray3d<varArray>>& du) {
  MSG_ASSERT(fi >= 0 && fi < static_cast<int>(solution_.size() - 1),
             "index for restriction out of range");
  return solution_[fi].Restriction(solution_[fi + 1], du);
}

double mgSolution::Relax(const int& ll, const int& sweeps, const physics& phys,
                       const input& inp, const int& rank,
                       const unique_ptr<linearSolver>& solver,
                       vector<blkMultiArray3d<varArray>>& du) const {
  return solver->Relax(solution_[ll], phys, inp, rank, sweeps, du);
}

double mgSolution::CycleAtLevel(const int& fl, const physics& phys,
                                const input& inp, const int& rank,
                                const unique_ptr<linearSolver>& solver,
                                vector<blkMultiArray3d<varArray>>& du) {
  // fl -- index for fine grid level
  MSG_ASSERT(fl >= 0 && fl < static_cast<int>(solution_.size()),
             "index for multigrid cycle out of range");
  auto matrixResid = 0.0;
  if (fl == this->NumGridLevels() - 1) {  // recursive base case
    matrixResid =
        this->Relax(fl, inp.MatrixSweeps(), phys, inp, rank, solver, du);
  } else {
    // pre-relaxation sweeps
    this->Relax(fl, 1, phys, inp, rank, solver, du);

    // coarse grid correction
    // restrict solution, residual, and implicit update to coarse grid
    auto cl = fl + 1;
    auto coarseDu = this->Restriction(fl, du);

    // recursive call to next coarse level
    for (auto ii = 0; ii < mgCycleIndex_; ++ii) {
      matrixResid = this->CycleAtLevel(cl, phys, inp, rank, solver, coarseDu);
    }
    for (auto ii = 0U; ii < du.size(); ++ii) {
      du[ii] -= coarseDu[ii];
    }

    // interpolate coarse level correction and add to solution
    this->Prolongation(cl, du);

    // post-relaxation sweeps
    matrixResid = this->Relax(fl, 1, phys, inp, rank, solver, du);
  }
  return matrixResid;
}

double mgSolution::MultigridCycle(const physics& phys, const input& inp,
                                  const int& rank,
                                  const unique_ptr<linearSolver>& solver) {
  auto du = this->Finest().InitializeMatrixUpdate(inp, phys);
  return this->CycleAtLevel(0, phys, inp, rank, solver, du);
}

double mgSolution::ImplicitUpdate(const input& inp, const physics& phys,
                                  const unique_ptr<linearSolver>& solver,
                                  const int& mm, residual& residL2,
                                  resid& residLinf, const int& rank) {
  // inp -- input variables
  // phys -- physics models
  // solver -- linear solver
  // mm -- nonlinear iteration
  // residL2 -- L2 residual
  // residLinf -- L infinity residual

  // initialize matrix error
  auto matrixError = 0.0;

  // add volume and time term and calculate inverse of main diagonal
  solution_[0].InvertDiagonal(inp);

  // initialize matrix update
  auto du = solution_[0].InitializeMatrixUpdate(inp, phys);

  // Solve Ax=b with supported solver and multigrid
  matrixError = this->CycleAtLevel(0, phys, inp, rank, solver, du);

  // Update blocks and reset main diagonal
  solution_[0].UpdateBlocks(inp, phys, mm, du, residL2, residLinf);

  return matrixError;
}

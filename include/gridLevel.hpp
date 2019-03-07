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

#ifndef GRIDLEVEL_HEADER_DEF
#define GRIDLEVEL_HEADER_DEF

#include <vector>
#include <array>
#include <string>
#include <memory>
#include "procBlock.hpp"
#include "boundaryConditions.hpp"
#include "vector3d.hpp"
#include "matMultiArray3d.hpp"
#include "linearSolver.hpp"
#include "mpi.h"

using std::string;
using std::unique_ptr;
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
  std::unique_ptr<linearSolver> solver_;

  // during restriction, traverse fine grid values in lexigraphical order,
  // applying volume weight factor, and adding to coarse grid, also in
  // lexigraphical order -- can use multArray3d<>
  // -----
  // during prolongation, traverse fine grid values in lexigraphical order,
  // use multiArray3d<> to get coarse grid index, then can use neighboring
  // coarse grid cells for trilinear interpolation
  vector<multiArray3d<vector3d<int>>> toCoarse_;
  vector<multiArray3d<double>> volWeightFactor_;
  vector<multiArray3d<std::array<double, 7>>> prolongCoeffs_;
  vector<blkMultiArray3d<varArray>> mgForcing_;

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
  void InvertDiagonal(const input &);
  void InitializeMatrixUpdate(const input&, const physics&);
  double Relax(const physics& phys, const input& inp, const int& rank,
             const int& sweeps) {
    return solver_->Relax(*this, phys, inp, rank, sweeps);
  }

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
  void UpdateBlocks(const input& inp, const physics& phys, const int& mm,
                    residual& residL2, resid& residLinf);

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
  gridLevel Coarsen(const decomposition& decomp, const input& inp,
                    const physics& phys);
  void Restriction(gridLevel& coarse) const;
  void Prolongation(gridLevel& fine) const;
  void SubtractFromUpdate(const vector<blkMultiArray3d<varArray>>& coarseDu);
  vector<blkMultiArray3d<varArray>> Update() const { return solver_->X(); }

  // Destructor
  ~gridLevel() noexcept {}
};

// function declarations
template <typename T>
void BlockProlongation(const T& coarse,
                       const multiArray3d<vector3d<int>>& toCoarse,
                       const multiArray3d<std::array<double, 7>>& coeffs,
                       T& fine) {
  // get coarse data at nodes
  const auto coarseNodes = ConvertCellToNode(coarse, true);
  // use trilinear interpolation
  for (auto kk = fine.StartK(); kk < fine.EndK(); ++kk) {
    for (auto jj = fine.StartJ(); jj < fine.EndJ(); ++jj) {
      for (auto ii = fine.StartI(); ii < fine.EndI(); ++ii) {
        const auto ci = toCoarse(ii, jj, kk);
        // data at coarse cell nodes
        const auto d0 = coarseNodes(ci.X()    , ci.Y()    , ci.Z());
        const auto d1 = coarseNodes(ci.X() + 1, ci.Y()    , ci.Z());
        const auto d2 = coarseNodes(ci.X()    , ci.Y() + 1, ci.Z());
        const auto d3 = coarseNodes(ci.X() + 1, ci.Y() + 1, ci.Z());
        const auto d4 = coarseNodes(ci.X()    , ci.Y()    , ci.Z() + 1);
        const auto d5 = coarseNodes(ci.X() + 1, ci.Y()    , ci.Z() + 1);
        const auto d6 = coarseNodes(ci.X()    , ci.Y() + 1, ci.Z() + 1);
        const auto d7 = coarseNodes(ci.X() + 1, ci.Y() + 1, ci.Z() + 1);
        const auto prolong =
            TrilinearInterp(coeffs(ii, jj, kk), d0, d1, d2, d3, d4, d5, d6, d7);
        fine.InsertBlock(ii, jj, kk, prolong);
      }
    }
  }
}

// member functions

#endif

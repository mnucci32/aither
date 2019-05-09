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
#include <string>
#include "gridLevel.hpp"
#include "utility.hpp"
#include "parallel.hpp"
#include "input.hpp"
#include "boundaryConditions.hpp"
#include "plot3d.hpp"
#include "physicsModels.hpp"
#include "output.hpp"
#include "resid.hpp"
#include "vector3d.hpp"
#include "kdtree.hpp"
#include "matMultiArray3d.hpp"
#include "linearSolver.hpp"
#include "macros.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// constructor
gridLevel::gridLevel(const vector<plot3dBlock>& mesh,
                     const vector<boundaryConditions>& bcs,
                     const decomposition& decomp, const physics& phys,
                     const vector<vector3d<int>>& origGridSizes,
                     const string& restartFile, input& inp, residual& first) {
  MSG_ASSERT(mesh.size() == bcs.size(), "block size mismatch");
  connections_ = GetConnectionBCs(bcs, mesh, decomp, inp);
  blocks_.reserve(mesh.size());
  mgForcing_.reserve(mesh.size());
  for (auto ll = 0U; ll < mesh.size(); ++ll) {
    blocks_.emplace_back(mesh[ll], decomp.ParentBlock(ll), bcs[ll], ll,
                         decomp.Rank(ll), decomp.LocalPosition(ll), inp);
    blocks_.back().InitializeStates(inp, phys);
    blocks_.back().AssignGhostCellsGeom();
    mgForcing_.emplace_back(
        blocks_.back().NumI(), blocks_.back().NumJ(), blocks_.back().NumK(), 0,
        blocks_.back().NumEquations(), blocks_.back().NumSpecies(), 0);
  }

  // if restart, get data from restart file
  if (inp.IsRestart()) {
    ReadRestart(*this, restartFile, decomp, inp, phys, first, origGridSizes);
  }

  // Swap geometry for interblock BCs
  for (auto& conn : connections_) {
    if (conn.IsInterblock()) {
      SwapGeomSlice(conn, blocks_[conn.BlockFirst()],
                    blocks_[conn.BlockSecond()]);
    }
  }
  // Get ghost cell edge data
  for (auto& block : blocks_) {
    block.AssignGhostCellsGeomEdge();
  }

  // Setup linear solver
  if (inp.IsImplicit()) {
    solver_ = inp.AssignLinearSolver(*this);
  }
}

/* Function to send procBlocks to their appropriate processor. This function is
called after the decomposition has been run. The procBlock data all
resides on the ROOT processor. In this function, the ROOT processor packs the
procBlocks and sends them to the appropriate processor. All the non-ROOT
processors receive and unpack the data from ROOT. This is used to send the
geometric block data from ROOT to all the processors at the beginning of the
simulation.
*/
gridLevel gridLevel::SendGridLevel(const int& rank, const int& numProcBlock,
                                   const MPI_Datatype& MPI_vec3d,
                                   const MPI_Datatype& MPI_vec3dMag,
                                   const MPI_Datatype& MPI_connection,
                                   const input& inp) const {
  // *this -- global gridLevel, only meaningful on ROOT processor
  // rank -- proc rank, used to determine if process should send or receive
  // numProcBlock -- number of procBlocks that the processor should have. (All
  //                 processors may give different values).
  // MPI_vec3d -- MPI_Datatype used for vector3d<double>  transmission
  // MPI_vec3dMag -- MPI_Datatype used for unitVec3dMag<double>  transmission
  // MPI_connection -- MPI_Datatype used for connection  transmission
  // input -- input variables

  gridLevel local(numProcBlock);
  //------------------------------------------------------------------------
  //                                  ROOT
  //------------------------------------------------------------------------
  if (rank == ROOTP) {  // may have to pack and send data
    // loop over ALL blocks
    for (const auto &global : blocks_) {
      // no need to send data because it is already on root processor
      if (global.Rank() == ROOTP) {
        const auto lp = global.LocalPosition();
        local.blocks_[lp] = global;
        local.mgForcing_[lp] = blkMultiArray3d<varArray>(
            global.NumI(), global.NumJ(), global.NumK(), 0,
            global.NumEquations(), global.NumSpecies());
      } else {  // send data to receiving processors
        // pack and send procBlock
        global.PackSendGeomMPI(MPI_vec3d, MPI_vec3dMag);
      }
    }
    //--------------------------------------------------------------------------
    //                                NON - ROOT
    //--------------------------------------------------------------------------
  } else {  // receive and unpack data (non-root)
    for (auto ii = 0; ii < numProcBlock; ++ii) {
      // recv and unpack procBlock
      procBlock tempBlock;
      tempBlock.RecvUnpackGeomMPI(MPI_vec3d, MPI_vec3dMag, inp);

      // add procBlock to output vector
      const auto lp = tempBlock.LocalPosition();
      local.blocks_[lp] = tempBlock;
      local.mgForcing_[lp] = blkMultiArray3d<varArray>(
          tempBlock.NumI(), tempBlock.NumJ(), tempBlock.NumK(), 0,
          tempBlock.NumEquations(), tempBlock.NumSpecies());
    }
  }

  // now do connections -------------------------------------------------------
  // first determine the number of connections and send that to all processors
  if (rank == ROOTP) {
    local.connections_ = connections_;
  }
  auto numCon = local.NumConnections();
  MPI_Bcast(&numCon, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  local.connections_.resize(numCon);  // allocate space to receive connections

  // broadcast all connections to all processors
  MPI_Bcast(&(*std::begin(local.connections_)), local.connections_.size(),
            MPI_connection, ROOTP, MPI_COMM_WORLD);

  // now allocate memory for linear solver
  local.solver_ = inp.AssignLinearSolver(local);

  return local;
}

/* Function to send procBlocks to the root processor. In this function, the
non-ROOT processors pack the procBlocks and send them to the ROOT processor.
The ROOT processor receives and unpacks the data from the non-ROOT processors.
This is used to get all the data on the ROOT processor to write out results.
*/

void gridLevel::GetGridLevel(const gridLevel& local, const int& rank,
                             const MPI_Datatype& MPI_uncoupledScalar,
                             const MPI_Datatype& MPI_vec3d,
                             const MPI_Datatype& MPI_tensorDouble,
                             const input& inp) {
  // *this -- global gridLevels - this is only meaningful on ROOT processor
  // local -- gridLevels local to each processor. These are sent to ROOT
  // rank -- processor rank, to determine if process should send or receive
  // MPI_uncoupledScalar -- MPI_Datatype used for uncoupledScalar transmission
  // MPI_vec3d -- MPI_Datatype used for vector3d<double> transmission
  // MPI_tensorDouble -- MPI_Datatype used for tensor<double> transmission
  // input -- input variables

  //--------------------------------------------------------------------------
  //                                      ROOT
  //--------------------------------------------------------------------------
  if (rank == ROOTP) {  // may have to recv and unpack data
    // loop over ALL blocks
    for (auto &global : blocks_) {
      if (global.Rank() == ROOTP) {  // data already on ROOT processor
        // assign local block to global block in order of local vector
        global = local.blocks_[global.LocalPosition()];
      } else {  // recv data from sending processors
        global.RecvUnpackSolMPI(MPI_uncoupledScalar, MPI_vec3d,
                                MPI_tensorDouble, inp);
      }
    }
    //-------------------------------------------------------------------------
    //                                   NON - ROOT
    //-------------------------------------------------------------------------
  } else {  // pack and send data (non-root)
    // get vector of local positions
    vector<std::pair<int, int>> localPos;
    localPos.reserve(local.NumBlocks());
    for (const auto &lb : local.blocks_) {
      localPos.push_back(std::make_pair(lb.LocalPosition(), lb.GlobalPos()));
    }
    // sort by global position
    // need to send data in order of global position, not local position to
    // prevent deadlock
    std::sort(
        std::begin(localPos), std::end(localPos),
        [](const auto &d1, const auto &d2) { return d1.second < d2.second; });

    for (auto &lp : localPos) {
      local.blocks_[lp.first].PackSendSolMPI(MPI_uncoupledScalar, MPI_vec3d,
                                             MPI_tensorDouble);
    }
  }
}

// function to calculate the distance to the nearest viscous wall of all
// cell centers
void gridLevel::CalcWallDistance(const kdtree &tree) {
  for (auto &block : blocks_) {
    block.CalcWallDistance(tree);
  }
}

void gridLevel::AssignSolToTimeN(const physics &phys) {
  for (auto &block : blocks_) {
    block.AssignSolToTimeN(phys);
  }
}

void gridLevel::AssignSolToTimeNm1() {
  for (auto &block : blocks_) {
    block.AssignSolToTimeNm1();
  }
}

void gridLevel::CalcTimeStep(const input &inp) {
  // states -- vector of all procBlocks on processor
  // inp -- input variables
  for (auto &block : blocks_) {
    // calculate time step
    block.CalcBlockTimeStep(inp);
  }
}


void gridLevel::ExplicitUpdate(const input& inp, const physics& phys,
                               const int& mm, residual& residL2,
                               resid& residLinf) {
  // create dummy update (not used in explicit update)
  blkMultiArray3d<varArray> du;
  // loop over all blocks and update
  for (auto &block : blocks_) {
    block.UpdateBlock(inp, phys, du, mm, residL2, residLinf);
  }
}

void gridLevel::SwapWallDist(const int& rank, const int& numGhosts) {
  // rank -- processor rank
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections_) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      blocks_[conn.LocalBlockFirst()].SwapWallDistSlice(
          conn, blocks_[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      blocks_[conn.LocalBlockFirst()].SwapWallDistSliceMPI(conn, rank);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      blocks_[conn.LocalBlockSecond()].SwapWallDistSliceMPI(conn, rank);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

/* Function to populate ghost cells with proper cell states for inviscid flow
calculation. This function operates on the entire grid and uses connection
boundaries to pass the correct data between grid blocks.
*/
void gridLevel::GetBoundaryConditions(const input& inp, const physics& phys,
                                      const int& rank) {
  // inp -- all input variables
  // phys -- physics models
  // rank -- processor rank

  // loop over all blocks and assign inviscid ghost cells
  for (auto &block : blocks_) {
    block.AssignInviscidGhostCells(inp, phys);
  }

  // loop over connections and swap ghost cells where needed
  for (auto &conn : connections_) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection on this processor, swap w/o mpi
      blocks_[conn.LocalBlockFirst()].SwapStateSlice(
          conn, blocks_[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      blocks_[conn.LocalBlockFirst()].SwapStateSliceMPI(conn, rank);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      blocks_[conn.LocalBlockSecond()].SwapStateSliceMPI(conn, rank);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }

  // loop over all blocks and get ghost cell edge data
  for (auto &block : blocks_) {
    block.AssignInviscidGhostCellsEdge(inp, phys);
  }
}

void gridLevel::SwapTurbVars(const int& rank, const int& numGhosts) {
  // rank -- processor rank
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections_) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      blocks_[conn.LocalBlockFirst()].SwapTurbSlice(
          conn, blocks_[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      blocks_[conn.LocalBlockFirst()].SwapTurbSliceMPI(conn, rank);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      blocks_[conn.LocalBlockSecond()].SwapTurbSliceMPI(conn, rank);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

void gridLevel::SwapEddyViscAndGradients(const int& rank,
                                         const MPI_Datatype& MPI_tensorDouble,
                                         const MPI_Datatype& MPI_vec3d,
                                         const int& numGhosts) {
  // rank -- processor rank
  // MPI_tensorDouble -- MPI datatype for tensor<double>
  // MPI_vec3d -- MPI datatype for vector3d<double>
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections_) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      blocks_[conn.LocalBlockFirst()].SwapEddyViscAndGradientSlice(
          conn, blocks_[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      blocks_[conn.LocalBlockFirst()].SwapEddyViscAndGradientSliceMPI(
          conn, rank, MPI_tensorDouble, MPI_vec3d);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      blocks_[conn.LocalBlockSecond()].SwapEddyViscAndGradientSliceMPI(
          conn, rank, MPI_tensorDouble, MPI_vec3d);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

void gridLevel::CalcResidual(const physics& phys, const input& inp,
                             const int& rank,
                             const MPI_Datatype& MPI_tensorDouble,
                             const MPI_Datatype& MPI_vec3d) {
  // phys -- physics models
  // inp -- input variables
  // rank -- processor rank
  // MPI_tensorDouble -- MPI datatype for tensor<double>
  // MPI_vec3d -- MPI datatype for vector3d<double>

  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    // calculate residual
    blocks_[bb].CalcResidualNoSource(phys, inp, solver_->A(bb));
  }
  // swap mut & gradients calculated during residual calculation
  this->SwapEddyViscAndGradients(rank, MPI_tensorDouble, MPI_vec3d,
                                 inp.NumberGhostLayers());

  if (inp.IsRANS()) {
    // swap turbulence variables calculated during residual calculation
    this->SwapTurbVars(rank, inp.NumberGhostLayers());
  }
  if (inp.IsRANS() || phys.Chemistry()->IsReacting()) {
    for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
      // calculate source terms for residual
      blocks_[bb].CalcSrcTerms(phys, inp, solver_->A(bb));
    }
  }
}

void gridLevel::InvertDiagonal(const input& inp) {
  // add volume and time term and calculate inverse of main diagonal
  solver_->InvertDiagonal(*this, inp);
}

void gridLevel::ResetDiagonal() {
  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    solver_->ZeroA(bb);
  }
}

void gridLevel::InitializeMatrixUpdate(const input& inp, const physics& phys) {
  solver_->InitializeMatrixUpdate(*this, inp, phys);
}

void gridLevel::UpdateBlocks(const input& inp, const physics& phys,
                             const int& mm,
                             residual& residL2, resid& residLinf) {
  // Update blocks
  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    // Update solution
    blocks_[bb].UpdateBlock(inp, phys, solver_->X(bb), mm, residL2, residLinf);

    // Assign time n to time n-1 at end of nonlinear iterations
    if (inp.IsMultilevelInTime() && mm == inp.NonlinearIterations() - 1) {
      blocks_[bb].AssignSolToTimeNm1();
    }
  }
}

void gridLevel::AuxillaryAndWidths(const physics& phys) {
  for (auto& block : blocks_) {
    block.UpdateAuxillaryVariables(phys, false);
    block.CalcCellWidths();
  }
}

gridLevel gridLevel::Coarsen(const decomposition& decomp, const input& inp,
                             const physics& phys) {
  // get plot3dBlocks and bcs for coarsened grid level
  vector<plot3dBlock> coarseMesh;
  coarseMesh.reserve(this->NumBlocks());
  vector<boundaryConditions> coarseBCs;
  coarseBCs.reserve(this->NumBlocks());
  toCoarse_.reserve(this->NumBlocks());
  volWeightFactor_.reserve(this->NumBlocks());
  for (const auto& blk : blocks_) {
    blk.GetCoarseMeshAndBCs(coarseMesh, coarseBCs, toCoarse_, volWeightFactor_);
  }

  gridLevel coarse;
  coarse.connections_ = GetConnectionBCs(coarseBCs, coarseMesh, decomp, inp);
  coarse.blocks_.reserve(coarseMesh.size());
  coarse.mgForcing_.reserve(coarseMesh.size());
  for (auto ll = 0U; ll < coarseMesh.size(); ++ll) {
    coarse.blocks_.emplace_back(coarseMesh[ll], blocks_[ll].ParentBlock(),
                                coarseBCs[ll], ll, blocks_[ll].Rank(),
                                blocks_[ll].LocalPosition(), inp);
    coarse.blocks_.back().InitializeStates(inp, phys);
    coarse.blocks_.back().AssignGhostCellsGeom();
    coarse.mgForcing_.emplace_back(
        coarse.blocks_.back().NumI(), coarse.blocks_.back().NumJ(),
        coarse.blocks_.back().NumK(), 0, coarse.blocks_.back().NumEquations(),
        coarse.blocks_.back().NumSpecies(), 0);
  }

  // Swap geometry for interblock BCs
  for (auto& conn : coarse.connections_) {
    if (conn.IsInterblock()) {
      SwapGeomSlice(conn, coarse.blocks_[conn.BlockFirst()],
                    coarse.blocks_[conn.BlockSecond()]);
    }
  }
  // Get ghost cell edge data
  for (auto& block : coarse.blocks_) {
    block.AssignGhostCellsGeomEdge();
  }

  // Calculate prolongation coefficients
  coarse.prolongCoeffs_.reserve(this->NumBlocks());
  for (auto ll = 0; ll < this->NumBlocks(); ++ll) {
    // size coeffs for fine grid
    coarse.prolongCoeffs_.emplace_back(blocks_[ll].NumI(), blocks_[ll].NumJ(),
                                       blocks_[ll].NumK(), 0);
    // loop over fine grid cells
    for (auto kk = blocks_[ll].StartK(); kk < blocks_[ll].EndK(); ++kk) {
      for (auto jj = blocks_[ll].StartJ(); jj < blocks_[ll].EndJ(); ++jj) {
        for (auto ii = blocks_[ll].StartI(); ii < blocks_[ll].EndI(); ++ii) {
          const auto ci = toCoarse_[ll](ii, jj, kk);
          // coordinates of fine cell center
          const auto fc = blocks_[ll].Center(ii, jj, kk);
          // nodal coordinates of bounding coarse cell
          const auto c0 = coarse.blocks_[ll].Node(ci.X(), ci.Y(), ci.Z());
          const auto c1 = coarse.blocks_[ll].Node(ci.X() + 1, ci.Y(), ci.Z());
          const auto c2 = coarse.blocks_[ll].Node(ci.X(), ci.Y() + 1, ci.Z());
          const auto c3 =
              coarse.blocks_[ll].Node(ci.X() + 1, ci.Y() + 1, ci.Z());
          const auto c4 = coarse.blocks_[ll].Node(ci.X(), ci.Y(), ci.Z() + 1);
          const auto c5 =
              coarse.blocks_[ll].Node(ci.X() + 1, ci.Y(), ci.Z() + 1);
          const auto c6 =
              coarse.blocks_[ll].Node(ci.X(), ci.Y() + 1, ci.Z() + 1);
          const auto c7 =
              coarse.blocks_[ll].Node(ci.X() + 1, ci.Y() + 1, ci.Z() + 1);
          coarse.prolongCoeffs_.back()(ii, jj, kk) =
              TrilinearInterpCoeff(c0, c1, c2, c3, c4, c5, c6, c7, fc);
        }
      }
    }
  }

  // Setup linear solver
  if (inp.IsImplicit()) {
    coarse.solver_ = inp.AssignLinearSolver(coarse);
  }

  return coarse;
}

void gridLevel::Restriction(gridLevel& coarse, const int &mm,
                            const vector<blkMultiArray3d<varArray>>& fineResid,
                            const input& inp, const physics& phys,
                            const int& rank,
                            const MPI_Datatype& MPI_tensorDouble,
                            const MPI_Datatype& MPI_vec3d) const {
  MSG_ASSERT(blocks_.size() == coarse.blocks_.size(),
             "gridLevel size mismatch");
  MSG_ASSERT(blocks_.size() == fineResid.size(), "residual size mismatch");
  MSG_ASSERT(inp.IsImplicit(), "calling gridLevel::Restriction for explicit");

  for (auto ii = 0; ii < this->NumBlocks(); ++ii) {
    // restrict solution
    coarse.blocks_[ii].RestrictState(blocks_[ii], toCoarse_[ii],
                                     volWeightFactor_[ii]);
    if (mm == 0) {  // need to store solution at time n for linear solvers
      coarse.AssignSolToTimeN(phys);
    }

    // restrict update
    coarse.solver_->X(ii).Zero();
    BlockRestriction(solver_->X(ii), toCoarse_[ii], volWeightFactor_[ii],
                     coarse.solver_->X(ii));
  }

  // DEBUG
  // should calculate solution on coarse grid level here - need to get Ax
  // should do everything so linear solver can call relax
  // Get boundary conditions for all blocks
  coarse.GetBoundaryConditions(inp, phys, rank);
  // Calculate residual (RHS)
  coarse.CalcResidual(phys, inp, rank, MPI_tensorDouble, MPI_vec3d);
  // Calculate time step
  coarse.CalcTimeStep(inp);
  // add volume and time term and calculate inverse of main diagonal
  coarse.InvertDiagonal(inp);
  // initialize matrix update
  // coarse.InitializeMatrixUpdate(inp, phys);
  // get Ax for coarse level
  const auto ax = coarse.AX(phys, inp);

  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    // DEBUG -- forcing term should be Ax - r
    // Ax is Ax of coarse state (now possible since update and state were
    //    restricted)
    // r is b - Ax for fine state, restricted down to coarse level
    // restrict matrix residuals
    auto tmpFactor = volWeightFactor_[bb];
    tmpFactor.Zero(1.0);
    coarse.mgForcing_[bb].Zero();
    //BlockRestriction(fineResid[bb], toCoarse_[bb], volWeightFactor_[bb],
    //                 coarse.mgForcing_[bb]);
    BlockRestriction(fineResid[bb], toCoarse_[bb], tmpFactor,
                     coarse.mgForcing_[bb]);
    //cout << "COARSE AX" << endl;
    //ax[bb].PrintPhysical(cout);
    //cout << endl;
    //cout << "RESTRICTED MATRIX RESID" << endl;
    //coarse.mgForcing_[bb].PrintPhysical(cout);
    //cout << endl;
    // DEBUG -- doing this because ax and mgForcing have different num ghosts
    for (auto kk = coarse.mgForcing_[bb].StartK(); kk < coarse.mgForcing_[bb].EndK(); ++kk) {
      for (auto jj = coarse.mgForcing_[bb].StartJ(); jj < coarse.mgForcing_[bb].EndJ(); ++jj) {
        for (auto ii = coarse.mgForcing_[bb].StartI(); ii < coarse.mgForcing_[bb].EndI(); ++ii) {
          coarse.mgForcing_[bb].InsertBlock(
              ii, jj, kk,
              ax[bb](ii, jj, kk) - coarse.mgForcing_[bb](ii, jj, kk));
        }
      }
    }
    //cout << "RESTRICTED MATRIX RESID AFTER OPERATION" << endl;
    //coarse.mgForcing_[bb].PrintPhysical(cout);
    //cout << endl;

    // coarse.mgForcing_[bb] = ax[bb] - coarse.mgForcing_[bb];
  }
}

void gridLevel::SubtractFromUpdate(
    const vector<blkMultiArray3d<varArray>>& coarseDu) {
  solver_->SubtractFromUpdate(coarseDu);
}

void gridLevel::Prolongation(gridLevel& fine) const {
  MSG_ASSERT(blocks_.size() == fine.blocks_.size(), "gridLevel size mismatch");
  vector<blkMultiArray3d<varArray>> fineCorrVec;
  fineCorrVec.reserve(this->NumBlocks());
  for (auto ii = 0; ii < this->NumBlocks(); ++ii) {
    blkMultiArray3d<varArray> fineCorrection(
        fine.blocks_[ii].NumI(), fine.blocks_[ii].NumJ(),
        fine.blocks_[ii].NumK(), fine.blocks_[ii].NumGhosts(),
        fine.blocks_[ii].NumEquations(), fine.blocks_[ii].NumSpecies());
  
    //if (blocks_[ii].NumI() < 3) {
    //  cout << "COARSE UPDATE BEFORE PROLONG" << endl;
    //  solver_->X(ii).PrintPhysical(cout);
    //}
  
    BlockProlongation(solver_->X(ii), fine.toCoarse_[ii], prolongCoeffs_[ii],
                      fineCorrection);
    
    //if (blocks_[ii].NumI() < 3) {
    //  cout << "PROLONGED UPDATE" << endl;
    //  fineCorrection.PrintPhysical(cout);
    //}
    
    fineCorrVec.push_back(fineCorrection);
  }
  
  //if (blocks_[0].NumI() < 3) {
  //  cout << "FINE UPDATE -- BEFORE ADD" << endl;
  //  fine.solver_->X(0).PrintPhysical(cout);
  //}
  
  fine.solver_->AddToUpdate(fineCorrVec);
  
  //if (blocks_[0].NumI() < 3) {
  //  cout << "PROLONGED UPDATE -- AFTER ADD" << endl;
  //  fine.solver_->X(0).PrintPhysical(cout);
  //}
  
}

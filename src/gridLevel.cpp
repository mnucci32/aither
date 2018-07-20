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
  auto connections_ = GetConnectionBCs(bcs, mesh, decomp, inp);
  blocks_.reserve(mesh.size());
  for (auto ll = 0U; ll < mesh.size(); ++ll) {
    blocks_.emplace_back(mesh[ll], decomp.ParentBlock(ll), bcs[ll], ll,
                         decomp.Rank(ll), decomp.LocalPosition(ll), inp);
    blocks_.back().InitializeStates(inp, phys);
    blocks_.back().AssignGhostCellsGeom();
  }

  // if restart, get data from restart file
  if (inp.IsRestart()) {
    ReadRestart(*this, restartFile, decomp, inp, phys, first, origGridSizes);
  }

  // Swap geometry for connection BCs
  for (auto& conn : connections_) {
    SwapGeomSlice(conn, blocks_[conn.BlockFirst()],
                  blocks_[conn.BlockSecond()]);
  }
  // Get ghost cell edge data
  for (auto& block : blocks_) {
    block.AssignGhostCellsGeomEdge();
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
        local.blocks_[global.LocalPosition()] = global;
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
      local.blocks_[tempBlock.LocalPosition()] = tempBlock;
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

void gridLevel::ResizeMatrix(const input& inp, const int &numProcBlock) {
  MSG_ASSERT(diagonal_.size() == 0,
             "only call when matrix diagonal hasn't been allocated");
  if (inp.IsImplicit()) {
    const auto fluxJac =
        inp.IsBlockMatrix()
            ? fluxJacobian(inp.NumFlowEquations(), inp.NumTurbEquations())
            : fluxJacobian(1, std::min(1, inp.NumTurbEquations()));

    diagonal_.reserve(numProcBlock);
    for (const auto& block : blocks_) {
      diagonal_.emplace_back(block.NumI(), block.NumJ(), block.NumK(), 0,
                             fluxJac);
    }
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
    blocks_[bb].CalcResidualNoSource(phys, inp, diagonal_[bb]);
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
      blocks_[bb].CalcSrcTerms(phys, inp, diagonal_[bb]);
    }
  }
}

double gridLevel::ImplicitUpdate(const input& inp, const physics& phys,
                                 const int& mm, residual& residL2,
                                 resid& residLinf, const int& rank) {
  // inp -- input variables
  // phys -- physics models
  // mm -- nonlinear iteration
  // residL2 -- L2 residual
  // residLinf -- L infinity residual

  // initialize matrix error
  auto matrixError = 0.0;

  const auto numG = this->Block(0).NumGhosts();

  // add volume and time term and calculate inverse of main diagonal
  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    blocks_[bb].InvertDiagonal(diagonal_[bb], inp);
  }

  // initialize matrix update
  vector<blkMultiArray3d<varArray>> du(this->NumBlocks());
  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    du[bb] = blocks_[bb].InitializeMatrixUpdate(inp, phys, diagonal_[bb]);
  }

  // Solve Ax=b with supported solver
  if (inp.MatrixSolver() == "lusgs" || inp.MatrixSolver() == "blusgs") {
    // calculate order by hyperplanes for each block
    vector<vector<vector3d<int>>> reorder(this->NumBlocks());
    for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
      reorder[bb] = HyperplaneReorder(blocks_[bb].NumI(), blocks_[bb].NumJ(),
                                      blocks_[bb].NumK());
    }

    // start sweeps through domain
    for (auto ii = 0; ii < inp.MatrixSweeps(); ii++) {
      // swap updates for ghost cells
      SwapImplicitUpdate(du, this->Connections(), rank, numG);

      // forward lu-sgs sweep
      for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
        blocks_[bb].LUSGS_Forward(reorder[bb], du[bb], phys, inp,
                                  diagonal_[bb], ii);
      }

      // swap updates for ghost cells
      SwapImplicitUpdate(du, this->Connections(), rank, numG);

      // backward lu-sgs sweep
      for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
        matrixError += blocks_[bb].LUSGS_Backward(reorder[bb], du[bb], phys,
                                                  inp, diagonal_[bb], ii);
      }
    }
  } else if (inp.MatrixSolver() == "dplur" || inp.MatrixSolver() == "bdplur") {
    for (auto ii = 0; ii < inp.MatrixSweeps(); ii++) {
      // swap updates for ghost cells
      SwapImplicitUpdate(du, this->Connections(), rank, numG);

      for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
        // Calculate correction (du)
        matrixError += blocks_[bb].DPLUR(du[bb], phys, inp, diagonal_[bb]);
      }
    }
  } else {
    cerr << "ERROR: Matrix solver " << inp.MatrixSolver() <<
        " is not recognized!" << endl;
    cerr << "Please choose lusgs, blusgs, dplur, or bdplur." << endl;
    exit(EXIT_FAILURE);
  }

  // Update blocks and reset main diagonal
  for (auto bb = 0; bb < this->NumBlocks(); ++bb) {
    // Update solution
    blocks_[bb].UpdateBlock(inp, phys, du[bb], mm, residL2, residLinf);

    // Assign time n to time n-1 at end of nonlinear iterations
    if (inp.IsMultilevelInTime() && mm == inp.NonlinearIterations() - 1) {
      blocks_[bb].AssignSolToTimeNm1();
    }

    // zero flux jacobians
    diagonal_[bb].Zero();
  }

  return matrixError;
}

void gridLevel::AuxillaryAndWidths(const physics& phys) {
  for (auto& block : blocks_) {
    block.UpdateAuxillaryVariables(phys, false);
    block.CalcCellWidths();
  }
}

/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

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

#include <iostream>               // cout, cerr, endl
#include <algorithm>              // max, min
#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include "utility.hpp"
#include "procBlock.hpp"
#include "eos.hpp"                 // equation of state
#include "transport.hpp"           // transport model
#include "thermodynamic.hpp"       // thermodynamic model
#include "input.hpp"               // inputVars
#include "genArray.hpp"
#include "turbulence.hpp"
#include "slices.hpp"
#include "fluxJacobian.hpp"
#include "kdtree.hpp"
#include "resid.hpp"
#include "primVars.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::min;
using std::unique_ptr;
using std::array;

/* Function to calculate the gradient of a vector at the cell center using the
Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a vector is calculated using the
Green-Gauss method. U is a vector. Aij is the area at face i (component j).
V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out
of the control volume.
 */
tensor<double> VectorGradGG(
    const vector3d<double> &vil, const vector3d<double> &viu,
    const vector3d<double> &vjl, const vector3d<double> &vju,
    const vector3d<double> &vkl, const vector3d<double> &vku,
    const vector3d<double> &ail, const vector3d<double> &aiu,
    const vector3d<double> &ajl, const vector3d<double> &aju,
    const vector3d<double> &akl, const vector3d<double> &aku,
    const double &vol) {
  // vil -- vector at the i-lower face of the cell at which the gradient
  //        is being calculated
  // viu -- vector at the i-upper face of the cell at which the gradient
  //        is being calculated
  // vjl -- vector at the j-lower face of the cell at which the gradient
  //        is being calculated
  // vju -- vector at the j-upper face of the cell at which the gradient
  //        is being calculated
  // vkl -- vector at the k-lower face of the cell at which the gradient
  //        is being calculated
  // vku -- vector at the k-upper face of the cell at which the gradient
  //        is being calculated

  // ail -- area vector at the lower i-face of the cell at which the gradient
  //        is being calculated
  // aiu -- area vector at the upper i-face of the cell at which the gradient
  //        is being calculated
  // ajl -- area vector at the lower j-face of the cell at which the gradient
  //        is being calculated
  // aju -- area vector at the upper j-face of the cell at which the gradient
  //        is being calculated
  // akl -- area vector at the lower k-face of the cell at which the gradient
  //        is being calculated
  // aku -- area vector at the upper k-face of the cell at which the gradient
  //        is being calculated

  // vol -- cell volume

  tensor<double> temp;
  const auto invVol = 1.0 / vol;

  // convention is for area vector to point out of cell, so lower values are
  // negative, upper are positive
  temp.SetXX(viu.X() * aiu.X() - vil.X() * ail.X() + vju.X() * aju.X() -
             vjl.X() * ajl.X() + vku.X() * aku.X() - vkl.X() * akl.X());
  temp.SetXY(viu.Y() * aiu.X() - vil.Y() * ail.X() + vju.Y() * aju.X() -
             vjl.Y() * ajl.X() + vku.Y() * aku.X() - vkl.Y() * akl.X());
  temp.SetXZ(viu.Z() * aiu.X() - vil.Z() * ail.X() + vju.Z() * aju.X() -
             vjl.Z() * ajl.X() + vku.Z() * aku.X() - vkl.Z() * akl.X());

  temp.SetYX(viu.X() * aiu.Y() - vil.X() * ail.Y() + vju.X() * aju.Y() -
             vjl.X() * ajl.Y() + vku.X() * aku.Y() - vkl.X() * akl.Y());
  temp.SetYY(viu.Y() * aiu.Y() - vil.Y() * ail.Y() + vju.Y() * aju.Y() -
             vjl.Y() * ajl.Y() + vku.Y() * aku.Y() - vkl.Y() * akl.Y());
  temp.SetYZ(viu.Z() * aiu.Y() - vil.Z() * ail.Y() + vju.Z() * aju.Y() -
             vjl.Z() * ajl.Y() + vku.Z() * aku.Y() - vkl.Z() * akl.Y());

  temp.SetZX(viu.X() * aiu.Z() - vil.X() * ail.Z() + vju.X() * aju.Z() -
             vjl.X() * ajl.Z() + vku.X() * aku.Z() - vkl.X() * akl.Z());
  temp.SetZY(viu.Y() * aiu.Z() - vil.Y() * ail.Z() + vju.Y() * aju.Z() -
             vjl.Y() * ajl.Z() + vku.Y() * aku.Z() - vkl.Y() * akl.Z());
  temp.SetZZ(viu.Z() * aiu.Z() - vil.Z() * ail.Z() + vju.Z() * aju.Z() -
             vjl.Z() * ajl.Z() + vku.Z() * aku.Z() - vkl.Z() * akl.Z());

  temp *= invVol;

  return temp;
}

/* Function to calculate the gradient of a scalar at the cell center using the
Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a scalar is calculated using the
Green-Gauss method. U is a scalar. Aij is the area at face i (component j).
V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out
of the control volume.
 */
vector3d<double> ScalarGradGG(
    const double &til, const double &tiu, const double &tjl, const double &tju,
    const double &tkl, const double &tku, const vector3d<double> &ail,
    const vector3d<double> &aiu, const vector3d<double> &ajl,
    const vector3d<double> &aju, const vector3d<double> &akl,
    const vector3d<double> &aku, const double &vol) {
  // til -- scalar value at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // tiu -- scalar value at the upper face of the cell at which the scalar
  //        gradient is being calculated
  // tjl -- scalar value at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // tju -- scalar value at the upper face of the cell at which the scalar
  //        gradient is being calculated
  // tkl -- scalar value at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // tku -- scalar value at the upper face of the cell at which the scalar
  //        gradient is being calculated

  // ail -- area vector at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // aiu -- area vector at the upper face of the cell at which the scalar
  //        gradient is being calculated
  // ajl -- area vector at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // aju -- area vector at the upper face of the cell at which the scalar
  //        gradient is being calculated
  // akl -- area vector at the lower face of the cell at which the scalar
  //        gradient is being calculated
  // aku -- area vector at the upper face of the cell at which the scalar
  //        gradient is being calculated

  // vol -- cell volume

  vector3d<double> temp;
  const auto invVol = 1.0 / vol;

  // define scalar gradient vector
  // convention is for area vector to point out of cell, so lower values are
  // negative, upper are positive
  temp.SetX(tiu * aiu.X() - til * ail.X() + tju * aju.X() -
            tjl * ajl.X() + tku * aku.X() - tkl * akl.X());
  temp.SetY(tiu * aiu.Y() - til * ail.Y() + tju * aju.Y() -
            tjl * ajl.Y() + tku * aku.Y() - tkl * akl.Y());
  temp.SetZ(tiu * aiu.Z() - til * ail.Z() + tju * aju.Z() -
            tjl * ajl.Z() + tku * aku.Z() - tkl * akl.Z());

  temp *= invVol;

  return temp;
}

/* Function to swap ghost cell geometry between two blocks at an connection
boundary. Slices are removed from the physical cells (extending into ghost cells
at the edges) of one block and inserted into the ghost cells of its partner
block. The reverse is also true. The slices are taken in the coordinate system
orientation of their parent block.

   Interior Cells    Ghost Cells               Ghost Cells   Interior Cells
   ________ ______|________ _________       _______________|_______ _________
Ui-3/2   Ui-1/2   |    Uj+1/2    Uj+3/2  Ui-3/2    Ui-1/2  |    Uj+1/2    Uj+3/2
  |        |      |        |         |     |        |      |       |         |
  | Ui-1   |  Ui  |  Uj    |  Uj+1   |     |  Ui-1  |   Ui |  Uj   |  Uj+1   |
  |        |      |        |         |     |        |      |       |         |
  |________|______|________|_________|     |________|______|_______|_________|
                  |                                        |

The above diagram shows the resulting values after the ghost cell swap. The
logic ensures that the ghost cells at the connection boundary exactly match
their partner block as if there were no separation in the grid.

Only 3 faces at each ghost cell need to be swapped (i.e. the lower face for Ui
is the upper face for Ui-1). At the end of a line (i-line, j-line or k-line),
both the upper and lower faces need to be swapped.
*/
void SwapGeomSlice(connection &inter, procBlock &blk1, procBlock &blk2) {
  // inter -- connection boundary information
  // blk1 -- first block involved in connection boundary
  // blk2 -- second block involved in connection boundary

  // Get indices for slice coming from first block to swap
  auto is1 = 0, ie1 = 0;
  auto js1 = 0, je1 = 0;
  auto ks1 = 0, ke1 = 0;

  inter.FirstSliceIndices(is1, ie1, js1, je1, ks1, ke1, blk1.NumGhosts());

  // Get indices for slice coming from second block to swap
  auto is2 = 0, ie2 = 0;
  auto js2 = 0, je2 = 0;
  auto ks2 = 0, ke2 = 0;

  inter.SecondSliceIndices(is2, ie2, js2, je2, ks2, ke2, blk2.NumGhosts());

  const auto geom1 = geomSlice(blk1, {is1, ie1}, {js1, je1}, {ks1, ke1});
  const auto geom2 = geomSlice(blk2, {is2, ie2}, {js2, je2}, {ks2, ke2});

  // change connections to work with slice and ghosts
  connection inter1 = inter;
  connection inter2 = inter;
  inter1.AdjustForSlice(false, blk1.NumGhosts());
  inter2.AdjustForSlice(true, blk2.NumGhosts());

  // put slices in proper blocks
  // return vector determining if any of the 4 edges of the connection need to
  // be updated for a "t" intersection
  const auto adjEdge1 = blk1.PutGeomSlice(geom2, inter2, blk2.NumGhosts());
  const auto adjEdge2 = blk2.PutGeomSlice(geom1, inter1, blk1.NumGhosts());

  // if an connection border needs to be updated, update
  for (auto ii = 0U; ii < adjEdge1.size(); ii++) {
    if (adjEdge1[ii]) {
      inter.UpdateBorderFirst(ii);
    }
    if (adjEdge2[ii]) {
      inter.UpdateBorderSecond(ii);
    }
  }
}


/* Function to populate ghost cells with proper cell states for inviscid flow
calculation. This function operates on the entire grid and uses connection
boundaries to pass the correct data between grid blocks.
*/
void GetBoundaryConditions(vector<procBlock> &states, const input &inp,
                           const unique_ptr<eos> &eqnState,
                           const unique_ptr<thermodynamic> &thermo,
                           const unique_ptr<transport> &trans,
                           const unique_ptr<turbModel> &turb,
                           vector<connection> &connections, const int &rank,
                           const MPI_Datatype &MPI_cellData) {
  // states -- vector of all procBlocks in the solution domain
  // inp -- all input variables
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // connections -- vector of connection boundary types
  // rank -- processor rank
  // MPI_cellData -- data type to pass primVars, genArray

  // loop over all blocks and assign inviscid ghost cells
  for (auto &state : states) {
    state.AssignInviscidGhostCells(inp, eqnState, thermo, trans, turb);
  }

  // loop over connections and swap ghost cells where needed
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection on this processor, swap w/o mpi
      states[conn.LocalBlockFirst()].SwapStateSlice(
          conn, states[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      states[conn.LocalBlockFirst()].SwapStateSliceMPI(conn, rank, MPI_cellData);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      states[conn.LocalBlockSecond()].SwapStateSliceMPI(conn, rank, MPI_cellData);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }

  // loop over all blocks and get ghost cell edge data
  for (auto &state : states) {
    state.AssignInviscidGhostCellsEdge(inp, eqnState, thermo, trans, turb);
  }
}


// function to get face centers of cells on viscous walls
vector<vector3d<double>> GetViscousFaceCenters(const vector<procBlock> &blks) {
  // blks -- vector of all procBlocks in simulation

  // get vector of BCs
  vector<boundaryConditions> bcs;
  bcs.reserve(blks.size());
  for (auto &blk : blks) {
    bcs.push_back(blk.BC());
  }

  // determine number of faces with viscous wall BC
  auto nFaces = 0;
  for (auto &bc : bcs) {
    nFaces += bc.NumViscousFaces();
  }

  // allocate vector for face centers
  vector<vector3d<double>> faceCenters;
  faceCenters.reserve(nFaces);

  // store viscous face centers
  auto aa = 0;
  for (auto &bc : bcs) {  // loop over BCs for each block
    for (auto bb = 0; bb < bc.NumSurfaces(); bb++) {  // loop over surfaces
      if (bc.GetBCTypes(bb) == "viscousWall") {
        // only store face center if surface is viscous wall

        if (bc.GetSurfaceType(bb) <= 2) {  // i-surface
          const auto ii = bc.GetIMin(bb);  // imin and imax are the same

          for (auto jj = bc.GetJMin(bb); jj < bc.GetJMax(bb); jj++) {
            for (auto kk = bc.GetKMin(bb); kk < bc.GetKMax(bb); kk++) {
              faceCenters.push_back(blks[aa].FCenterI(ii, jj, kk));
            }
          }
        } else if (bc.GetSurfaceType(bb) <= 4) {  // j-surface
          const auto jj = bc.GetJMin(bb);  // jmin and jmax are the same

          for (auto ii = bc.GetIMin(bb); ii < bc.GetIMax(bb); ii++) {
            for (auto kk = bc.GetKMin(bb); kk < bc.GetKMax(bb); kk++) {
              faceCenters.push_back(blks[aa].FCenterJ(ii, jj, kk));
            }
          }
        } else {  // k-surface
          const auto kk = bc.GetKMin(bb);  // kmin and kmax are the same

          for (auto ii = bc.GetIMin(bb); ii < bc.GetIMax(bb); ii++) {
            for (auto jj = bc.GetJMin(bb); jj < bc.GetJMax(bb); jj++) {
              faceCenters.push_back(blks[aa].FCenterK(ii, jj, kk));
            }
          }
        }
      }
    }
    aa++;
  }
  return faceCenters;
}


// function to calculate the distance to the nearest viscous wall of all
// cell centers
void CalcWallDistance(vector<procBlock> &localBlocks, const kdtree &tree) {
  for (auto &block : localBlocks) {
    block.CalcWallDistance(tree);
  }
}

void AssignSolToTimeN(vector<procBlock> &blocks,
                      const unique_ptr<eos> &eqnState,
                      const unique_ptr<thermodynamic> &thermo) {
  for (auto &block : blocks) {
    block.AssignSolToTimeN(eqnState, thermo);
  }
}

void AssignSolToTimeNm1(vector<procBlock> &blocks) {
  for (auto &block : blocks) {
    block.AssignSolToTimeNm1();
  }
}

void ExplicitUpdate(vector<procBlock> &blocks, const input &inp,
                    const unique_ptr<eos> &eqnState,
                    const unique_ptr<thermodynamic> &thermo,
                    const unique_ptr<transport> &trans,
                    const unique_ptr<turbModel> &turb, const int &mm,
                    genArray &residL2, resid &residLinf) {
  // create dummy update (not used in explicit update)
  multiArray3d<genArray> du(1, 1, 1, 0);
  // loop over all blocks and update
  for (auto &block : blocks) {
    block.UpdateBlock(inp, eqnState, thermo, trans, du, turb, mm, residL2,
                      residLinf);
  }
}

double ImplicitUpdate(vector<procBlock> &blocks,
                      vector<multiArray3d<fluxJacobian>> &mainDiagonal,
                      const input &inp, const unique_ptr<eos> &eqnState,
                      const unique_ptr<thermodynamic> &thermo,
                      const unique_ptr<transport> &trans,
                      const unique_ptr<turbModel> &turb, const int &mm,
                      genArray &residL2, resid &residLinf,
                      const vector<connection> &connections, const int &rank,
                      const MPI_Datatype &MPI_cellData) {
  // blocks -- vector of procBlocks on current processor
  // mainDiagonal -- main diagonal of A matrix for all blocks on processor
  // inp -- input variables
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // trans -- viscous transport model
  // turb -- turbulence model
  // mm -- nonlinear iteration
  // residL2 -- L2 residual
  // residLinf -- L infinity residual

  // initialize matrix error
  auto matrixError = 0.0;

  const auto numG = blocks[0].NumGhosts();

  // add volume and time term and calculate inverse of main diagonal
  for (auto bb = 0U; bb < blocks.size(); bb++) {
    blocks[bb].InvertDiagonal(mainDiagonal[bb], inp);
  }

  // initialize matrix update
  vector<multiArray3d<genArray>> du(blocks.size());
  for (auto bb = 0U; bb < blocks.size(); bb++) {
    du[bb] = blocks[bb].InitializeMatrixUpdate(inp, eqnState, thermo,
                                               mainDiagonal[bb]);
  }

  // Solve Ax=b with supported solver
  if (inp.MatrixSolver() == "lusgs" || inp.MatrixSolver() == "blusgs") {
    // calculate order by hyperplanes for each block
    vector<vector<vector3d<int>>> reorder(blocks.size());
    for (auto bb = 0U; bb < blocks.size(); bb++) {
      reorder[bb] = HyperplaneReorder(blocks[bb].NumI(), blocks[bb].NumJ(),
                                      blocks[bb].NumK());
    }

    // start sweeps through domain
    for (auto ii = 0; ii < inp.MatrixSweeps(); ii++) {
      // swap updates for ghost cells
      SwapImplicitUpdate(du, connections, rank, MPI_cellData, numG);

      // forward lu-sgs sweep
      for (auto bb = 0U; bb < blocks.size(); bb++) {
        blocks[bb].LUSGS_Forward(reorder[bb], du[bb], eqnState, inp, thermo,
                                 trans, turb, mainDiagonal[bb], ii);
      }

      // swap updates for ghost cells
      SwapImplicitUpdate(du, connections, rank, MPI_cellData, numG);

      // backward lu-sgs sweep
      for (auto bb = 0U; bb < blocks.size(); bb++) {
        matrixError += blocks[bb].LUSGS_Backward(reorder[bb], du[bb], eqnState,
                                                 inp, thermo, trans, turb,
                                                 mainDiagonal[bb], ii);
      }
    }
  } else if (inp.MatrixSolver() == "dplur" || inp.MatrixSolver() == "bdplur") {
    for (auto ii = 0; ii < inp.MatrixSweeps(); ii++) {
      // swap updates for ghost cells
      SwapImplicitUpdate(du, connections, rank, MPI_cellData, numG);

      for (auto bb = 0U; bb < blocks.size(); bb++) {
        // Calculate correction (du)
        matrixError += blocks[bb].DPLUR(du[bb], eqnState, inp, thermo, trans,
                                        turb, mainDiagonal[bb]);
      }
    }
  } else {
    cerr << "ERROR: Matrix solver " << inp.MatrixSolver() <<
        " is not recognized!" << endl;
    cerr << "Please choose lusgs, blusgs, dplur, or bdplur." << endl;
    exit(EXIT_FAILURE);
  }

  // Update blocks and reset main diagonal
  for (auto bb = 0U; bb < blocks.size(); bb++) {
    // Update solution
    blocks[bb].UpdateBlock(inp, eqnState, thermo, trans, du[bb], turb, mm,
                           residL2, residLinf);

    // Assign time n to time n-1 at end of nonlinear iterations
    if (inp.IsMultilevelInTime() && mm == inp.NonlinearIterations() - 1) {
      blocks[bb].AssignSolToTimeNm1();
    }

    // zero flux jacobians
    mainDiagonal[bb].Zero();
  }

  return matrixError;
}

void SwapImplicitUpdate(vector<multiArray3d<genArray>> &du,
                        const vector<connection> &connections, const int &rank,
                        const MPI_Datatype &MPI_cellData,
                        const int &numGhosts) {
  // du -- implicit update in conservative variables
  // conn -- connection boundary conditions
  // rank -- processor rank
  // MPI_cellData -- datatype to pass primVars or genArray
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      du[conn.LocalBlockFirst()].SwapSlice(conn, du[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      du[conn.LocalBlockFirst()].SwapSliceMPI(conn, rank, MPI_cellData);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      du[conn.LocalBlockSecond()].SwapSliceMPI(conn, rank, MPI_cellData);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}


void SwapTurbVars(vector<procBlock> &states,
                  const vector<connection> &connections, const int &rank,
                  const int &numGhosts) {
  // states -- vector of all procBlocks in the solution domain
  // conn -- connection boundary conditions
  // rank -- processor rank
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      states[conn.LocalBlockFirst()].SwapTurbSlice(
          conn, states[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      states[conn.LocalBlockFirst()].SwapTurbSliceMPI(conn, rank);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      states[conn.LocalBlockSecond()].SwapTurbSliceMPI(conn, rank);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

void SwapEddyViscAndGradients(vector<procBlock> &states,
                              const vector<connection> &connections,
                              const int &rank,
                              const MPI_Datatype &MPI_tensorDouble,
                              const MPI_Datatype &MPI_vec3d,
                              const int &numGhosts) {
  // states -- vector of all procBlocks in the solution domain
  // conn -- connection boundary conditions
  // rank -- processor rank
  // MPI_tensorDouble -- MPI datatype for tensor<double>
  // MPI_vec3d -- MPI datatype for vector3d<double>
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      states[conn.LocalBlockFirst()].SwapEddyViscAndGradientSlice(
          conn, states[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      states[conn.LocalBlockFirst()].SwapEddyViscAndGradientSliceMPI(
          conn, rank, MPI_tensorDouble, MPI_vec3d);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      states[conn.LocalBlockSecond()].SwapEddyViscAndGradientSliceMPI(
          conn, rank, MPI_tensorDouble, MPI_vec3d);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

void SwapWallDist(vector<procBlock> &states,
                  const vector<connection> &connections, const int &rank,
                  const int &numGhosts) {
  // states -- vector of all procBlocks in the solution domain
  // conn -- connection boundary conditions
  // rank -- processor rank
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      states[conn.LocalBlockFirst()].SwapWallDistSlice(
          conn, states[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      states[conn.LocalBlockFirst()].SwapWallDistSliceMPI(conn, rank);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      states[conn.LocalBlockSecond()].SwapWallDistSliceMPI(conn, rank);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
  }
}

void CalcResidual(vector<procBlock> &states,
                  vector<multiArray3d<fluxJacobian>> &mainDiagonal,
                  const unique_ptr<transport> &trans,
                  const unique_ptr<thermodynamic> &thermo,
                  const unique_ptr<eos> &eqnState, const input &inp,
                  const unique_ptr<turbModel> &turb,
                  const vector<connection> &connections, const int &rank,
                  const MPI_Datatype &MPI_tensorDouble,
                  const MPI_Datatype &MPI_vec3d) {
  // states -- vector of all procBlocks on processor
  // mainDiagonal -- main diagonal of A matrix for implicit solve
  // trans -- viscous transport model
  // thermo -- thermodynamic model
  // eqnState -- equation of state
  // inp -- input variables
  // turb -- turbulence model
  // connections -- connection boundary conditions
  // rank -- processor rank
  // MPI_tensorDouble -- MPI datatype for tensor<double>
  // MPI_vec3d -- MPI datatype for vector3d<double>

  for (auto bb = 0U; bb < states.size(); bb++) {
    // calculate residual
    states[bb].CalcResidualNoSource(trans, thermo, eqnState, inp, turb,
                                    mainDiagonal[bb]);
  }
  // swap mut & gradients calculated during residual calculation
  SwapEddyViscAndGradients(states, connections, rank, MPI_tensorDouble,
                           MPI_vec3d, inp.NumberGhostLayers());

  if (inp.IsRANS()) {
    // swap turbulence variables calculated during residual calculation
    SwapTurbVars(states, connections, rank, inp.NumberGhostLayers());

    for (auto bb = 0U; bb < states.size(); bb++) {
      // calculate source terms for residual
      states[bb].CalcSrcTerms(trans, turb, inp, mainDiagonal[bb]);
    }
  }
}

void CalcTimeStep(vector<procBlock> &states, const input &inp) {
  // states -- vector of all procBlocks on processor
  // inp -- input variables

  for (auto &state : states) {
    // calculate time step
    state.CalcBlockTimeStep(inp);
  }
}


// function to reorder the block by hyperplanes
/* A hyperplane is a plane of i+j+k=constant within an individual block. The
LUSGS solver must sweep along these hyperplanes to avoid
calculating a flux jacobian. Ex. The solver must visit all points on hyperplane
1 before visiting any points on hyperplane 2.
*/
vector<vector3d<int>> HyperplaneReorder(const int &imax, const int &jmax,
                                        const int &kmax) {
  // total number of hyperplanes in a given block
  const auto numPlanes = imax + jmax + kmax - 2;
  vector<vector3d<int>> reorder;
  reorder.reserve(imax * jmax * kmax);

  for (auto pp = 0; pp < numPlanes; pp++) {
    for (auto kk = 0; kk < kmax; kk++) {
      for (auto jj = 0; jj < jmax; jj++) {
        for (auto ii = 0; ii < imax; ii++) {
          if (ii + jj + kk == pp) {  // if sum of ii, jj, and kk equals pp than
                                     // point is on hyperplane pp
            reorder.push_back(vector3d<int>(ii, jj, kk));
          }
        }
      }
    }
  }

  return reorder;
}

void ResizeArrays(const vector<procBlock> &states, const input &inp,
                  vector<multiArray3d<fluxJacobian>> &jac) {
  // states -- all states on processor
  // sol -- vector of solutions to be resized
  // jac -- vector of flux jacobians to be resized

  const auto fluxJac = inp.IsBlockMatrix() ?
      fluxJacobian(inp.NumFlowEquations(), inp.NumTurbEquations()) :
      fluxJacobian(1, 1);

  for (auto bb = 0U; bb < states.size(); bb++) {
    jac[bb].ClearResize(states[bb].NumI(), states[bb].NumJ(), states[bb].NumK(),
                        0, fluxJac);
  }
}

vector3d<double> TauNormal(const tensor<double> &velGrad,
                           const vector3d<double> &area, const double &mu,
                           const double &mut,
                           const unique_ptr<transport> &trans) {
  // get 2nd coefficient of viscosity assuming bulk viscosity is 0 (Stoke's)
  const auto lambda = trans->Lambda(mu + mut);

  // wall shear stress
  return lambda * velGrad.Trace() * area + (mu + mut) *
      (velGrad + velGrad.Transpose()).MatMult(area);
}

vector3d<double> TauShear(const tensor<double> &velGrad,
                          const vector3d<double> &area, const double &mu,
                          const double &mut,
                          const unique_ptr<transport> &trans) {
  auto tauN = TauNormal(velGrad, area, mu, mut, trans);
  return tauN - tauN.DotProd(area) * area;
}

// This function calculates the coefficients for three third order accurate
// reconstructions used in a 5th order WENO scheme. This follows the formula
// 2.20 from ICASE report No 97-65 by Shu.
vector<double> LagrangeCoeff(const vector<double> &cellWidth,
                             const unsigned int &degree,
                             const int & rr, const int &ii) {
  // cellWidth -- vector of cell widths in stencil
  // degree -- degree of polynomial to use in reconstruction
  // rr -- number of cells to left of reconstruction location - 1
  // ii -- location of the first upwind cell width in the cellWidth vector

  vector<double> coeffs(degree + 1, 0.0);

  for (auto jj = 0U; jj < coeffs.size(); ++jj) {
    for (auto mm = jj + 1; mm <= degree + 1; ++mm) {
      auto numer = 0.0;
      auto denom = 1.0;
      for (auto ll = 0U; ll <= degree + 1; ++ll) {
        // calculate numerator
        if (ll != mm) {
          auto numProd = 1.0;
          for (auto qq = 0U; qq <= degree + 1; ++qq) {
            if (qq != mm && qq != ll) {
              numProd *= StencilWidth(cellWidth, ii - rr + qq, ii + 1);
            }
          }
          numer += numProd;

          // calculate denominator
          denom *= StencilWidth(cellWidth, ii - rr + ll, ii - rr + mm);
        }
      }
      coeffs[jj] += numer / denom;
    }
    coeffs[jj] *= cellWidth[ii - rr + jj];
  }
  return coeffs;
}

template <typename T>
double StencilWidth(const T &cellWidth, const int &start, const int &end) {
  auto width = 0.0;
  if (end > start) {
    width = std::accumulate(std::begin(cellWidth) + start,
                            std::begin(cellWidth) + end, 0.0);
  } else if (start > end) {  // width is negative
    width = -1.0 * std::accumulate(std::begin(cellWidth) + end,
                                   std::begin(cellWidth) + start, 0.0);
  }
  return width;
}

primVars BetaIntegral(const primVars &deriv1, const primVars &deriv2,
                      const double &dx, const double &x) {
  return (deriv1.Squared() * x + deriv1 * deriv2 * x * x +
          deriv2.Squared() * pow(x, 3.0) / 3.0) * dx +
      deriv2.Squared() * x * pow(dx, 3.0);
}

primVars BetaIntegral(const primVars &deriv1, const primVars &deriv2,
                      const double &dx, const double &xl, const double &xh) {
  return BetaIntegral(deriv1, deriv2, dx, xh) -
      BetaIntegral(deriv1, deriv2, dx, xl);
}

primVars Beta0(const double &x_0, const double &x_1, const double &x_2,
               const primVars &y_0, const primVars &y_1, const primVars &y_2) {
  const auto deriv2nd = Derivative2nd(x_0, x_1, x_2, y_0, y_1, y_2);
  const auto deriv1st = (y_2 - y_1) / (0.5 * (x_2 + x_1)) + 0.5 * x_2 * deriv2nd;

  return BetaIntegral(deriv1st, deriv2nd, x_2, -0.5 * x_2, 0.5 * x_2);
}

primVars Beta1(const double &x_0, const double &x_1, const double &x_2,
               const primVars &y_0, const primVars &y_1, const primVars &y_2) {
  const auto deriv2nd = Derivative2nd(x_0, x_1, x_2, y_0, y_1, y_2);
  const auto deriv1st = (y_2 - y_1) / (0.5 * (x_2 + x_1)) - 0.5 * x_1 * deriv2nd;

  return BetaIntegral(deriv1st, deriv2nd, x_1, -0.5 * x_1, 0.5 * x_1);
}

primVars Beta2(const double &x_0, const double &x_1, const double &x_2,
               const primVars &y_0, const primVars &y_1, const primVars &y_2) {
  const auto deriv2nd = Derivative2nd(x_0, x_1, x_2, y_0, y_1, y_2);
  const auto deriv1st = (y_1 - y_0) / (0.5 * (x_1 + x_0)) - 0.5 * x_0 * deriv2nd;

  return BetaIntegral(deriv1st, deriv2nd, x_0, -0.5 * x_0, 0.5 * x_0);
}

// function to calculate the velocity gradients at a cell face using the Thin
// Shear Layer approximation
tensor<double> CalcVelGradTSL(const primVars &left, const primVars &right,
                              const vector3d<double> &normArea,
                              const double &dist) {
  // left -- left state (primative)
  // right -- right state (primative)
  // normArea -- unit area vector of face
  // dist -- distance between centroid of left cell and right cell

  // calculate velocity derivatives
  const auto velDeriv = (right.Velocity() - left.Velocity()) / dist;

  // populate velocity gradient tensor
  tensor<double> velGrad(
      velDeriv.X() * normArea.X(),
      velDeriv.Y() * normArea.X(),
      velDeriv.Z() * normArea.X(),
      velDeriv.X() * normArea.Y(),
      velDeriv.Y() * normArea.Y(),
      velDeriv.Z() * normArea.Y(),
      velDeriv.X() * normArea.Z(),
      velDeriv.Y() * normArea.Z(),
      velDeriv.Z() * normArea.Z());

  return velGrad;
}

// function to read in cloud of points from file and create a kdtree for 
// nearest neighbor search. File format is space delimited as follows.
//
// numberOfPoints
// species1 species2 ...
// x y z rho u v w p tke omega mf1 mf2 ...
// ...
//
kdtree CalcTreeFromCloud(const string &fname, const input &inp,
                         const unique_ptr<transport> &trans,
                         vector<primVars> &states, vector<string> &species) {
  // fname -- name of file to open
  // inp -- input variables
  // trans -- transport model
  // states -- vector of states read from file
  // species -- species present in file

  // open file
  ifstream inFile(fname, ios::in);
  if (inFile.fail()) {
    cerr << "ERROR: Error in CalcTreeFromCloud(). Input file " << fname
         << " did not open correctly!" << endl;
    exit(EXIT_FAILURE);
  }

  vector<vector3d<double>> points;
  auto count = 0;
  string line = "";
  while (getline(inFile, line)) {
    // remove leading and trailing whitespace and ignore comments
    line = Trim(line);
    if (line.length() > 0) {  // only proceed if line has data
      // split line at variable separator
      auto tokens = Tokenize(line, " ");

      if (count == 0) {  // first line has number of points
        auto numPts = std::stoi(tokens[0]);
        points.resize(numPts);
        states.resize(numPts);
      } else if (count == 1) {  // second line has species
        species = tokens;
        if (species.size() != 1) {
          cerr << "ERROR in CalcTreeFromCloud(), only single species currently "
                  "supported"
               << endl;
          exit(EXIT_FAILURE);
        }
      } else if (tokens.size() != 10 + species.size()) {
        cerr << "ERROR in CalcTreeFromCloud(). Expecting "
             << 10 + species.size() << " data points on line " << count
             << " but only found " << tokens.size() << endl;
        exit(EXIT_FAILURE);
      } else {
        vector3d<double> point;
        point[0] = std::stod(tokens[0]) / inp.LRef();
        point[1] = std::stod(tokens[1]) / inp.LRef();
        point[2] = std::stod(tokens[2]) / inp.LRef();
        points[count - 2] = point;
        auto rho = std::stod(tokens[3]) / inp.RRef();
        auto uVel = std::stod(tokens[4]) / inp.ARef();
        auto vVel = std::stod(tokens[5]) / inp.ARef();
        auto wVel = std::stod(tokens[6]) / inp.ARef();
        auto pressure =
            std::stod(tokens[7]) / (inp.RRef() * inp.ARef() * inp.ARef());
        auto tke = std::stod(tokens[8]) / (inp.ARef() * inp.ARef());
        auto omega = std::stod(tokens[9]) * trans->MuRef() /
                     (inp.RRef() * inp.ARef() * inp.ARef());
        vector<double> massFractions(species.size(), 0.0);
        for (auto ii = 0U; ii < massFractions.size(); ++ii) {
          massFractions[ii] = std::stod(tokens[ii + 10]);
        }
        primVars state(rho, uVel, vVel, wVel, pressure, tke, omega);
        states[count - 2] = state;
      }
    }
    count++;
  }

  // create kd tree
  return kdtree(points);
}
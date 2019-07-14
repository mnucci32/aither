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

#include <iostream>               // cout, cerr, endl
#include <algorithm>              // max, min
#include <vector>
#include <string>
#include <memory>
#include <numeric>
#include "utility.hpp"
#include "procBlock.hpp"
#include "eos.hpp"                 // equation of state
#include "physicsModels.hpp"       // physics models
#include "input.hpp"               // inputVars
#include "varArray.hpp"
#include "turbulence.hpp"
#include "slices.hpp"
#include "matMultiArray3d.hpp"
#include "kdtree.hpp"
#include "resid.hpp"
#include "primitive.hpp"
#include "macros.hpp"

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

  // if a connection border needs to be updated, update
  for (auto ii = 0U; ii < adjEdge1.size(); ii++) {
    if (adjEdge1[ii]) {
      inter.UpdateBorderFirst(ii);
    }
    if (adjEdge2[ii]) {
      inter.UpdateBorderSecond(ii);
    }
  }
}

void SwapGeomSliceMPI(connection &inter, procBlock &blk, const int &tag,
                      const MPI_Datatype &MPI_vec3d,
                      const MPI_Datatype &MPI_vec3dMag) {
  // inter -- connection boundary information
  // blk -- first block involved in connection boundary
  // Get indices for slice coming from block to swap
  auto is = 0, ie = 0;
  auto js = 0, je = 0;
  auto ks = 0, ke = 0;

  const auto rank = blk.Rank();
  if (rank == inter.RankFirst()) {  // local block first in connection
    inter.FirstSliceIndices(is, ie, js, je, ks, ke, blk.NumGhosts());
  } else if (rank == inter.RankSecond()) {  // local block second in connection
    inter.SecondSliceIndices(is, ie, js, je, ks, ke, blk.NumGhosts());
  } else {
    cerr << "ERROR: Error in SwapGeomSliceMPI(). Processor rank does "
            "not match either of connection ranks!" << endl;
    exit(EXIT_FAILURE);
  }

  // get local geomslice to swap
  auto slice = geomSlice(blk, {is, ie}, {js, je}, {ks, ke});

  // swap geomSlices with partner block
  slice.PackSwapUnpackMPI(inter, MPI_vec3d, MPI_vec3dMag, rank, tag);

  // change connections to work with slice and ghosts
  auto interAdj = inter;
  // block to insert into is first in connection
  if (rank == inter.RankFirst()) {
    interAdj.AdjustForSlice(true, blk.NumGhosts());
  } else {  // block to insert into is second in connection, so pass swapped
            // version
    interAdj.AdjustForSlice(false, blk.NumGhosts());
  }

  // insert geomSlice into procBlock
  // return vector determining if any of the 4 edges of the connection need to
  // be updated for a "t" intersection
  const auto adjEdge = blk.PutGeomSlice(slice, interAdj, blk.NumGhosts());

  // if a connection border needs to be updated, update
  for (auto ii = 0U; ii < adjEdge.size(); ++ii) {
    if (adjEdge[ii]) {
      inter.UpdateBorderFirst(ii);
    }
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
            reorder.emplace_back(ii, jj, kk);
          }
        }
      }
    }
  }

  return reorder;
}

void SwapImplicitUpdate(vector<blkMultiArray3d<varArray>> &du,
                        const vector<connection> &connections, const int &rank,
                        const int &numGhosts) {
  // du -- implicit update in conservative variables
  // conn -- connection boundary conditions
  // rank -- processor rank
  // numGhosts -- number of ghost cells

  // loop over all connections and swap connection updates when necessary
  for (auto &conn : connections) {
    if (conn.RankFirst() == rank && conn.RankSecond() == rank) {
      // both sides of connection are on this processor, swap w/o mpi
      du[conn.LocalBlockFirst()].SwapSlice(conn, du[conn.LocalBlockSecond()]);
    } else if (conn.RankFirst() == rank) {
      // rank matches rank of first side of connection, swap over mpi
      du[conn.LocalBlockFirst()].SwapSliceMPI(conn, rank, MPI_DOUBLE);
    } else if (conn.RankSecond() == rank) {
      // rank matches rank of second side of connection, swap over mpi
      du[conn.LocalBlockSecond()].SwapSliceMPI(conn, rank, MPI_DOUBLE);
    }
    // if rank doesn't match either side of connection, then do nothing and
    // move on to the next connection
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

// function to calculate the velocity gradients at a cell face using the Thin
// Shear Layer approximation
tensor<double> CalcVelGradTSL(const primitive &left, const primitive &right,
                              const vector3d<double> &normArea,
                              const double &dist) {
  // left -- left state (primitive)
  // right -- right state (primitive)
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
                         vector<primitive> &states, vector<string> &species) {
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
  std::map<int, int> speciesMap;
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
        states.resize(numPts, {inp.NumEquations(), inp.NumSpecies()});
      } else if (count == 1) {  // second line has species
        species = tokens;
        // check that species are defined
        inp.CheckSpecies(species);
        // create map from species order in file to order in simulation
        for (auto ii = 0U; ii < species.size(); ++ii) {
          speciesMap.insert(std::make_pair(ii, inp.SpeciesIndex(species[ii])));
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
        primitive state(inp.NumEquations(), species.size());
        for (auto ii = 0U; ii < massFractions.size(); ++ii) {
          state[speciesMap[ii]] = rho * massFractions[ii];
        }
        state[state.MomentumXIndex()] = uVel;
        state[state.MomentumYIndex()] = vVel;
        state[state.MomentumZIndex()] = wVel;
        state[state.EnergyIndex()] = pressure;
        if (state.HasTurbulenceData()) {
          state[state.TurbulenceIndex()] = tke;
          state[state.TurbulenceIndex() + 1] = omega;
        }
        states[count - 2] = state;
      }
    }
    count++;
  }

  // create kd tree
  return kdtree(points);
}

void AssertWithMessage(const char *exprStr, bool expr, const char *file, 
                       int line, const char *msg) {
  if (!expr) {
    cerr << "Assert failed: " << msg << endl;
    cerr << "Condition: " << exprStr << endl;
    cerr << "At: " << file << ":" << line << endl;
    exit(EXIT_FAILURE);
  }
}

string GetEnvironmentVariable(const string &var) {
  auto val = getenv(var.c_str());
  return val == NULL ? "" : string(val);
}

double Kronecker(const int &ii, const int &jj) {
  return (ii == jj) ? 1.0 : 0.0;
}

double LinearInterpCoeff(const vector3d<double> &x0, const vector3d<double> &x1,
                         const vector3d<double> &x) {
  const auto dir = (x1 - x0).Normalize();
  const auto dist = x0.Distance(x1);
  return (x - x0).DotProd(dir) / dist;
}

std::array<double, 7> TrilinearInterpCoeff(
    const vector3d<double> &x0, const vector3d<double> &x1,
    const vector3d<double> &x2, const vector3d<double> &x3,
    const vector3d<double> &x4, const vector3d<double> &x5,
    const vector3d<double> &x6, const vector3d<double> &x7,
    const vector3d<double> &x) {
  std::array<double, 7> coeffs;
  // 4 linear interpolations to convert to 2D
  coeffs[0] = LinearInterpCoeff(x0, x4, x);
  const auto x04 = LinearInterp(x0, x4, coeffs[0]);
  coeffs[1] = LinearInterpCoeff(x1, x5, x);
  const auto x15 = LinearInterp(x1, x5, coeffs[1]);
  coeffs[2] = LinearInterpCoeff(x2, x6, x);
  const auto x26 = LinearInterp(x2, x6, coeffs[2]);
  coeffs[3] = LinearInterpCoeff(x3, x7, x);
  const auto x37 = LinearInterp(x3, x7, coeffs[3]);

  // 2 linear interpolations to convert to 1D
  coeffs[4] = LinearInterpCoeff(x04, x15, x);
  const auto x0415 = LinearInterp(x04, x15, coeffs[4]);
  coeffs[5] = LinearInterpCoeff(x26, x37, x);
  const auto x2637 = LinearInterp(x26, x37, coeffs[5]);

  // 1 linear interpolation to complete trilinear interpolation
  coeffs[6] = LinearInterpCoeff(x0415, x2637, x);
  return coeffs;
}

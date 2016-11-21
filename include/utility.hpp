/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef UTILITYHEADERDEF  // only if the macro UTILITYHEADERDEF is not
                            // defined execute these lines of code
#define UTILITYHEADERDEF  // define the macro

#include <vector>                  // vector
#include <memory>
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "multiArray3d.hpp"        // multiArray3d
#include "tensor.hpp"              // tensor
#include "macros.hpp"

using std::vector;
using std::unique_ptr;

// forward class declarations
class procBlock;
class idealGas;
class sutherland;
class input;
class genArray;
class turbModel;
class fluxJacobian;
class kdtree;
class resid;


// function definitions
template <typename T>
inline T FaceReconCentral(const T &, const T &, const vector3d<double> &,
                          const vector3d<double> &, const vector3d<double> &);

tensor<double> VectorGradGG(const vector3d<double> &, const vector3d<double> &,
                            const vector3d<double> &, const vector3d<double> &,
                            const vector3d<double> &, const vector3d<double> &,
                            const vector3d<double> &, const vector3d<double> &,
                            const vector3d<double> &, const vector3d<double> &,
                            const vector3d<double> &, const vector3d<double> &,
                            const double &);

vector3d<double> ScalarGradGG(
    const double &, const double &, const double &, const double &,
    const double &, const double &, const vector3d<double> &,
    const vector3d<double> &, const vector3d<double> &,
    const vector3d<double> &, const vector3d<double> &,
    const vector3d<double> &, const double &);

void SwapGeomSlice(interblock &, procBlock &, procBlock &);

void GetBoundaryConditions(vector<procBlock> &, const input &, const idealGas &,
                           const sutherland &, const unique_ptr<turbModel> &,
                           vector<interblock> &, const int &,
                           const MPI_Datatype &);

vector<vector3d<double>> GetViscousFaceCenters(const vector<procBlock> &);
void CalcWallDistance(vector<procBlock> &, const kdtree &);

vector<multiArray3d<genArray>> GetCopyConsVars(const vector<procBlock> &,
                                               const idealGas &);

void ExplicitUpdate(vector<procBlock> &, const input &, const idealGas &,
                    const double &, const sutherland &,
                    const vector<multiArray3d<genArray>> &,
                    const unique_ptr<turbModel> &, const int &, genArray &,
                    resid &);
double ImplicitUpdate(vector<procBlock> &, vector<multiArray3d<fluxJacobian>> &,
                      const input &, const idealGas &, const double &,
                      const sutherland &,
                      const vector<multiArray3d<genArray>> &,
                      const vector<multiArray3d<genArray>> &,
                      vector<multiArray3d<genArray>> &,
                      const unique_ptr<turbModel> &,
                      const int &, genArray &, resid &,
                      const vector<interblock> &, const int &,
                      const MPI_Datatype &);

void SwapImplicitUpdate(vector<multiArray3d<genArray>> &,
                        const vector<interblock> &, const int &,
                        const MPI_Datatype &, const int &);
void SwapTurbVars(vector<procBlock> &, const vector<interblock> &, const int &,
                  const int &);

void CalcResidual(vector<procBlock> &,
                  vector<multiArray3d<fluxJacobian>> &,
                  const sutherland &, const idealGas &, const input &,
                  const unique_ptr<turbModel> &,
                  const vector<interblock> &, const int &);

void CalcTimeStep(vector<procBlock> &, const input &, const double &);

void GetSolMMinusN(vector<multiArray3d<genArray>> &, const vector<procBlock> &,
                   const vector<multiArray3d<genArray>> &,
                   const idealGas &, const input &, const int &);

// function to reorder block by hyperplanes
vector<vector3d<int>> HyperplaneReorder(const int &, const int &, const int &);

void ResizeArrays(const vector<procBlock> &, const input &,
                  vector<multiArray3d<genArray>> &,
                  vector<multiArray3d<fluxJacobian>> &);

vector3d<double> TauNormal(const tensor<double> &, const vector3d<double> &,
                           const double &, const double &, const sutherland &);
// ---------------------------------------------------------------------------
// inline function definitions

// function to reconstruct cell variables to the face using central
// differences
template <typename T>
T FaceReconCentral(const T &varU, const T &varD, const vector3d<double> &pU,
                   const vector3d<double> &pD, const vector3d<double> &pF) {
  // varU -- variable at the cell center of the upwind cell
  // varD -- variable at the cell center of the downwind cell
  // pU -- position of the cell center of the upwind cell
  // pD -- position of the cell center of the downwind cell
  // pF -- position of the face center of the face on which the reconstruction
  // is happening

  // distance from cell center to cell center
  const auto cen2cen = pU.Distance(pD);
  // distance from upwind cell center to cell face
  const auto up2face = pU.Distance(pF);
  // ratio of distance from upwind cell center to cell face to center to center
  const auto upRatio = up2face / cen2cen;

  // reconstruct with central difference
  return varD * upRatio + varU * (1.0 - upRatio);
}


#endif


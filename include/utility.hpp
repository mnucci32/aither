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
class primVars;

// function definitions
template <typename T>
inline T FaceReconCentral(const T &, const T &, const vector<double> &);
template <typename T>
inline T FaceReconCentral4th(const T &, const T &, const T &, const T&,
                             const vector<double> &);

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
void SwapGradients(vector<procBlock> &, const vector<interblock> &, const int &,
                   const MPI_Datatype &, const MPI_Datatype &, const int &);

void CalcResidual(vector<procBlock> &,
                  vector<multiArray3d<fluxJacobian>> &,
                  const sutherland &, const idealGas &, const input &,
                  const unique_ptr<turbModel> &,
                  const vector<interblock> &, const int &, const MPI_Datatype &,
                  const MPI_Datatype &);

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

vector<double> LagrangeCoeff(const vector<double> &, const unsigned int &,
                             const int &, const int &);
template <typename T>
double StencilWidth(const T &, const int &, const int &);

template <typename T>
T Derivative2nd(const double &, const double &, const double &,
                const T &, const T &, const T &);

primVars Beta0(const double &, const double &, const double &,
               const primVars &, const primVars &, const primVars &);
primVars Beta1(const double &, const double &, const double &,
               const primVars &, const primVars &, const primVars &);
primVars Beta2(const double &, const double &, const double &,
               const primVars &, const primVars &, const primVars &);
primVars BetaIntegral(const primVars &, const primVars &, const double &,
                      const double &);
primVars BetaIntegral(const primVars &, const primVars &, const double &,
                      const double &, const double &);

// ---------------------------------------------------------------------------
// inline function definitions

// function to reconstruct cell variables to the face using central
// differences
template <typename T>
T FaceReconCentral(const T &varU, const T &varD,
                   const vector<double> &cellWidth) {
  // varU -- variable at the cell center of the upwind cell
  // varD -- variable at the cell center of the downwind cell
  // cellWidth -- width of cells in stencil

  // get coefficients
  const auto coeffs = LagrangeCoeff(cellWidth, 1, 0, 0);

  // reconstruct with central difference
  return coeffs[0] * varD + coeffs[1] * varU;
}

// function to reconstruct cell variables to the face using central
// differences (4th order)
template <typename T>
T FaceReconCentral4th(const T &varU2, const T &varU1, const T &varD1,
                      const T &varD2, const vector<double> &cellWidth) {
  // varU2 -- variable at the cell center of the second upwind cell
  // varU1 -- variable at the cell center of the first upwind cell
  // varD1 -- variable at the cell center of the first downwind cell
  // varD2 -- variable at the cell center of the second downwind cell
  // cellWidth -- width of cells in stencil

  // get coefficients
  const auto coeffs = LagrangeCoeff(cellWidth, 3, 1, 1);

  // reconstruct with central difference
  return coeffs[0] * varU2 + coeffs[1] * varU1 + coeffs[2] * varD1 +
      coeffs[3] * varD2;
}

template <typename T>
T Derivative2nd(const double &x_0, const double &x_1, const double &x_2,
                const T &y_0, const T &y_1, const T &y_2) {
  const auto fwdDiff1stOrder = (y_2 - y_1) / (0.5 * (x_2 + x_1));
  const auto bckDiff1stOrder = (y_1 - y_0) / (0.5 * (x_1 + x_0));
  return (fwdDiff1stOrder - bckDiff1stOrder) / (0.25 * (x_2 + x_0) + 0.5 * x_1);
}

#endif


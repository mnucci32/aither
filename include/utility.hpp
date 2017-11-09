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
#include <cmath>
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "multiArray3d.hpp"        // multiArray3d
#include "blkMultiArray3d.hpp"     // blkMultiArray3dd
#include "tensor.hpp"              // tensor
#include "macros.hpp"

using std::vector;
using std::unique_ptr;

// forward class declarations
class procBlock;
class eos;
class transport;
class thermodynamic;
class input;
class residual;
class turbModel;
class matMultiArray3d;
class kdtree;
class resid;
class primitive;
class varArray;

// function definitions
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

void SwapGeomSlice(connection &, procBlock &, procBlock &);

void GetBoundaryConditions(vector<procBlock> &, const input &,
                           const unique_ptr<eos> &,
                           const unique_ptr<thermodynamic> &,
                           const unique_ptr<transport> &,
                           const unique_ptr<turbModel> &, vector<connection> &,
                           const int &);

vector<vector3d<double>> GetViscousFaceCenters(const vector<procBlock> &);
void CalcWallDistance(vector<procBlock> &, const kdtree &);

void AssignSolToTimeN(vector<procBlock> &, const unique_ptr<eos> &,
                      const unique_ptr<thermodynamic> &);
void AssignSolToTimeNm1(vector<procBlock> &);

void ExplicitUpdate(vector<procBlock> &, const input &, const unique_ptr<eos> &,
                    const unique_ptr<thermodynamic> &,
                    const unique_ptr<transport> &,
                    const unique_ptr<turbModel> &, const int &, residual &,
                    resid &);
double ImplicitUpdate(vector<procBlock> &, vector<matMultiArray3d> &,
                      const input &, const unique_ptr<eos> &,
                      const unique_ptr<thermodynamic> &,
                      const unique_ptr<transport> &,
                      const unique_ptr<turbModel> &, const int &, residual &,
                      resid &, const vector<connection> &, const int &);

void SwapImplicitUpdate(vector<blkMultiArray3d<varArray>> &,
                        const vector<connection> &, const int &, const int &);
void SwapTurbVars(vector<procBlock> &, const vector<connection> &, const int &,
                  const int &);
void SwapWallDist(vector<procBlock> &, const vector<connection> &, const int &,
                  const int &);
void SwapEddyViscAndGradients(vector<procBlock> &, const vector<connection> &,
                              const int &, const MPI_Datatype &,
                              const MPI_Datatype &, const int &);

void CalcResidual(vector<procBlock> &, vector<matMultiArray3d> &,
                  const unique_ptr<transport> &,
                  const unique_ptr<thermodynamic> &, const unique_ptr<eos> &,
                  const input &, const unique_ptr<turbModel> &,
                  const vector<connection> &, const int &, const MPI_Datatype &,
                  const MPI_Datatype &);

void CalcTimeStep(vector<procBlock> &, const input &);

// function to reorder block by hyperplanes
vector<vector3d<int>> HyperplaneReorder(const int &, const int &, const int &);

void ResizeArrays(const vector<procBlock> &, const input &,
                  vector<matMultiArray3d> &);

vector3d<double> TauNormal(const tensor<double> &, const vector3d<double> &,
                           const double &, const double &,
                           const unique_ptr<transport> &);
vector3d<double> TauShear(const tensor<double> &, const vector3d<double> &,
                          const double &, const double &,
                          const unique_ptr<transport> &);

vector<double> LagrangeCoeff(const vector<double> &, const unsigned int &,
                             const int &, const int &);
template <typename T>
double StencilWidth(const T &, const int &, const int &);

template <typename T>
auto Derivative2nd(const double &, const double &, const double &,
                const T &, const T &, const T &);

tensor<double> CalcVelGradTSL(const primitive&, const primitive&,
                              const vector3d<double>&, const double&);

kdtree CalcTreeFromCloud(const string &, const input &,
                         const unique_ptr<transport> &, vector<primitive> &,
                         vector<string> &);

string GetEnvironmentVariable(const string &);

double Kronecker(const int &, const int &);

// ---------------------------------------------------------------------------
// inline function definitions
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

template <typename T>
auto Derivative2nd(const double &x_0, const double &x_1, const double &x_2,
                const T &y_0, const T &y_1, const T &y_2) {
  const auto fwdDiff1stOrder = (y_2 - y_1) / (0.5 * (x_2 + x_1));
  const auto bckDiff1stOrder = (y_1 - y_0) / (0.5 * (x_1 + x_0));
  return (fwdDiff1stOrder - bckDiff1stOrder) / (0.25 * (x_2 + x_0) + 0.5 * x_1);
}

// function to implement mathematical sign function
template <typename T>
int Sign (const T &val) {
  return (T(0) < val) - (val < T(0));
}

// function to find the root of a function using Ridder's method
template <typename T1, typename T2>
T1 FindRoot(const T2 &func, T1 x1, T1 x2, const T1 &tol,
            const int &maxIter = 100) {
  // check that x1 and x2 bracket root
  auto f1 = func(x1);
  auto f2 = func(x2);
  if (Sign(f1) == Sign(f2) && Sign(f1) != 0.0) {
    cerr << "ERROR. Root is not within specified interval!" << endl;
    return 0.5 * (x1 + x2);
  }
  
  // intialize variable outside of loop so it can be accessed if routine does
  // not converge
  auto x4 = x1;
  for (auto ii = 0; ii < maxIter; ++ii) {
    auto x3 = 0.5 * (x1 + x2);
    auto f3 = func(x3);
    if (f3 == 0.0) {  // root found
      return x3;
    }
    
    auto denom = sqrt(std::abs(f3 * f3 - f1 * f2));
    if (denom == 0.0) {
      return x3;
    }
    auto fac = Sign(f1 - f2);
    x4 = x3 + (x3 - x1) * (fac * f3) / denom;
    auto f4 = func(x4);

    if (f4 == 0.0) {  // root found
      return x4;
    }

    if (Sign(f4) != Sign(f3)) {  // f4/f3 opposite sign
      x1 = x3;
      f1 = f3;
      x2 = x4;
      f2 = f4;
    } else if (Sign(f4) != Sign(f1)) {  // f4/f1 opposite sign
      x2 = x4;
      f2 = f4;
    } else {
      x1 = x4;
      f1 = f4;
    }

    if (std::abs(x2 - x1) <= tol) {
      return x4;
    }
  }

  cerr << "ERROR: FindRoot did not converge!" << endl;
  return x4;
}

#endif


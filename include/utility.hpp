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
class physics;
class transport;
class input;
class residual;
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
vector<vector3d<double>> GetViscousFaceCenters(const vector<procBlock> &);
void SwapImplicitUpdate(vector<blkMultiArray3d<varArray>> &,
                        const vector<connection> &, const int &, const int &);

// function to reorder block by hyperplanes
vector<vector3d<int>> HyperplaneReorder(const int &, const int &, const int &);

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

template <typename T>
T ConvertCellToNode(const T &cellData, const bool &ignoreEdge = false) {
  T nodeData(cellData.NumINoGhosts() + 1, cellData.NumJNoGhosts() + 1,
             cellData.NumKNoGhosts() + 1, 0, cellData.BlockInfo());
  const auto haveGhosts = cellData.GhostLayers() > 0;
  if (haveGhosts) {
    string dir = "";
    for (auto kk = cellData.PhysStartK() - 1; kk <= cellData.PhysEndK(); ++kk) {
      for (auto jj = cellData.PhysStartJ() - 1; jj <= cellData.PhysEndJ();
           ++jj) {
        for (auto ii = cellData.PhysStartI() - 1; ii <= cellData.PhysEndI();
             ++ii) {
          if (cellData.IsPhysical(ii, jj, kk)) {
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj, kk, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii + 1, jj, kk, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii + 1, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii + 1, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
            }
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii + 1, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
            }
          } else if (!(ignoreEdge && (cellData.AtEdge(ii, jj, kk, dir) ||
                                      cellData.AtCorner(ii, jj, kk)))) {
            // cell data is ghost cell - or -
            // we are ignoring edge and corner ghost cells, and not at 
            // edge/corner
            if (nodeData.IsInRange(ii, jj, kk)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii, jj, kk, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii, jj + 1, kk)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii, jj, kk + 1)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii, jj + 1, kk + 1)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii + 1, jj, kk)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii + 1, jj, kk, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii + 1, jj + 1, kk)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii + 1, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii + 1, jj, kk + 1)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii + 1, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
              }
            }
            if (nodeData.IsInRange(ii + 1, jj + 1, kk + 1)) {
              for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
                nodeData(ii + 1, jj + 1, kk + 1, bb) +=
                    cellData(ii, jj, kk, bb);
              }
            }
          }
        }
      }
    }
  } else {  // no ghost layers in cell data
    for (auto kk = cellData.PhysStartK(); kk < cellData.PhysEndK(); ++kk) {
      for (auto jj = cellData.PhysStartJ(); jj < cellData.PhysEndJ(); ++jj) {
        for (auto ii = cellData.PhysStartI(); ii < cellData.PhysEndI(); ++ii) {
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii, jj, kk, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii + 1, jj, kk, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii + 1, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii + 1, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
          }
          for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
            nodeData(ii + 1, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
          }
        }
      }
    }
  }
  constexpr auto eighth = 1.0 / 8.0;
  if (ignoreEdge) {
    const auto edgeFactor = haveGhosts ? 1.0 / 6.0 : 1.0 / 2.0;
    const auto cornerFactor = haveGhosts ? 1.0 / 4.0 : 1.0;
    string edge = "";
    for (auto kk = nodeData.PhysStartK(); kk < nodeData.PhysEndK(); ++kk) {
      for (auto jj = nodeData.PhysStartJ(); jj < nodeData.PhysEndJ(); ++jj) {
        for (auto ii = nodeData.PhysStartI(); ii < nodeData.PhysEndI(); ++ii) {
          if (nodeData.AtInteriorCorner(ii, jj, kk)) {
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj, kk, bb) *= cornerFactor;
            }
          } else if (nodeData.AtInteriorEdge(ii, jj, kk, edge)) {
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj, kk, bb) *= edgeFactor;
            }
          } else {
            for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
              nodeData(ii, jj, kk, bb) *= eighth;
            }
          }
        }
      }
    }
  } else {
    nodeData *= eighth;
  }
  return nodeData;
}


template <typename T>
T ConvertGradCellToNode(const T &cellData) {
  T nodeData(cellData.NumINoGhosts() + 1, cellData.NumJNoGhosts() + 1,
             cellData.NumKNoGhosts() + 1, 0, cellData.BlockInfo());

  // convert interior data with standard averaging
  for (auto kk = cellData.PhysStartK(); kk < cellData.PhysEndK(); ++kk) {
    for (auto jj = cellData.PhysStartJ(); jj < cellData.PhysEndJ(); ++jj) {
      for (auto ii = cellData.PhysStartI(); ii < cellData.PhysEndI(); ++ii) {
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii, jj, kk, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii + 1, jj, kk, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii + 1, jj + 1, kk, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii + 1, jj + 1, kk + 1, bb) += cellData(ii, jj, kk, bb);
        }
        for (auto bb = 0; bb < nodeData.BlockSize(); ++bb) {
          nodeData(ii + 1, jj, kk + 1, bb) += cellData(ii, jj, kk, bb);
        }
      }
    }
  }
  constexpr auto eighth = 1.0 / 8.0;
  nodeData *= eighth;

  // recalculate gradients on the boundary and reassign

  return nodeData;
}

template <typename T>
T TrilinearInterpolation(const vector3d<double> &il, const vector3d<double> &iu,
                         const vector3d<double> &jl, const vector3d<double> &ju,
                         const vector3d<double> &kl, const vector3d<double> &ku,
                         const vector3d<double> &xc,
                         const T &d1, const T &d2, const T &d3, const T &d4,
                         const T &d5, const T &d6, const T &d7, const T &d8,
                         const vector3d<double> &x) {
  // 1 - bottom, left, aft
  // 2 - bottom, left, fore
  // 3 - bottom, right, aft
  // 4 - bottom, right, fore
  // 5 - top, left, aft
  // 6 - top, left, fore
  // 7 - top, right, aft
  // 8 - top, right, fore

  constexpr auto eighth = 1.0 / 8.0;
  const auto intPt1 = -1.0 / sqrt(3.0);
  const auto intPt2 = -intPt1;

  // convert x to isoparametric coordinates
  auto iNorm = 2.0 * (iu - il).Normalize() - 1.0;
  auto jNorm = 2.0 * (ju - jl).Normalize() - 1.0;
  auto kNorm = 2.0 * (ku - kl).Normalize() - 1.0;
  // x = xc + iNorm * zeta + jNorm * eta + kNorm * mu
  // zeta = (x - xc) / iNorm
  auto zeta = -1.0;
  auto eta = -1.0;
  auto mu = -1.0;

  auto shape =
      (1.0 + zeta * intPt1) * (1.0 + eta * intPt1) * (1.0 + mu * intPt1) * d1 +
      (1.0 + zeta * intPt1) * (1.0 + eta * intPt1) * (1.0 + mu * intPt2) * d2 +
      (1.0 + zeta * intPt1) * (1.0 + eta * intPt2) * (1.0 + mu * intPt1) * d3 +
      (1.0 + zeta * intPt1) * (1.0 + eta * intPt2) * (1.0 + mu * intPt2) * d4 +
      (1.0 + zeta * intPt2) * (1.0 + eta * intPt1) * (1.0 + mu * intPt1) * d5 +
      (1.0 + zeta * intPt2) * (1.0 + eta * intPt1) * (1.0 + mu * intPt2) * d6 +
      (1.0 + zeta * intPt2) * (1.0 + eta * intPt2) * (1.0 + mu * intPt1) * d7 +
      (1.0 + zeta * intPt2) * (1.0 + eta * intPt2) * (1.0 + mu * intPt2) * d8;
  shape *= eighth;

  return shape;
}

#endif


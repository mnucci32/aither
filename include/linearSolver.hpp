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

#ifndef LINEAR_SOLVER_HEADER_DEF
#define LINEAR_SOLVER_HEADER_DEF

#include <vector>                  // vector
#include <string>                  // string
#include "matMultiArray3d.hpp"
#include "macros.hpp"

using std::string;
using std::vector;

// forward class declarations
class procBlock;
class input;
class gridLevel;

// class to solve system of linear equations Ax=b
class linearSolver {
  string solverType_;

 public:
  // constructors
  linearSolver(const string &type) : solverType_(type) {
    MSG_ASSERT(type == "lusgs" || type == "blusgs" || type == "dplur" ||
                   type == "bdplur",
               "linear solver type not supported");
  }

  // move constructor and assignment operator
  linearSolver(linearSolver&&) noexcept = default;
  linearSolver& operator=(linearSolver&&) noexcept = default;

  // copy constructor and assignment operator
  linearSolver(const linearSolver&) = default;
  linearSolver& operator=(const linearSolver&) = default;

  // member functions
  void InvertDiagonal(const procBlock &, const input &,
                      matMultiArray3d &) const;
  blkMultiArray3d<varArray> InitializeMatrixUpdate(
      const procBlock &, const input &, const physics &,
      const matMultiArray3d &) const;

  void LUSGS_Forward(const procBlock &, const vector<vector3d<int>> &,
                     const physics &, const input &, const matMultiArray3d &,
                     const int &, blkMultiArray3d<varArray> &) const;
  double LUSGS_Backward(const procBlock &, const vector<vector3d<int>> &,
                        const physics &, const input &, const matMultiArray3d &,
                        const int &, blkMultiArray3d<varArray> &) const;
  double LUSGS_Relax(const gridLevel &, const physics &, const input &,
                     const int &, const int &,
                     vector<blkMultiArray3d<varArray>> &) const;

  double DPLUR(const procBlock &, const physics &, const input &,
               const matMultiArray3d &, blkMultiArray3d<varArray> &) const;
  double DPLUR_Relax(const gridLevel &, const physics &, const input &,
                     const int &, const int &,
                     vector<blkMultiArray3d<varArray>> &) const;

  // destructor
  ~linearSolver() noexcept {}
};

// ----------------------------------------------------------------------------
// function definitions

#endif









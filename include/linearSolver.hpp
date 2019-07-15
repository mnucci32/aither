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

#ifndef LINEAR_SOLVER_HEADER_DEF
#define LINEAR_SOLVER_HEADER_DEF

#include <vector>                  // vector
#include <string>                  // string
#include "matMultiArray3d.hpp"
#include "blkMultiArray3d.hpp"
#include "macros.hpp"

using std::string;
using std::vector;

// forward class declarations
class procBlock;
class input;
class gridLevel;

// classes to solve system of linear equations Ax=b

// abstract base class
class linearSolver {
  string solverType_;
  vector<matMultiArray3d> a_;
  vector<matMultiArray3d> aInv_;
 protected:
  vector<blkMultiArray3d<varArray>> x_;

 public:
  // constructors
  linearSolver(const input &inp, const gridLevel &level);

  // move constructor and assignment operator
  linearSolver(linearSolver&&) noexcept = default;
  linearSolver& operator=(linearSolver&&) noexcept = default;

  // copy constructor and assignment operator
  linearSolver(const linearSolver&) = default;
  linearSolver& operator=(const linearSolver&) = default;

  // member functions
  int NumBlocks() const { return a_.size(); }
  const vector<blkMultiArray3d<varArray>> &X() const { return x_; }
  const blkMultiArray3d<varArray> &X(const int &bb) const { return x_[bb]; }

  const vector<matMultiArray3d> &A() const { return a_; }
  const matMultiArray3d &A(const int &bb) const { return a_[bb]; }
  matMultiArray3d &A(const int &bb) { return a_[bb]; }

  const matMultiArray3d &AInv(const int &bb) const { return aInv_[bb]; }

  vector<blkMultiArray3d<varArray>> AXmB(const gridLevel &, const physics &,
                                       const input &) const;
  vector<blkMultiArray3d<varArray>> Residual(const gridLevel &, const physics &,
                                             const input &) const;

  void AddDiagonalTerms(const gridLevel &, const input &);
  void Invert();
  void InitializeMatrixUpdate(const gridLevel &, const input &,
                              const physics &);
  void SwapUpdate(const vector<connection> &, const int &, const int &);
  void SubtractFromUpdate(const vector<blkMultiArray3d<varArray>>& coarseDu);
  void AddToUpdate(const vector<blkMultiArray3d<varArray>>& correction);
  void ZeroA(const int &bb) { a_[bb].Zero(); }
  void Restriction(unique_ptr<linearSolver> &coarse,
                   const vector<connection> &conn,
                   const vector<multiArray3d<vector3d<int>>> &toCoarse,
                   const vector<multiArray3d<double>> &volWeightFactor,
                   const int &rank) const;

  virtual vector<blkMultiArray3d<varArray>> Relax(const gridLevel &,
                                                  const physics &,
                                                  const input &, const int &,
                                                  const int &) = 0;

  // destructor
  virtual ~linearSolver() noexcept {}
};

// --------------------------------------------------------------------------
class lusgs : public linearSolver {
  vector<vector<vector3d<int>>> reorder_;

  // private member functions
  void LUSGS_Forward(const procBlock &, const vector<vector3d<int>> &,
                     const physics &, const input &, const matMultiArray3d &,
                     const int &, const blkMultiArray3d<varArray> &,
                     blkMultiArray3d<varArray> &) const;
  void LUSGS_Backward(const procBlock &, const vector<vector3d<int>> &,
                      const physics &, const input &, const matMultiArray3d &,
                      const matMultiArray3d &, const int &,
                      const blkMultiArray3d<varArray> &,
                      blkMultiArray3d<varArray> &) const;

 public:
  // constructors
  lusgs(const input &inp, const gridLevel &level);

  // move constructor and assignment operator
  lusgs(lusgs &&solver) noexcept : linearSolver(std::move(solver)) {}
  lusgs &operator=(lusgs &&) noexcept = default;

  // copy constructor and assignment operator
  lusgs(const lusgs &solver) : linearSolver(solver) {}
  lusgs &operator=(const lusgs &) = default;

  // member functions
  vector<blkMultiArray3d<varArray>> Relax(const gridLevel &, const physics &,
                                          const input &, const int &,
                                          const int &) override;

  // destructor
  virtual ~lusgs() noexcept {}
};

// --------------------------------------------------------------------------
class dplur : public linearSolver {

  // private member functions
  void DPLUR(const procBlock &, const physics &, const input &,
             const matMultiArray3d &, const matMultiArray3d &,
             const blkMultiArray3d<varArray> &,
             blkMultiArray3d<varArray> &) const;

 public:
  // constructors
  dplur(const input &inp, const gridLevel &level) : linearSolver(inp, level) {}

  // move constructor and assignment operator
  dplur(dplur &&solver) noexcept : linearSolver(std::move(solver)) {}
  dplur &operator=(dplur &&) noexcept = default;

  // copy constructor and assignment operator
  dplur(const dplur &solver) : linearSolver(solver) {}
  dplur &operator=(const dplur &) = default;

  // member functions
  vector<blkMultiArray3d<varArray>> Relax(const gridLevel &, const physics &,
                                          const input &, const int &,
                                          const int &) override;

  // destructor
  virtual ~dplur() noexcept {}
};

// ----------------------------------------------------------------------------
// function definitions

#endif


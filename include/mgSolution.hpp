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

#ifndef MGSOLUTION_HEADER_DEF
#define MGSOLUTION_HEADER_DEF

#include <vector>
#include "gridLevel.hpp"
#include "vector3d.hpp"
#include "multiArray3d.hpp"

using std::vector;

// forward class declaration
class plot3dBlock;
class boundaryConditions;
class input;
class decomposition;
class input;
class physics;
class residual;
class kdtree;

class mgSolution {
  vector<gridLevel> solution_;
  int mgCycleIndex_;

  // private member functions
  void Restriction(const int&);
  template<typename T>
  void Prolongation(const int&, const vector<T>&);
  template <typename T>
  void CycleAtLevel(const int&, const input &, T &);
  template <typename T>
  void Relax(const int&, T&) const;

 public:
  // Constructor
  mgSolution(const int& numLevels, const int& cycle);
  mgSolution() : mgSolution(1, 1) {}
  mgSolution(const input& inp)
      : mgSolution(inp.MultigridLevels(), inp.MultigridCycleIndex()) {}

  // move constructor and assignment operator
  mgSolution(mgSolution&&) noexcept = default;
  mgSolution& operator=(mgSolution&&) noexcept = default;

  // copy constructor and assignment operator
  mgSolution(const mgSolution&) = default;
  mgSolution& operator=(const mgSolution&) = default;

  // Member functions
  int NumGridLevels() const { return solution_.size(); }

  const gridLevel& operator[](const int &a) const { return solution_[a]; }
  gridLevel& operator[](const int &a) { return solution_[a]; }
  const gridLevel& Finest() const { return solution_.front(); }
  const gridLevel& Coarsest() const { return solution_.back(); }

  void ConstructFinestLevel(const vector<plot3dBlock>& mesh,
                            const vector<boundaryConditions>& bcs,
                            const decomposition& decomp, const physics& phys,
                            const string& restartFile, input& inp,
                            residual& first);
  void ConstructMultigrids(const decomposition& decomp, const input& inp,
                           const physics& phys);
  mgSolution SendFinestGridLevel(const int& rank, const int& numProcBlock,
                                 const MPI_Datatype& MPI_vec3d,
                                 const MPI_Datatype& MPI_vec3dMag,
                                 const MPI_Datatype& MPI_connection,
                                 const input& inp) const;
  void GetFinestGridLevel(const mgSolution& local, const int& rank,
                          const MPI_Datatype& MPI_uncoupledScalar,
                          const MPI_Datatype& MPI_vec3d,
                          const MPI_Datatype& MPI_tensorDouble,
                          const input& inp);
  void AuxillaryAndWidths(const physics& phys);
  void CalcWallDistance(const kdtree& tree);
  void SwapWallDist(const int& rank, const int& numGhosts);
  void ResizeMatrix(const input& inp, const int& numProcBlock);
  void FullMultigridCycle();

  // Destructor
  ~mgSolution() noexcept {}
};

// multigrid prolongation - coarse grid to fine grid operator
template <typename T>
void mgSolution::Prolongation(const int &ci, const vector<T> &correction) {
  MSG_ASSERT(ci > 0 && ci < static_cast<int>(solution_.size()),
             "index for prolongation out of range");
  solution_[ci].Prolongation(correction, solution_[ci - 1]);
}

// relax solution
template <typename T>
void mgSolution::Relax(const int& sweeps, T& data) const {
  
}

template <typename T>
void mgSolution::CycleAtLevel(const int &fl, const input &inp, T &sol) {
  // fl -- index for fine grid level
  MSG_ASSERT(fl >= 0 && fl < static_cast<int>(solution_.size()),
             "index for multigrid cycle out of range");

  if (fl == this->NumGridLevels() - 1) {  // recursive base case

  } else {
    // pre-relaxation sweeps


    // coarse grid correction
    auto cl = fl + 1;
    // calc residual, restrict to forcing term of coarse grid
    this->Restriction(fl);

    // recursive call to next coarse level
    vector<blkMultiArray3d<varArray>> correction;
    correction.reserve(solution_[cl].NumBlocks());
    for (const auto &blk : solution_[cl].Blocks()) {
      correction.emplace_back(blk.NumI(), blk.NumJ(), blk.NumK(), 0,
                              blk.NumEquations(), blk.NumSpecies());
    }
    for (auto ii = 0; ii < mgCycleIndex_; ++ii) {
      this->CycleAtLevel(cl, inp, correction);
    }

    // interpolate coarse level correction and add to solution
    this->Prolongation(cl, correction);

    // post-relaxation sweeps

  }

}


#endif

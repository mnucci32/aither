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

#ifndef MGSOLUTION_HEADER_DEF
#define MGSOLUTION_HEADER_DEF

#include <vector>
#include "gridLevel.hpp"
#include "vector3d.hpp"
#include "multiArray3d.hpp"
#include "blkMultiArray3d.hpp"
#include "physicsModels.hpp"
#include "linearSolver.hpp"

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
  double ImplicitUpdate(const input& inp, const physics& phys, const int& mm,
                        const int& rank, const MPI_Datatype& MPI_tensorDouble,
                        const MPI_Datatype& MPI_vec3d, residual& residL2,
                        resid& residLinf);
  void Restriction(const int&, const int&,
                   const vector<blkMultiArray3d<varArray>>&, const input& inp,
                   const physics& phys, const int& rank,
                   const MPI_Datatype& MPI_tensorDouble,
                   const MPI_Datatype& MPI_vec3d);
  void Prolongation(const int&);
  double CycleAtLevel(const int&, const int&, const physics&, const input&,
                      const int&, const MPI_Datatype&, const MPI_Datatype&);
  vector<blkMultiArray3d<varArray>> Relax(const int&, const int&,
                                          const physics&, const input&,
                                          const int&);

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
  int NumBlocks() const { return this->Finest().NumBlocks(); }

  const gridLevel& operator[](const int &a) const { return solution_[a]; }
  gridLevel& operator[](const int &a) { return solution_[a]; }
  const gridLevel& Finest() const { return solution_.front(); }
  int FinestIndex() const { return 0; }
  const gridLevel& Coarsest() const { return solution_.back(); }

  void ConstructFinestLevel(const vector<plot3dBlock>& mesh,
                            const vector<boundaryConditions>& bcs,
                            const decomposition& decomp, const physics& phys,
                            const string& restartFile, input& inp,
                            residual& first);
  void ConstructMultigrids(const decomposition& decomp, const input& inp,
                           const physics& phys, const int &rank,
                           const MPI_Datatype& MPI_connection,
                           const MPI_Datatype& MPI_vec3d,
                           const MPI_Datatype& MPI_vec3dMag);
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
  void StoreOldSolution(const input& inp, const physics& phys, const int &iter);
  void CalcWallDistance(const kdtree& tree);
  void SwapWallDist(const int& rank, const int& numGhosts);
  void SubtractFromUpdate(const int& ll,
                          const vector<blkMultiArray3d<varArray>>& coarseDu);
  double Iterate(const input& inp, const physics& phys,
                 const MPI_Datatype& MPI_tensorDouble,
                 const MPI_Datatype& MPI_vec3d, const int& mm, const int& rank,
                 residual& residL2, resid& residLinf);

  // Destructor
  ~mgSolution() noexcept {}
};


#endif

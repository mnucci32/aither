/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef PROCBLOCKHEADERDEF  // only if the macro PROCBLOCKHEADERDEF is not
                            // defined execute these lines of code
#define PROCBLOCKHEADERDEF  // define the macro

#include <vector>                  // vector
#include <string>                  // string
#include <fstream>
#include <iostream>
#include <memory>
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "multiArray3d.hpp"        // multiArray3d
#include "tensor.hpp"              // tensor
#include "primVars.hpp"            // primVars
#include "genArray.hpp"            // genArray
#include "boundaryConditions.hpp"  // interblock, patch
#include "macros.hpp"
#include "kdtree.hpp"              // kdtree
#include "uncoupledScalar.hpp"     // uncoupledScalar

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;

// forward class declarations
class idealGas;
class sutherland;
class inviscidFlux;
class viscousFlux;
class input;
class gradients;
class geomSlice;
class source;
class turbModel;
class plot3dBlock;
class resid;
class fluxJacobian;

class procBlock {
  multiArray3d<primVars> state_;  // primative variables at cell center

  multiArray3d<genArray> residual_;  // cell residual

  multiArray3d<unitVec3dMag<double>> fAreaI_;  // face area vector for i-faces
  multiArray3d<unitVec3dMag<double>> fAreaJ_;  // face area vector for j-faces
  multiArray3d<unitVec3dMag<double>> fAreaK_;  // face area vector for k-faces

  multiArray3d<vector3d<double>> center_;  // coordinates of cell center
  multiArray3d<vector3d<double>> fCenterI_;  // coordinates of i-face centers
  multiArray3d<vector3d<double>> fCenterJ_;  // coordinates of j-face centers
  multiArray3d<vector3d<double>> fCenterK_;  // coordinates of k-face centers

  multiArray3d<uncoupledScalar> specRadius_;  // maximum wave speed for cell
  multiArray3d<double> vol_;  // cell volume
  multiArray3d<double> dt_;  // cell time step
  multiArray3d<double> wallDist_;  // distance to nearest viscous wall

  boundaryConditions bc_;  // boundary conditions for block

  int numGhosts_;  // number of layers of ghost cells surrounding block
  int parBlock_;  // parent block number
  int rank_;  // processor rank
  int localPos_;  // position on local processor
  int globalPos_;  // global position of procBlock in decomposed vector of
                   // procBlocks

  // private member functions
  void CalcInvFluxI(const idealGas &, const input &,
                    const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);
  void CalcInvFluxJ(const idealGas &, const input &,
                    const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);
  void CalcInvFluxK(const idealGas &, const input &,
                    const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);

  void CalcViscFluxI(const sutherland &, const idealGas &, const input &,
                     const gradients &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);
  void CalcViscFluxJ(const sutherland &, const idealGas &, const input &,
                     const gradients &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);
  void CalcViscFluxK(const sutherland &, const idealGas &, const input &,
                     const gradients &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);

  void CalcCellDt(const int &, const int &, const int &, const double &);

  void ExplicitEulerTimeAdvance(const idealGas &, const unique_ptr<turbModel> &,
                                const int &, const int &, const int &,
                                const int &, const int &, const int &);
  void ImplicitTimeAdvance(const genArray &, const idealGas &,
                           const unique_ptr<turbModel> &, const int &,
                           const int &, const int &);
  void RK4TimeAdvance(const genArray &, const idealGas &,
                      const unique_ptr<turbModel> &, const int &, const int &,
                      const int &, const int &, const int &, const int &,
                      const int &);

  void CalcSrcTerms(const gradients &, const sutherland &, const idealGas &,
                    const unique_ptr<turbModel> &, const input &,
                    multiArray3d<fluxJacobian> &);

  void AddToResidual(const inviscidFlux &, const int &, const int &,
                     const int &);
  void AddToResidual(const viscousFlux &, const int &, const int &,
                     const int &);
  void SubtractFromResidual(const inviscidFlux &, const int &, const int &,
                            const int &);
  void SubtractFromResidual(const viscousFlux &, const int &, const int &,
                            const int &);
  void SubtractFromResidual(const source &, const int &, const int &,
                            const int &);

 public:
  // constructors
  procBlock(const primVars &, const plot3dBlock &, const int &, const int &,
            const boundaryConditions &, const int &, const int &, const int &);
  procBlock(const int &, const int &, const int &, const int &);
  procBlock() : procBlock(1, 1, 1, 0) {}

  // move constructor and assignment operator
  procBlock(procBlock&&) noexcept = default;
  procBlock& operator=(procBlock&&) noexcept = default;

  // copy constructor and assignment operator
  procBlock(const procBlock&) = default;
  procBlock& operator=(const procBlock&) = default;

  // member functions
  int NumCells() const { return residual_.Size(); }
  int NumCellsGhosts() const { return state_.Size(); }
  int NumI() const { return residual_.NumI(); }
  int NumJ() const { return residual_.NumJ(); }
  int NumK() const { return residual_.NumK(); }
  int NumIG() const { return state_.NumI(); }
  int NumJG() const { return state_.NumJ(); }
  int NumKG() const { return state_.NumK(); }
  int NumGhosts() const { return numGhosts_; }
  int ParentBlock() const { return parBlock_; }
  int LocalPosition() const { return localPos_; }
  int Rank() const { return rank_; }
  int GlobalPos() const { return globalPos_; }

  boundaryConditions BC() const { return bc_; }

  primVars State(const int &ii, const int &jj, const int &kk) const {
    return state_(ii, jj, kk);
  }

  multiArray3d<primVars> SliceState(const int &, const int &, const int &,
                                    const int &, const int &,
                                    const int &) const;

  multiArray3d<genArray> GetCopyConsVars(const idealGas &) const;

  double Vol(const int &ii, const int &jj, const int &kk) const {
    return vol_(ii, jj, kk);
  }
  vector3d<double> Center(const int &ii, const int &jj, const int &kk) const {
    return center_(ii, jj, kk);
  }
  vector3d<double> FAreaUnitI(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaI_(ii, jj, kk).UnitVector();
  }
  vector3d<double> FAreaUnitJ(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaJ_(ii, jj, kk).UnitVector();
  }
  vector3d<double> FAreaUnitK(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaK_(ii, jj, kk).UnitVector();
  }
  double FAreaMagI(const int &ii, const int &jj, const int &kk) const {
    return fAreaI_(ii, jj, kk).Mag();
  }
  double FAreaMagJ(const int &ii, const int &jj, const int &kk) const {
    return fAreaJ_(ii, jj, kk).Mag();
  }
  double FAreaMagK(const int &ii, const int &jj, const int &kk) const {
    return fAreaK_(ii, jj, kk).Mag();
  }
  unitVec3dMag<double> FAreaI(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaI_(ii, jj, kk);
  }
  unitVec3dMag<double> FAreaJ(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaJ_(ii, jj, kk);
  }
  unitVec3dMag<double> FAreaK(const int &ii, const int &jj,
                              const int &kk) const {
    return fAreaK_(ii, jj, kk);
  }

  vector3d<double> FCenterI(const int &ii, const int &jj, const int &kk) const {
    return fCenterI_(ii, jj, kk);
  }
  vector3d<double> FCenterJ(const int &ii, const int &jj, const int &kk) const {
    return fCenterJ_(ii, jj, kk);
  }
  vector3d<double> FCenterK(const int &ii, const int &jj, const int &kk) const {
    return fCenterK_(ii, jj, kk);
  }

  uncoupledScalar SpectralRadius(const int &ii, const int &jj, const int &kk) const {
    return specRadius_(ii, jj, kk);
  }
  double Dt(const int &ii, const int &jj, const int &kk) const {
    return dt_(ii, jj, kk);
  }
  double WallDist(const int &ii, const int &jj, const int &kk) const {
    return wallDist_(ii, jj, kk);
  }

  genArray Residual(const int &ii, const int &jj, const int &kk) const {
    return residual_(ii, jj, kk);
  }
  double Residual(const int &ii, const int &jj, const int &kk,
                  const int &a) const {
    return residual_(ii, jj, kk)[a];
  }

  void CalcBlockTimeStep(const input &, const double &);
  void UpdateBlock(const input &, const idealGas &, const double &,
                   const sutherland &, const multiArray3d<genArray> &,
                   const multiArray3d<genArray> &,
                   const unique_ptr<turbModel> &, const int &, genArray &,
                   resid &);

  void CalcResidual(const sutherland &, const idealGas &, const input &,
                    const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);

  void ResetResidWS();
  void CleanResizeVecs(const int &, const int &, const int &);

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input &, const idealGas &,
                                const sutherland &,
                                const unique_ptr<turbModel> &);
  void AssignInviscidGhostCellsEdge(const input &, const idealGas &,
                                    const sutherland &,
                                    const unique_ptr<turbModel> &);

  void AssignViscousGhostCells(const input &, const idealGas &,
                               const sutherland &,
                               const unique_ptr<turbModel> &);
  void AssignViscousGhostCellsEdge(const input &, const idealGas &,
                                   const sutherland &,
                                   const unique_ptr<turbModel> &);

  void CalcGradsI(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsJ(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsK(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;

  void CalcWallDistance(const kdtree &);

  multiArray3d<genArray> DeltaNMinusOne(const multiArray3d<genArray> &,
                                        const idealGas &, const double &,
                                        const double &) const;
  multiArray3d<genArray> SolTimeMMinusN(const multiArray3d<genArray> &,
                                        const idealGas &, const input &,
                                        const int &) const;
  void InvertDiagonal(multiArray3d<fluxJacobian> &, const input &) const;

  multiArray3d<genArray> InitializeMatrixUpdate(const input &,
						const multiArray3d<genArray> &,
						const multiArray3d<genArray> &,
						const multiArray3d<fluxJacobian> &) const;
  void LUSGS_Forward(const vector<vector3d<int>> &, multiArray3d<genArray> &,
		     const multiArray3d<genArray> &, const multiArray3d<genArray> &,
		     const idealGas &, const input &, const sutherland &,
		     const unique_ptr<turbModel> &,
		     const multiArray3d<fluxJacobian> &, const int &) const;
  double LUSGS_Backward(const vector<vector3d<int>> &, multiArray3d<genArray> &,
			const multiArray3d<genArray> &, const multiArray3d<genArray> &,
			const idealGas &, const input &, const sutherland &,
			const unique_ptr<turbModel> &,
			const multiArray3d<fluxJacobian> &, const int &) const;

  double DPLUR(multiArray3d<genArray> &,
               const multiArray3d<genArray> &, const multiArray3d<genArray> &,
               const idealGas &, const input &, const sutherland &,
               const unique_ptr<turbModel> &,
               const multiArray3d<fluxJacobian> &) const;

  bool IsPhysical(const int &, const int &, const int &, const bool &) const;
  bool AtCorner(const int &, const int &, const int &, const bool &) const;
  bool AtEdge(const int &, const int &, const int &, const bool &,
              string &) const;

  vector<bool> PutGeomSlice(const geomSlice &, interblock &, const int &,
                            const int &);
  void PutStateSlice(const multiArray3d<primVars> &, const interblock &,
                     const int &, const int &);

  procBlock Split(const string &, const int &, const int &,
                  vector<boundarySurface> &);
  void Join(const procBlock &, const string &, vector<boundarySurface> &);

  void SwapStateSlice(const interblock &, procBlock &);
  void SwapStateSliceMPI(const interblock &, const int &, const MPI_Datatype &);

  void PackSendGeomMPI(const MPI_Datatype &, const MPI_Datatype &,
                       const MPI_Datatype &) const;
  void RecvUnpackGeomMPI(const MPI_Datatype &, const MPI_Datatype &,
                         const MPI_Datatype &);
  void PackSendSolMPI(const MPI_Datatype &, const MPI_Datatype &) const;
  void RecvUnpackSolMPI(const MPI_Datatype &, const MPI_Datatype &);

  // destructor
  ~procBlock() noexcept {}
};

// function definitions
template <typename T>
T FaceReconCentral(const T &, const T &, const vector3d<double> &,
                   const vector3d<double> &, const vector3d<double> &);

template <typename T>
multiArray3d<T> PadWithGhosts(const multiArray3d<T> &, const int &);

tensor<double> CalcVelGradGG(const vector3d<double> &, const vector3d<double> &,
                             const vector3d<double> &, const vector3d<double> &,
                             const vector3d<double> &, const vector3d<double> &,
                             const vector3d<double> &, const vector3d<double> &,
                             const vector3d<double> &, const vector3d<double> &,
                             const vector3d<double> &, const vector3d<double> &,
                             const double &);

vector3d<double> CalcScalarGradGG(
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
                      const input &, const idealGas &, const double &, const sutherland &,
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
#endif









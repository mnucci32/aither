/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
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
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "tensor.hpp"              // tensor
#include "primVars.hpp"            // primVars
#include "matrix.hpp"              // genArray
#include "boundaryConditions.hpp"  // interblock, patch
#include "macros.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// forward class declarations
class idealGas;
class sutherland;
class inviscidFlux;
class viscousFlux;
class input;
class gradients;
class geomSlice;
class stateSlice;
class source;
class turbModel;
class plot3dBlock;
class resid;

class procBlock {
  vector<primVars> state_;  // primative variables at cell center

  vector<vector3d<double> > center_;  // coordinates of cell center
  vector<vector3d<double> > fAreaI_;  // face area vector for i-faces
  vector<vector3d<double> > fAreaJ_;  // face area vector for j-faces
  vector<vector3d<double> > fAreaK_;  // face area vector for k-faces
  vector<vector3d<double> > fCenterI_;  // coordinates of i-face centers
  vector<vector3d<double> > fCenterJ_;  // coordinates of j-face centers
  vector<vector3d<double> > fCenterK_;  // coordinates of k-face centers

  vector<genArray> residual_;  // cell residual

  vector<double> vol_;  // cell volume
  vector<double> avgWaveSpeed_;  // maximum wave speed for cell
  vector<double> dt_;  // cell time step

  boundaryConditions bc_;  // boundary conditions for block

  int numCells_;  // number of cells in block
  int numVars_;  // number of variables stored at cell
  int numI_;  // i-dimension of block (cells)
  int numJ_;  // j-dimension of block (cells)
  int numK_;  // k-dimension of block (cells)
  int numGhosts_;  // number of layers of ghost cells surrounding block
  int parBlock_;  // parent block number
  int rank_;  // processor rank_
  int localPos_;  // position on local processor
  int globalPos_;  // global position of procBlock in decomposed vector of
                   // procBlocks

 public:
  // constructors
  procBlock();
  procBlock(const plot3dBlock &, const int &, const int &);
  procBlock(const double, const double, const vector3d<double>,
            const plot3dBlock &, const int &, const int &,
            const boundaryConditions &);
  procBlock(const primVars &, const plot3dBlock &, const int &, const int &,
            const boundaryConditions &, const int &, const int &, const int &);
  procBlock(const int &, const int &, const int &, const int &);

  // member functions
  int NumCells() const { return numCells_; }
  int NumVars() const { return numVars_; }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }
  int NumGhosts() const { return numGhosts_; }
  int ParentBlock() const { return parBlock_; }
  int LocalPosition() const { return localPos_; }
  int Rank() const { return rank_; }
  int GlobalPos() const { return globalPos_; }

  boundaryConditions BC() const { return bc_; }

  primVars State(const int &ind) const { return state_[ind]; }
  vector<genArray> GetCopyConsVars(const idealGas &) const;

  double Vol(const int &ind) const { return vol_[ind]; }
  vector3d<double> Center(const int &ind) const { return center_[ind]; }
  vector3d<double> FAreaI(const int &ind) const { return fAreaI_[ind]; }
  vector3d<double> FAreaJ(const int &ind) const { return fAreaJ_[ind]; }
  vector3d<double> FAreaK(const int &ind) const { return fAreaK_[ind]; }
  vector3d<double> FCenterI(const int &ind) const { return fCenterI_[ind]; }
  vector3d<double> FCenterJ(const int &ind) const { return fCenterJ_[ind]; }
  vector3d<double> FCenterK(const int &ind) const { return fCenterK_[ind]; }

  double AvgWaveSpeed(const int &ind) const { return avgWaveSpeed_[ind]; }
  double Dt(const int &ind) const { return dt_[ind]; }

  void AddToResidual(const inviscidFlux &, const int &);
  void AddToResidual(const viscousFlux &, const int &);
  void AddToResidual(const source &, const int &);
  genArray Residual(const int &ind) const { return residual_[ind]; }
  double Residual(const int &ind, const int &a) const {
    return residual_[ind][a];
  }

  void CalcCellDt(const int &, const int &, const int &, const double &);

  void CalcInvFluxI(const idealGas &, const input &);
  void CalcInvFluxJ(const idealGas &, const input &);
  void CalcInvFluxK(const idealGas &, const input &);

  void CalcBlockTimeStep(const input &, const double &);
  void UpdateBlock(const input &, const int &, const idealGas &, const double &,
                   const vector<genArray> &, genArray &, resid &);

  void ExplicitEulerTimeAdvance(const idealGas &, const int &, const int &);
  void ImplicitTimeAdvance(const genArray &, const idealGas &, const int &);
  void RK4TimeAdvance(const primVars &, const idealGas &, const double &,
                      const int &, const int &, const int &);

  void ResetResidWS();
  void CleanResizeVecs();

  vector<genArray> AddVolTime(const vector<genArray> &,
                              const vector<genArray> &, const double &,
                              const double &) const;
  void DeltaNMinusOne(vector<genArray> &, const vector<genArray> &,
                      const idealGas &, const double &, const double &);

  double LUSGS(const vector<vector3d<int> > &, vector<genArray> &,
               const vector<genArray> &, const vector<genArray> &,
               const idealGas &, const input &, const sutherland &,
               const turbModel *) const;

  void CalcViscFluxI(const sutherland &, const idealGas &, const input &,
                     const gradients &, const turbModel *);
  void CalcViscFluxJ(const sutherland &, const idealGas &, const input &,
                     const gradients &, const turbModel *);
  void CalcViscFluxK(const sutherland &, const idealGas &, const input &,
                     const gradients &, const turbModel *);

  void CalcGradsI(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsJ(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsK(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;

  void CalcSrcTerms(const gradients &, const sutherland &, const turbModel *);

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input &, const idealGas &,
                                const sutherland &);
  void AssignInviscidGhostCellsEdge(const input &, const idealGas &,
                                    const sutherland &);

  void AssignViscousGhostCells(const input &, const idealGas &,
                               const sutherland &);
  void AssignViscousGhostCellsEdge(const input &, const idealGas &,
                                   const sutherland &);

  bool IsPhysical(const int &, const int &, const int &, const bool &) const;
  bool AtCorner(const int &, const int &, const int &, const bool &) const;
  bool AtEdge(const int &, const int &, const int &, const bool &,
              string &) const;

  vector<bool> PutGeomSlice(const geomSlice &, interblock &, const int &,
                            const int &);
  void PutStateSlice(const stateSlice &, const interblock &, const int &,
                     const int &);

  procBlock Split(const string &, const int &, const int &,
                  vector<boundarySurface> &);
  void Join(const procBlock &, const string &, vector<boundarySurface> &);

  void SwapSliceMPI(const interblock &, const int &, const MPI_Datatype &);
  void PackSendGeomMPI(const MPI_Datatype &, const MPI_Datatype &) const;
  void RecvUnpackGeomMPI(const MPI_Datatype &, const MPI_Datatype &);
  void PackSendSolMPI(const MPI_Datatype &) const;
  void RecvUnpackSolMPI(const MPI_Datatype &);

  // destructor
  ~procBlock() {}
};

// function definitions
double CellSpectralRadius(const vector3d<double> &, const vector3d<double> &,
                          const primVars &, const idealGas &);
double ViscCellSpectralRadius(const vector3d<double> &,
                              const vector3d<double> &, const primVars &,
                              const idealGas &, const sutherland &,
                              const double &, const double &);

template <class T>
T FaceReconCentral(const T &, const T &, const vector3d<double> &,
                   const vector3d<double> &, const vector3d<double> &);

template <class T>
vector<T> PadWithGhosts(const vector<T> &, const int &, const int &,
                        const int &, const int &);

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

vector3d<int> GetSwapLoc(const int &, const int &, const int &,
                         const interblock &, const bool &);
void SwapSlice(interblock &, procBlock &, procBlock &, const bool &);

void GetBoundaryConditions(vector<procBlock> &, const input &, const idealGas &,
                           const sutherland &, vector<interblock> &,
                           const int &, const MPI_Datatype &);

#endif

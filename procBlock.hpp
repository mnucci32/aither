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
#include "plot3d.hpp"              // plot3d
#include "eos.hpp"                 // idealGas
#include "primVars.hpp"            // primVars
#include "inviscidFlux.hpp"        // inviscidFlux
#include "viscousFlux.hpp"         // viscousFlux
#include "input.hpp"               // inputVars
#include "matrix.hpp"              // squareMatrix, matrixDiagonal
#include "boundaryConditions.hpp"  // interblock, patch
#include "macros.hpp"
#include "turbulence.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

// forward declarations
class geomSlice;
class stateSlice;
class gradients;

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
               const idealGas &, const input &, const sutherland &) const;

  void CalcViscFluxI(const sutherland &, const idealGas &, const input &,
                     const gradients &);
  void CalcViscFluxJ(const sutherland &, const idealGas &, const input &,
                     const gradients &);
  void CalcViscFluxK(const sutherland &, const idealGas &, const input &,
                     const gradients &);

  void CalcGradsI(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsJ(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsK(const int &, const int &, const int &, const idealGas &,
                  const bool &, tensor<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;

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

  bool IsPhysical(const int &, const int &, const int &) const;
  bool AtCorner(const int &, const int &, const int &) const;
  bool AtEdge(const int &, const int &, const int &, string &) const;

  geomSlice GetGeomSlice(const int &, const int &, const int &, const int &,
                         const int &, const int &, const bool = false,
                         const bool = false, const bool = false) const;
  vector<bool> PutGeomSlice(const geomSlice &, interblock &, const int &,
                            const int &);

  stateSlice GetStateSlice(const int &, const int &, const int &, const int &,
                           const int &, const int &, const bool = false,
                           const bool = false, const bool = false) const;
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

class geomSlice {
  vector<vector3d<double> > center_;  // coordinates of cell center_
  vector<vector3d<double> > fAreaI_;  // face area vector for i-faces
  vector<vector3d<double> > fAreaJ_;  // face area vector for j-faces
  vector<vector3d<double> > fAreaK_;  // face area vector for k-faces
  vector<vector3d<double> > fCenterI_;  // coordinates of i-face centers
  vector<vector3d<double> > fCenterJ_;  // coordinates of j-face centers
  vector<vector3d<double> > fCenterK_;  // coordinates of k-face centers

  vector<double> vol_;  // cell volume

  int numCells_;  // number of cells in block
  int numI_;  // i-dimension of block (cells)
  int numJ_;  // j-dimension of block (cells)
  int numK_;  // k-dimension of block (cells)
  int parBlock_;  // parent block number

 public:
  // constructors
  geomSlice();
  geomSlice(const int &, const int &, const int &, const int &);

  friend geomSlice procBlock::GetGeomSlice(const int &, const int &,
                                           const int &, const int &,
                                           const int &, const int &,
                                           const bool = false,
                                           const bool = false,
                                           const bool = false) const;

  // member functions
  int NumCells() const { return numCells_; }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }
  int ParentBlock() const { return parBlock_; }

  double Vol(const int &ind) const { return vol_[ind]; }
  vector3d<double> Center(const int &ind) const { return center_[ind]; }
  vector3d<double> FAreaI(const int &ind) const { return fAreaI_[ind]; }
  vector3d<double> FAreaJ(const int &ind) const { return fAreaJ_[ind]; }
  vector3d<double> FAreaK(const int &ind) const { return fAreaK_[ind]; }
  vector3d<double> FCenterI(const int &ind) const { return fCenterI_[ind]; }
  vector3d<double> FCenterJ(const int &ind) const { return fCenterJ_[ind]; }
  vector3d<double> FCenterK(const int &ind) const { return fCenterK_[ind]; }

  // destructor
  ~geomSlice() {}
};

class stateSlice {
  vector<primVars> state_;  // cell states

  int numCells_;  // number of cells in block
  int numI_;  // i-dimension of block (cells)
  int numJ_;  // j-dimension of block (cells)
  int numK_;  // k-dimension of block (cells)
  int parBlock_;  // parent block number

 public:
  // constructors
  stateSlice();
  stateSlice(const int &, const int &, const int &, const int &);

  friend stateSlice procBlock::GetStateSlice(const int &, const int &,
                                             const int &, const int &,
                                             const int &, const int &,
                                             const bool = false,
                                             const bool = false,
                                             const bool = false) const;

  // member functions
  int NumCells() const { return numCells_; }
  int NumI() const { return numI_; }
  int NumJ() const { return numJ_; }
  int NumK() const { return numK_; }
  int ParentBlock() const { return parBlock_; }

  primVars State(const int &ind) const { return state_[ind]; }

  void PackSwapUnpackMPI(const interblock &, const MPI_Datatype &, const int &);

  // destructor
  ~stateSlice() {}
};

class gradients {
  vector<tensor<double> > velocityI_;  // velocity gradients at cell i-face
  vector<tensor<double> > velocityJ_;  // velocity gradients at cell j-face
  vector<tensor<double> > velocityK_;  // velocity gradients at cell k-face
  vector<vector3d<double> > temperatureI_;  // temperature gradients at cell
                                            // i-face
  vector<vector3d<double> > temperatureJ_;  // temperature gradients at cell
                                            // j-face
  vector<vector3d<double> > temperatureK_;  // temperature gradients at cell
                                            // k-face
  vector<vector3d<double> > tkeI_;  // tke gradients at cell i-face
  vector<vector3d<double> > tkeJ_;  // tke gradients at cell j-face
  vector<vector3d<double> > tkeK_;  // tke gradients at cell k-face
  vector<vector3d<double> > omegaI_;  // omega gradients at cell i-face
  vector<vector3d<double> > omegaJ_;  // omega gradients at cell j-face
  vector<vector3d<double> > omegaK_;  // omega gradients at cell k-face

  int imax_;  // number of cells in i-direction
  int jmax_;  // number of cells in j-direction
  int kmax_;  // number of cells in k-direction

 public:
  // constructors
  gradients();
  gradients(const bool &, const procBlock &, const idealGas &);

  // member functions
  tensor<double> VelGradI(const int &a) const { return velocityI_[a]; }
  tensor<double> VelGradJ(const int &a) const { return velocityJ_[a]; }
  tensor<double> VelGradK(const int &a) const { return velocityK_[a]; }
  vector3d<double> TempGradI(const int &a) const { return temperatureI_[a]; }
  vector3d<double> TempGradJ(const int &a) const { return temperatureJ_[a]; }
  vector3d<double> TempGradK(const int &a) const { return temperatureK_[a]; }
  vector3d<double> TkeGradI(const int &a) const { return tkeI_[a]; }
  vector3d<double> TkeGradJ(const int &a) const { return tkeJ_[a]; }
  vector3d<double> TkeGradK(const int &a) const { return tkeK_[a]; }
  vector3d<double> OmegaGradI(const int &a) const { return omegaI_[a]; }
  vector3d<double> OmegaGradJ(const int &a) const { return omegaJ_[a]; }
  vector3d<double> OmegaGradK(const int &a) const { return omegaK_[a]; }

  int NumI() const { return imax_; }
  int NumJ() const { return jmax_; }
  int NumK() const { return kmax_; }

  tensor<double> VelGradCell(const int &, const int &, const int &) const;
  vector3d<double> TempGradCell(const int &, const int &, const int &) const;
  vector3d<double> TkeGradCell(const int &, const int &, const int &) const;
  vector3d<double> OmegaGradCell(const int &, const int &, const int &) const;

  // destructor
  ~gradients() {}
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

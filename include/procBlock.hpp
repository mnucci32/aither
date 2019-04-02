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

#ifndef PROCBLOCKHEADERDEF  // only if the macro PROCBLOCKHEADERDEF is not
                            // defined execute these lines of code
#define PROCBLOCKHEADERDEF  // define the macro

#include <vector>                  // vector
#include <string>                  // string
#include <fstream>
#include <iostream>
#include <memory>
#include <algorithm>
#include <cstdlib>
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "multiArray3d.hpp"        // multiArray3d
#include "plot3d.hpp"              // plot3dBlk
#include "blkMultiArray3d.hpp"     // blkMultiArray3d
#include "tensor.hpp"              // tensor
#include "primitive.hpp"           // primitive
#include "varArray.hpp"            // varArray
#include "arrayView.hpp"           // primitiveView
#include "boundaryConditions.hpp"  // connection, patch
#include "macros.hpp"
#include "uncoupledScalar.hpp"     // uncoupledScalar
#include "wallData.hpp"
#include "utility.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::unique_ptr;

// forward class declarations
class inviscidFlux;
class viscousFlux;
class input;
class geomSlice;
class source;
class resid;
class kdtree;
class conserved;
class matMultiArray3d;
class physics;
class turbModel;
class eos;

class procBlock {
  blkMultiArray3d<primitive> state_;  // primitive vars at cell center
  blkMultiArray3d<conserved> consVarsN_;  // conserved vars at t=n
  blkMultiArray3d<conserved> consVarsNm1_;  // conserved vars at t=n-1

  blkMultiArray3d<residual> residual_;  // cell residual

  multiArray3d<unitVec3dMag<double>> fAreaI_;  // face area vector for i-faces
  multiArray3d<unitVec3dMag<double>> fAreaJ_;  // face area vector for j-faces
  multiArray3d<unitVec3dMag<double>> fAreaK_;  // face area vector for k-faces

  plot3dBlock nodes_;  // coordinates of nodes
  multiArray3d<vector3d<double>> center_;  // coordinates of cell center
  multiArray3d<vector3d<double>> fCenterI_;  // coordinates of i-face centers
  multiArray3d<vector3d<double>> fCenterJ_;  // coordinates of j-face centers
  multiArray3d<vector3d<double>> fCenterK_;  // coordinates of k-face centers

  multiArray3d<double> cellWidthI_;  // i-width of cell
  multiArray3d<double> cellWidthJ_;  // j-width of cell
  multiArray3d<double> cellWidthK_;  // k-width of cell

  multiArray3d<uncoupledScalar> specRadius_;  // maximum wave speed for cell
  multiArray3d<double> vol_;  // cell volume
  multiArray3d<double> dt_;  // cell time step
  multiArray3d<double> wallDist_;  // distance to nearest viscous wall

  // gradients
  multiArray3d<tensor<double>> velocityGrad_;
  multiArray3d<vector3d<double>> temperatureGrad_;
  multiArray3d<vector3d<double>> densityGrad_;
  multiArray3d<vector3d<double>> pressureGrad_;
  multiArray3d<vector3d<double>> tkeGrad_;
  multiArray3d<vector3d<double>> omegaGrad_;
  multiArray3d<vector3d<double>> mixtureGrad_;

  // auxillary variables
  multiArray3d<double> temperature_;
  multiArray3d<double> viscosity_;
  multiArray3d<double> eddyViscosity_;
  multiArray3d<double> f1_;
  multiArray3d<double> f2_;

  boundaryConditions bc_;  // boundary conditions for block

  vector<wallData> wallData_;  // wall variables at viscous walls

  int numGhosts_;  // number of layers of ghost cells surrounding block
  int parBlock_;  // parent block number
  int rank_;  // processor rank
  int localPos_;  // position on local processor
  int globalPos_;  // global position of procBlock in decomposed vector of
                   // procBlocks
  bool isViscous_;
  bool isTurbulent_;
  bool isRANS_;
  bool storeTimeN_;
  bool isMultiLevelTime_;
  bool isMultiSpecies_;

  // private member functions
  void CalcInvFluxI(const physics &, const input &, matMultiArray3d &);
  void CalcInvFluxJ(const physics &, const input &, matMultiArray3d &);
  void CalcInvFluxK(const physics &, const input &, matMultiArray3d &);

  void CalcViscFluxI(const physics &, const input &, matMultiArray3d &);
  void CalcViscFluxJ(const physics &, const input &, matMultiArray3d &);
  void CalcViscFluxK(const physics &, const input &, matMultiArray3d &);

  void CalcCellDt(const int &, const int &, const int &, const double &);

  void ExplicitEulerTimeAdvance(const physics &, const int &, const int &,
                                const int &);
  void ImplicitTimeAdvance(const varArrayView &, const physics &, const int &,
                           const int &, const int &);
  void RK4TimeAdvance(const conservedView &, const physics &, const int &,
                      const int &, const int &, const int &);
  template <typename T>
  void AddToResidual(const int &, const int &, const int &, const T &);
  template <typename T>
  void SubtractFromResidual(const int &, const int &, const int &, const T &);
  vector<wallData> SplitWallData(const string &, const int &);
  void JoinWallData(const vector<wallData> &, const string &);

  void CalcGradsI();
  void CalcGradsJ();
  void CalcGradsK();

 public:
  // constructors
  procBlock(const plot3dBlock &, const int &, const boundaryConditions &,
            const int &, const int &, const int &, const input &);
  procBlock(const int &, const int &, const int &, const int &, const int &,
            const int &, const bool &, const bool &, const bool &, const bool &,
            const bool &, const bool &);
  procBlock()
      : procBlock(0, 0, 0, 0, 0, 0, false, false, false, false, false, false) {}

  // move constructor and assignment operator
  procBlock(procBlock&&) noexcept = default;
  procBlock& operator=(procBlock&&) noexcept = default;

  // copy constructor and assignment operator
  procBlock(const procBlock&) = default;
  procBlock& operator=(const procBlock&) = default;

  // member functions
  int NumCells() const { return residual_.NumBlocks(); }
  int NumCellsGhosts() const { return state_.NumBlocks(); }
  int NumEquations() const { return residual_(0, 0, 0).Size(); }
  int NumSpecies() const { return residual_(0, 0, 0).NumSpecies(); }
  int NumI() const { return residual_.NumI(); }
  int NumJ() const { return residual_.NumJ(); }
  int NumK() const { return residual_.NumK(); }
  int NumIG() const { return state_.NumI(); }
  int NumJG() const { return state_.NumJ(); }
  int NumKG() const { return state_.NumK(); }

  int StartI() const { return residual_.StartI(); }
  int StartJ() const { return residual_.StartJ(); }
  int StartK() const { return residual_.StartK(); }
  int StartIG() const { return state_.StartI(); }
  int StartJG() const { return state_.StartJ(); }
  int StartKG() const { return state_.StartK(); }
  int EndI() const { return residual_.EndI(); }
  int EndJ() const { return residual_.EndJ(); }
  int EndK() const { return residual_.EndK(); }
  int EndIG() const { return state_.EndI(); }
  int EndJG() const { return state_.EndJ(); }
  int EndKG() const { return state_.EndK(); }

  int Start(const string &dir) const {
    if (dir == "i") return this->StartI();
    else if (dir == "j") return this->StartJ();
    else if (dir == "k") return this->StartK();
    else {
      cerr << "ERROR: direction " << dir << " not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  int End(const string &dir) const {
    if (dir == "i") return this->EndI();
    else if (dir == "j") return this->EndJ();
    else if (dir == "k") return this->EndK();
    else {
      cerr << "ERROR: direction " << dir << " not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  const int & NumGhosts() const {return numGhosts_;}
  const int & ParentBlock() const {return parBlock_;}
  const int & LocalPosition() const {return localPos_;}
  const int & Rank() const {return rank_;}
  const int & GlobalPos() const {return globalPos_;}
  const bool & IsViscous() const {return isViscous_;}
  const bool & IsTurbulent() const {return isTurbulent_;}
  const bool & IsRANS() const {return isRANS_;}

  const boundaryConditions & BC() const {return bc_;}

  primitiveView State(const int &ii, const int &jj, const int &kk) const {
    return state_(ii, jj, kk);
  }
  void RestrictState(const procBlock &fine,
                     const multiArray3d<vector3d<int>> &toCoarse,
                     const multiArray3d<double> &volWeightFactor);
  conservedView ConsVarsN(const int &ii, const int &jj, const int &kk) const {
    return consVarsN_(ii, jj, kk);
  }
  conservedView ConsVarsNm1(const int &ii, const int &jj, const int &kk) const {
    return consVarsNm1_(ii, jj, kk);
  }

  blkMultiArray3d<primitive> SliceState(const int &, const int &, const int &,
                                        const int &, const int &,
                                        const int &) const;
  multiArray3d<vector3d<double>> SliceBoundaryCenters(const int &) const;

  void AssignSolToTimeN(const physics &);
  void AssignSolToTimeNm1();
  double SolDeltaNCoeff(const int &, const int &, const int &,
                        const input &) const;
  double SolDeltaNm1Coeff(const int &, const int &, const int &,
                          const input &) const;
  varArray SolDeltaMmN(const int &, const int &, const int &, const input &,
                       const physics &) const;
  varArray SolDeltaNm1(const int &, const int &, const int &,
                       const input &) const;

  const double &Vol(const int &ii, const int &jj, const int &kk) const {
    return vol_(ii, jj, kk);
  }
  const vector3d<double> &Node(const int &ii, const int &jj,
                               const int &kk) const {
    return nodes_.Coords(ii, jj, kk);
  }
  const plot3dBlock &Nodes() const { return nodes_; }
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
  const unitVec3dMag<double> &FAreaI(const int &ii, const int &jj,
                                     const int &kk) const {
    return fAreaI_(ii, jj, kk);
  }
  const unitVec3dMag<double> &FAreaJ(const int &ii, const int &jj,
                                     const int &kk) const {
    return fAreaJ_(ii, jj, kk);
  }
  const unitVec3dMag<double> &FAreaK(const int &ii, const int &jj,
                                     const int &kk) const {
    return fAreaK_(ii, jj, kk);
  }

  const vector3d<double> &FCenterI(const int &ii, const int &jj,
                                   const int &kk) const {
    return fCenterI_(ii, jj, kk);
  }
  const vector3d<double> &FCenterJ(const int &ii, const int &jj,
                                   const int &kk) const {
    return fCenterJ_(ii, jj, kk);
  }
  const vector3d<double> &FCenterK(const int &ii, const int &jj,
                                   const int &kk) const {
    return fCenterK_(ii, jj, kk);
  }

  const double &CellWidthI(const int &ii, const int &jj, const int &kk) const {
    return cellWidthI_(ii, jj, kk);
  }
  const double &CellWidthJ(const int &ii, const int &jj, const int &kk) const {
    return cellWidthJ_(ii, jj, kk);
  }
  const double &CellWidthK(const int &ii, const int &jj, const int &kk) const {
    return cellWidthK_(ii, jj, kk);
  }
  double MaxCellWidth(const int &ii, const int &jj, const int &kk) const {
    return std::max(std::max(cellWidthI_(ii, jj, kk), cellWidthJ_(ii, jj, kk)),
                    cellWidthK_(ii, jj, kk));
  }
  double MinCellWidth(const int &ii, const int &jj, const int &kk) const {
    return std::min(std::min(cellWidthI_(ii, jj, kk), cellWidthJ_(ii, jj, kk)),
                    cellWidthK_(ii, jj, kk));
  }

  const uncoupledScalar &SpectralRadius(const int &ii, const int &jj,
                                        const int &kk) const {
    return specRadius_(ii, jj, kk);
  }
  const double &Dt(const int &ii, const int &jj, const int &kk) const {
    return dt_(ii, jj, kk);
  }
  const double &WallDist(const int &ii, const int &jj, const int &kk) const {
    return wallDist_(ii, jj, kk);
  }

  residualView Residual(const int &ii, const int &jj, const int &kk) const {
    return residual_(ii, jj, kk);
  }
  const double &Residual(const int &ii, const int &jj, const int &kk,
                         const int &a) const {
    return residual_(ii, jj, kk)[a];
  }
  const blkMultiArray3d<residual> &Residuals() const { return residual_; }

  const tensor<double> &VelGrad(const int &ii, const int &jj,
                                const int &kk) const {
    return velocityGrad_(ii, jj, kk);
  }
  const vector3d<double> &TempGrad(const int &ii, const int &jj,
                                   const int &kk) const {
    return temperatureGrad_(ii, jj, kk);
  }
  const vector3d<double> &DensityGrad(const int &ii, const int &jj,
                                      const int &kk) const {
    return densityGrad_(ii, jj, kk);
  }
  const vector3d<double> &PressureGrad(const int &ii, const int &jj,
                                       const int &kk) const {
    return pressureGrad_(ii, jj, kk);
  }
  vector3d<double> TkeGrad(const int &ii, const int &jj, const int &kk) const {
    return isRANS_ ? tkeGrad_(ii, jj, kk) : vector3d<double>();
  }
  vector3d<double> OmegaGrad(const int &ii, const int &jj,
                             const int &kk) const {
    return isRANS_ ? omegaGrad_(ii, jj, kk) : vector3d<double>();
  }
  vector3d<double> MixtureGrad(const int &ii, const int &jj, const int &kk,
                               const int &ll) const {
    return isMultiSpecies_ ? mixtureGrad_(ii, jj, kk, ll) : vector3d<double>();
  }

  const double &Temperature(const int &ii, const int &jj, const int &kk) const {
    return temperature_(ii, jj, kk);
  }
  double Viscosity(const int &ii, const int &jj, const int &kk) const {
    return isViscous_ ? viscosity_(ii, jj, kk) : 0.0;
  }
  double EddyViscosity(const int &ii, const int &jj, const int &kk) const {
    return isTurbulent_ ? eddyViscosity_(ii, jj, kk) : 0.0;
  }
  double F1(const int &ii, const int &jj, const int &kk) const {
    return isRANS_ ? f1_(ii, jj, kk) : 0.0;
  }
  double F2(const int &ii, const int &jj, const int &kk) const {
    return isRANS_ ? f2_(ii, jj, kk) : 0.0;
  }

  void CalcBlockTimeStep(const input &);
  void UpdateBlock(const input &, const physics &,
                   const blkMultiArray3d<varArray> &, const int &, residual &,
                   resid &);

  void CalcResidualNoSource(const physics &, const input &, matMultiArray3d &);
  void CalcSrcTerms(const physics &, const input &, matMultiArray3d &);

  void ResetResidWS();
  void ResetGradients();
  void ResetTurbVars();
  void CleanResizeVecs(const int &, const int &, const int &, const int &,
                       const int &, const int &);

  void InitializeStates(const input &, const physics &);

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input &, const physics &);
  void AssignInviscidGhostCellsEdge(const input &, const physics &);
  void AssignCornerGhostCells();
  void AssignViscousGhostCells(const input &, const physics &);
  void AssignViscousGhostCellsEdge(const input &, const physics &);
  blkMultiArray3d<primitive> GetGhostStates(
      const blkMultiArray3d<primitive> &, const string &,
      const multiArray3d<unitVec3dMag<double>> &, const multiArray3d<double> &,
      const boundarySurface &, const input &, const physics &, const int &,
      const multiArray3d<double> & = {},
      const multiArray3d<double> & = {},
      const blkMultiArray3d<conserved> & = {},
      const multiArray3d<vector3d<double>> & = {},
      const multiArray3d<tensor<double>> & = {});

  void CalcGradsI(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &,
                  vector<vector3d<double>> &) const;
  void CalcGradsJ(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &,
                  vector<vector3d<double>> &) const;
  void CalcGradsK(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &,
                  vector<vector3d<double>> &) const;

  void CalcWallDistance(const kdtree &);

  varArray ImplicitLower(const int &, const int &, const int &,
                         const blkMultiArray3d<varArray> &, const physics &,
                         const input &) const;
  varArray ImplicitUpper(const int &, const int &, const int &,
                         const blkMultiArray3d<varArray> &, const physics &,
                         const input &) const;

  bool IsPhysical(const int &ii, const int &jj, const int &kk) const {
    return state_.IsPhysical(ii, jj, kk);
  }
  bool AtCorner(const int &ii, const int &jj, const int &kk) const {
    return state_.AtCorner(ii, jj, kk);
  }
  bool AtEdge(const int &ii, const int &jj, const int &kk, string &dir) const {
    return state_.AtEdge(ii, jj, kk, dir);
  }
  bool AtEdgeInclusive(const int &ii, const int &jj, const int &kk,
                       string &dir) const {
    return state_.AtEdgeInclusive(ii, jj, kk, dir);
  }
  bool AtGhostNonEdge(const int &ii, const int &jj, const int &kk, string &dir,
                      int &type) const {
    return state_.AtGhostNonEdge(ii, jj, kk, dir, type);
  }

  vector<bool> PutGeomSlice(const geomSlice &, connection &, const int &);
  void PutStateSlice(const blkMultiArray3d<primitive> &, const connection &,
                     const int &, const int &);

  procBlock Split(const string &, const int &, const int &,
                  vector<boundarySurface> &);
  void Join(const procBlock &, const string &, vector<boundarySurface> &);

  void SwapStateSlice(const connection &, procBlock &);
  void SwapStateSliceMPI(const connection &, const int &);
  void SwapTurbSlice(const connection &, procBlock &);
  void SwapTurbSliceMPI(const connection &, const int &);
  void SwapWallDistSlice(const connection &, procBlock &);
  void SwapWallDistSliceMPI(const connection &, const int &);
  void SwapEddyViscAndGradientSlice(const connection &, procBlock &);
  void SwapEddyViscAndGradientSliceMPI(const connection &, const int &,
                                       const MPI_Datatype &,
                                       const MPI_Datatype &);

  void PackSendGeomMPI(const MPI_Datatype &, const MPI_Datatype &) const;
  void RecvUnpackGeomMPI(const MPI_Datatype &, const MPI_Datatype &,
                         const input &);
  void PackSendSolMPI(const MPI_Datatype &, const MPI_Datatype &,
                      const MPI_Datatype &) const;
  void RecvUnpackSolMPI(const MPI_Datatype &, const MPI_Datatype &,
                        const MPI_Datatype &, const input &);

  void UpdateAuxillaryVariables(const physics &, const bool = true);
  void UpdateUnlimTurbEddyVisc(const unique_ptr<turbModel> &, const bool &);

  double ProjC2CDist(const int &, const int &, const int &,
                     const string &) const;

  void DumpToFile(const string &, const string &) const;
  void CalcCellWidths();
  void GetStatesFromRestart(const blkMultiArray3d<primitive> &);
  void GetSolNm1FromRestart(const blkMultiArray3d<conserved> &);

  int WallDataIndex(const boundarySurface &) const;
  int WallDataSize() const {return wallData_.size();}
  bool HasWallData() const {return this->WallDataSize() > 0;}
  boundarySurface WallSurface(const int &ii) const {
    return wallData_[ii].Surface();
  }
  int NumWallSurfI(const int &ii) const { return wallData_[ii].NumI(); }
  int NumWallSurfJ(const int &ii) const { return wallData_[ii].NumJ(); }
  int NumWallSurfK(const int &ii) const { return wallData_[ii].NumK(); }
  double WallYplus(const int &ss, const int &ii, const int &jj,
                   const int &kk) const {
    return wallData_[ss].Yplus(ii, jj, kk);
  }
  double WallHeatFlux(const int &ss, const int &ii, const int &jj,
                      const int &kk) const {
    return wallData_[ss].WallHeatFlux(ii, jj, kk);
  }
  vector3d<double> WallShearStress(const int &ss, const int &ii, const int &jj,
                                   const int &kk) const {
    return wallData_[ss].WallShearStress(ii, jj, kk);
  }
  double WallTemperature(const int &ss, const int &ii, const int &jj,
                         const int &kk) const {
    return wallData_[ss].WallTemperature(ii, jj, kk);
  }
  double WallEddyVisc(const int &ss, const int &ii, const int &jj,
                      const int &kk) const {
    return wallData_[ss].WallEddyViscosity(ii, jj, kk);
  }
  double WallViscosity(const int &ss, const int &ii, const int &jj,
                       const int &kk) const {
    return wallData_[ss].WallViscosity(ii, jj, kk);
  }
  double WallDensity(const int &ss, const int &ii, const int &jj,
                     const int &kk) const {
    return wallData_[ss].WallDensity(ii, jj, kk);
  }
  double WallFrictionVelocity(const int &ss, const int &ii, const int &jj,
                              const int &kk) const {
    return wallData_[ss].WallFrictionVelocity(ii, jj, kk);
  }
  double WallPressure(const int &ss, const int &ii, const int &jj,
                      const int &kk, const unique_ptr<eos> &eqnState) const {
    return wallData_[ss].WallPressure(ii, jj, kk, eqnState);
  }
  double WallTke(const int &ss, const int &ii, const int &jj,
                 const int &kk) const {
    return wallData_[ss].WallTke(ii, jj, kk);
  }
  double WallSdr(const int &ss, const int &ii, const int &jj,
                 const int &kk) const {
    return wallData_[ss].WallSdr(ii, jj, kk);
  }
  void GetCoarseMeshAndBCs(vector<plot3dBlock> &mesh,
                           vector<boundaryConditions> &bcs,
                           vector<multiArray3d<vector3d<int>>> &toCoarse,
                           vector<multiArray3d<double>> &volFac) const;
  procBlock CellToNode() const;
  void AddCoarseGridCorrection(const blkMultiArray3d<varArray> &correction) {
    state_ += correction;
  }

  // destructor
  ~procBlock() noexcept {}
};

// ----------------------------------------------------------------------------
// function definitions
template <typename T>
void procBlock::AddToResidual(const int &ii, const int &jj, const int &kk,
                              const T &arr) {
  static_assert(std::is_base_of<varArray, T>::value, "T must be varArray type");
  MSG_ASSERT(arr.Size() == residual_.BlockSize(),
             "array block size must match residual block size");
  for (auto bb = 0; bb < residual_.BlockSize(); ++bb) {
    residual_(ii, jj, kk, bb) += arr[bb];
  }
}

template <typename T>
void procBlock::SubtractFromResidual(const int &ii, const int &jj,
                                     const int &kk, const T &arr) {
  static_assert(std::is_base_of<varArray, T>::value, "T must be varArray type");
  MSG_ASSERT(arr.Size() == residual_.BlockSize(),
             "array block size must match residual block size");
  for (auto bb = 0; bb < residual_.BlockSize(); ++bb) {
    residual_(ii, jj, kk, bb) -= arr[bb];
  }
}

/* Function to pad a multiArray3d with a specified number of ghost cells
           ___ ___ ___ ___ ___ ___ ___ ___
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|

In the above diagram, the cells marked with an "X" represent physical cells. The
entire diagram represents the block (in 2D) padded with 2 layers of ghost cells.
The cells marked with "G" are regualar ghost cells. The cells marked with "E" are
ghost cells located along one of the 12 edges that form a plot3d block. In 3D
there are also "corner" cells located at the 8 corners that form the plot3d block.
These cells are not used though. There is a place in the vector for them to make
accessing the padded vector of cells the same as for a plot3d block without ghost
cells.
*/
template <typename T>
T PadWithGhosts(const T &var, const int &numGhosts) {
  // var -- vector of variables to pad (no ghost cells included)
  // numGhosts -- number of layers of ghost cells to pad var with
  // T should be multiArray3d or blkMultiArray3d type

  // initialize added array
  T padBlk(var.NumI(), var.NumJ(), var.NumK(), numGhosts, var.BlockInfo());

  padBlk.Insert(var.RangeI(), var.RangeJ(), var.RangeK(), var);
  return padBlk;
}


template <typename T1, typename T2>
void BlockRestriction(const blkMultiArray3d<T1>& fine,
                      const multiArray3d<vector3d<int>>& toCoarse,
                      const multiArray3d<double>& volFac,
                      blkMultiArray3d<T2>& coarse) {
  // use volume weighted average
  for (auto kk = fine.PhysStartK(); kk < fine.PhysEndK(); ++kk) {
    for (auto jj = fine.PhysStartJ(); jj < fine.PhysEndJ(); ++jj) {
      for (auto ii = fine.PhysStartI(); ii < fine.PhysEndI(); ++ii) {
        const auto ci = toCoarse(ii, jj, kk);
        T2 restricted =
            coarse(ci[0], ci[1], ci[2]) + volFac(ii, jj, kk) * fine(ii, jj, kk);
        coarse.InsertBlock(ci[0], ci[1], ci[2], restricted);
      }
    }
  }
}

#endif









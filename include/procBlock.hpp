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

#ifndef PROCBLOCKHEADERDEF  // only if the macro PROCBLOCKHEADERDEF is not
                            // defined execute these lines of code
#define PROCBLOCKHEADERDEF  // define the macro

#include <vector>                  // vector
#include <string>                  // string
#include <fstream>
#include <iostream>
#include <memory>
#include <algorithm>
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"            // vector3d
#include "multiArray3d.hpp"        // multiArray3d
#include "tensor.hpp"              // tensor
#include "primVars.hpp"            // primVars
#include "genArray.hpp"            // genArray
#include "boundaryConditions.hpp"  // connection, patch
#include "macros.hpp"
#include "uncoupledScalar.hpp"     // uncoupledScalar
#include "wallData.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::unique_ptr;

// forward class declarations
class eos;
class transport;
class thermodynamic;
class inviscidFlux;
class viscousFlux;
class input;
class geomSlice;
class source;
class turbModel;
class plot3dBlock;
class resid;
class fluxJacobian;
class kdtree;

class procBlock {
  multiArray3d<primVars> state_;  // primative variables at cell center
  multiArray3d<genArray> consVarsN_;  // conserved variables at time n
  multiArray3d<genArray> consVarsNm1_;  // conserved variables at time n-1

  multiArray3d<genArray> residual_;  // cell residual

  multiArray3d<unitVec3dMag<double>> fAreaI_;  // face area vector for i-faces
  multiArray3d<unitVec3dMag<double>> fAreaJ_;  // face area vector for j-faces
  multiArray3d<unitVec3dMag<double>> fAreaK_;  // face area vector for k-faces

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

  // private member functions
  void CalcInvFluxI(const unique_ptr<eos> &, const unique_ptr<thermodynamic> &,
                    const input &, const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);
  void CalcInvFluxJ(const unique_ptr<eos> &, const unique_ptr<thermodynamic> &,
                    const input &, const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);
  void CalcInvFluxK(const unique_ptr<eos> &, const unique_ptr<thermodynamic> &,
                    const input &, const unique_ptr<turbModel> &,
                    multiArray3d<fluxJacobian> &);

  void CalcViscFluxI(const unique_ptr<transport> &,
                     const unique_ptr<thermodynamic> &, const unique_ptr<eos> &,
                     const input &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);
  void CalcViscFluxJ(const unique_ptr<transport> &,
                     const unique_ptr<thermodynamic> &, const unique_ptr<eos> &,
                     const input &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);
  void CalcViscFluxK(const unique_ptr<transport> &,
                     const unique_ptr<thermodynamic> &, const unique_ptr<eos> &,
                     const input &, const unique_ptr<turbModel> &,
                     multiArray3d<fluxJacobian> &);

  void CalcCellDt(const int &, const int &, const int &, const double &);

  void ExplicitEulerTimeAdvance(const unique_ptr<eos> &,
                                const unique_ptr<thermodynamic> &,
                                const unique_ptr<turbModel> &, const int &,
                                const int &, const int &);
  void ImplicitTimeAdvance(const genArray &, const unique_ptr<eos> &,
                           const unique_ptr<thermodynamic> &,
                           const unique_ptr<turbModel> &, const int &,
                           const int &, const int &);
  void RK4TimeAdvance(const genArray &, const unique_ptr<eos> &,
                      const unique_ptr<thermodynamic> &,
                      const unique_ptr<turbModel> &, const int &, const int &,
                      const int &, const int &);

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
  vector<wallData> SplitWallData(const string &, const int &);
  void JoinWallData(const vector<wallData> &, const string &);

  void CalcGradsI();
  void CalcGradsJ();
  void CalcGradsK();

 public:
  // constructors
  procBlock(const plot3dBlock &, const int &, const boundaryConditions &,
            const int &, const int &, const int &, const input &);
  procBlock(const int &, const int &, const int &, const int &, const bool &,
            const bool &, const bool &, const bool &, const bool &);
  procBlock() : procBlock(1, 1, 1, 0, false, false, false, false, false) {}

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

  int NumGhosts() const {return numGhosts_;}
  int ParentBlock() const {return parBlock_;}
  int LocalPosition() const {return localPos_;}
  int Rank() const {return rank_;}
  int GlobalPos() const {return globalPos_;}
  bool IsViscous() const {return isViscous_;}
  bool IsTurbulent() const {return isTurbulent_;}
  bool IsRANS() const {return isRANS_;}

  boundaryConditions BC() const {return bc_;}

  primVars State(const int &ii, const int &jj, const int &kk) const {
    return state_(ii, jj, kk);
  }
  genArray ConsVarsN(const int &ii, const int &jj, const int &kk) const {
    return consVarsN_(ii, jj, kk);
  }
  genArray ConsVarsNm1(const int &ii, const int &jj, const int &kk) const {
    return consVarsNm1_(ii, jj, kk);
  }

  multiArray3d<primVars> SliceState(const int &, const int &, const int &,
                                    const int &, const int &,
                                    const int &) const;
  multiArray3d<vector3d<double>> SliceBoundaryCenters(const int &) const;

  void AssignSolToTimeN(const unique_ptr<eos> &,
                        const unique_ptr<thermodynamic> &);
  void AssignSolToTimeNm1();
  double SolDeltaNCoeff(const int &, const int &, const int &,
                        const input &) const;
  double SolDeltaNm1Coeff(const int &, const int &, const int &,
                          const input &) const;
  genArray SolDeltaMmN(const int &, const int &, const int &, const input &,
                       const unique_ptr<eos> &,
                       const unique_ptr<thermodynamic> &) const;
  genArray SolDeltaNm1(const int &, const int &, const int &,
                       const input &) const;

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

  double CellWidthI(const int &ii, const int &jj, const int &kk) const {
    return cellWidthI_(ii, jj, kk);
  }
  double CellWidthJ(const int &ii, const int &jj, const int &kk) const {
    return cellWidthJ_(ii, jj, kk);
  }
  double CellWidthK(const int &ii, const int &jj, const int &kk) const {
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

  uncoupledScalar SpectralRadius(const int &ii, const int &jj,
                                 const int &kk) const {
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

  tensor<double> VelGrad(const int &ii, const int &jj, const int &kk) const {
    return velocityGrad_(ii, jj, kk);
  }
  vector3d<double> TempGrad(const int &ii, const int &jj, const int &kk) const {
    return temperatureGrad_(ii, jj, kk);
  }
  vector3d<double> DensityGrad(const int &ii, const int &jj,
                               const int &kk) const {
    return densityGrad_(ii, jj, kk);
  }
  vector3d<double> PressureGrad(const int &ii, const int &jj,
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

  double Temperature(const int &ii, const int &jj, const int &kk) const {
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
  void UpdateBlock(const input &, const unique_ptr<eos> &,
                   const unique_ptr<thermodynamic> &,
                   const unique_ptr<transport> &,
                   const multiArray3d<genArray> &,
                   const unique_ptr<turbModel> &, const int &, genArray &,
                   resid &);

  void CalcResidualNoSource(const unique_ptr<transport> &,
                            const unique_ptr<thermodynamic> &,
                            const unique_ptr<eos> &, const input &,
                            const unique_ptr<turbModel> &,
                            multiArray3d<fluxJacobian> &);
  void CalcSrcTerms(const unique_ptr<transport> &,
                    const unique_ptr<turbModel> &, const input &,
                    multiArray3d<fluxJacobian> &);

  void ResetResidWS();
  void ResetGradients();
  void ResetTurbVars();
  void CleanResizeVecs(const int &, const int &, const int &, const int &);

  void InitializeStates(const input &, const unique_ptr<eos> &,
                        const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &);

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input &, const unique_ptr<eos> &,
                                const unique_ptr<thermodynamic> &,
                                const unique_ptr<transport> &,
                                const unique_ptr<turbModel> &);
  void AssignInviscidGhostCellsEdge(const input &, const unique_ptr<eos> &,
                                    const unique_ptr<thermodynamic> &thermo,
                                    const unique_ptr<transport> &,
                                    const unique_ptr<turbModel> &);

  void AssignViscousGhostCells(const input &, const unique_ptr<eos> &,
                               const unique_ptr<thermodynamic> &,
                               const unique_ptr<transport> &,
                               const unique_ptr<turbModel> &);
  void AssignViscousGhostCellsEdge(const input &, const unique_ptr<eos> &,
                                   const unique_ptr<thermodynamic> &,
                                   const unique_ptr<transport> &,
                                   const unique_ptr<turbModel> &);
  multiArray3d<primVars> GetGhostStates(
      const multiArray3d<primVars> &, const string &,
      const multiArray3d<unitVec3dMag<double>> &, const multiArray3d<double> &,
      const boundarySurface &, const input &, const unique_ptr<eos> &,
      const unique_ptr<thermodynamic> &, const unique_ptr<transport> &,
      const unique_ptr<turbModel> &, const int &,
      const multiArray3d<double> & = {}, const multiArray3d<genArray> & = {},
      const multiArray3d<vector3d<double>> & = {},
      const multiArray3d<tensor<double>> & = {});

  void CalcGradsI(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsJ(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;
  void CalcGradsK(const int &, const int &, const int &, tensor<double> &,
                  vector3d<double> &, vector3d<double> &, vector3d<double> &,
                  vector3d<double> &, vector3d<double> &) const;

  void CalcWallDistance(const kdtree &);

  multiArray3d<genArray> DeltaNMinusOne(const multiArray3d<genArray> &,
                                        const unique_ptr<eos> &, const double &,
                                        const double &) const;
  multiArray3d<genArray> SolTimeMMinusN(const multiArray3d<genArray> &,
                                        const unique_ptr<eos> &, const input &,
                                        const int &) const;
  void InvertDiagonal(multiArray3d<fluxJacobian> &, const input &) const;

  multiArray3d<genArray> InitializeMatrixUpdate(
      const input &, const unique_ptr<eos> &eos,
      const unique_ptr<thermodynamic> &,
      const multiArray3d<fluxJacobian> &) const;
  void LUSGS_Forward(const vector<vector3d<int>> &, multiArray3d<genArray> &,
                     const unique_ptr<eos> &, const input &,
                     const unique_ptr<thermodynamic> &,
                     const unique_ptr<transport> &,
                     const unique_ptr<turbModel> &,
                     const multiArray3d<fluxJacobian> &, const int &) const;
  double LUSGS_Backward(const vector<vector3d<int>> &, multiArray3d<genArray> &,
                        const unique_ptr<eos> &, const input &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &,
                        const multiArray3d<fluxJacobian> &, const int &) const;

  double DPLUR(multiArray3d<genArray> &, const unique_ptr<eos> &, const input &,
               const unique_ptr<thermodynamic> &, const unique_ptr<transport> &,
               const unique_ptr<turbModel> &,
               const multiArray3d<fluxJacobian> &) const;

  bool IsPhysical(const int &, const int &, const int &) const;
  bool AtCorner(const int &, const int &, const int &) const;
  bool AtEdge(const int &, const int &, const int &, string &) const;
  bool AtEdgeInclusive(const int &, const int &, const int &, string &) const;
  bool AtGhostNonEdge(const int &, const int &, const int &, string &,
                      int &) const;

  vector<bool> PutGeomSlice(const geomSlice &, connection &, const int &);
  void PutStateSlice(const multiArray3d<primVars> &, const connection &,
                     const int &, const int &);

  procBlock Split(const string &, const int &, const int &,
                  vector<boundarySurface> &);
  void Join(const procBlock &, const string &, vector<boundarySurface> &);

  void SwapStateSlice(const connection &, procBlock &);
  void SwapStateSliceMPI(const connection &, const int &, const MPI_Datatype &);
  void SwapTurbSlice(const connection &, procBlock &);
  void SwapTurbSliceMPI(const connection &, const int &);
  void SwapWallDistSlice(const connection &, procBlock &);
  void SwapWallDistSliceMPI(const connection &, const int &);
  void SwapEddyViscAndGradientSlice(const connection &, procBlock &);
  void SwapEddyViscAndGradientSliceMPI(const connection &, const int &,
                                       const MPI_Datatype &,
                                       const MPI_Datatype &);

  void PackSendGeomMPI(const MPI_Datatype &, const MPI_Datatype &,
                       const MPI_Datatype &, const MPI_Datatype &) const;
  void RecvUnpackGeomMPI(const MPI_Datatype &, const MPI_Datatype &,
                         const MPI_Datatype &, const MPI_Datatype &,
                         const input &);
  void PackSendSolMPI(const MPI_Datatype &, const MPI_Datatype &,
                      const MPI_Datatype &, const MPI_Datatype &,
                      const MPI_Datatype &) const;
  void RecvUnpackSolMPI(const MPI_Datatype &, const MPI_Datatype &,
                        const MPI_Datatype &, const MPI_Datatype &,
                        const MPI_Datatype &, const input &);

  void UpdateAuxillaryVariables(const unique_ptr<eos> &,
                                const unique_ptr<transport> &,
                                const bool = true);
  void UpdateUnlimTurbEddyVisc(const unique_ptr<turbModel> &, const bool &);

  double ProjC2CDist(const int &, const int &, const int &,
                     const string &) const;

  void DumpToFile(const string &, const string &) const;
  void CalcCellWidths();
  void GetStatesFromRestart(const multiArray3d<primVars> &);
  void GetSolNm1FromRestart(const multiArray3d<genArray> &);

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

  // destructor
  ~procBlock() noexcept {}
};

// function definitions
template <typename T>
multiArray3d<T> PadWithGhosts(const multiArray3d<T> &, const int &);


#endif









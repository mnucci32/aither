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

#ifndef PROCBLOCKHEADERDEF             //only if the macro PROCBLOCKHEADERDEF is not defined execute these lines of code
#define PROCBLOCKHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.hpp" //vector3d
#include "tensor.hpp" //tensor
#include "plot3d.hpp" //plot3d
#include "eos.hpp" //idealGas
#include "primVars.hpp" //primVars
#include "inviscidFlux.hpp" //inviscidFlux
#include "viscousFlux.hpp" //viscousFlux
#include "input.hpp" //inputVars
#include "matrix.hpp" //squareMatrix, matrixDiagonal
#include "boundaryConditions.hpp" //interblock, patch
#include "mpi.h" //parallelism
#include <fstream>
#include <iostream>
#include "macros.hpp"
#include "turbulence.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

//forward declarations
class geomSlice;
class stateSlice;
class gradients;

class procBlock {
  vector<primVars> state;            //primative variables at cell center

  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;        //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;        //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;        //coordinates of k-face centers

  vector<genArray> residual;                 //cell residual

  vector<double> vol ;                      //cell volume
  vector<double> avgWaveSpeed;             //maximum wave speed for cell
  vector<double> dt;                        //cell time step

  boundaryConditions bc;                   //boundary conditions for block

  int numCells;                            //number of cells in block
  int numVars;                             //number of variables stored at cell
  int numI;                                //i-dimension of block (cells)
  int numJ;                                //j-dimension of block (cells)
  int numK;                                //k-dimension of block (cells)
  int numGhosts;                           //number of layers of ghost cells surrounding block
  int parBlock;                            //parent block number
  int rank;                                //processor rank
  int localPos;                            //position on local processor
  int globalPos;                           //global position of procBlock in decomposed vector of procBlocks

 public:
  //constructors
  procBlock();
  procBlock( const plot3dBlock& , const int&, const int& );
  procBlock( const double, const double, const vector3d<double>, const plot3dBlock&, const int&, const int&, const boundaryConditions& );
  procBlock( const primVars&, const plot3dBlock&, const int &, const int&, const boundaryConditions&, const int&, const int&, const int& );
  procBlock( const int&, const int&, const int&, const int& );

  //member functions
  int NumCells() const {return numCells;}
  int NumVars() const {return numVars;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}
  int NumGhosts() const {return numGhosts;}
  int ParentBlock() const {return parBlock;}
  int LocalPosition() const {return localPos;}
  int Rank() const {return rank;}
  int GlobalPos() const {return globalPos;}

  boundaryConditions BC() const {return bc;}

  primVars State(const int &ind) const {return state[ind];}
  vector<genArray> GetCopyConsVars(const idealGas &) const; 

  double Vol(const int &ind) const {return vol[ind];}
  vector3d<double> Center(const int &ind) const {return center[ind];}
  vector3d<double> FAreaI(const int &ind) const {return fAreaI[ind];}
  vector3d<double> FAreaJ(const int &ind) const {return fAreaJ[ind];}
  vector3d<double> FAreaK(const int &ind) const {return fAreaK[ind];}
  vector3d<double> FCenterI(const int &ind) const {return fCenterI[ind];}
  vector3d<double> FCenterJ(const int &ind) const {return fCenterJ[ind];}
  vector3d<double> FCenterK(const int &ind) const {return fCenterK[ind];}

  double AvgWaveSpeed(const int &ind) const {return avgWaveSpeed[ind];}
  double Dt(const int &ind) const {return dt[ind];}

  void AddToResidual( const inviscidFlux &, const int &);
  void AddToResidual( const viscousFlux &, const int &);
  genArray Residual(const int &ind) const {return residual[ind];}
  double Residual(const int &ind, const int &a) const {return residual[ind][a];}

  void CalcCellDt(const int&, const int&, const int&, const double&);

  void CalcInvFluxI(const idealGas&, const input&);
  void CalcInvFluxJ(const idealGas&, const input&);
  void CalcInvFluxK(const idealGas&, const input&);

  void CalcBlockTimeStep(const input&, const double&);
  void UpdateBlock(const input&, const int&, const idealGas&, const double&, const vector<genArray> &, genArray &, resid &);

  void ExplicitEulerTimeAdvance(const idealGas&, const int&, const int&);
  void ImplicitTimeAdvance(const genArray&, const idealGas&, const int&);
  void RK4TimeAdvance(const primVars&, const idealGas&, const double&, const int&, const int&, const int&);

  void ResetResidWS();
  void CleanResizeVecs();

  vector<genArray> AddVolTime(const vector<genArray>&, const vector<genArray>&, const double &, const double &)const;
  void DeltaNMinusOne( vector<genArray> &, const vector<genArray> &, const idealGas &, const double &, const double &);

  double LUSGS( const vector<vector3d<int> > &, vector<genArray> &, const vector<genArray> &, const vector<genArray> &, 
		const idealGas&, const input&, const sutherland&)const;

  void CalcViscFluxI(const sutherland&, const idealGas&, const input&);
  void CalcViscFluxJ(const sutherland&, const idealGas&, const input&);
  void CalcViscFluxK(const sutherland&, const idealGas&, const input&);

  gradients CalcGradients();

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input&, const idealGas&);
  void AssignInviscidGhostCellsEdge(const input&, const idealGas&);

  void AssignViscousGhostCells(const input&, const idealGas&);
  void AssignViscousGhostCellsEdge(const input&, const idealGas&);

  bool IsPhysical(const int&, const int&, const int&)const;
  bool AtCorner(const int&, const int&, const int&)const;
  bool AtEdge(const int&, const int&, const int&, string&)const;

  geomSlice GetGeomSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;
  vector<bool> PutGeomSlice(const geomSlice&, interblock&, const int&, const int&);

  stateSlice GetStateSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;
  void PutStateSlice(const stateSlice&, const interblock&, const int&, const int&);

  procBlock Split(const string&, const int&, const int&, vector<boundarySurface>&);
  void Join(const procBlock&, const string&, vector<boundarySurface>&);

  void SwapSliceMPI(const interblock&, const int&, const MPI_Datatype& );
  void PackSendGeomMPI(const MPI_Datatype&, const MPI_Datatype&)const;
  void RecvUnpackGeomMPI(const MPI_Datatype&, const MPI_Datatype&);
  void PackSendSolMPI(const MPI_Datatype&)const;
  void RecvUnpackSolMPI(const MPI_Datatype&);

  //destructor
  ~procBlock() {}

};

class geomSlice {

  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;      //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;      //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;      //coordinates of k-face centers

  vector<double> vol ;                     //cell volume

  int numCells;                            //number of cells in block
  int numI;                                //i-dimension of block (cells)
  int numJ;                                //j-dimension of block (cells)
  int numK;                                //k-dimension of block (cells)
  int parBlock;                            //parent block number

 public:
  //constructors
  geomSlice();
  geomSlice( const int&, const int&, const int&, const int& );

  friend geomSlice procBlock::GetGeomSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;

  //member functions
  int NumCells() const {return numCells;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}
  int ParentBlock() const {return parBlock;}

  double Vol(const int &ind) const {return vol[ind];}
  vector3d<double> Center(const int &ind) const {return center[ind];}
  vector3d<double> FAreaI(const int &ind) const {return fAreaI[ind];}
  vector3d<double> FAreaJ(const int &ind) const {return fAreaJ[ind];}
  vector3d<double> FAreaK(const int &ind) const {return fAreaK[ind];}
  vector3d<double> FCenterI(const int &ind) const {return fCenterI[ind];}
  vector3d<double> FCenterJ(const int &ind) const {return fCenterJ[ind];}
  vector3d<double> FCenterK(const int &ind) const {return fCenterK[ind];}

  //destructor
  ~geomSlice() {}

};

class stateSlice {
  vector<primVars> state ;                     //cell states

  int numCells;                            //number of cells in block
  int numI;                                //i-dimension of block (cells)
  int numJ;                                //j-dimension of block (cells)
  int numK;                                //k-dimension of block (cells)
  int parBlock;                            //parent block number

 public:
  //constructors
  stateSlice();
  stateSlice( const int&, const int&, const int&, const int& );

  friend stateSlice procBlock::GetStateSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;

  //member functions
  int NumCells() const {return numCells;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}
  int ParentBlock() const {return parBlock;}

  primVars State(const int &ind) const {return state[ind];}

  void PackSwapUnpackMPI( const interblock&, const MPI_Datatype&, const int& );

  //destructor
  ~stateSlice() {}

};


//function definitions
double CellSpectralRadius(const vector3d<double> &, const vector3d<double> &, const primVars&, const idealGas&);
double ViscCellSpectralRadius(const vector3d<double>&, const vector3d<double>&, const primVars&, const idealGas&, const sutherland&, const double&, const double&);

template<class T>
T FaceReconCentral(const T&, const T&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&);

template<class T>
vector<T> PadWithGhosts(const vector<T>&, const int&, const int&, const int&, const int&);

tensor<double> CalcVelGradGG(const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
			     const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
			     const vector3d<double>&, const vector3d<double>&, const double&);

vector3d<double> CalcScalarGradGG(const double&, const double&, const double&, const double&, const double&, const double&, const vector3d<double>&, const vector3d<double>&,
				const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const double&);

vector3d<int> GetSwapLoc( const int&, const int&, const int&, const interblock&, const bool&);
void SwapSlice(interblock&, procBlock&, procBlock&, const bool&);

void GetBoundaryConditions(vector<procBlock>&, const input&, const idealGas&, vector<interblock>&, const int &rank, const MPI_Datatype &MPI_cellData);

#endif

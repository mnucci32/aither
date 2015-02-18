#ifndef PROCBLOCKHEADERDEF             //only if the macro PROCBLOCKHEADERDEF is not defined execute these lines of code
#define PROCBLOCKHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "plot3d.h" //plot3d
#include "eos.h" //idealGas
#include "primVars.h" //primVars
#include "inviscidFlux.h" //inviscidFlux
#include "viscousFlux.h" //viscousFlux
#include "input.h" //inputVars
#include "matrix.h" //squareMatrix, matrixDiagonal
#include "boundaryConditions.h" //interblock, patch
#include "mpi.h" //parallelism
#include <fstream>
#include <iostream>
#include "macros.h"

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
  int parBlockStartI;                      //parent block starting index for i
  int parBlockEndI;                        //parent block ending index for i
  int parBlockStartJ;                      //parent block starting index for j
  int parBlockEndJ;                        //parent block ending index for j
  int parBlockStartK;                      //parent block starting index for k
  int parBlockEndK;                        //parent block ending index for k
  int rank;                                //processor rank
  int globalPos;                           //global position of procBlock in decomposed vector of procBlocks

 public:
  //constructors
  procBlock();
  procBlock( const plot3dBlock& , const int&, const int&, const string& );
  procBlock( const double, const double, const vector3d<double>, const plot3dBlock&, const int&, const int&, const string&, const boundaryConditions& );
  procBlock( const primVars&, const plot3dBlock&, const int &, const int&, const string&, const boundaryConditions& );
  procBlock( const int&, const int&, const int&, const int& );

  //member functions
  int NumCells() const {return numCells;}
  int NumVars() const {return numVars;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}
  int NumGhosts() const {return numGhosts;}
  int ParentBlock() const {return parBlock;}
  int ParentBlockStartI() const {return parBlockStartI;}
  int ParentBlockEndI() const {return parBlockEndI;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  int ParentBlockStartK() const {return parBlockStartK;}
  int ParentBlockEndK() const {return parBlockEndK;}
  void SetRank( const int &a){rank = a;}                             //setter needed for parallel decomposition
  int Rank() const {return rank;}
  void SetGlobalPos( const int &a){globalPos = a;}                   //setter needed for parallel decomposition
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

  void AssignGhostCellsGeom(const input&);
  void AssignGhostCellsGeomEdge(const input&);

  void AssignInviscidGhostCells(const input&, const idealGas&);
  void AssignInviscidGhostCellsEdge(const input&, const idealGas&);

  void AssignViscousGhostCells(const input&, const idealGas&);
  void AssignViscousGhostCellsEdge(const input&, const idealGas&);

  bool AtCorner(const int&, const int&, const int&)const;
  bool AtEdge(const int&, const int&, const int&)const;

  geomSlice GetGeomSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;
  void PutGeomSlice(const geomSlice&, const interblock&, const int&);

  stateSlice GetStateSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;
  void PutStateSlice(const stateSlice&, const interblock&, const int&);

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
  int parBlockStartI;                      //parent block starting index for i
  int parBlockEndI;                        //parent block ending index for i
  int parBlockStartJ;                      //parent block starting index for j
  int parBlockEndJ;                        //parent block ending index for j
  int parBlockStartK;                      //parent block starting index for k
  int parBlockEndK;                        //parent block ending index for k

 public:
  //constructors
  geomSlice();
  geomSlice( const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&,
	     const int&, const int& );

  friend geomSlice procBlock::GetGeomSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;

  //member functions
  int NumCells() const {return numCells;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}

  int ParentBlock() const {return parBlock;}
  int ParentBlockStartI() const {return parBlockStartI;}
  int ParentBlockEndI() const {return parBlockEndI;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  int ParentBlockStartK() const {return parBlockStartK;}
  int ParentBlockEndK() const {return parBlockEndK;}

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
  int parBlockStartI;                      //parent block starting index for i
  int parBlockEndI;                        //parent block ending index for i
  int parBlockStartJ;                      //parent block starting index for j
  int parBlockEndJ;                        //parent block ending index for j
  int parBlockStartK;                      //parent block starting index for k
  int parBlockEndK;                        //parent block ending index for k

 public:
  //constructors
  stateSlice();
  stateSlice( const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&,
	     const int&, const int& );

  friend stateSlice procBlock::GetStateSlice(const int&, const int&, const int&, const int&, const int&, const int&, const bool=false, const bool=false, const bool=false)const;

  //member functions
  int NumCells() const {return numCells;}
  int NumI() const {return numI;}
  int NumJ() const {return numJ;}
  int NumK() const {return numK;}

  int ParentBlock() const {return parBlock;}
  int ParentBlockStartI() const {return parBlockStartI;}
  int ParentBlockEndI() const {return parBlockEndI;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  int ParentBlockStartK() const {return parBlockStartK;}
  int ParentBlockEndK() const {return parBlockEndK;}

  primVars State(const int &ind) const {return state[ind];}

  void PackSwapUnpackMPI( const interblock&, const MPI_Datatype&, const int& );

  //destructor
  ~stateSlice() {}

};


//function definitions
double CellSpectralRadius(const vector3d<double> &, const vector3d<double> &, const primVars&, const idealGas&);
double ViscCellSpectralRadius(const vector3d<double>&, const vector3d<double>&, const primVars&, const idealGas&, const sutherland&, const double&);

template<class T>
T FaceReconCentral(const T&, const T&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&);

template<class T>
vector<T> PadWithGhosts(const vector<T>&, const int&, const int&, const int&, const int&);

tensor<double> CalcVelGradGG(const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
			     const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
			     const vector3d<double>&, const vector3d<double>&, const double&);

vector3d<double> CalcTempGradGG(const double&, const double&, const double&, const double&, const double&, const double&, const vector3d<double>&, const vector3d<double>&,
				const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const double&);

vector<int> GetSwapLoc( const int&, const int&, const int&, const interblock&, const bool&);
void SwapSlice(const interblock&, procBlock&, procBlock&, const bool&);

void GetBoundaryConditions(vector<procBlock>&, const input&, const idealGas&, const vector<interblock>&, const int &rank, const MPI_Datatype &MPI_cellData);

#endif

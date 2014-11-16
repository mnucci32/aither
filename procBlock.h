#ifndef PROCBLOCKHEADERDEF             //only if the macro PROCBLOCKHEADERDEF is not defined execute these lines of code
#define PROCBLOCKHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "plot3d.h" //plot3d
#include "eos.h"
#include "primVars.h" //primVars
#include "inviscidFlux.h" //inviscidFlux
#include "viscousFlux.h" //viscousFlux
#include "input.h" //inputVars
#include "matrix.h" //squareMatrix, matrixDiagonal
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

class procBlock {
  vector<primVars> state;            //primative variables at cell center

  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;        //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;        //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;        //coordinates of k-face centers

  vector<colMatrix> residual;               //cell residual

  vector<double> vol ;                      //cell volume
  vector<double> avgWaveSpeed;             //maximum wave speed normal to i-faces
  vector<double> dt;                        //cell time step

  vector<tensor<double> > velGrad;          //velocity gradient on i-faces
  vector<vector3d<double> > tempGrad;       //temperature gradient on i-faces

  int numCells;
  int numVars;
  int numI;
  int numJ;
  int numK;
  int numGhosts;
  int parBlock;
  int parBlockStartI;
  int parBlockEndI;
  int parBlockStartJ;
  int parBlockEndJ;
  int parBlockStartK;
  int parBlockEndK;

 public:
  //constructors
  procBlock();
  procBlock( const plot3dBlock& , const int&, const int&, const string& );
  procBlock( const double, const double, const vector3d<double>, const plot3dBlock&, const int&, const int&, const string& );
  procBlock( const primVars&, const plot3dBlock&, const int &, const int&, const string& );

  //member functions
  void SetNumCells( const int &a){numCells = a;}
  int NumCells() const {return numCells;}
  void SetI( const int &a){numI = a;}
  int NumI() const {return numI;}
  void SetJ( const int &a){numJ = a;}
  int NumJ() const {return numJ;}
  void SetK( const int &a){numK = a;}
  int NumK() const {return numK;}
  void SetNumGhosts( const int &a){numGhosts = a;}
  int NumGhosts() const {return numGhosts;}
  void SetParentBlock( const int &a){parBlock = a;}
  int ParentBlock() const {return parBlock;}
  void SetParentBlockStartI( const int &a){parBlockStartI = a;}
  int ParentBlockStartI() const {return parBlockStartI;}
  void SetParentBlockEndI( const int &a){parBlockEndI = a;}
  int ParentBlockEndI() const {return parBlockEndI;}
  void SetParentBlockStartJ( const int &a){parBlockStartJ = a;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  void SetParentBlockEndJ( const int &a){parBlockEndJ = a;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  void SetParentBlockStartK( const int &a){parBlockStartK = a;}
  int ParentBlockStartK() const {return parBlockStartK;}
  void SetParentBlockEndK( const int &a){parBlockEndK = a;}
  int ParentBlockEndK() const {return parBlockEndK;}

  void SetState( const primVars &a, const int &ind){state[ind] = a;}
  primVars State(const int &ind) const {return state[ind];}
  vector<primVars> GetCopyState() const {return state;}
  vector<colMatrix> GetCopyConsVars(const idealGas &) const; 
  const vector<primVars> & GetRefState() const {return state;}

  void SetVol( const double &a, const int &ind){vol[ind] = a;}
  double Vol(const int &ind) const {return vol[ind];}
  void SetCenter( const vector3d<double> &a, const int &ind){center[ind] = a;}
  vector3d<double> Center(const int &ind) const {return center[ind];}

  void SetFAreaI( const vector3d<double> &a, const int &ind){fAreaI[ind] = a;}
  vector3d<double> FAreaI(const int &ind) const {return fAreaI[ind];}
  void SetFAreaJ( const vector3d<double> &a, const int &ind){fAreaJ[ind] = a;}
  vector3d<double> FAreaJ(const int &ind) const {return fAreaJ[ind];}
  void SetFAreaK( const vector3d<double> &a, const int &ind){fAreaK[ind] = a;}
  vector3d<double> FAreaK(const int &ind) const {return fAreaK[ind];}
  void SetFCenterI( const vector3d<double> &a, const int &ind){fCenterI[ind] = a;}
  vector3d<double> FCenterI(const int &ind) const {return fCenterI[ind];}
  void SetFCenterJ( const vector3d<double> &a, const int &ind){fCenterJ[ind] = a;}
  vector3d<double> FCenterJ(const int &ind) const {return fCenterJ[ind];}
  void SetFCenterK( const vector3d<double> &a, const int &ind){fCenterK[ind] = a;}
  vector3d<double> FCenterK(const int &ind) const {return fCenterK[ind];}

  void SetAvgWaveSpeed( const double &a, const int &ind){avgWaveSpeed[ind] = a;}
  double AvgWaveSpeed(const int &ind) const {return avgWaveSpeed[ind];}

  void SetDt( const double &a, const int &ind){dt[ind] = a;}
  double Dt(const int &ind) const {return dt[ind];}

  void SetResidual( const colMatrix &a, const int &ind){residual[ind] = a;}
  void AddToResidual( const inviscidFlux &, const int &);
  void AddToResidual( const viscousFlux &, const int &);

  const vector<colMatrix> & Residual() const {return residual;}
  colMatrix Residual(const int &ind) const {return residual[ind];}
  double Residual(const int &ind, const int &a) const {return residual[ind].Data(a);}

  void CalcCellDt(const int&, const int&, const int&, const double&);

  void CalcInvFluxI(const idealGas&, const input&);
  void CalcInvFluxJ(const idealGas&, const input&);
  void CalcInvFluxK(const idealGas&, const input&);

  void CalcBlockTimeStep(const input&, const double&);
  void UpdateBlock(const input&, const int&, const idealGas&, const double&, const int&, const vector<colMatrix> &, colMatrix &, colMatrix &, int &);

  void ExplicitEulerTimeAdvance(const idealGas&, const int&, const int&);
  void ImplicitTimeAdvance(const colMatrix&, const idealGas&, const int&);
  void RK4TimeAdvance(const primVars&, const idealGas&, const double&, const int&, const int&, const int&);

  void ResetResidWS();

  vector<colMatrix> AddVolTime(const vector<colMatrix>&, const vector<colMatrix>&, const double &, const double &)const;
  void DeltaNMinusOne( vector<colMatrix> &, const vector<colMatrix> &, const idealGas &, const double &, const double &);

  double LUSGS( const vector<vector3d<int> > &, vector<colMatrix> &, const vector<colMatrix> &, const vector<colMatrix> &, const idealGas&, const input&, const sutherland&)const;

  void SetVelGrad( const tensor<double> &a, const int &ind){velGrad[ind] = a;}
  tensor<double> VelGrad(const int &ind) const {return velGrad[ind];}
  void SetTempGrad( const vector3d<double> &a, const int &ind){tempGrad[ind] = a;}
  vector3d<double> TempGrad(const int &ind) const {return tempGrad[ind];}

  void CalcVelGradGG(const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const double&, const int&);
  void CalcTempGradGG(const double&, const double&, const vector3d<double>&, const vector3d<double>&, const double&, const int&);

  void InitializeGrads();

  void CalcCellGradsI(const idealGas&, const sutherland&, const input&);
  void CalcCellGradsJ(const idealGas&, const sutherland&, const input&);
  void CalcCellGradsK(const idealGas&, const sutherland&, const input&);

  void CalcViscFluxI(const sutherland&, const idealGas&, const input&);
  void CalcViscFluxJ(const sutherland&, const idealGas&, const input&);
  void CalcViscFluxK(const sutherland&, const idealGas&, const input&);

  void AssignGhostCellsGeom();
  void AssignGhostCellsGeomEdge();

  void AssignInviscidGhostCells(const input&, const idealGas&);
  void AssignInviscidGhostCellsEdge(const input&, const idealGas&);

  void AssignViscousGhostCells(const input&, const idealGas&);
  void AssignViscousGradGhostCells(const input&);

  //destructor
  ~procBlock() {}

};

//function definitions
double CellSpectralRadius(const vector3d<double> &, const vector3d<double> &, const primVars&, const idealGas&);
double ViscCellSpectralRadius(const vector3d<double>&, const vector3d<double>&, const primVars&, const idealGas&, const sutherland&, const double&);

template<class T>
T FaceReconCentral(const T&, const T&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&);

template<class T>
vector<T> PadWithGhosts(const vector<T>&, const int&, const int&, const int&, const int&);

vector<primVars> PadStateWithGhosts(const primVars&, const int&, const int&, const int&, const int&);

#endif

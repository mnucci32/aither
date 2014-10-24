#ifndef BLOCKVARSHEADERDEF             //only if the macro BLOCKVARSHEADERDEF is not defined execute these lines of code
#define BLOCKVARSHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
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

class blockVars {
  int length;
  int numI;
  int numJ;
  int numK;
  vector<primVars> state;            //conserved variables at cell center

  vector<double> vol ;                      //cell volume
  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;        //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;        //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;        //coordinates of k-face centers

  vector<double> avgWaveSpeed;             //maximum wave speed normal to i-faces

  vector<double> dt;                        //cell time step

  vector<colMatrix> residual;               //cell residual

 public:
  //constructors
  blockVars();
  blockVars( const plot3dBlock& );
  blockVars( const double, const double, const vector3d<double>, const plot3dBlock& );
  blockVars( const primVars&, const plot3dBlock& );

  //member functions
  void SetLen( const int &a){length = a;}
  int Len() const {return length;}
  void SetI( const int &a){numI = a;}
  int NumI() const {return numI;}
  void SetJ( const int &a){numJ = a;}
  int NumJ() const {return numJ;}
  void SetK( const int &a){numK = a;}
  int NumK() const {return numK;}

  void SetState( const primVars &a, const int &ind){state[ind] = a;}
  primVars State(const int &ind) const {return state[ind];}
  vector<primVars> GetCopyState() const {return state;}
  vector<colMatrix> GetCopyConsVars(const idealGas &) const; 
  const vector<primVars> & GetRefState() const {return state;}

  void SetVol( const double &a, const int &ind){vol[ind] = a;}
  double Vol(const int &ind) const {return vol[ind];}
  void SetCenter( const vector<vector3d<double> > &a){center = a;}
  vector3d<double> Center(const int &ind) const {return center[ind];}

  void SetFAreaI( const vector<vector3d<double> > &a){fAreaI = a;}
  vector3d<double> FAreaI(const int &ind) const {return fAreaI[ind];}
  void SetFAreaJ( const vector<vector3d<double> > &a){fAreaJ = a;}
  vector3d<double> FAreaJ(const int &ind) const {return fAreaJ[ind];}
  void SetFAreaK( const vector<vector3d<double> > &a){fAreaK = a;}
  vector3d<double> FAreaK(const int &ind) const {return fAreaK[ind];}
  void SetFCenterI( const vector<vector3d<double> > &a){fCenterI = a;}
  vector3d<double> FCenterI(const int &ind) const {return fCenterI[ind];}
  void SetFCenterJ( const vector<vector3d<double> > &a){fCenterJ = a;}
  vector3d<double> FCenterJ(const int &ind) const {return fCenterJ[ind];}
  void SetFCenterK( const vector<vector3d<double> > &a){fCenterK = a;}
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

  void CalcInvFluxI(const idealGas&, const input&, const int&);
  void CalcInvFluxJ(const idealGas&, const input&, const int&);
  void CalcInvFluxK(const idealGas&, const input&, const int&);
  void CalcBlockTimeStep(const input&, const double&);
  void UpdateBlock(const input&, const int&, const idealGas&, const double&, const int&, const vector<colMatrix> &, vector<double> &, vector<double> &, int &);

  void ExplicitEulerTimeAdvance(const idealGas&, const int&);
  void ImplicitTimeAdvance(const colMatrix&, const idealGas&, const int&);
  void RK4TimeAdvance(const primVars&, const idealGas&, const double&, const int&, const int&);

  //void TotalResidual( vector<double> &, vector<double> &, int &, const int & );

  /* void CalcInvFluxJacI(const idealGas&, const input&, const int&, colMatrix&)const; */
  /* void CalcInvFluxJacJ(const idealGas&, const input&, const int&, colMatrix&)const; */
  /* void CalcInvFluxJacK(const idealGas&, const input&, const int&, colMatrix&)const; */

  void ResetResidWS();

  void AddVolTime(colMatrix&, const double &, const double &, const double &)const;
  vector<colMatrix> AddVolTime(const vector<colMatrix>&, const vector<colMatrix>&, const double &, const double &)const;

  void DeltaNMinusOne( vector<colMatrix> &, const vector<colMatrix> &, const idealGas &, const double &, const double &);

  double LUSGS( const vector<vector3d<int> > &, vector<colMatrix> &, const vector<colMatrix> &, const vector<colMatrix> &, const idealGas&, const input&, const sutherland&)const;



  //destructor
  ~blockVars() {}

};

//function definitions
//double ViscFaceSpecRadTSL(const primVars&, const idealGas&, const sutherland&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&);
double CellSpectralRadius(const vector3d<double> &, const vector3d<double> &, const primVars&, const idealGas&);
double ViscCellSpectralRadius(const vector3d<double>&, const vector3d<double>&, const primVars&, const idealGas&, const sutherland&, const double&);

#endif

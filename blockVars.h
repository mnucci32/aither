#ifndef BLOCKVARSHEADERDEF             //only if the macro BLOCKVARSHEADERDEF is not defined execute these lines of code
#define BLOCKVARSHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "plot3d.h" //plot3d
#include "eos.h"
#include "primVars.h" //primVars
#include "inviscidFlux.h" //inviscidFlux
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
  vector<inviscidFlux> invFluxI;            //inviscid flux on i-faces
  vector<inviscidFlux> invFluxJ;            //inviscid flux on j-faces
  vector<inviscidFlux> invFluxK;            //inviscid flux on k-faces

  vector<double> vol ;                      //cell volume
  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;        //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;        //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;        //coordinates of k-face centers

  vector<double> maxWaveSpeedI;             //maximum wave speed normal to i-faces
  vector<double> maxWaveSpeedJ;             //maximum wave speed normal to i-faces
  vector<double> maxWaveSpeedK;             //maximum wave speed normal to i-faces

  vector<double> dt;                        //cell time step

  vector<vector<double> > residual;         //cell residual

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
  void SetInvFluxI( const inviscidFlux &a, const int &ind){invFluxI[ind] = a;}
  inviscidFlux InvFluxI(const int &ind) const {return invFluxI[ind];}
  void SetInvFluxJ( const inviscidFlux &a, const int &ind){invFluxJ[ind] = a;}
  inviscidFlux InvFluxJ(const int &ind) const {return invFluxJ[ind];}
  void SetInvFluxK( const inviscidFlux &a, const int &ind){invFluxK[ind] = a;}
  inviscidFlux InvFluxK(const int &ind) const {return invFluxK[ind];}

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

  void SetMaxWaveSpeedI( const double &a, const int &ind){maxWaveSpeedI[ind] = a;}
  double MaxWaveSpeedI(const int &ind) const {return maxWaveSpeedI[ind];}
  void SetMaxWaveSpeedJ( const double &a, const int &ind){maxWaveSpeedJ[ind] = a;}
  double MaxWaveSpeedJ(const int &ind) const {return maxWaveSpeedJ[ind];}
  void SetMaxWaveSpeedK( const double &a, const int &ind){maxWaveSpeedK[ind] = a;}
  double MaxWaveSpeedK(const int &ind) const {return maxWaveSpeedK[ind];}

  void SetDt( const double &a, const int &ind){dt[ind] = a;}
  double Dt(const int &ind) const {return dt[ind];}

  void SetResidual( const vector<double> &a, const int &ind){residual[ind] = a;}
  vector<double> Residual(const int &ind) const {return residual[ind];}

  void CalcCellDt(const int&, const int&, const int&, const double&);
  void CalcCellResidual(const int&, const int&, const int&, const int&, const int&);

  void CalcInvFluxI(const idealGas&, const input&, const int&);
  void CalcInvFluxJ(const idealGas&, const input&, const int&);
  void CalcInvFluxK(const idealGas&, const input&, const int&);
  void CalcBlockResidDT(const input&, const double&);
  void UpdateBlock(const input&, const idealGas&, const double&, const int&, vector<double> &, vector<double> &, int &);

  void ExplicitEulerTimeAdvance(const idealGas&, const int&);
  void RK4TimeAdvance(const primVars&, const idealGas&, const double&, const int&, const int&);

  //void TotalResidual( vector<double> &, vector<double> &, int &, const int & );

  void CalcInvFluxJacI(const idealGas&, const input&, const int&, matrixDiagonal&, matrixDiagonal&, matrixDiagonal&)const;
  void CalcInvFluxJacJ(const idealGas&, const input&, const int&, matrixDiagonal&, matrixDiagonal&, matrixDiagonal&)const;
  void CalcInvFluxJacK(const idealGas&, const input&, const int&, matrixDiagonal&, matrixDiagonal&, matrixDiagonal&)const;



  //destructor
  ~blockVars() {}

};

//function definitions


#endif

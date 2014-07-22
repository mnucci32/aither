#ifndef VISCBLOCKVARSHEADERDEF             //only if the macro VISCBLOCKVARSHEADERDEF is not defined execute these lines of code
#define VISCBLOCKVARSHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "plot3d.h" //plot3d
#include "viscousFlux.h" //viscousFlux
#include "blockVars.h" //blockVars
#include "eos.h" //idealGas
#include "input.h" //inputVars
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

class viscBlockVars {
  int length;
  int numI;
  int numJ;
  int numK;

  vector<viscousFlux> viscFluxI;            //viscous flux on i-faces
  vector<viscousFlux> viscFluxJ;            //viscous flux on j-faces
  vector<viscousFlux> viscFluxK;            //viscous flux on k-faces

  vector<tensor<double> > velGrad;          //velocity gradient on i-faces

  vector<vector3d<double> > tempGrad;       //temperature gradient on i-faces

 public:
  //constructors
  viscBlockVars();
  viscBlockVars( const plot3dBlock& );

  //member functions
  void SetLen( const int &a){length = a;}
  int Len() const {return length;}
  void SetI( const int &a){numI = a;}
  int NumI() const {return numI;}
  void SetJ( const int &a){numJ = a;}
  int NumJ() const {return numJ;}
  void SetK( const int &a){numK = a;}
  int NumK() const {return numK;}

  void SetViscFluxI( const viscousFlux &a, const int &ind){viscFluxI[ind] = a;}
  viscousFlux ViscFluxI(const int &ind) const {return viscFluxI[ind];}
  void SetViscFluxJ( const viscousFlux &a, const int &ind){viscFluxJ[ind] = a;}
  viscousFlux ViscFluxJ(const int &ind) const {return viscFluxJ[ind];}
  void SetViscFluxK( const viscousFlux &a, const int &ind){viscFluxK[ind] = a;}
  viscousFlux ViscFluxK(const int &ind) const {return viscFluxK[ind];}

  void SetVelGrad( const tensor<double> &a, const int &ind){velGrad[ind] = a;}
  tensor<double> VelGrad(const int &ind) const {return velGrad[ind];}

  void SetTempGrad( const vector3d<double> &a, const int &ind){tempGrad[ind] = a;}
  vector3d<double> TempGrad(const int &ind) const {return tempGrad[ind];}

  void CalcVelGradGG(const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, 
		     const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
		     const double&, const int&);

  void CalcTempGradGG(const double&, const double&, const double&, const double&, const double&, const double&, 
		      const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&,
		      const double&, const int&);

  template<class T>
  T FaceReconCentral(const T&, const T&, const vector3d<double>&, const vector3d<double>&, const vector3d<double>&)const;

  void CalcCellGrads(const blockVars&, const idealGas&, const input&, const int&);

  void CalcViscFluxI(const blockVars&, const sutherland&, const idealGas&, const input&, const int&);
  void CalcViscFluxJ(const blockVars&, const sutherland&, const idealGas&, const input&, const int&);
  void CalcViscFluxK(const blockVars&, const sutherland&, const idealGas&, const input&, const int&);


  void CalcCellResidual(blockVars&, const int&, const int&, const int&, const int&, const int&)const;

  void CalcBlockResidDT(blockVars&, const input&, const double&);

  //destructor
  ~viscBlockVars() {}

};

//function definitions


#endif

#ifndef VISCFLUXHEADERDEF             //only if the macro VISCFLUXHEADERDEF is not defined execute these lines of code
#define VISCFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "eos.h"  //idealGas
#include "primVars.h" //primVars
#include "matrix.h" //squareMatrix
#include "input.h" //input
#include <iostream> //cout

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class viscousFlux {
  double data[4];           //viscous flux for x-momentum equation
                            //viscous flux for y-momentum equation
                            //viscous flux for z-momentum equation
                            //viscous flux for energy equation

 public:
  //constructors
  viscousFlux() : data{0.0, 0.0, 0.0, 0.0} {}
  viscousFlux( const tensor<double>&, const vector3d<double>&, const double&, const sutherland&, const idealGas&, const vector3d<double>&, const vector3d<double>& );

  //member functions
  void SetMomX( const double &a){data[0] = a;}
  double MomX() const {return data[0];}
  void SetMomY( const double &a){data[1] = a;}
  double MomY() const {return data[1];}
  void SetMomZ( const double &a){data[2] = a;}
  double MomZ() const {return data[2];}
  void SetEngy( const double &a){data[3] = a;}
  double Engy() const {return data[3];}

  void SetFlux( const tensor<double>&, const vector3d<double>&, const double&, const sutherland&, const idealGas&, const vector3d<double>&, const vector3d<double>& );

  viscousFlux operator * (const double&);
  viscousFlux operator / (const double&);

  friend ostream & operator<< (ostream &os, viscousFlux&);
  friend viscousFlux operator * (const double&, const viscousFlux&);
  friend viscousFlux operator / (const double&, const viscousFlux&);

  //destructor
  ~viscousFlux() {}

};

//function definitions
void CalcTSLFluxJac( const double&, const idealGas&, const vector3d<double>&, const primVars&, const primVars&, const double&, squareMatrix&, squareMatrix&, const sutherland&);

tensor<double> CalcVelGradTSL(const primVars&, const primVars&, const vector3d<double>&, const double&);

#endif

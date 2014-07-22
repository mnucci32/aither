#ifndef VISCFLUXHEADERDEF             //only if the macro VISCFLUXHEADERDEF is not defined execute these lines of code
#define VISCFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "eos.h"  //idealGas
#include "primVars.h" //primVars
#include "input.h" //input
#include "boundaryConditions.h" //boundaryConditions
#include <iostream> //cout

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class viscousFlux {
  double momX;           //viscous flux for x-momentum equation
  double momY;           //viscous flux for y-momentum equation
  double momZ;           //viscous flux for z-momentum equation
  double engy;           //viscous flux for energy equation

 public:
  //constructors
  viscousFlux() : momX(0.0), momY(0.0), momZ(0.0), engy(0.0) {}
  viscousFlux( const tensor<double>&, const vector3d<double>&, const double&, const sutherland&, const idealGas&, const vector3d<double>&, const vector3d<double>& );

  //member functions
  void SetMomX( const double &a){momX = a;}
  double MomX() const {return momX;}
  void SetMomY( const double &a){momY = a;}
  double MomY() const {return momY;}
  void SetMomZ( const double &a){momZ = a;}
  double MomZ() const {return momZ;}
  void SetEngy( const double &a){engy = a;}
  double Engy() const {return engy;}

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


#endif

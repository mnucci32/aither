#ifndef INVFLUXHEADERDEF             //only if the macro INVFLUXHEADERDEF is not defined execute these lines of code
#define INVFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "eos.h"  //idealGas
#include "primVars.h" //primVars
#include "input.h" //input
//#include "boundaryConditions.h" //boundaryConditions
#include <iostream> //cout
#include "matrix.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class inviscidFlux {
  double data[5];
  /* double rhoVel;            //rho dot velocity vector */
  /* double rhoVelU;           //rho dot velocity vector * u-velocity + pressure * i-dir-vector */
  /* double rhoVelV;           //rho dot velocity vector * v-velocity + pressure * j-dir-vector */
  /* double rhoVelW;           //rho dot velocity vector * w-velocity + pressure * k-dir-vector */
  /* double rhoVelH;           //rho dot velocity vector * enthalpy */

 public:
  //constructors
  inviscidFlux() : data{0.0, 0.0, 0.0, 0.0, 0.0} {}
  inviscidFlux( const primVars&, const idealGas&, const vector3d<double>& );
  inviscidFlux( const colMatrix&, const idealGas&, const vector3d<double>& );

  //member functions
  void SetFlux( const primVars&, const idealGas&, const vector3d<double>&);

  void SetRhoVel( const double &a){data[0] = a;}
  double RhoVel() const {return data[0];}
  void SetRhoVelU( const double &a){data[1] = a;}
  double RhoVelU() const {return data[1];}
  void SetRhoVelV( const double &a){data[2] = a;}
  double RhoVelV() const {return data[2];}
  void SetRhoVelW( const double &a){data[3] = a;}
  double RhoVelW() const {return data[3];}
  void SetRhoVelH( const double &a){data[4] = a;}
  double RhoVelH() const {return data[4];}

  inviscidFlux operator * (const double&);
  inviscidFlux operator / (const double&);

  inviscidFlux operator + (const inviscidFlux&)const;
  inviscidFlux operator - (const inviscidFlux&)const;

  friend ostream & operator<< (ostream &os, inviscidFlux&);
  friend inviscidFlux operator * (const double&, const inviscidFlux&);
  friend inviscidFlux operator / (const double&, const inviscidFlux&);

  colMatrix ConvertToColMatrix()const;

  //destructor
  ~inviscidFlux() {}

};

//function definitions
//inviscidFlux BoundaryFlux( const string&, const vector3d<double>&, const primVars&, const primVars&, const idealGas&, const input&, const string&, double&, const double=0.5, const double=1.0 );
//function to calculate Roe flux with entropy fix
inviscidFlux RoeFlux( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double& );                  
//inviscidFlux LaxFriedrichsFlux( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double& );                  
//function to calculate Roe flux with entropy fix for implicit methods
void ApproxRoeFluxJacobian( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double&, squareMatrix&, squareMatrix& );                  
//void LaxFriedrichsFluxJacobian( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double&, double&, squareMatrix&, squareMatrix& );                  

//double BoundaryInvSpecRad( const string&, const vector3d<double>&, const primVars&, const idealGas&, const string&, const input&);

colMatrix ConvectiveFluxUpdate( const primVars&, const idealGas &, const vector3d<double> &, const colMatrix&);
//double ConvSpecRad( const vector3d<double> &, const primVars&, const primVars&, const idealGas&);


#endif

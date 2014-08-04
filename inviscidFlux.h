#ifndef INVFLUXHEADERDEF             //only if the macro INVFLUXHEADERDEF is not defined execute these lines of code
#define INVFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "eos.h"  //idealGas
#include "primVars.h" //primVars
#include "input.h" //input
#include "boundaryConditions.h" //boundaryConditions
#include <iostream> //cout
#include "matrix.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class inviscidFlux {
  double rhoVel;            //rho dot velocity vector
  double rhoVelU;           //rho dot velocity vector * u-velocity + pressure * i-dir-vector
  double rhoVelV;           //rho dot velocity vector * v-velocity + pressure * j-dir-vector
  double rhoVelW;           //rho dot velocity vector * w-velocity + pressure * k-dir-vector
  double rhoVelH;           //rho dot velocity vector * enthalpy

 public:
  //constructors
  inviscidFlux() : rhoVel(0.0), rhoVelU(0.0), rhoVelV(0.0), rhoVelW(0.0), rhoVelH(0.0) {}
  inviscidFlux( const primVars&, const idealGas&, const vector3d<double>& );

  //member functions
  void SetRhoVel( const double &a){rhoVel = a;}
  double RhoVel() const {return rhoVel;}
  void SetRhoVelU( const double &a){rhoVelU = a;}
  double RhoVelU() const {return rhoVelU;}
  void SetRhoVelV( const double &a){rhoVelV = a;}
  double RhoVelV() const {return rhoVelV;}
  void SetRhoVelW( const double &a){rhoVelW = a;}
  double RhoVelW() const {return rhoVelW;}
  void SetRhoVelH( const double &a){rhoVelH = a;}
  double RhoVelH() const {return rhoVelH;}

  inviscidFlux operator * (const double&);
  inviscidFlux operator / (const double&);

  friend ostream & operator<< (ostream &os, inviscidFlux&);
  friend inviscidFlux operator * (const double&, const inviscidFlux&);
  friend inviscidFlux operator / (const double&, const inviscidFlux&);

  //destructor
  ~inviscidFlux() {}

};

//function definitions
//function to update the state variables given the flux and area vector magnitude at each face
/* vector<double> CalculateResidual(inviscidFlux const &, inviscidFlux const&, inviscidFlux const&, inviscidFlux const &, inviscidFlux const&, inviscidFlux const&,                  */
/*         	                 double const&, double const&, double const&, double const&, double const&, double const& ); */
inviscidFlux BoundaryFlux( const string&, const vector3d<double>&, const primVars&, const primVars&, const idealGas&, const input&, const string&, double&, const double=0.5, const double=1.0 );
//function to calculate Roe flux with entropy fix
inviscidFlux RoeFlux( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double& );                  
//function to calculate Roe flux with entropy fix for implicit methods
void RoeFluxJacobian( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double&, squareMatrix&, squareMatrix& );                  

#endif

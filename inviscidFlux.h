#ifndef INVFLUXHEADERDEF             //only if the macro INVFLUXHEADERDEF is not defined execute these lines of code
#define INVFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "eos.h"  //idealGas
#include "primVars.h" //primVars
#include "input.h" //input
#include <iostream> //cout
#include "matrix.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class inviscidFlux {
  double data[NUMVARS];           //rho dot velocity vector
                                  //rho dot velocity vector * u-velocity + pressure * i-dir-vector
                                  //rho dot velocity vector * v-velocity + pressure * j-dir-vector
                                  //rho dot velocity vector * w-velocity + pressure * k-dir-vector
                                  //rho dot velocity vector * enthalpy

 public:
  //constructors
  inviscidFlux() : data{0.0, 0.0, 0.0, 0.0, 0.0} {}
  inviscidFlux( const primVars&, const idealGas&, const vector3d<double>& );
  inviscidFlux( const genArray&, const idealGas&, const vector3d<double>& );

  //member functions
  double RhoVel() const {return data[0];}
  double RhoVelU() const {return data[1];}
  double RhoVelV() const {return data[2];}
  double RhoVelW() const {return data[3];}
  double RhoVelH() const {return data[4];}

  void RoeFlux( const inviscidFlux&, const double(&) [NUMVARS] );

  inviscidFlux operator * (const double&);
  inviscidFlux operator / (const double&);

  inviscidFlux operator + (const inviscidFlux&)const;
  inviscidFlux operator - (const inviscidFlux&)const;

  friend ostream & operator<< (ostream &os, inviscidFlux&);
  friend inviscidFlux operator * (const double&, const inviscidFlux&);
  friend inviscidFlux operator / (const double&, const inviscidFlux&);

  genArray ConvertToGenArray()const;

  //destructor
  ~inviscidFlux() {}

};

//function definitions
//function to calculate Roe flux with entropy fix
inviscidFlux RoeFlux( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double& );                  

//function to calculate Roe flux with entropy fix for implicit methods
void ApproxRoeFluxJacobian( const primVars&, const primVars&, const idealGas&, const vector3d<double>&, double&, squareMatrix&, squareMatrix& );                  

genArray ConvectiveFluxUpdate( const primVars&, const idealGas &, const vector3d<double> &, const genArray&);


#endif

#include "viscousFlux.h"
#include <cmath> //sqrt

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;

//constructors
// viscousFlux::viscousFlux(){
//   momX = 0.0;     
//   momY = 0.0;     
//   momZ = 0.0;
//   engy = 0.0;     
// }
//constructor -- initialize flux from velocity gradient
viscousFlux::viscousFlux( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec){
  vector3d<double> normArea = areaVec / areaVec.Mag();

  double lambda = suth.GetLambda(mu);

  tensor<double> sumVelGrad = velGrad + velGrad.Transpose();
  double velGradTrace = velGrad.Trace();

  tensor<double> tau = mu * sumVelGrad + lambda * velGradTrace ;
  vector3d<double> tauX(tau.XX(), tau.XY(), tau.XZ());
  vector3d<double> tauY(tau.YX(), tau.YY(), tau.YZ());
  vector3d<double> tauZ(tau.ZX(), tau.ZY(), tau.ZZ());

  momX = -1.0 * tauX.DotProd(normArea);     
  momY = -1.0 * tauY.DotProd(normArea);     
  momZ = -1.0 * tauZ.DotProd(normArea);     
  engy = -1.0 * (tauX.DotProd(normArea) * vel.X() + tauY.DotProd(normArea) * vel.Y() + tauZ.DotProd(normArea) * vel.Z()) - (mu/((eqnState.Gamma() - 1.0) * eqnState.GetPrandtl() )) * tGrad.DotProd(normArea);

}

//non-member functions -----------------------------------------------------------------------------------------------------------//
//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, viscousFlux &flux){

  os << "0.0   " << flux.momX << "   " << flux.momY << "   " << flux.momZ << "   " << flux.engy << endl;

  return os;
}

//member function for scalar multiplication
viscousFlux  viscousFlux::operator * (const double &scalar){
  viscousFlux temp = *this;
  temp.momX *= scalar;
  temp.momY *= scalar;
  temp.momZ *= scalar;
  temp.engy *= scalar;
  return temp;
}

//friend function to allow multiplication from either direction
viscousFlux operator* (const double &scalar, const viscousFlux &flux){
  viscousFlux temp;
  temp.SetMomX(flux.MomX() * scalar);
  temp.SetMomY(flux.MomY() * scalar);
  temp.SetMomZ(flux.MomZ() * scalar);
  temp.SetEngy(flux.Engy() * scalar);
  return temp;
}

//member function for scalar division
viscousFlux  viscousFlux::operator / (const double &scalar){
  viscousFlux temp = *this;
  temp.momX /= scalar;
  temp.momY /= scalar;
  temp.momZ /= scalar;
  temp.engy /= scalar;
  return temp;
}

//friend function to allow division from either direction
viscousFlux operator/ (const double &scalar, const viscousFlux &flux){
  viscousFlux temp;
  temp.SetMomX(scalar / flux.MomX());
  temp.SetMomY(scalar / flux.MomY());
  temp.SetMomZ(scalar / flux.MomZ());
  temp.SetEngy(scalar / flux.Engy());
  return temp;
}


//member function to set the viscous flux
void viscousFlux::SetFlux( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec){
  vector3d<double> normArea = areaVec / areaVec.Mag();

  double lambda = suth.GetLambda(mu);

  tensor<double> sumVelGrad = velGrad + velGrad.Transpose();
  double velGradTrace = velGrad.Trace();

  tensor<double> tau = mu * sumVelGrad + lambda * velGradTrace ;
  vector3d<double> tauX(tau.XX(), tau.XY(), tau.XZ());
  vector3d<double> tauY(tau.YX(), tau.YY(), tau.YZ());
  vector3d<double> tauZ(tau.ZX(), tau.ZY(), tau.ZZ());

  momX = -1.0 * tauX.DotProd(normArea);     
  momY = -1.0 * tauY.DotProd(normArea);     
  momZ = -1.0 * tauZ.DotProd(normArea);     
  engy = -1.0 * (tauX.DotProd(normArea) * vel.X() + tauY.DotProd(normArea) * vel.Y() + tauZ.DotProd(normArea) * vel.Z()) - (mu/((eqnState.Gamma() - 1.0) * eqnState.GetPrandtl() )) * tGrad.DotProd(normArea);

}

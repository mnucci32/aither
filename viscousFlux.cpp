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

  //tensor<double> sumVelGrad = velGrad + velGrad.Transpose();
  double velGradTrace = velGrad.Trace();

  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

  //tensor<double> tau = mu * sumVelGrad + lambda * velGradTrace ;
  // tensor<double> tau;
  // tau.SetXX(2.0 * mu * velGrad.XX() + lambda * velGradTrace);
  // tau.SetYY(2.0 * mu * velGrad.YY() + lambda * velGradTrace);
  // tau.SetZZ(2.0 * mu * velGrad.ZZ() + lambda * velGradTrace);

  // tau.SetXY(mu * ( velGrad.XY() + velGrad.YX() ));
  // tau.SetYX(mu * ( velGrad.XY() + velGrad.YX() )); //symmetric tensor
  // tau.SetXZ(mu * ( velGrad.XZ() + velGrad.ZX() ));
  // tau.SetZX(mu * ( velGrad.XZ() + velGrad.ZX() )); //symmetric tensor
  // tau.SetYZ(mu * ( velGrad.YZ() + velGrad.ZY() ));
  // tau.SetZY(mu * ( velGrad.YZ() + velGrad.ZY() )); //symmetric tensor

  // vector3d<double> tauX(tau.XX(), tau.XY(), tau.XZ());
  // vector3d<double> tauY(tau.YX(), tau.YY(), tau.YZ());
  // vector3d<double> tauZ(tau.ZX(), tau.ZY(), tau.ZZ());

  // momX = tauX.DotProd(normArea);     
  // momY = tauY.DotProd(normArea);     
  // momZ = tauZ.DotProd(normArea);     
  // engy = (tauX.DotProd(normArea) * vel.X() + tauY.DotProd(normArea) * vel.Y() + tauZ.DotProd(normArea) * vel.Z()) + eqnState.GetConductivity(mu) * tGrad.DotProd(normArea);


  momX = tau.X();
  momY = tau.Y();
  momZ = tau.Z();
  engy = tau.DotProd(vel) + eqnState.GetConductivity(mu) * tGrad.DotProd(normArea);

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

//function to calculate the viscous flux jacobian using the thin shear layer approximation
squareMatrix CalcTSLFluxJac(const double &mu, const idealGas &eqnState, const double & vol, const vector3d<double> &areaVec, const primVars &left, const primVars &right){

  primVars diff = (right - left) ;
  double d1_rho = 1.0/right.Rho() - 1.0/left.Rho() ;

  double du_rho = (right.U()/right.Rho() - left.U()/left.Rho() ) ;
  double dv_rho = (right.V()/right.Rho() - left.V()/left.Rho() ) ;
  double dw_rho = (right.W()/right.Rho() - left.W()/left.Rho() ) ;

  double duu_rho = (right.U()*right.U()/right.Rho() - left.U()*left.U()/left.Rho() ) ;
  double dvv_rho = (right.V()*right.V()/right.Rho() - left.V()*left.V()/left.Rho() ) ;
  double dww_rho = (right.W()*right.W()/right.Rho() - left.W()*left.V()/left.Rho() ) ;

  double duv_rho = (right.U()*right.V()/right.Rho() - left.U()*left.V()/left.Rho() ) ;
  double duw_rho = (right.U()*right.W()/right.Rho() - left.U()*left.W()/left.Rho() ) ;
  double dvw_rho = (right.V()*right.W()/right.Rho() - left.V()*left.W()/left.Rho() ) ;

  double dE_rho = ( eqnState.GetEnergy( eqnState.GetSpecEnergy(right.P(), right.Rho() ), right.Velocity().Mag() )/ right.Rho() - eqnState.GetEnergy( eqnState.GetSpecEnergy(left.P(), left.Rho() ), left.Velocity().Mag() )/ left.Rho()) ;

  //calculate coefficients
  double a1 = (4.0/3.0) * areaVec.DotProd(areaVec);
  double a2 = (1.0/3.0) * areaVec.X() * areaVec.Y();
  double a3 = (1.0/3.0) * areaVec.X() * areaVec.Z();
  double a5 = (1.0/3.0) * areaVec.Y() * areaVec.Z();

  double a4 = areaVec.X() * areaVec.X() + (4.0/3.0) * areaVec.Y() * areaVec.Y() + areaVec.Z() * areaVec.Z();
  double a6 = areaVec.X() * areaVec.X() + areaVec.Y() * areaVec.Y() + (4.0/3.0) * areaVec.Z() * areaVec.Z();
  double a7 = eqnState.Gamma()/eqnState.GetPrandtl() * areaVec.DotProd(areaVec);

  double b21 = -a1 * du_rho -a2 * dv_rho -a3 * dw_rho;
  double b31 = -a2 * du_rho -a4 * dv_rho -a5 * dw_rho;
  double b41 = -a3 * du_rho -a5 * dv_rho -a6 * dw_rho;

  double b51 = a7 * (duu_rho * dvv_rho * dww_rho - dE_rho) -a1 * duu_rho -a4 * dvv_rho -a6 * dww_rho - 2.0 * a2 * duv_rho - 2.0 * a3 * duw_rho - 2.0 * a5 * dvw_rho;

  double b52 = -a7 * du_rho - b21;
  double b53 = -a7 * dv_rho - b31;
  double b54 = -a7 * dw_rho - b41;

  //calculate viscous flux jacobian
  squareMatrix vFluxJac(5);

  //column 0
  vFluxJac.SetData(0, 0, 0.0);
  vFluxJac.SetData(1, 0, b21);
  vFluxJac.SetData(2, 0, b31);
  vFluxJac.SetData(3, 0, b41);
  vFluxJac.SetData(4, 0, b51);

  //column 1
  vFluxJac.SetData(0, 1, 0.0);
  vFluxJac.SetData(1, 1, a1 * d1_rho);
  vFluxJac.SetData(2, 1, a2 * d1_rho);
  vFluxJac.SetData(3, 1, a3 * d1_rho);
  vFluxJac.SetData(4, 1, b52);

  //column 2
  vFluxJac.SetData(0, 2, 0.0);
  vFluxJac.SetData(1, 2, a2 * d1_rho);
  vFluxJac.SetData(2, 2, a4 * d1_rho);
  vFluxJac.SetData(3, 2, a5 * d1_rho);
  vFluxJac.SetData(4, 2, b53);

  //column 3
  vFluxJac.SetData(0, 1, 0.0);
  vFluxJac.SetData(1, 1, a3 * d1_rho);
  vFluxJac.SetData(2, 1, a5 * d1_rho);
  vFluxJac.SetData(3, 1, a6 * d1_rho);
  vFluxJac.SetData(4, 1, b54);

  //column 4
  vFluxJac.SetData(0, 1, 0.0);
  vFluxJac.SetData(1, 1, 0.0);
  vFluxJac.SetData(2, 1, 0.0);
  vFluxJac.SetData(3, 1, 0.0);
  vFluxJac.SetData(4, 1, a7 * d1_rho);

  return (mu/vol) * vFluxJac;
}

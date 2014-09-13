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

//member function to calculate the viscous flux jacobian
void viscousFlux::CalcFluxJac( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec, squareMatrix &dF_dUl, squareMatrix &dF_dUr, const primVars &left, const primVars &right){

  //check to see that input matricies are the correct size
  if ( dF_dUl.Size() != 5 || dF_dUr.Size() != 5 ){
    cerr << "ERROR: The matrices for the viscous flux jacobian are not the correct size!" << endl;
    exit(0);
  }

  //assemble dW_dUl -- derivative of primative variables at left side with respect to conserved variables at left side
  squareMatrix dW_dUl(dF_dUl.Size());
  //first row
  dW_dUl.SetData(0, 0, 1.0);
  dW_dUl.SetData(0, 1, 0.0);
  dW_dUl.SetData(0, 2, 0.0);
  dW_dUl.SetData(0, 3, 0.0);
  dW_dUl.SetData(0, 4, 0.0);
  //second row
  dW_dUl.SetData(1, 0, -1.0 * left.U() / left.Rho());
  dW_dUl.SetData(1, 1, 1.0 / left.Rho() );
  dW_dUl.SetData(1, 2, 0.0);
  dW_dUl.SetData(1, 3, 0.0);
  dW_dUl.SetData(1, 4, 0.0);
  //third row
  dW_dUl.SetData(2, 0, -1.0 * left.V() / left.Rho());
  dW_dUl.SetData(2, 1, 0.0);
  dW_dUl.SetData(2, 2, 1.0 / left.Rho() );
  dW_dUl.SetData(2, 3, 0.0);
  dW_dUl.SetData(2, 4, 0.0);
  //fourth row
  dW_dUl.SetData(3, 0, -1.0 * left.W() / left.Rho());
  dW_dUl.SetData(3, 1, 0.0);
  dW_dUl.SetData(3, 2, 0.0);
  dW_dUl.SetData(3, 3, 1.0 / left.Rho() );
  dW_dUl.SetData(3, 4, 0.0);
  //fifth row
  dW_dUl.SetData(4, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() );
  dW_dUl.SetData(4, 1, -1.0 * (eqnState.Gamma() - 1.0) * left.U() );
  dW_dUl.SetData(4, 2, -1.0 * (eqnState.Gamma() - 1.0) * left.V());
  dW_dUl.SetData(4, 3, -1.0 * (eqnState.Gamma() - 1.0) * left.W() );
  dW_dUl.SetData(4, 4, eqnState.Gamma() - 1.0 );

  //calculate dF_dWl -- derivative of flux at face with respect to primative variables at left side
  squareMatrix dF_dWl(dF_dUl.Size());
  //first row
  dF_dWl.SetData(0, 0, 0.0);
  dF_dWl.SetData(0, 1, 0.0);
  dF_dWl.SetData(0, 2, 0.0);
  dF_dWl.SetData(0, 3, 0.0);
  dF_dWl.SetData(0, 4, 0.0);




  //assemble dW_dUr -- derivative of primative variables at right side with respect to conservative variables at right side
  squareMatrix dW_dUr(dF_dUr.Size());
  //first row
  dW_dUr.SetData(0, 0, 1.0);
  dW_dUr.SetData(0, 1, 0.0);
  dW_dUr.SetData(0, 2, 0.0);
  dW_dUr.SetData(0, 3, 0.0);
  dW_dUr.SetData(0, 4, 0.0);
  //second row
  dW_dUr.SetData(1, 0, -1.0 * right.U() / right.Rho());
  dW_dUr.SetData(1, 1, 1.0 / right.Rho() );
  dW_dUr.SetData(1, 2, 0.0);
  dW_dUr.SetData(1, 3, 0.0);
  dW_dUr.SetData(1, 4, 0.0);
  //third row
  dW_dUr.SetData(2, 0, -1.0 * right.V() / right.Rho());
  dW_dUr.SetData(2, 1, 0.0);
  dW_dUr.SetData(2, 2, 1.0 / right.Rho() );
  dW_dUr.SetData(2, 3, 0.0);
  dW_dUr.SetData(2, 4, 0.0);
  //fourth row
  dW_dUr.SetData(3, 0, -1.0 * right.W() / right.Rho());
  dW_dUr.SetData(3, 1, 0.0);
  dW_dUr.SetData(3, 2, 0.0);
  dW_dUr.SetData(3, 3, 1.0 / right.Rho() );
  dW_dUr.SetData(3, 4, 0.0);
  //fifth row
  dW_dUr.SetData(4, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() );
  dW_dUr.SetData(4, 1, -1.0 * (eqnState.Gamma() - 1.0) * right.U() );
  dW_dUr.SetData(4, 2, -1.0 * (eqnState.Gamma() - 1.0) * right.V());
  dW_dUr.SetData(4, 3, -1.0 * (eqnState.Gamma() - 1.0) * right.W() );
  dW_dUr.SetData(4, 4, eqnState.Gamma() - 1.0 );







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

#include "viscousFlux.h"
#include <cmath> //sqrt

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;

//constructor -- initialize flux from velocity gradient
viscousFlux::viscousFlux( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec){
  vector3d<double> normArea = areaVec / areaVec.Mag();

  double lambda = suth.GetLambda(mu);

  double velGradTrace = velGrad.Trace();
  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

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

  double velGradTrace = velGrad.Trace();
  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

  momX = tau.X();
  momY = tau.Y();
  momZ = tau.Z();
  engy = tau.DotProd(vel) + eqnState.GetConductivity(mu) * tGrad.DotProd(normArea);

}

void CalcTSLFluxJac(const double &mu, const idealGas &eqnState, const vector3d<double> &areaVec, const primVars &left, const primVars &right, const double &dist, squareMatrix &dFv_dUl, squareMatrix &dFv_dUr, const sutherland &suth){

  //check to make sure square matrices are of correct size
  if (dFv_dUl.Size() != 5 || dFv_dUr.Size() != 5){
    cerr << "ERROR: Error in thin shear layer viscous jacobian calculation. The input jacobian matrices are not of the correct size!" << endl;
    exit(0);
  }

  //normalize area vector
  vector3d<double> normArea = areaVec / areaVec.Mag();

  //get velocity at face
  vector3d<double> vel = 0.5 * (right.Velocity() + left.Velocity());

  //calculate thin shear layer velocity gradients
  tensor<double> velGradTSL = CalcVelGradTSL(left, right, normArea, dist);

  //calculate bulk viscosity
  double lambda = suth.GetLambda(mu);

  //calculate shear stress at face
  double velGradTrace = velGradTSL.Trace();
  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGradTSL.MatMult(normArea) + velGradTSL.Transpose().MatMult(normArea));

  //calculate coefficients
  double theta = normArea.MagSq();
  double thetaX = (4.0/3.0) * normArea.X() * normArea.X() + normArea.Y() * normArea.Y() + normArea.Z() * normArea.Z();
  double thetaY = normArea.X() * normArea.X() + (4.0/3.0) * normArea.Y() * normArea.Y() + normArea.Z() * normArea.Z();
  double thetaZ = normArea.X() * normArea.X() + normArea.Y() * normArea.Y() + (4.0/3.0) * normArea.Z() * normArea.Z();

  double etaX = (1.0/3.0) * normArea.Y() * normArea.Z();
  double etaY = (1.0/3.0) * normArea.X() * normArea.Z();
  double etaZ = (1.0/3.0) * normArea.X() * normArea.Y();

  double piX = vel.X() * thetaX + vel.Y() * etaZ   + vel.Z() * etaY;
  double piY = vel.X() * etaZ   + vel.Y() * thetaY + vel.Z() * etaX;
  double piZ = vel.X() * etaY   + vel.Y() * etaX   + vel.Z() * thetaZ;

  double phiRhoL = -1.0 * eqnState.GetConductivity(mu) * left.Temperature(eqnState) / (mu * left.Rho());
  double phiRhoR = -1.0 * eqnState.GetConductivity(mu) * right.Temperature(eqnState) / (mu * right.Rho());

  double phiPressL = eqnState.GetConductivity(mu) / (mu * left.Rho());
  double phiPressR = eqnState.GetConductivity(mu) / (mu * right.Rho());

  //calculate matrix - derivative of left primative vars wrt left conservative vars
  squareMatrix dWl_dUl(5);
  dWl_dUl.Zero();

  //column 0
  dWl_dUl.SetData(0,0, 1.0);
  dWl_dUl.SetData(1,0, -1.0 * left.U() / left.Rho());
  dWl_dUl.SetData(2,0, -1.0 * left.V() / left.Rho());
  dWl_dUl.SetData(3,0, -1.0 * left.W() / left.Rho());
  dWl_dUl.SetData(4,0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq());

  //column 1
  dWl_dUl.SetData(1,1, 1.0/left.Rho());
  dWl_dUl.SetData(4,1, -1.0 * (eqnState.Gamma() - 1.0) * left.U());

  //column 2
  dWl_dUl.SetData(2,2, 1.0/left.Rho());
  dWl_dUl.SetData(4,2, -1.0 * (eqnState.Gamma() - 1.0) * left.V());

  //column 3
  dWl_dUl.SetData(3,3, 1.0/left.Rho());
  dWl_dUl.SetData(4,3, -1.0 * (eqnState.Gamma() - 1.0) * left.W());

  //column 4
  dWl_dUl.SetData(4,4, eqnState.Gamma() - 1.0);


  //calculate matrix - derivative of right primative vars wrt right conservative vars
  squareMatrix dWr_dUr(5);
  dWr_dUr.Zero();

  //column 0
  dWr_dUr.SetData(0,0, 1.0);
  dWr_dUr.SetData(1,0, -1.0 * right.U() / right.Rho());
  dWr_dUr.SetData(2,0, -1.0 * right.V() / right.Rho());
  dWr_dUr.SetData(3,0, -1.0 * right.W() / right.Rho());
  dWr_dUr.SetData(4,0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq());

  //column 1
  dWr_dUr.SetData(1,1, 1.0/right.Rho());
  dWr_dUr.SetData(4,1, -1.0 * (eqnState.Gamma() - 1.0) * right.U());

  //column 2
  dWr_dUr.SetData(2,2, 1.0/right.Rho());
  dWr_dUr.SetData(4,2, -1.0 * (eqnState.Gamma() - 1.0) * right.V());

  //column 3
  dWr_dUr.SetData(3,3, 1.0/right.Rho());
  dWr_dUr.SetData(4,3, -1.0 * (eqnState.Gamma() - 1.0) * right.W());

  //column 4
  dWr_dUr.SetData(4,4, eqnState.Gamma() - 1.0);


  //calculate matrix - derivative of viscous flux wrt left primative vars

  //column 0
  dFv_dUl.SetData(0,0, 0.0);
  dFv_dUl.SetData(1,0, 0.0);
  dFv_dUl.SetData(2,0, 0.0);
  dFv_dUl.SetData(3,0, 0.0);
  dFv_dUl.SetData(4,0, phiRhoL * theta);

  //column 1
  dFv_dUl.SetData(0,1, 0.0);
  dFv_dUl.SetData(1,1, thetaX);
  dFv_dUl.SetData(2,1, etaZ);
  dFv_dUl.SetData(3,1, etaY);
  dFv_dUl.SetData(4,1, -0.5 * (dist/mu) * tau.X() + piX);

  //column 2
  dFv_dUl.SetData(0,2, 0.0);
  dFv_dUl.SetData(1,2, etaZ);
  dFv_dUl.SetData(2,2, thetaY);
  dFv_dUl.SetData(3,2, etaX);
  dFv_dUl.SetData(4,2, -0.5 * (dist/mu) * tau.Y() + piY);

  //column 3
  dFv_dUl.SetData(0,3, 0.0);
  dFv_dUl.SetData(1,3, etaY);
  dFv_dUl.SetData(2,3, etaX);
  dFv_dUl.SetData(3,3, thetaZ);
  dFv_dUl.SetData(4,3, -0.5 * (dist/mu) * tau.Z() + piZ);

  //column 4
  dFv_dUl.SetData(0,4, 0.0);
  dFv_dUl.SetData(1,4, 0.0);
  dFv_dUl.SetData(2,4, 0.0);
  dFv_dUl.SetData(3,4, 0.0);
  dFv_dUl.SetData(4,4, phiPressL * theta);

  dFv_dUl = -1.0 * (mu/dist) * dFv_dUl;

  //calculate matrix - derivative of viscous flux wrt right primative vars

  //column 0
  dFv_dUr.SetData(0,0, 0.0);
  dFv_dUr.SetData(1,0, 0.0);
  dFv_dUr.SetData(2,0, 0.0);
  dFv_dUr.SetData(3,0, 0.0);
  dFv_dUr.SetData(4,0, phiRhoR * theta);

  //column 1
  dFv_dUr.SetData(0,1, 0.0);
  dFv_dUr.SetData(1,1, thetaX);
  dFv_dUr.SetData(2,1, etaZ);
  dFv_dUr.SetData(3,1, etaY);
  dFv_dUr.SetData(4,1, 0.5 * (dist/mu) * tau.X() + piX);

  //column 2
  dFv_dUr.SetData(0,2, 0.0);
  dFv_dUr.SetData(1,2, etaZ);
  dFv_dUr.SetData(2,2, thetaY);
  dFv_dUr.SetData(3,2, etaX);
  dFv_dUr.SetData(4,2, 0.5 * (dist/mu) * tau.Y() + piY);

  //column 3
  dFv_dUr.SetData(0,3, 0.0);
  dFv_dUr.SetData(1,3, etaY);
  dFv_dUr.SetData(2,3, etaX);
  dFv_dUr.SetData(3,3, thetaZ);
  dFv_dUr.SetData(4,3, 0.5 * (dist/mu) * tau.Z() + piZ);

  //column 4
  dFv_dUr.SetData(0,4, 0.0);
  dFv_dUr.SetData(1,4, 0.0);
  dFv_dUr.SetData(2,4, 0.0);
  dFv_dUr.SetData(3,4, 0.0);
  dFv_dUr.SetData(4,4, phiPressR * theta);

  dFv_dUr = (mu/dist) * dFv_dUr;

  //multiply by dW_dU to get flux jacobian derivative wrt conservative variables
  //cout << "dFv_dUl" << endl << dFv_dUl << endl;
  //cout << "dFv_dUr" << endl << dFv_dUr << endl;

  dFv_dUl = dFv_dUl * dWl_dUl;
  dFv_dUr = dFv_dUr * dWr_dUr;

  primVars faceState = 0.5 * (left + right);
  dFv_dUl.Identity();
  dFv_dUr.Identity();
  double specRad = mu * eqnState.Gamma() / (eqnState.GetPrandtl() * faceState.Rho() * dist) ;


  dFv_dUl = -1.0 * specRad * dFv_dUl;
  dFv_dUr = specRad * dFv_dUr;

}


//function to calculate the velocity gradients at a cell face using the Thin Shear Layer approximation
tensor<double> CalcVelGradTSL(const primVars &left, const primVars &right, const vector3d<double> &areaVec, const double &dist){

  //normalize area vector
  vector3d<double> normArea = areaVec / areaVec.Mag();

  //initialize velocity gradient tensor
  tensor<double> velGrad;

  //calculate velocity derivatives
  vector3d<double> velDeriv = (right.Velocity() - left.Velocity()) / dist;

  //populate velocity gradient tensor
  velGrad.SetXX( velDeriv.X() * normArea.X() );
  velGrad.SetXY( velDeriv.Y() * normArea.X() );
  velGrad.SetXZ( velDeriv.Z() * normArea.X() );

  velGrad.SetYX( velDeriv.X() * normArea.Y() );
  velGrad.SetYY( velDeriv.Y() * normArea.Y() );
  velGrad.SetYZ( velDeriv.Z() * normArea.Y() );

  velGrad.SetZX( velDeriv.X() * normArea.Z() );
  velGrad.SetZY( velDeriv.Y() * normArea.Z() );
  velGrad.SetZZ( velDeriv.Z() * normArea.Z() );

  return velGrad;
}

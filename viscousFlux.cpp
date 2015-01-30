#include "viscousFlux.h"
#include <cmath> //sqrt

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;

//constructor -- initialize flux from velocity gradient
/*
Viscous flux normal to face:
F = [ 0,
      taux,
      tauy,
      tauz,
      tau (dot) vel + K * tGrad (dot) area ]

In the above equation tau is the wall shear stress. Taux, tauy, and tauz are the rows of the wall shear stress tensor i.e. 
taux = tauxx + tauxy + tauxz. K is the thermal conductivity, tGrad is the temperature gradient, and area is the normalized
face area.

Wall shear stress:
tau = lambda * velGradTrace * area + mu * ( velGrad * area + velGrad' * area)

In the above equation lambda is the bulk viscosity, velGradTrace is the trace of the velocity gradient, area is the normalized
face area, mu is the dynamic viscosity, and velGrad is the velocity gradient tensor.
*/
viscousFlux::viscousFlux( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, 
			  const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec){
  // velGrad -- velocity gradient tensor
  // vel -- velocity vector
  // mu -- dynamic viscosity
  // suth -- method to get viscosity as a function of temperature (Sutherland's law)
  // eqnState -- equation of state
  // tGrad -- temperature gradient
  // areaVec -- area vector of face

  vector3d<double> normArea = areaVec / areaVec.Mag(); //normalize face area

  double lambda = suth.GetLambda(mu); //get 2nd coefficient of viscosity assuming bulk viscosity is 0 (Stoke's hypothesis)

  double velGradTrace = velGrad.Trace(); //trace of velocity gradient
  //wall shear stress
  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

  data[0] = tau.X();
  data[1] = tau.Y();
  data[2] = tau.Z();
  data[3] = tau.DotProd(vel) + eqnState.GetConductivity(mu) * tGrad.DotProd(normArea);

}

//non-member functions -----------------------------------------------------------------------------------------------------------//
//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, viscousFlux &flux){
  os << "0.0, " << flux.data[0] << ", " << flux.data[1] << ", " << flux.data[2] << ", " << flux.data[3] << endl;
  return os;
}

//member function for scalar multiplication
viscousFlux  viscousFlux::operator * (const double &scalar){
  viscousFlux temp = *this;
  temp.data[0] *= scalar;
  temp.data[1] *= scalar;
  temp.data[2] *= scalar;
  temp.data[3] *= scalar;
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
  temp.data[0] /= scalar;
  temp.data[1] /= scalar;
  temp.data[2] /= scalar;
  temp.data[3] /= scalar;
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
void viscousFlux::SetFlux( const tensor<double> &velGrad, const vector3d<double> &vel, const double &mu, const sutherland &suth, 
			   const idealGas &eqnState, const vector3d<double> &tGrad, const vector3d<double> &areaVec){
  // velGrad -- velocity gradient tensor
  // vel -- velocity vector
  // mu -- dynamic viscosity
  // suth -- method to get viscosity as a function of temperature (Sutherland's law)
  // eqnState -- equation of state
  // tGrad -- temperature gradient
  // areaVec -- area vector of face

  vector3d<double> normArea = areaVec / areaVec.Mag(); //normalize face area

  double lambda = suth.GetLambda(mu); //get 2nd coefficient of viscosity assuming bulk viscosity is 0 (Stoke's hypothesis)

  double velGradTrace = velGrad.Trace(); //get trace of velocity tensor
  //calculate wall shear stress
  vector3d<double> tau = lambda * velGradTrace * normArea + mu * (velGrad.MatMult(normArea) + velGrad.Transpose().MatMult(normArea));

  data[0] = tau.X();
  data[1] = tau.Y();
  data[2] = tau.Z();
  data[3] = tau.DotProd(vel) + eqnState.GetConductivity(mu) * tGrad.DotProd(normArea);
}

//function to calculate the thin shear layer flux jacobian -- NOT USED in LUSGS formulation
void CalcTSLFluxJac(const double &mu, const idealGas &eqnState, const vector3d<double> &areaVec, const primVars &left, const primVars &right, 
		    const double &dist, squareMatrix &dFv_dUl, squareMatrix &dFv_dUr, const sutherland &suth){
  // mu -- dynamic viscosity
  // eqnState -- equation of state
  // areaVec -- area vector of face
  // left -- left state (primative)
  // right -- right state (primative)
  // dist -- distance from centroid of left cell to centroid of right cell
  // dFv_dUl -- flux jacobian of viscous flux with respect to left state
  // dFV_dUr -- flux jacobian of viscous flux with respect to right state
  // suth -- method to get viscosity as a function of temperature

  //check to make sure square matrices are of correct size
  if (dFv_dUl.Size() != 5 || dFv_dUr.Size() != 5){
    cerr << "ERROR: Error in viscousFlux.cpp:CalcTSLFluxJac. Problem with thin shear layer viscous jacobian calculation. The input jacobian matrices are not of the correct size!" << endl;
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

  //calculate coefficients (from Blazek)
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

  //------------------------------------------------------------------------------------------------------------
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

  //--------------------------------------------------------------------------------------------------------------
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

  //--------------------------------------------------------------------------------------------------------
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
  dFv_dUl = dFv_dUl * dWl_dUl;
  dFv_dUr = dFv_dUr * dWr_dUr;

  //calculate spectral radius
  primVars faceState = 0.5 * (left + right);
  dFv_dUl.Identity();
  dFv_dUr.Identity();
  double specRad = mu * eqnState.Gamma() / (eqnState.GetPrandtl() * faceState.Rho() * dist) ;

  //add or subtract spectral radius to flux jacobian
  dFv_dUl = -1.0 * specRad * dFv_dUl;
  dFv_dUr = specRad * dFv_dUr;

}


//function to calculate the velocity gradients at a cell face using the Thin Shear Layer approximation
//NOT USED in LUSGS formulation
tensor<double> CalcVelGradTSL(const primVars &left, const primVars &right, const vector3d<double> &areaVec, const double &dist){
  // left -- left state (primative)
  // right -- right state (primative)
  // areaVec -- area vector of face
  // dist -- distance between centroid of left cell and right cell

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

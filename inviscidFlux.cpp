#include "inviscidFlux.h"
#include <cmath> //sqrt

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;
using std::copysign;

//constructor -- initialize flux from state vector
inviscidFlux::inviscidFlux( const primVars &state, const idealGas &eqnState, const vector3d<double>& areaVec){
  vector3d<double> normArea = areaVec / areaVec.Mag();
  vector3d<double> vel = state.Velocity();

  rhoVel  = state.Rho() * vel.DotProd(normArea);      
  rhoVelU = state.Rho() * vel.DotProd(normArea) * vel.X() + state.P() * normArea.X();     
  rhoVelV = state.Rho() * vel.DotProd(normArea) * vel.Y() + state.P() * normArea.Y();     
  rhoVelW = state.Rho() * vel.DotProd(normArea) * vel.Z() + state.P() * normArea.Z();     
  rhoVelH = state.Rho() * vel.DotProd(normArea) * state.Enthalpy(eqnState);
}

//constructor -- initialize flux from state vector using conservative variables
inviscidFlux::inviscidFlux( const colMatrix &cons, const idealGas &eqnState, const vector3d<double>& areaVec){

  //check to see that colMatrix is correct size
  if (cons.Size() != 5){
    cerr << "ERROR: Error in inviscidFlux::inviscidFlux. Column matrix of conservative variables is not the correct size!" << endl;
    exit(0);
  }

  primVars state;
  state.SetRho(cons.Data(0));
  state.SetU(cons.Data(1)/cons.Data(0));
  state.SetV(cons.Data(2)/cons.Data(0));
  state.SetW(cons.Data(3)/cons.Data(0));
  double energy = cons.Data(4)/cons.Data(0);
  state.SetP(eqnState.GetPressFromEnergy(state.Rho(), energy, state.Velocity().Mag() ));

  vector3d<double> normArea = areaVec / areaVec.Mag();
  vector3d<double> vel = state.Velocity();

  rhoVel  = state.Rho() * vel.DotProd(normArea);      
  rhoVelU = state.Rho() * vel.DotProd(normArea) * vel.X() + state.P() * normArea.X();     
  rhoVelV = state.Rho() * vel.DotProd(normArea) * vel.Y() + state.P() * normArea.Y();     
  rhoVelW = state.Rho() * vel.DotProd(normArea) * vel.Z() + state.P() * normArea.Z();     
  rhoVelH = state.Rho() * vel.DotProd(normArea) * state.Enthalpy(eqnState);
}

void inviscidFlux::SetFlux( const primVars &state, const idealGas &eqnState, const vector3d<double>& areaVec){
  vector3d<double> normArea = areaVec / areaVec.Mag();
  vector3d<double> vel = state.Velocity();

  rhoVel  = state.Rho() * vel.DotProd(normArea);      
  rhoVelU = state.Rho() * vel.DotProd(normArea) * vel.X() + state.P() * normArea.X();     
  rhoVelV = state.Rho() * vel.DotProd(normArea) * vel.Y() + state.P() * normArea.Y();     
  rhoVelW = state.Rho() * vel.DotProd(normArea) * vel.Z() + state.P() * normArea.Z();     
  rhoVelH = state.Rho() * vel.DotProd(normArea) * state.Enthalpy(eqnState);
}


//function to calculate pressure from conserved variables and equation of state
inviscidFlux RoeFlux( const primVars &left, const primVars &right, const idealGas &eqnState, const vector3d<double>& areaVec, double &maxWS){

  //compute Rho averaged quantities
  double denRatio = sqrt(right.Rho()/left.Rho());
  double rhoR = left.Rho() * denRatio;  //Roe averaged density
  double uR = (left.U() + denRatio * right.U()) / (1.0 + denRatio);  //Roe averaged u-velocity
  double vR = (left.V() + denRatio * right.V()) / (1.0 + denRatio);  //Roe averaged v-velocity
  double wR = (left.W() + denRatio * right.W()) / (1.0 + denRatio);  //Roe averaged w-velocity
  double hR = (left.Enthalpy(eqnState) + denRatio * right.Enthalpy(eqnState)) / (1.0 + denRatio);  //Roe averaged total enthalpy
  double aR = sqrt( (eqnState.Gamma() - 1.0) * (hR - 0.5 * (uR*uR + vR*vR + wR*wR)) );  //Roe averaged speed of sound
  //Roe averaged face normal velocity
  vector3d<double> velR(uR,vR,wR);

  //cout << "Roe velocity " << velR << endl;

  vector3d<double> areaNorm = areaVec / areaVec.Mag();
  double velRSum = velR.DotProd(areaNorm);

  //calculate wave strengths
  // double denDiff = right.Rho() - left.Rho();
  // double pressDiff = right.Pressure(eqnState) - left.Pressure(eqnState);

  double normVelDiff = right.Velocity().DotProd(areaNorm) - left.Velocity().DotProd(areaNorm);

  //vector<double> waveStrength(4);
  double waveStrength[4] = {((right.P() - left.P()) - rhoR * aR * normVelDiff) / (2.0 * aR * aR), 
			    (right.Rho() - left.Rho()) - (right.P() - left.P()) / (aR * aR), 
			    ((right.P() - left.P()) + rhoR * aR * normVelDiff) / (2.0 * aR * aR), 
			    rhoR};

  // waveStrength[0] = (pressDiff - rhoR * aR * normVelDiff) / (2.0 * aR * aR);       //left moving acoustic wave strength
  // waveStrength[1] = denDiff - pressDiff / (aR * aR);                         //entropy wave strength
  // waveStrength[2] = (pressDiff + rhoR * aR * normVelDiff) / (2.0 * aR * aR);        //right moving acoustic wave strength
  // waveStrength[3] = rhoR;                                                    // shear waves get combined into one factor


  //calculate absolute value of wave speeds
  //vector<double> waveSpeed(4);
  double waveSpeed[4] = {fabs(velRSum - aR), 
			 fabs(velRSum), 
			 fabs(velRSum + aR), 
			 fabs(velRSum)};

  // waveSpeed[0] = fabs(velRSum - aR);                                          //left moving acoustic wave speed
  // waveSpeed[1] = fabs(velRSum);                                               //entropy wave speed
  // waveSpeed[2] = fabs(velRSum + aR);                                          //right moving acoustic wave speed
  // waveSpeed[3] = fabs(velRSum);                                               //shear wave speed

  //calculate entropy fix
  double entropyFix = 0.1;                                                            // default setting for entropy fix to kick in

  if ( waveSpeed[0] < entropyFix ){
    waveSpeed[0] = 0.5 * (waveSpeed[0] * waveSpeed[0] / entropyFix + entropyFix);
  }
  if ( waveSpeed[2] < entropyFix ){
    waveSpeed[2] = 0.5 * (waveSpeed[2] * waveSpeed[2] / entropyFix + entropyFix);
  }

  maxWS = fabs(velRSum) + aR;

  //calculate eigenvectors
  // vector<double> lAcousticEigV(5);
  // lAcousticEigV[0] = 1.0;
  // lAcousticEigV[1] = uR - aR * areaNorm.X();
  // lAcousticEigV[2] = vR - aR * areaNorm.Y();
  // lAcousticEigV[3] = wR - aR * areaNorm.Z();
  // lAcousticEigV[4] = hR - aR * velRSum;

  double lAcousticEigV[5] = {1.0, 
			     uR - aR * areaNorm.X(), 
			     vR - aR * areaNorm.Y(), 
			     wR - aR * areaNorm.Z(), 
			     hR - aR * velRSum};

  // vector<double> entropyEigV(5);
  // entropyEigV[0] = 1.0;
  // entropyEigV[1] = uR;
  // entropyEigV[2] = vR;
  // entropyEigV[3] = wR;
  // entropyEigV[4] = 0.5 * ( uR * uR + vR * vR + wR * wR);

  double entropyEigV[5] = {1.0, 
			   uR, 
			   vR, 
			   wR, 
			   0.5 * ( uR * uR + vR * vR + wR * wR)};

  // vector<double> rAcousticEigV(5);
  // rAcousticEigV[0] = 1.0;
  // rAcousticEigV[1] = uR + aR * areaNorm.X();
  // rAcousticEigV[2] = vR + aR * areaNorm.Y();
  // rAcousticEigV[3] = wR + aR * areaNorm.Z();
  // rAcousticEigV[4] = hR + aR * velRSum;

  double rAcousticEigV[5] = {1.0, 
			     uR + aR * areaNorm.X(), 
			     vR + aR * areaNorm.Y(), 
			     wR + aR * areaNorm.Z(), 
			     hR + aR * velRSum};

  // vector<double> shearEigV(5);
  // shearEigV[0] = 0.0;
  // shearEigV[1] = (right.Vx() - left.Vx()) - normVelDiff * areaNorm.X();
  // shearEigV[2] = (right.Vy() - left.Vy()) - normVelDiff * areaNorm.Y();
  // shearEigV[3] = (right.Vz() - left.Vz()) - normVelDiff * areaNorm.Z();
  // shearEigV[4] = uR * (right.Vx() - left.Vx()) + vR * (right.Vy() - left.Vy()) + 
  //                wR * (right.Vz() - left.Vz()) - velRSum * normVelDiff;

  double shearEigV[5] = {0.0, 
			 (right.U() - left.U()) - normVelDiff * areaNorm.X(), 
			 (right.V() - left.V()) - normVelDiff * areaNorm.Y(), 
			 (right.W() - left.W()) - normVelDiff * areaNorm.Z(), 
			 uR * (right.U() - left.U()) + vR * (right.V() - left.V()) + wR * (right.W() - left.W()) - velRSum * normVelDiff};

  //calculate dissipation term ( eigenvector * wave speed * wave strength)
  //vector<double> dissipation(5,0.0);

  double dissipation[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  unsigned int ii = 0;
  for ( ii=0; ii < 5; ii++ ) {                                    
    dissipation[ii] = waveSpeed[0] * waveStrength[0] * lAcousticEigV[ii]   //contribution from left acoustic wave
                    + waveSpeed[1] * waveStrength[1] * entropyEigV[ii]     //contribution from entropy wave
                    + waveSpeed[2] * waveStrength[2] * rAcousticEigV[ii]   //contribution from right acoustic wave
                    + waveSpeed[3] * waveStrength[3] * shearEigV[ii];      //contribution from shear wave
  }


  //cout << "Roe Dissipation term " << dissipation[0] << ", " << dissipation[1] << ", " << dissipation[2] << ", " << dissipation[3] << ", " << dissipation[4] << endl;

  //calculate left/right physical flux
  inviscidFlux leftFlux(left, eqnState, areaNorm);
  inviscidFlux rightFlux(right, eqnState, areaNorm);

  //calculate numerical flux
  inviscidFlux flux;
  flux.SetRhoVel(  0.5 * (leftFlux.RhoVel()  + rightFlux.RhoVel()  - dissipation[0]) );
  flux.SetRhoVelU( 0.5 * (leftFlux.RhoVelU() + rightFlux.RhoVelU() - dissipation[1]) );
  flux.SetRhoVelV( 0.5 * (leftFlux.RhoVelV() + rightFlux.RhoVelV() - dissipation[2]) );
  flux.SetRhoVelW( 0.5 * (leftFlux.RhoVelW() + rightFlux.RhoVelW() - dissipation[3]) );
  flux.SetRhoVelH( 0.5 * (leftFlux.RhoVelH() + rightFlux.RhoVelH() - dissipation[4]) );

  return flux;

}

//function to calculate exact Roe flux jacobians
void ApproxRoeFluxJacobian( const primVars &left, const primVars &right, const idealGas &eqnState, const vector3d<double>& areaVec, double &maxWS, squareMatrix &dF_dUl, squareMatrix &dF_dUr){

  //left --> primative variables from left side
  //right --> primative variables from right side
  //eqnStat --> ideal gas equation of state
  //areaVec --> face area vector
  //maxWS --> maximum wave speed
  //dF_dUl --> dF/dUl, derivative of the Roe flux wrt the left state (conservative variables)
  //dF_dUr --> dF/dUlr, derivative of the Roe flux wrt the right state (conservative variables)


  //check to see that output matricies are correct size
  if( (dF_dUl.Size() != 5) || (dF_dUr.Size() != 5)){
    cerr << "ERROR: Input matricies to RoeFLuxJacobian function are not the correct size!" << endl;
  }

  //compute Rho averaged quantities
  double denRatio = sqrt(right.Rho()/left.Rho());
  //double rhoR = left.Rho() * denRatio;  //Roe averaged density
  double uR = (left.U() + denRatio * right.U()) / (1.0 + denRatio);  //Roe averaged u-velocity
  double vR = (left.V() + denRatio * right.V()) / (1.0 + denRatio);  //Roe averaged v-velocity
  double wR = (left.W() + denRatio * right.W()) / (1.0 + denRatio);  //Roe averaged w-velocity
  double hR = (left.Enthalpy(eqnState) + denRatio * right.Enthalpy(eqnState)) / (1.0 + denRatio);  //Roe averaged total enthalpy
  double aR = sqrt( (eqnState.Gamma() - 1.0) * (hR - 0.5 * (uR*uR + vR*vR + wR*wR)) );  //Roe averaged speed of sound
  //Roe averaged face normal velocity
  vector3d<double> velR(uR,vR,wR);

  vector3d<double> areaNorm = areaVec / areaVec.Mag();  //normalize area vector to unit vector

  //dot product of velocities (Roe, left, right) with unit area vector
  double velRNorm = velR.DotProd(areaNorm);
  double velLeftNorm = left.Velocity().DotProd(areaNorm);
  double velRightNorm = right.Velocity().DotProd(areaNorm);

  //calculate diagonal eigenvalue matrix |lambda|
  squareMatrix lambda(5);
  lambda.Zero();
  lambda.SetData(0, 0, fabs(velRNorm - aR) );
  lambda.SetData(1, 1, fabs(velRNorm) );
  lambda.SetData(2, 2, fabs(velRNorm + aR) );
  lambda.SetData(3, 3, fabs(velRNorm) );
  lambda.SetData(4, 4, fabs(velRNorm) );

  //calculate Roe jacobian matrix A
  //contribution due to normal velocity eigenvalues
  squareMatrix A(5);

  //column zero
  A.SetData(0, 0, 1.0 - 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) );
  A.SetData(1, 0, -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * uR + velRNorm * areaNorm.X() );
  A.SetData(2, 0, -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * vR + velRNorm * areaNorm.Y() );
  A.SetData(3, 0, -(1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * wR + velRNorm * areaNorm.Z() );
  A.SetData(4, 0, velRNorm * velRNorm - 0.5 * velR.MagSq() * (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) ); 

  //column one
  A.SetData(0, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * uR );
  A.SetData(1, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * uR + 1.0 - areaNorm.X() * areaNorm.X() );
  A.SetData(2, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * uR       - areaNorm.Y() * areaNorm.X() );
  A.SetData(3, 1, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * uR       - areaNorm.Z() * areaNorm.X() );
  A.SetData(4, 1, (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * uR - velRNorm * areaNorm.X() );

  //column two
  A.SetData(0, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * vR );
  A.SetData(1, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * vR       - areaNorm.X() * areaNorm.Y() );
  A.SetData(2, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * vR + 1.0 - areaNorm.Y() * areaNorm.Y() );
  A.SetData(3, 2, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * vR       - areaNorm.Z() * areaNorm.Y() );
  A.SetData(4, 2, (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * vR - velRNorm * areaNorm.Y() );

  //column three
  A.SetData(0, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * wR );
  A.SetData(1, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * uR * wR       - areaNorm.X() * areaNorm.Z() );
  A.SetData(2, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * vR * wR       - areaNorm.Y() * areaNorm.Z() );
  A.SetData(3, 3, (eqnState.Gamma() - 1.0) / (aR * aR) * wR * wR + 1.0 - areaNorm.Z() * areaNorm.Z() );
  A.SetData(4, 3, (1.0 + 0.5 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR)) * wR - velRNorm * areaNorm.Z() );

  //column four
  A.SetData(0, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) );
  A.SetData(1, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * uR );
  A.SetData(2, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * vR );
  A.SetData(3, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * wR );
  A.SetData(4, 4, -(eqnState.Gamma() - 1.0) / (aR * aR) * velR.MagSq() / (aR * aR) );

  A = fabs(velRNorm) * A;

  //contribution due to u - c wave
  squareMatrix temp(5);

  //column zero
  temp.SetData(0, 0, 0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) + 0.5 * velRNorm / aR );
  temp.SetData(1, 0, (uR - aR * areaNorm.X()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) + 0.5 * velRNorm / aR) );
  temp.SetData(2, 0, (vR - aR * areaNorm.Y()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) + 0.5 * velRNorm / aR) );
  temp.SetData(3, 0, (wR - aR * areaNorm.Z()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) + 0.5 * velRNorm / aR) );
  temp.SetData(4, 0, (hR - velRNorm * aR * areaNorm.X()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) + 0.5 * velRNorm / aR) );

  //column one
  temp.SetData(0, 1, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR - areaNorm.X() / (2.0 * aR ) );
  temp.SetData(1, 1,  (uR - aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR - areaNorm.X() / (2.0 * aR)) );
  temp.SetData(2, 1,  (vR - aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR - areaNorm.X() / (2.0 * aR)) );
  temp.SetData(3, 1,  (wR - aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR - areaNorm.X() / (2.0 * aR)) );
  temp.SetData(4, 1,  (hR - aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR - areaNorm.X() / (2.0 * aR)) );

  //column two
  temp.SetData(0, 2, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR - areaNorm.Y() / (2.0 * aR ) );
  temp.SetData(1, 2,  (uR - aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR - areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(2, 2,  (vR - aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR - areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(3, 2,  (wR - aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR - areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(4, 2,  (hR - aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR - areaNorm.Y() / (2.0 * aR)) );

  //column three
  temp.SetData(0, 3, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR - areaNorm.Z() / (2.0 * aR ) );
  temp.SetData(1, 3,  (uR - aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR - areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(2, 3,  (vR - aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR - areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(3, 3,  (wR - aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR - areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(4, 3,  (hR - aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR - areaNorm.Z() / (2.0 * aR)) );

  //column four
  temp.SetData(0, 4, (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(1, 4, (uR - aR * areaNorm.X()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(2, 4, (vR - aR * areaNorm.Y()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(3, 4, (wR - aR * areaNorm.Z()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(4, 4, (hR - aR * velRNorm)     * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );

  A = A + (fabs(velRNorm - aR) * temp);

  //contribution due to u + c wave

  //column zero
  temp.SetData(0, 0, 0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) - 0.5 * velRNorm / aR );
  temp.SetData(1, 0, (uR + aR * areaNorm.X()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) - 0.5 * velRNorm / aR) );
  temp.SetData(2, 0, (vR + aR * areaNorm.Y()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) - 0.5 * velRNorm / aR) );
  temp.SetData(3, 0, (wR + aR * areaNorm.Z()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) - 0.5 * velRNorm / aR) );
  temp.SetData(4, 0, (hR + velRNorm * aR * areaNorm.X()) * (0.25 * (eqnState.Gamma() - 1.0) * velR.MagSq() / (aR * aR) - 0.5 * velRNorm / aR) );

  //column one
  temp.SetData(0, 1, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR + areaNorm.X() / (2.0 * aR ) );
  temp.SetData(1, 1,  (uR + aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR + areaNorm.X() / (2.0 * aR)) );
  temp.SetData(2, 1,  (vR + aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR + areaNorm.X() / (2.0 * aR)) );
  temp.SetData(3, 1,  (wR + aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR + areaNorm.X() / (2.0 * aR)) );
  temp.SetData(4, 1,  (hR + aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * uR + areaNorm.X() / (2.0 * aR)) );

  //column two
  temp.SetData(0, 2, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR + areaNorm.Y() / (2.0 * aR ) );
  temp.SetData(1, 2,  (uR + aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR + areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(2, 2,  (vR + aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR + areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(3, 2,  (wR + aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR + areaNorm.Y() / (2.0 * aR)) );
  temp.SetData(4, 2,  (hR + aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * vR + areaNorm.Y() / (2.0 * aR)) );

  //column three
  temp.SetData(0, 3, -(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR + areaNorm.Z() / (2.0 * aR ) );
  temp.SetData(1, 3,  (uR + aR * areaNorm.X()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR + areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(2, 3,  (vR + aR * areaNorm.Y()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR + areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(3, 3,  (wR + aR * areaNorm.Z()) * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR + areaNorm.Z() / (2.0 * aR)) );
  temp.SetData(4, 3,  (hR + aR * velRNorm)     * (-(eqnState.Gamma() - 1.0) / (2.0 * aR * aR) * wR + areaNorm.Z() / (2.0 * aR)) );

  //column four
  temp.SetData(0, 4, (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(1, 4, (uR + aR * areaNorm.X()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(2, 4, (vR + aR * areaNorm.Y()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(3, 4, (wR + aR * areaNorm.Z()) * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );
  temp.SetData(4, 4, (hR + aR * velRNorm)     * (eqnState.Gamma() - 1.0) / (2.0 * aR * aR) );

  A = A + (fabs(velRNorm + aR) * temp);

  //begin jacobian calculation ////////////////////////////////////////////////////////////////////////////////////////

  //derivative of Roe flux wrt left conservative variables
  dF_dUl.Zero();

  //column zero
  dF_dUl.SetData(0, 0, 0.0);
  dF_dUl.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.X() - left.U() * velLeftNorm);
  dF_dUl.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.Y() - left.V() * velLeftNorm);
  dF_dUl.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.Z() - left.W() * velLeftNorm);
  dF_dUl.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() - left.Enthalpy(eqnState)) * velLeftNorm); 
		       
  //column one
  dF_dUl.SetData(1, 0, areaNorm.X());
  dF_dUl.SetData(1, 1, left.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.X() + velLeftNorm);
  dF_dUl.SetData(1, 2, left.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Y());
  dF_dUl.SetData(1, 3, left.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Z());
  dF_dUl.SetData(1, 4, left.Enthalpy(eqnState) * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * velLeftNorm);

  //column two
  dF_dUl.SetData(2, 0, areaNorm.Y());
  dF_dUl.SetData(2, 1, left.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.X());
  dF_dUl.SetData(2, 2, left.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Y() + velLeftNorm);
  dF_dUl.SetData(2, 3, left.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Z());
  dF_dUl.SetData(2, 4, left.Enthalpy(eqnState) * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * velLeftNorm);

  //column three
  dF_dUl.SetData(3, 0, areaNorm.Z());
  dF_dUl.SetData(3, 1, left.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.X());
  dF_dUl.SetData(3, 2, left.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Y());
  dF_dUl.SetData(3, 3, left.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Z() + velLeftNorm);
  dF_dUl.SetData(3, 4, left.Enthalpy(eqnState) * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * velLeftNorm);

  //column four
  dF_dUl.SetData(4, 0, 0.0);
  dF_dUl.SetData(4, 1, (eqnState.Gamma() - 1.0) * areaNorm.X());
  dF_dUl.SetData(4, 2, (eqnState.Gamma() - 1.0) * areaNorm.Y());
  dF_dUl.SetData(4, 3, (eqnState.Gamma() - 1.0) * areaNorm.Z());
  dF_dUl.SetData(4, 4, eqnState.Gamma() * velLeftNorm);

  dF_dUl = 0.5 * (dF_dUl + A);


  //Compute derivative of flux wrt right conservative variables //////////////////////////////////////////////////////////////
  dF_dUr.Zero();

  //calculate flux derivatives
  //column zero
  dF_dUr.SetData(0, 0, 0.0);
  dF_dUr.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.X() - right.U() * velRightNorm);
  dF_dUr.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.Y() - right.V() * velRightNorm);
  dF_dUr.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.Z() - right.W() * velRightNorm);
  dF_dUr.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() - right.Enthalpy(eqnState)) * velRightNorm); 
		       
  //column one
  dF_dUr.SetData(1, 0, areaNorm.X());
  dF_dUr.SetData(1, 1, right.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.X() + velRightNorm);
  dF_dUr.SetData(1, 2, right.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Y());
  dF_dUr.SetData(1, 3, right.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Z());
  dF_dUr.SetData(1, 4, right.Enthalpy(eqnState) * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * velRightNorm);

  //column two
  dF_dUr.SetData(2, 0, areaNorm.Y());
  dF_dUr.SetData(2, 1, right.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.X());
  dF_dUr.SetData(2, 2, right.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Y() + velRightNorm);
  dF_dUr.SetData(2, 3, right.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Z());
  dF_dUr.SetData(2, 4, right.Enthalpy(eqnState) * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * velRightNorm);

  //column three
  dF_dUr.SetData(3, 0, areaNorm.Z());
  dF_dUr.SetData(3, 1, right.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.X());
  dF_dUr.SetData(3, 2, right.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Y());
  dF_dUr.SetData(3, 3, right.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Z() + velRightNorm);
  dF_dUr.SetData(3, 4, right.Enthalpy(eqnState) * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * velRightNorm);

  //column four
  dF_dUr.SetData(4, 0, 0.0);
  dF_dUr.SetData(4, 1, (eqnState.Gamma() - 1.0) * areaNorm.X());
  dF_dUr.SetData(4, 2, (eqnState.Gamma() - 1.0) * areaNorm.Y());
  dF_dUr.SetData(4, 3, (eqnState.Gamma() - 1.0) * areaNorm.Z());
  dF_dUr.SetData(4, 4, eqnState.Gamma() * velRightNorm);

  dF_dUr = 0.5 * (dF_dUr - A);

}

// //function to calculate Lax-Friedrichs flux jacobians
// void LaxFriedrichsFluxJacobian( const primVars &left, const primVars &right, const idealGas &eqnState, const vector3d<double>& areaVec, double &specRadL, double &specRadR, squareMatrix &dF_dUl, squareMatrix &dF_dUr){

//   //left --> primative variables from left side
//   //right --> primative variables from right side
//   //eqnStat --> ideal gas equation of state
//   //areaVec --> face area vector
//   //maxWS --> maximum wave speed
//   //dF_dUl --> dF/dUl, derivative of the Roe flux wrt the left state (conservative variables)
//   //dF_dUr --> dF/dUlr, derivative of the Roe flux wrt the right state (conservative variables)


//   //check to see that output matricies are correct size
//   if( (dF_dUl.Size() != 5) || (dF_dUr.Size() != 5)){
//     cerr << "ERROR: Input matricies to LaxFreidrichsFLuxJacobian function are not the correct size!" << endl;
//   }

//   vector3d<double> areaNorm = areaVec / areaVec.Mag();  //normalize area vector to unit vector

//   //dot product of velocities with unit area vector
//   double maxWS = ConvSpecRad(areaNorm, left, right, eqnState);

//   double velLeftNorm = left.Velocity().DotProd(areaNorm);
//   double velRightNorm = right.Velocity().DotProd(areaNorm);

//   //calculate spectral radii
//   specRadL = maxWS;
//   specRadR = maxWS;

//   //form spectral radii identity matrices
//   squareMatrix dissLeft(5);
//   dissLeft.Identity();
//   dissLeft = specRadL * dissLeft;

//   squareMatrix dissRight(5);
//   dissRight.Identity();
//   dissRight = specRadR * dissRight;

//   //begin jacobian calculation ////////////////////////////////////////////////////////////////////////////////////////

//   dF_dUl.Zero();
//   dF_dUr.Zero();

//   //calculate flux derivatives wrt left state
//   //column zero
//   dF_dUl.SetData(0, 0, 0.0);
//   dF_dUl.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.X() - left.U() * velLeftNorm);
//   dF_dUl.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.Y() - left.V() * velLeftNorm);
//   dF_dUl.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() * areaNorm.Z() - left.W() * velLeftNorm);
//   dF_dUl.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) * left.Velocity().MagSq() - left.Enthalpy(eqnState)) * velLeftNorm); 
		       
//   //column one
//   dF_dUl.SetData(0, 1, areaNorm.X());
//   dF_dUl.SetData(1, 1, left.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.X() + velLeftNorm);
//   dF_dUl.SetData(2, 1, left.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Y());
//   dF_dUl.SetData(3, 1, left.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * areaNorm.Z());
//   dF_dUl.SetData(4, 1, left.Enthalpy(eqnState) * areaNorm.X() - (eqnState.Gamma() - 1.0) * left.U() * velLeftNorm);

//   //column two
//   dF_dUl.SetData(0, 2, areaNorm.Y());
//   dF_dUl.SetData(1, 2, left.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.X());
//   dF_dUl.SetData(2, 2, left.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Y() + velLeftNorm);
//   dF_dUl.SetData(3, 2, left.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * areaNorm.Z());
//   dF_dUl.SetData(4, 2, left.Enthalpy(eqnState) * areaNorm.Y() - (eqnState.Gamma() - 1.0) * left.V() * velLeftNorm);

//   //column three
//   dF_dUl.SetData(0, 3, areaNorm.Z());
//   dF_dUl.SetData(1, 3, left.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.X());
//   dF_dUl.SetData(2, 3, left.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Y());
//   dF_dUl.SetData(3, 3, left.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * areaNorm.Z() + velLeftNorm);
//   dF_dUl.SetData(4, 3, left.Enthalpy(eqnState) * areaNorm.Z() - (eqnState.Gamma() - 1.0) * left.W() * velLeftNorm);

//   //column four
//   dF_dUl.SetData(0, 4, 0.0);
//   dF_dUl.SetData(1, 4, (eqnState.Gamma() - 1.0) * areaNorm.X());
//   dF_dUl.SetData(2, 4, (eqnState.Gamma() - 1.0) * areaNorm.Y());
//   dF_dUl.SetData(3, 4, (eqnState.Gamma() - 1.0) * areaNorm.Z());
//   dF_dUl.SetData(4, 4, eqnState.Gamma() * velLeftNorm);

//   dF_dUl = 0.5 * (dF_dUl + dissLeft);


//   //calculate flux derivatives wrt right state
//   //column zero
//   dF_dUr.SetData(0, 0, 0.0);
//   dF_dUr.SetData(1, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.X() - right.U() * velRightNorm);
//   dF_dUr.SetData(2, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.Y() - right.V() * velRightNorm);
//   dF_dUr.SetData(3, 0, 0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() * areaNorm.Z() - right.W() * velRightNorm);
//   dF_dUr.SetData(4, 0, (0.5 * (eqnState.Gamma() - 1.0) * right.Velocity().MagSq() - right.Enthalpy(eqnState)) * velRightNorm); 
		       
//   //column one
//   dF_dUr.SetData(0, 1, areaNorm.X());
//   dF_dUr.SetData(1, 1, right.U() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.X() + velRightNorm);
//   dF_dUr.SetData(2, 1, right.V() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Y());
//   dF_dUr.SetData(3, 1, right.W() * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * areaNorm.Z());
//   dF_dUr.SetData(4, 1, right.Enthalpy(eqnState) * areaNorm.X() - (eqnState.Gamma() - 1.0) * right.U() * velRightNorm);

//   //column two
//   dF_dUr.SetData(0, 2, areaNorm.Y());
//   dF_dUr.SetData(1, 2, right.U() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.X());
//   dF_dUr.SetData(2, 2, right.V() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Y() + velRightNorm);
//   dF_dUr.SetData(3, 2, right.W() * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * areaNorm.Z());
//   dF_dUr.SetData(4, 2, right.Enthalpy(eqnState) * areaNorm.Y() - (eqnState.Gamma() - 1.0) * right.V() * velRightNorm);

//   //column three
//   dF_dUr.SetData(0, 3, areaNorm.Z());
//   dF_dUr.SetData(1, 3, right.U() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.X());
//   dF_dUr.SetData(2, 3, right.V() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Y());
//   dF_dUr.SetData(3, 3, right.W() * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * areaNorm.Z() + velRightNorm);
//   dF_dUr.SetData(4, 3, right.Enthalpy(eqnState) * areaNorm.Z() - (eqnState.Gamma() - 1.0) * right.W() * velRightNorm);

//   //column four
//   dF_dUr.SetData(0, 4, 0.0);
//   dF_dUr.SetData(1, 4, (eqnState.Gamma() - 1.0) * areaNorm.X());
//   dF_dUr.SetData(2, 4, (eqnState.Gamma() - 1.0) * areaNorm.Y());
//   dF_dUr.SetData(3, 4, (eqnState.Gamma() - 1.0) * areaNorm.Z());
//   dF_dUr.SetData(4, 4, eqnState.Gamma() * velRightNorm);

//   dF_dUr = 0.5 * (dF_dUr - dissRight);


// }


//member function to return flux on boundaries
inviscidFlux BoundaryFlux( const string &bcName, const vector3d<double>& areaVec, const primVars &state1, const primVars &state2, const idealGas& eqnState, const input& inputVars, const string &surf, double &maxWS, const double upwind, const double upwind2 ){

  inviscidFlux flux;

  // vector3d<double> vel;
  //vector3d<double> velFace;

  vector3d<double> normArea = areaVec / areaVec.Mag();

  primVars state;

  if (bcName == "slipWall" || bcName == "viscousWall"){
    //state = (2.0 * state1) - state2;
    state = state1;
  }
  else{
    state = state1;
    //state = (2.0 * state1) - state2;
  }

  //Apply correct flux based on boundary condition to be applied 
  if ( bcName == "subsonicInflow" || bcName == "subsonicOutflow" || bcName == "supersonicInflow" || bcName == "supersonicOutflow" ||
       bcName == "characteristic" || bcName == "stagnationInlet" || bcName == "pressureOutlet"){

    if (inputVars.Kappa() == -2.0){ //first order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	lState = ghostState1.FaceReconConst();
	rState = state.FaceReconConst();
      }
      else {
	rState = ghostState1.FaceReconConst();
	lState = state.FaceReconConst();
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);
    }
    else{ //second order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars ghostState2 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState, 2 );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	//upwind cells are ghost cells, same distance as boundary cell
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind, upwind ); 
	//downwind cell is ghost cell, same distance as boundary cell
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind2, upwind ); 
      }
      else {
	//upwind cells are ghost cells, same distance as boundary cell
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind, upwind );
	//downwind cell is ghost cell, same distance as boundary cell
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind2, upwind );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

    }


  }
  else if ( (bcName == "slipWall") || (bcName == "viscousWall") ){

    if (inputVars.Kappa() == -2.0 ){ //first order
      primVars ghostState1 = state1.GetGhostState( "slipWall", normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state1.FaceReconConst();
	lState = ghostState1.FaceReconConst();
      }
      else {
	lState = state1.FaceReconConst();
	rState = ghostState1.FaceReconConst();
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);
    }
    else{ //second order
      primVars ghostState1 = state1.GetGhostState( "slipWall", normArea, surf, inputVars, eqnState );
      primVars ghostState2 = state2.GetGhostState( "slipWall", normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	//upwind cells are ghost cells, same distance as boundary cell
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind, upwind );
	//downwind cell is ghost cell, same distance as boundary cell
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind2, upwind );
      }
      else {
	//upwind cells are ghost cells, same distance as boundary cell
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind, upwind );
	//downwind cell is ghost cell, same distance as boundary cell
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), upwind, upwind2, upwind );

      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

    }

  }
  else{
    cerr << "ERROR: Boundary condition " << bcName << " is not recognized!" << endl;
  }

  return flux;

}

// double BoundaryInvSpecRad( const string &bcName, const vector3d<double>& areaVec, const primVars &state, const idealGas& eqnState, const string &surf, const input &inputVars){

//   double maxWS = 0.0;

//   vector3d<double> normArea = areaVec / areaVec.Mag();

//   //Apply correct flux based on boundary condition to be applied 
//   if ( bcName == "subsonicInflow" || bcName == "subsonicOutflow" || bcName == "supersonicInflow" || bcName == "supersonicOutflow" || bcName == "characteristic"){

//     primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
//     maxWS = ConvSpecRad(normArea, state, ghostState1, eqnState);

//   }
//   else if ( bcName == "slipWall" || "viscousWall" ){

//     maxWS = state.SoS(eqnState);

//   }
//   else{
//     cerr << "ERROR: Boundary condition " << bcName << " is not recognized!" << endl;
//   }

//   return maxWS;

// }


//non-member functions -----------------------------------------------------------------------------------------------------------//

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, inviscidFlux &flux){

  os << flux.rhoVel << "   " << flux.rhoVelU << "   " << flux.rhoVelV << "   " << flux.rhoVelW << "   " << flux.rhoVelH << endl;

  return os;
}

//member function for scalar multiplication
inviscidFlux  inviscidFlux::operator * (const double &scalar){
  inviscidFlux temp = *this;
  temp.rhoVel *= scalar;
  temp.rhoVelU *= scalar;
  temp.rhoVelV *= scalar;
  temp.rhoVelW *= scalar;
  temp.rhoVelH *= scalar;
  return temp;
}

//friend function to allow multiplication from either direction
inviscidFlux operator* (const double &scalar, const inviscidFlux &flux){
  inviscidFlux temp;
  temp.SetRhoVel(flux.RhoVel() * scalar);
  temp.SetRhoVelU(flux.RhoVelU() * scalar);
  temp.SetRhoVelV(flux.RhoVelV() * scalar);
  temp.SetRhoVelW(flux.RhoVelW() * scalar);
  temp.SetRhoVelH(flux.RhoVelH() * scalar);
  return temp;
}

//operator overload for addition
inviscidFlux inviscidFlux::operator + (const inviscidFlux& invf2)const{
  inviscidFlux invf1 = *this;
  invf1.rhoVel += invf2.rhoVel;
  invf1.rhoVelU += invf2.rhoVelU;
  invf1.rhoVelV += invf2.rhoVelV;
  invf1.rhoVelW += invf2.rhoVelW;
  invf1.rhoVelH += invf2.rhoVelH;
  return invf1;
}

//operator overload for subtraction
inviscidFlux inviscidFlux::operator - (const inviscidFlux& invf2)const{
  inviscidFlux invf1 = *this;
  invf1.rhoVel -= invf2.rhoVel;
  invf1.rhoVelU -= invf2.rhoVelU;
  invf1.rhoVelV -= invf2.rhoVelV;
  invf1.rhoVelW -= invf2.rhoVelW;
  invf1.rhoVelH -= invf2.rhoVelH;
  return invf1;
}


//member function for scalar division
inviscidFlux  inviscidFlux::operator / (const double &scalar){
  inviscidFlux temp = *this;
  temp.rhoVel /= scalar;
  temp.rhoVelU /= scalar;
  temp.rhoVelV /= scalar;
  temp.rhoVelW /= scalar;
  temp.rhoVelH /= scalar;
  return temp;
}

//friend function to allow division from either direction
inviscidFlux operator/ (const double &scalar, const inviscidFlux &flux){
  inviscidFlux temp;
  temp.SetRhoVel(scalar / flux.RhoVel());
  temp.SetRhoVelU(scalar / flux.RhoVelU());
  temp.SetRhoVelV(scalar / flux.RhoVelV());
  temp.SetRhoVelW(scalar / flux.RhoVelW());
  temp.SetRhoVelH(scalar / flux.RhoVelH());
  return temp;
}


//convert the inviscid flux to a column matrix
colMatrix inviscidFlux::ConvertToColMatrix()const{
  colMatrix temp(5);
  temp.SetData(0, (*this).RhoVel());
  temp.SetData(1, (*this).RhoVelU());
  temp.SetData(2, (*this).RhoVelV());
  temp.SetData(3, (*this).RhoVelW());
  temp.SetData(4, (*this).RhoVelH());
  return temp;
}

//function to take in the primative variables, equation of state, face area vector, and conservative variable update and calculate the change in the convective flux
colMatrix ConvectiveFluxUpdate( const primVars &state, const idealGas &eqnState, const vector3d<double> &fArea, const colMatrix &du){

  inviscidFlux oldFlux(state, eqnState, fArea);

  primVars stateUpdate = state.UpdateWithConsVars(eqnState, du);
  inviscidFlux newFlux(stateUpdate, eqnState, fArea);

  inviscidFlux dFlux = newFlux - oldFlux;
    
  return dFlux.ConvertToColMatrix();
}


// inviscidFlux LaxFriedrichsFlux( const primVars &left, const primVars &right, const idealGas &eqnState, const vector3d<double> &fArea, double &maxWS ){

//   inviscidFlux lFlux(left, eqnState, fArea);
//   inviscidFlux rFlux(right, eqnState, fArea);

//   maxWS = ConvSpecRad(fArea, left, right, eqnState);

//   colMatrix lCons = left.ConsVars(eqnState);
//   colMatrix rCons = right.ConsVars(eqnState);

//   inviscidFlux lfFlux;
//   lfFlux.SetRhoVel(  0.5 * (rFlux.RhoVel()  + lFlux.RhoVel()  - maxWS * (rCons.Data(0) - lCons.Data(0)) ) );
//   lfFlux.SetRhoVelU( 0.5 * (rFlux.RhoVelU() + lFlux.RhoVelU() - maxWS * (rCons.Data(1) - lCons.Data(1)) ) );
//   lfFlux.SetRhoVelV( 0.5 * (rFlux.RhoVelV() + lFlux.RhoVelV() - maxWS * (rCons.Data(2) - lCons.Data(2)) ) );
//   lfFlux.SetRhoVelW( 0.5 * (rFlux.RhoVelW() + lFlux.RhoVelW() - maxWS * (rCons.Data(3) - lCons.Data(3)) ) );
//   lfFlux.SetRhoVelH( 0.5 * (rFlux.RhoVelH() + lFlux.RhoVelH() - maxWS * (rCons.Data(4) - lCons.Data(4)) ) );

//   return lfFlux;

// }

// //function to return the convective spectral radius given a face state, equation of state, and area vector
// double ConvSpecRad(const vector3d<double> &fArea, const primVars &state1, const primVars &state2, const idealGas &eqnState){

//   vector3d<double> normArea = fArea / fArea.Mag();

//   double a1 = state1.SoS(eqnState);
//   double u1 = state1.Velocity().DotProd(normArea);

//   double a2 = state2.SoS(eqnState);
//   double u2 = state2.Velocity().DotProd(normArea);

//   return max(fabs(u1) + a1, fabs(u2) + a2);

// }


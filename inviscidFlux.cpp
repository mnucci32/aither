#include "inviscidFlux.h"
#include <cmath> //sqrt

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::max;

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

//function to calculate pressure from conserved variables and equation of state
inviscidFlux RoeFluxJacobian( const primVars &left, const primVars &right, const idealGas &eqnState, const vector3d<double>& areaVec, double &maxWS){

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

  double normVelDiff = right.Velocity().DotProd(areaNorm) - left.Velocity().DotProd(areaNorm);

  double waveStrength[4] = {((right.P() - left.P()) - rhoR * aR * normVelDiff) / (2.0 * aR * aR), 
			    (right.Rho() - left.Rho()) - (right.P() - left.P()) / (aR * aR), 
			    ((right.P() - left.P()) + rhoR * aR * normVelDiff) / (2.0 * aR * aR), 
			    rhoR};

  //left moving acoustic wave speed, entropy wave speed, right moving acoustic wave speed, shear wave speed
  double waveSpeed[4] = {fabs(velRSum - aR), 
			 fabs(velRSum), 
			 fabs(velRSum + aR), 
			 fabs(velRSum)};

  //calculate entropy fix
  double entropyFix = 0.1;                                                            // default setting for entropy fix to kick in

  if ( waveSpeed[0] < entropyFix ){
    waveSpeed[0] = 0.5 * (waveSpeed[0] * waveSpeed[0] / entropyFix + entropyFix);
  }
  if ( waveSpeed[2] < entropyFix ){
    waveSpeed[2] = 0.5 * (waveSpeed[2] * waveSpeed[2] / entropyFix + entropyFix);
  }

  maxWS = fabs(velRSum) + aR;

  double lAcousticEigV[5] = {1.0, 
			     uR - aR * areaNorm.X(), 
			     vR - aR * areaNorm.Y(), 
			     wR - aR * areaNorm.Z(), 
			     hR - aR * velRSum};

  double entropyEigV[5] = {1.0, 
			   uR, 
			   vR, 
			   wR, 
			   0.5 * ( uR * uR + vR * vR + wR * wR)};

  double rAcousticEigV[5] = {1.0, 
			     uR + aR * areaNorm.X(), 
			     vR + aR * areaNorm.Y(), 
			     wR + aR * areaNorm.Z(), 
			     hR + aR * velRSum};

  double shearEigV[5] = {0.0, 
			 (right.U() - left.U()) - normVelDiff * areaNorm.X(), 
			 (right.V() - left.V()) - normVelDiff * areaNorm.Y(), 
			 (right.W() - left.W()) - normVelDiff * areaNorm.Z(), 
			 uR * (right.U() - left.U()) + vR * (right.V() - left.V()) + wR * (right.W() - left.W()) - velRSum * normVelDiff};

  //begin jacobian calculation



  //calculate dissipation term ( eigenvector * wave speed * wave strength)
  double dissipation[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  unsigned int ii = 0;
  for ( ii=0; ii < 5; ii++ ) {                                    
    dissipation[ii] = waveSpeed[0] * waveStrength[0] * lAcousticEigV[ii]   //contribution from left acoustic wave
                    + waveSpeed[1] * waveStrength[1] * entropyEigV[ii]     //contribution from entropy wave
                    + waveSpeed[2] * waveStrength[2] * rAcousticEigV[ii]   //contribution from right acoustic wave
                    + waveSpeed[3] * waveStrength[3] * shearEigV[ii];      //contribution from shear wave
  }


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

//member function to return flux on boundaries
inviscidFlux BoundaryFlux( const string &bcName, const vector3d<double>& areaVec, const primVars &state1, const primVars &state2, const idealGas& eqnState, const input& inputVars, const string &surf, double &maxWS, const double up2face, const double upwind ){

  inviscidFlux flux;

  // vector3d<double> vel;
  //vector3d<double> velFace;

  vector3d<double> normArea = areaVec / areaVec.Mag();

  primVars state;

  if (bcName == "slipWall"){
    //state = (2.0 * state1) - (1.0 * state2);
    state = state1;
  }
  else{
    state = state1;
    //state = (2.0 * state1) - (1.0 * state2);
  }

  //Apply correct flux based on boundary condition to be applied 
  if ( bcName == "subsonicInflow" ){

    if (inputVars.Kappa() == -2.0 ){ //first order
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
      primVars ghostState2 = ghostState1;
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

    }


  }
  else if ( bcName == "subsonicOutflow" ){

    if (inputVars.Kappa() == -2.0 ){ //first order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state.FaceReconConst();
	lState = ghostState1.FaceReconConst();
      }
      else {
	lState = state.FaceReconConst();
	rState = ghostState1.FaceReconConst();
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);
    }
    else{ //second order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars ghostState2 = ghostState1;
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS); 

    }

  }
  else if ( bcName == "supersonicInflow" ){

    if (inputVars.Kappa() == -2.0 ){ //first order
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
      primVars ghostState2 = ghostState1;
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

    }

    // //physical boundary conditions 
    // pressure = 1.0 / eqnState.Gamma(); //set pressure, velocity, and energy equal to nondimensional freestream value
    // rho = 1.0;
    // sos = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef());
    // specEnergy = 1.0 / (eqnState.Gamma() * (eqnState.Gamma() - 1.0));
    // vel = inputVars.VelRef() / sos;                       //nondimensional velocity
    // energy = eqnState.GetEnergy(specEnergy,vel.Mag());
    // enthalpy = eqnState.GetEnthalpy(energy, pressure, rho);
    // //numerical boundary conditions - none

    // //calculate flux
    // flux.SetRhoVel(rho * vel.DotProd(normArea));
    // flux.SetRhoVelU(rho * vel.DotProd(normArea) * vel.X() + pressure * normArea.X());
    // flux.SetRhoVelV(rho * vel.DotProd(normArea) * vel.Y() + pressure * normArea.Y());
    // flux.SetRhoVelW(rho * vel.DotProd(normArea) * vel.Z() + pressure * normArea.Z());
    // flux.SetRhoVelH(rho * vel.DotProd(normArea) * enthalpy);

    // //Calculate max wave speed
    // sos = eqnState.GetSoS(pressure, rho);
    // maxWS = fabs(vel.DotProd(normArea)) + sos;

  }
  else if ( bcName == "supersonicOutflow" ){

    if (inputVars.Kappa() == -2.0 ){ //first order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state.FaceReconConst();
	lState = ghostState1.FaceReconConst();
      }
      else {
	lState = state.FaceReconConst();
	rState = ghostState1.FaceReconConst();
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);
    }
    else{ //second order
      primVars ghostState1 = state.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars ghostState2 = ghostState1;
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS); 

    }

    // //physical boundary conditions - none

    // //numerical boundary conditions - set pressure, rho, and velocity at face equal to values at cell center of cell on boundary
    // rho = state.Rho();
    // pressure = state.Pressure(eqnState);
    // vel.SetX(state.Vx());
    // vel.SetY(state.Vy());
    // vel.SetZ(state.Vz());
    // enthalpy = state.Enthalpy(eqnState);

    // //calculate flux
    // flux.SetRhoVel(rho * vel.DotProd(normArea));
    // flux.SetRhoVelU(rho * vel.DotProd(normArea) * vel.X() + pressure * normArea.X());
    // flux.SetRhoVelV(rho * vel.DotProd(normArea) * vel.Y() + pressure * normArea.Y());
    // flux.SetRhoVelW(rho * vel.DotProd(normArea) * vel.Z() + pressure * normArea.Z());
    // flux.SetRhoVelH(rho * vel.DotProd(normArea) * enthalpy);

    // //calculate max wave speed
    // sos = eqnState.GetSoS(pressure, rho);
    // maxWS = fabs(vel.DotProd(normArea)) + sos;

  }
  else if ( bcName == "slipWall" || "viscousWall" ){

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
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      // cout << "at slip wall " << endl;
      // cout << "ghost state one " << ghostState1;
      // cout << "ghost state two " << ghostState2;
      // cout << "left state " << lState;
      // cout << "right state " << rState;

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

      // cout << flux;

    }

  }
  else if ( bcName == "characteristic" ){

    if (inputVars.Kappa() == -2.0 ){ //first order
      primVars ghostState1 = state1.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
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
      primVars ghostState1 = state1.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars ghostState2 = state2.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
      primVars lState, rState;

      if (surf == "il" || surf == "jl" || surf == "kl"){
	rState = state1.FaceReconMUSCL( state2, ghostState1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	lState = ghostState1.FaceReconMUSCL( ghostState2, state1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }
      else {
	lState = state1.FaceReconMUSCL( state2, ghostState1, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
	rState = ghostState1.FaceReconMUSCL( ghostState2, state1, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
      }

      flux = RoeFlux( lState, rState, eqnState, normArea, maxWS);

    }

  //   if (inputVars.Kappa() == -2.0 ){  //first order
  //     //primVars ghostState1 = state1.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
  //     primVars lState, rState;

  //     rState = state1.FaceReconConst();
  //     vector3d<double> boundaryNormArea;
  //     if (surf == "iu" || surf == "ju" || surf == "ku"){
  // 	boundaryNormArea = normArea * -1.0;
  //     }
  //     else{
  // 	boundaryNormArea = normArea;
  //     }

  //     //calculate ui, u0, ci, co, mi
  //     vector3d<double> freeVel = inputVars.VelRef();
  //     double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
  //     freeVel = freeVel / freeSoS;
  //     freeSoS = freeSoS / freeSoS;

  //     double freeSpecEnergy = eqnState.GetSpecEnergy(1.0/eqnState.Gamma(), 1.0);
  //     double freeEnergy = eqnState.GetEnergy(freeSpecEnergy, freeVel.Mag());

  //     primVars freeState(1.0, freeEnergy, freeVel);
  //     double rhoFree = freeState.Rho();
  //     double velFreeNorm = freeVel.DotProd(boundaryNormArea); //CHANGED from boundaryNormArea
  //     double SoSFree = 1.0;

  //     vector3d<double> freeVelTan;                                  //freestream tangent velocity
  //     freeVelTan.SetX(freeVel.X() - boundaryNormArea.X()*velFreeNorm);
  //     freeVelTan.SetY(freeVel.Y() - boundaryNormArea.Y()*velFreeNorm);
  //     freeVelTan.SetZ(freeVel.Z() - boundaryNormArea.Z()*velFreeNorm);

  //     vector3d<double> velInt;
  //     velInt.SetX(state1.Vx());
  //     velInt.SetY(state1.Vy());
  //     velInt.SetZ(state1.Vz());
  //     double velIntNorm;

  //     //if (surf == "iu" || surf == "ju" || surf == "ku"){
  //     //velIntNorm = velInt.DotProd(velIntDir*-1.0); //CHANGED from boundaryNormArea
  // 	//}
  // 	//else{
  // 	velIntNorm = velInt.DotProd(boundaryNormArea); //CHANGED from boundaryNormArea
  // 	//}
      
  //     double rhoInt = state1.Rho();
  //     double SoSInt = eqnState.GetSoS(state1.Pressure(eqnState), rhoInt);

  //     vector3d<double> velIntTan;                                  //internal tangent velocity
  //     velIntTan.SetX(velInt.X() - boundaryNormArea.X()*velIntNorm);
  //     velIntTan.SetY(velInt.Y() - boundaryNormArea.Y()*velIntNorm);
  //     velIntTan.SetZ(velInt.Z() - boundaryNormArea.Z()*velIntNorm);

  //     double machInt = fabs(velIntNorm)/SoSInt;
  //     //cout << "mach interior " << machInt << endl;

  //     double riemannInvarPlus, riemannInvarMinus;

  //     if ( machInt >= 1.0 && velIntNorm > 0.0 ){ //supersonic inflow

  // 	//cout << "supersonic inflow bc applied" << endl;

  // 	//characteristics all go into the domain, so use freestream values for both riemann invariants
  // 	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarMinus = velFreeNorm - (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	lState = freeState;
  //     }
  //     else if ( machInt >= 1.0 && velIntNorm < 0.0 ){ //supersonic outflow

  // 	//cout << "supersonic outflow bc applied" << endl;

  // 	//characteristics all leave the domain, so use interior values for both riemann invariants
  // 	riemannInvarPlus = velIntNorm + (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	lState = state1;
  //     }
  //     else if ( machInt < 1.0 && velIntNorm > 0.0 ){ //subsonic inflow

  // 	//cout << "subsonic inflow bc applied" << endl;

  // 	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus); //this is the normal velocity at the boundary
  // 	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);
  // 	//double SoSBound = sqrt(0.5 * (eqnState.Gamma() - 1.0) * (freeState.Enthalpy(eqnState) - velB*velB) );

  // 	double entropyBound = (SoSFree * SoSFree)/(eqnState.Gamma() * pow(rhoFree , eqnState.Gamma()-1.0));
  // 	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
  // 	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
  // 	vector3d<double> velBound;
  // 	// velBound.SetX( freeVel.X() + freeVelDir.X() * (velB - velFreeNorm) );  //changed boundaryNormArea.X()
  // 	// velBound.SetY( freeVel.Y() + freeVelDir.Y() * (velB - velFreeNorm) );
  // 	// velBound.SetZ( freeVel.Z() + freeVelDir.Z() * (velB - velFreeNorm) );

  // 	velBound.SetX( (velB * boundaryNormArea.X() + freeVelTan.X()) );  
  // 	velBound.SetY( (velB * boundaryNormArea.Y() + freeVelTan.Y()) );
  // 	velBound.SetZ( (velB * boundaryNormArea.Z() + freeVelTan.Z()) );

  // 	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
  // 	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

  // 	lState.SetRho(rhoBound);
  // 	lState.SetRhoU(rhoBound * velBound.X());
  // 	lState.SetRhoV(rhoBound * velBound.Y());
  // 	lState.SetRhoW(rhoBound * velBound.Z());
  // 	lState.SetRhoE(rhoBound * energyBound);

  // 	// vector3d<double> dummyVec;
  // 	// dummyVec.SetX(freeVelDir.X() * boundaryNormArea.X() * velB);
  // 	// dummyVec.SetY(freeVelDir.Y() * boundaryNormArea.Y() * velB);
  // 	// dummyVec.SetZ(freeVelDir.Z() * boundaryNormArea.Z() * velB);

  // 	// cout << "boundary velocity " << velBound << endl;
  // 	// cout << "vel norm " << velB*boundaryNormArea << endl;

  // 	// dummyVec.SetX(freeVelTan.X() * freeVelDir.X());
  // 	// dummyVec.SetY(freeVelTan.Y() * freeVelDir.Y());
  // 	// dummyVec.SetZ(freeVelTan.Z() * freeVelDir.Z());

  // 	// cout << "free vel tan " << freeVelTan << endl;
  // 	// cout << "vel norm dot vel tan " << freeVelTan.DotProd(velFreeNorm*boundaryNormArea) << endl;
  // 	// cout << "boundary density " << rhoBound << endl;
  // 	// cout << "riemann vel " << velB << endl;
  // 	// cout << "riemann sos " << SoSBound << endl;
  // 	// cout << "entropy bound " << entropyBound << endl;
  // 	// cout << "riemann invars " << riemannInvarPlus << ", " << riemannInvarMinus << endl;
  // 	// cout << "norm vel free, int " << velFreeNorm << ", " << velIntNorm << endl;

  //     }
  //     else if ( machInt < 1.0 && velIntNorm < 0.0 ){ //subsonic outflow

  // 	//cout << "subsonic outflow bc applied" << endl;

  // 	//velIntNorm = fabs(velIntNorm) * -1.0;
  // 	//velIntDir = velIntDir * -1.0;

  // 	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarPlus = -1.0 * fabs(velFreeNorm) + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus);
  // 	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);

  // 	double entropyBound = (SoSInt * SoSInt)/(eqnState.Gamma() * pow(rhoInt , eqnState.Gamma()-1.0));
  // 	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
  // 	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
  // 	vector3d<double> velBound;
  // 	// velBound.SetX( velInt.X() + velIntDir.X() * (velB - velIntNorm) );  //changed from boundaryNormArea.X()
  // 	// velBound.SetY( velInt.Y() + velIntDir.Y() * (velB - velIntNorm) );
  // 	// velBound.SetZ( velInt.Z() + velIntDir.Z() * (velB - velIntNorm) );

  // 	velBound.SetX( (velB * boundaryNormArea.X() + velIntTan.X()) );
  // 	velBound.SetY( (velB * boundaryNormArea.Y() + velIntTan.Y()) );
  // 	velBound.SetZ( (velB * boundaryNormArea.Z() + velIntTan.Z()) );

  // 	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
  // 	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

  // 	lState.SetRho(rhoBound);
  // 	lState.SetRhoU(rhoBound * velBound.X());
  // 	lState.SetRhoV(rhoBound * velBound.Y());
  // 	lState.SetRhoW(rhoBound * velBound.Z());
  // 	lState.SetRhoE(rhoBound * energyBound);

  // 	// vector3d<double> dummyVec;
  // 	// dummyVec.SetX(velIntDir.X() * boundaryNormArea.X() * velB);
  // 	// dummyVec.SetY(velIntDir.Y() * boundaryNormArea.Y() * velB);
  // 	// dummyVec.SetZ(velIntDir.Z() * boundaryNormArea.Z() * velB);

  // 	// cout << "boundary velocity " << velBound << endl;
  // 	// cout << "int vel norm " << velB*boundaryNormArea << endl;

  // 	// dummyVec.SetX(velIntTan.X() * velIntDir.X());
  // 	// dummyVec.SetY(velIntTan.Y() * velIntDir.Y());
  // 	// dummyVec.SetZ(velIntTan.Z() * velIntDir.Z());

  // 	// cout << "int vel tan " << velIntTan << endl;
  // 	//cout << "vel norm dot vel tan " << velIntTan.DotProd(velIntNorm*boundaryNormArea) << endl;


  // 	// cout << "boundary velocity " << velBound << endl;
  // 	// cout << "boundary density " << rhoBound << endl;


  //     }
  //     else {
  // 	cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
  // 	exit(0);
  //     }

  //     flux = RoeFlux( lState, rState, eqnState, normArea, maxWS); //roe flux needs original area
  //   }

  //   else{  //2nd order
  //     //primVars ghostState1 = state1.GetGhostState( bcName, normArea, surf, inputVars, eqnState );
  //     primVars lState, rState;

  //     //rState = state1.FaceReconConst();

  //     vector3d<double> boundaryNormArea;
  //     if (surf == "iu" || surf == "ju" || surf == "ku"){
  // 	boundaryNormArea = normArea * -1.0;
  //     }
  //     else{
  // 	boundaryNormArea = normArea;
  //     }

  //     //calculate ui, u0, ci, co, mi
  //     vector3d<double> freeVel = inputVars.VelRef();
  //     double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
  //     freeVel = freeVel / freeSoS;
  //     freeSoS = freeSoS / freeSoS;

  //     double freeSpecEnergy = eqnState.GetSpecEnergy(1.0/eqnState.Gamma(), 1.0);
  //     double freeEnergy = eqnState.GetEnergy(freeSpecEnergy, freeVel.Mag());

  //     primVars freeState(1.0, freeEnergy, freeVel);
  //     double rhoFree = freeState.Rho();
  //     double velFreeNorm = freeVel.DotProd(boundaryNormArea);
  //     double SoSFree = 1.0;

  //     vector3d<double> freeVelTan;                                  //freestream tangent velocity
  //     freeVelTan.SetX(freeVel.X() - boundaryNormArea.X()*velFreeNorm);
  //     freeVelTan.SetY(freeVel.Y() - boundaryNormArea.Y()*velFreeNorm);
  //     freeVelTan.SetZ(freeVel.Z() - boundaryNormArea.Z()*velFreeNorm);

  //     vector3d<double> velInt;
  //     velInt.SetX(state.Vx());
  //     velInt.SetY(state.Vy());
  //     velInt.SetZ(state.Vz());
  //     double velIntNorm;

  //     //if (surf == "iu" || surf == "ju" || surf == "ku"){
  //     //velIntNorm = velInt.DotProd(velIntDir*-1.0); //CHANGED from boundaryNormArea
  // 	//}
  // 	//else{
  // 	velIntNorm = velInt.DotProd(boundaryNormArea); //CHANGED from boundaryNormArea
  // 	//}

  //     double rhoInt = state.Rho();
  //     double SoSInt = eqnState.GetSoS(state.Pressure(eqnState), rhoInt);

  //     vector3d<double> velIntTan;                                  //freestream tangent velocity
  //     velIntTan.SetX(velInt.X() - boundaryNormArea.X()*velIntNorm);
  //     velIntTan.SetY(velInt.Y() - boundaryNormArea.Y()*velIntNorm);
  //     velIntTan.SetZ(velInt.Z() - boundaryNormArea.Z()*velIntNorm);

  //     double machInt = fabs(velIntNorm)/SoSInt;
  //     //cout << "mach interior " << machInt << endl;

  //     double riemannInvarPlus, riemannInvarMinus;

  //     if ( machInt >= 1.0 && velIntNorm > 0.0 ){ //supersonic inflow

  // 	//cout << "supersonic inflow bc applied" << endl;

  // 	//characteristics all go into the domain, so use freestream values for both riemann invariants
  // 	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarMinus = velFreeNorm - (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	lState = freeState;
  //     }
  //     else if ( machInt >= 1.0 && velIntNorm < 0.0 ){ //supersonic outflow

  // 	//cout << "supersonic outflow bc applied" << endl;

  // 	//characteristics all leave the domain, so use interior values for both riemann invariants
  // 	riemannInvarPlus = velIntNorm + (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	lState = state1;
  //     }
  //     else if ( machInt < 1.0 && velIntNorm > 0.0 ){ //subsonic inflow

  // 	//cout << "subsonic inflow bc applied" << endl;


  // 	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus);
  // 	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);

  // 	double entropyBound = (SoSFree * SoSFree)/(eqnState.Gamma() * pow(rhoFree , eqnState.Gamma()-1.0));
  // 	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
  // 	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
  // 	vector3d<double> velBound;
  // 	// velBound.SetX( freeVel.X() + freeVelDir.X() * (velB - velFreeNorm) );  //changed boundaryNormArea.X()
  // 	// velBound.SetY( freeVel.Y() + freeVelDir.Y() * (velB - velFreeNorm) );
  // 	// velBound.SetZ( freeVel.Z() + freeVelDir.Z() * (velB - velFreeNorm) );

  // 	velBound.SetX( (velB * boundaryNormArea.X() + freeVelTan.X()) );  
  // 	velBound.SetY( (velB * boundaryNormArea.Y() + freeVelTan.Y()) );
  // 	velBound.SetZ( (velB * boundaryNormArea.Z() + freeVelTan.Z()) );

  // 	//cout << "free velocity direction, mag " << freeVelDir << ", " << freeVel.Mag() << endl;
  // 	//cout << "boundary velocity direction, mag " << velBound / velBound.Mag() << ", " << velBound.Mag() << endl;

  // 	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
  // 	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

  // 	lState.SetRho(rhoBound);
  // 	lState.SetRhoU(rhoBound * velBound.X());
  // 	lState.SetRhoV(rhoBound * velBound.Y());
  // 	lState.SetRhoW(rhoBound * velBound.Z());
  // 	lState.SetRhoE(rhoBound * energyBound);

  // 	// cout << "boundary velocity " << velBound << endl;
  // 	// cout << "boundary density " << rhoBound << endl;
  // 	// cout << "riemann vel " << velB << endl;
  // 	// cout << "riemann sos " << SoSBound << endl;
  // 	// cout << "entropy bound " << entropyBound << endl;
  // 	// cout << "riemann invars " << riemannInvarPlus << ", " << riemannInvarMinus << endl;
  // 	// cout << "norm vel free, int " << velFreeNorm << ", " << velIntNorm << endl;

  //     }
  //     else if ( machInt < 1.0 && velIntNorm < 0.0 ){ //subsonic outflow

  // 	// cout << "subsonic outflow bc applied" << endl;

  // 	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
  // 	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
  // 	riemannInvarPlus = -1.0 * fabs(velFreeNorm) + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
  // 	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus);
  // 	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);

  // 	double entropyBound = (SoSInt * SoSInt)/(eqnState.Gamma() * pow(rhoInt , eqnState.Gamma()-1.0));
  // 	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
  // 	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
  // 	vector3d<double> velBound;
  // 	// velBound.SetX( velInt.X() + velIntDir.X() * (velB - velIntNorm)  );
  // 	// velBound.SetY( velInt.Y() + velIntDir.Y() * (velB - velIntNorm)  );
  // 	// velBound.SetZ( velInt.Z() + velIntDir.Z() * (velB - velIntNorm)  );

  // 	velBound.SetX( (velB * boundaryNormArea.X() + velIntTan.X()) );
  // 	velBound.SetY( (velB * boundaryNormArea.Y() + velIntTan.Y()) );
  // 	velBound.SetZ( (velB * boundaryNormArea.Z() + velIntTan.Z()) );


  // 	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
  // 	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

  // 	lState.SetRho(rhoBound);
  // 	lState.SetRhoU(rhoBound * velBound.X());
  // 	lState.SetRhoV(rhoBound * velBound.Y());
  // 	lState.SetRhoW(rhoBound * velBound.Z());
  // 	lState.SetRhoE(rhoBound * energyBound);

  // 	// cout << "boundary velocity " << velBound << endl;
  // 	// cout << "boundary density " << rhoBound << endl;


  //     }
  //     else {
  // 	cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
  // 	exit(0);
  //     }

  //     if (surf == "il" || surf == "jl" || surf == "kl"){
  // 	rState = state1.FaceReconMUSCL( state2, lState, "right", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
  //     }
  //     else {
  // 	//on upper surface "left" state is determined by characteristics above, and "right" state actually comes from lower indices and the left
  // 	rState = state1.FaceReconMUSCL( state2, lState, "left", inputVars.Kappa(), inputVars.Limiter(), up2face, upwind, up2face*2.0 );
  //     }


  //     flux = RoeFlux( lState, rState, eqnState, normArea, maxWS); //roe flux needs original area

  //   }

  }
  
  else{
    cerr << "ERROR: Boundary condition " << bcName << " is not recognized!" << endl;
  }

  return flux;

}


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

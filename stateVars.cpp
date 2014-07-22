#include "stateVars.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::max;
using std::min;

//constructors
// stateVars::stateVars(){
//   rho = 0.0;      
//   rhoU = 0.0;     
//   rhoV = 0.0;     
//   rhoW = 0.0;
//   rhoE = 0.0;     
// }

//constructor -- assign passed variables to initialize state vector
// stateVars::stateVars( double density, double energy, vector3d<double> vel ){

//   rho = density;      
//   rhoU = density * vel.X();     
//   rhoV = density * vel.Y();     
//   rhoW = density * vel.Z();
//   rhoE = density * energy;     

// }

//member function to advance the state vector to time n+1 using explicit Euler method
// stateVars stateVars::ExplicitEulerTimeAdvance( double dt, double vol, vector<double>& resid ){

//   stateVars newState;

//   newState.SetRho( rho + dt / vol * resid[0] );
//   newState.SetRhoU( rhoU + dt / vol * resid[1] );
//   newState.SetRhoV( rhoV + dt / vol * resid[2] );
//   newState.SetRhoW( rhoW + dt / vol * resid[3] );
//   newState.SetRhoE( rhoE + dt / vol * resid[4] );

//   return newState;
// }

// //member function to advance the state vector to time n+1 using 4th order Runge-Kutta method
// stateVars stateVars::RK4TimeAdvance( double dt, double vol, vector<double>& resid, int rk ){

//   stateVars newState;

//   vector<double> alpha(4,0.0);  //runge-kutta coefficients
//   alpha[0] = 0.25;
//   alpha[1] = 1.0/3.0;
//   alpha[2] = 0.5;
//   alpha[3] = 1.0;

//   newState.SetRho(  rho  + dt / vol * alpha[rk] * resid[0] );
//   newState.SetRhoU( rhoU + dt / vol * alpha[rk] * resid[1] );
//   newState.SetRhoV( rhoV + dt / vol * alpha[rk] * resid[2] );
//   newState.SetRhoW( rhoW + dt / vol * alpha[rk] * resid[3] );
//   newState.SetRhoE( rhoE + dt / vol * alpha[rk] * resid[4] );

//   return newState;
// }


//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const stateVars &st){

  os << st.rho << "   " << st.rhoU << "   " << st.rhoV << "   " << st.rhoW << "   " << st.rhoE << endl;

  return os;
}

//operator overload for addition
stateVars stateVars::operator + (const stateVars& s2)const{
  stateVars s1 = *this;
  s1.rho += s2.rho;
  s1.rhoU += s2.rhoU;
  s1.rhoV += s2.rhoV;
  s1.rhoW += s2.rhoW;
  s1.rhoE += s2.rhoE;
  return s1;
}

//operator overload for addition with a scalar
stateVars operator+ (const double &scalar, const stateVars &s2){
  stateVars s1;
  s1.rho = s2.rho + scalar;
  s1.rhoU = s2.rhoU + scalar;
  s1.rhoV = s2.rhoV + scalar;
  s1.rhoW = s2.rhoW + scalar;
  s1.rhoE = s2.rhoE + scalar;
  return s1;
}

//operator overload for subtraction
stateVars stateVars::operator - (const stateVars& s2)const{
  stateVars s1 = *this;
  s1.rho -= s2.rho;
  s1.rhoU -= s2.rhoU;
  s1.rhoV -= s2.rhoV;
  s1.rhoW -= s2.rhoW;
  s1.rhoE -= s2.rhoE;
  return s1;
}

//operator overload for subtraction with a scalar
stateVars operator- (const double &scalar, const stateVars &s2){
  stateVars s1;
  s1.rho = s2.rho - scalar;
  s1.rhoU = s2.rhoU - scalar;
  s1.rhoV = s2.rhoV - scalar;
  s1.rhoW = s2.rhoW - scalar;
  s1.rhoE = s2.rhoE - scalar;
  return s1;
}

//operator overload for elementwise multiplication
stateVars stateVars::operator * (const stateVars& s2)const{
  stateVars s1 = *this;
  s1.rho *= s2.rho;
  s1.rhoU *= s2.rhoU;
  s1.rhoV *= s2.rhoV;
  s1.rhoW *= s2.rhoW;
  s1.rhoE *= s2.rhoE;
  return s1;
}

//operator overload for multiplication with a scalar
stateVars operator* (const double &scalar, const stateVars &s2){
  stateVars s1;
  s1.rho = s2.rho * scalar;
  s1.rhoU = s2.rhoU * scalar;
  s1.rhoV = s2.rhoV * scalar;
  s1.rhoW = s2.rhoW * scalar;
  s1.rhoE = s2.rhoE * scalar;
  return s1;
}


//operator overload for elementwise division
stateVars stateVars::operator / (const stateVars& s2)const{
  stateVars s1 = *this;
  s1.rho /= s2.rho;
  s1.rhoU /= s2.rhoU;
  s1.rhoV /= s2.rhoV;
  s1.rhoW /= s2.rhoW;
  s1.rhoE /= s2.rhoE;
  return s1;
}

//operator overload for division with a scalar
stateVars operator/ (const double &scalar, const stateVars &s2){
  stateVars s1;
  s1.rho = scalar / s2.rho ;
  s1.rhoU = scalar / s2.rhoU ;
  s1.rhoV = scalar / s2.rhoV ;
  s1.rhoW = scalar / s2.rhoW ;
  s1.rhoE = scalar / s2.rhoE ;
  return s1;
}


//member function to calculate reconstruction of state variables from cell center to cell face
//this function uses muscle extrapolation resulting in higher order accuracy
stateVars stateVars::FaceReconMUSCL( const stateVars &stateUW2, const stateVars &stateDW1, const string &side, const double &kappa, const string &lim, double upwind2face, double upwind, double central)const{

  //stateUW2 is the upwind cell furthest from the face at which the state is being reconstructed. If the state is being reconstructed from the left it would be at p-1
  //stateUW1 is the upwind cell nearest to the face at which the state is being reconstructed. If the state is being reconstructed from the left it would be at p
  //         it is the object being operated on
  //stateDW1 is the downwind cell. If the state is being reconstructed from the left it would be at p+1
  //side is a string the says whether the reconstruction is from the left or right
  //kappa is the parameter that determines which scheme is implemented

  stateVars faceState;
  stateVars r;
  stateVars limiter;
  stateVars invLimiter;
  double eps = 1.0e-20;

  stateVars stateUW1 = *this;

  stateVars num = (1.0 / central) * (eps + (stateDW1 - stateUW1));
  stateVars denom = (1.0 / upwind) * (eps + (stateUW1 - stateUW2));

  if (side == "left" || "right"){
    r =  num / denom;  //divided differences to base ratio on; eps must be listed to left of stateVars
  }
  // else if (side == "right"){
  //   r =  (eps + (stateUW1 - stateDW1)) / (eps + (stateUW2 - stateUW1)) ;  //divided differences to base ratio on; eps must be listed to left of stateVars
  // }
  else{
    cerr << "ERROR: Face reconstruction from side " << side << " is not recognized!" << endl;
  }


  if (lim == "none"){
    limiter = LimiterNone();    //limiter defined on inverse of r
    invLimiter = LimiterNone();
  }
  else if (lim == "vanAlbada"){
    limiter = LimiterVanAlbada(r);  //limiter defined on inverse of r
    invLimiter = LimiterVanAlbada(1.0/r);
  }
  else if (lim == "minmod"){
    limiter = LimiterMinmod(stateUW1 - stateUW2, stateDW1 - stateUW1, kappa);
  }
  else{
    cerr << "ERROR: Limiter " << lim << " is not recognized!" << endl;
  }


  // cout << "numerator " << num << endl;
  // cout << "denominator " << denom << endl;
  // cout << stateUW1 << endl;
  // cout << stateUW2 << endl;
  // cout << "r: " << r << endl;
  // cout << "limiter " << limiter << endl;
  // cout << "invLimiter " << invLimiter << endl;

  //CHANGED last term from (stateUW1 - stateUW2)
  if (side == "left" || "right"){
    //faceState = stateUW1 + 0.25 * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter) * (stateDW1 - stateUW1);
    //faceState = stateUW1 + 0.25 * ( (1.0 - kappa) * invLimiter * (stateUW1 - stateUW2) + (1.0 + kappa) * limiter * (stateDW1 - stateUW1) );
    //faceState = stateUW1 + 0.25 * (stateUW1 - stateUW2) * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter );
    //faceState = stateUW1 + 0.5 * (stateUW1 - stateUW2) * limiter ;
    faceState = stateUW1 + 0.25 * (stateUW1 - stateUW2) * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter );
  }
  // else if (side == "right"){
  //   //faceState = stateUW1 + 0.25 * ( (1.0 + kappa) * limiter + (1.0 - kappa) * r * invLimiter) * (stateDW1 - stateUW1);
  //   //faceState = stateUW1 - 0.25 * ( (1.0 + kappa) * invLimiter * (stateUW1 - stateUW2) + (1.0 - kappa) * limiter * (stateDW1 - stateUW1) );
  //   //faceState = stateUW1 - 0.25 * (stateUW2 - stateUW1) * ( (1.0 + kappa) * r * invLimiter + (1.0 - kappa) * limiter );
  //   //faceState = stateUW1 - 0.25 * (1.0 - kappa) * (stateUW2 - stateUW1) ;
  //   faceState = stateUW1 - 0.5 * (stateUW2 - stateUW1) * limiter ;
  // }
  else {
    cerr << "ERROR: Face reconstruction from side " << side << " is not recognized!" << endl;
  }
  //cout << "face state " << side << " = " << faceState << endl;

  // if ( ((stateDW1.Rho() <= faceState.Rho()) && (faceState.Rho() <= stateUW1.Rho())) || ((stateDW1.Rho() >= faceState.Rho()) && (faceState.Rho() >= stateUW1.Rho())) ){
  //   //do nothing
  // }
  // else{
  //   cout << "ERROR extrapolated state not inbetween upwind and downwind states!" << endl;
  //   cout << "side " << side << endl;
  //   cout << "downwind state " << stateDW1 ;
  //   cout << "upwind state " << stateUW1 ;
  //   cout << "upwind state two " << stateUW2 ;
  //   cout << "extrapolated state " << faceState ;
  //   cout << "limiter " << limiter ;
  //   cout << "inverse limiter " << invLimiter ;
  //   cout << "extrapolation factor left " << 0.25 * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter ) ;
  //   cout << "extrapolation factor right " << 0.25 * ( (1.0 + kappa) * r * invLimiter + (1.0 - kappa) * limiter ) ;
  //   cout << "state diff left" << stateUW1 - stateUW2 ;
  //   cout << "state diff right" << stateUW2 - stateUW1 ;
  //   cout << "r " << r << endl ;
  // }

  return faceState;

}

//member function to calculate minmod limiter
stateVars stateVars::LimiterMinmod( const stateVars &upwind, const stateVars &downwind, const double kap)const{
  stateVars limiter, sign;

  double beta = (3.0 - kap) / (1.0 - kap);

  if (upwind.Rho() > 0.0 ){
    sign.SetRho(1.0);
  }
  else if (upwind.Rho() < 0.0 ){
    sign.SetRho(-1.);
  }
  else{
    sign.SetRho(0.0);
  }

  if (upwind.RhoU() > 0.0 ){
    sign.SetRhoU(1.0);
  }
  else if (upwind.RhoU() < 0.0 ){
    sign.SetRhoU(-1.);
  }
  else{
    sign.SetRhoU(0.0);
  }

  if (upwind.RhoV() > 0.0 ){
    sign.SetRhoV(1.0);
  }
  else if (upwind.RhoV() < 0.0 ){
    sign.SetRhoV(-1.);
  }
  else{
    sign.SetRhoV(0.0);
  }

  if (upwind.RhoW() > 0.0 ){
    sign.SetRhoW(1.0);
  }
  else if (upwind.RhoW() < 0.0 ){
    sign.SetRhoW(-1.);
  }
  else{
    sign.SetRhoW(0.0);
  }

  if (upwind.RhoE() > 0.0 ){
    sign.SetRhoE(1.0);
  }
  else if (upwind.RhoE() < 0.0 ){
    sign.SetRhoE(-1.);
  }
  else{
    sign.SetRhoE(0.0);
  }


  limiter.SetRho( sign.Rho() * max(0.0, min( fabs(upwind.Rho()), sign.Rho()*downwind.Rho()*beta ) ) );
  limiter.SetRhoU( sign.RhoU() * max(0.0, min( fabs(upwind.RhoU()), sign.RhoU()*downwind.RhoU()*beta ) ) );
  limiter.SetRhoV( sign.RhoV() * max(0.0, min( fabs(upwind.RhoV()), sign.RhoV()*downwind.RhoV()*beta ) ) );
  limiter.SetRhoW( sign.RhoW() * max(0.0, min( fabs(upwind.RhoW()), sign.RhoW()*downwind.RhoW()*beta ) ) );
  limiter.SetRhoE( sign.RhoE() * max(0.0, min( fabs(upwind.RhoE()), sign.RhoE()*downwind.RhoE()*beta ) ) );

  return limiter;
}

//member function to calculate Van Albada limiter
stateVars stateVars::LimiterVanAlbada( const stateVars &r)const{
  stateVars limiter;
  limiter = (r + r*r)/(1 + r*r);
  //if value is negative, return zero
  limiter.SetRho( max(0.0, limiter.Rho()) );
  limiter.SetRhoU( max(0.0, limiter.RhoU()) );
  limiter.SetRhoV( max(0.0, limiter.RhoV()) );
  limiter.SetRhoW( max(0.0, limiter.RhoW()) );
  limiter.SetRhoE( max(0.0, limiter.RhoE()) );
  return limiter;
}

//member function to return no limiter
stateVars stateVars::LimiterNone()const{
  stateVars limiter;
  limiter.SetRho(1.0);
  limiter.SetRhoU(1.0);
  limiter.SetRhoV(1.0);
  limiter.SetRhoW(1.0);
  limiter.SetRhoE(1.0);
  return limiter;
}

//member function to return the state of the appropriate ghost cell
stateVars stateVars::GetGhostState( const string &bcType, const vector3d<double> &areaVec, const string &surf, const input &inputVars, const idealGas &eqnState )const{

  stateVars ghostState = *this;

  vector3d<double> normArea;

  if (surf == "il" || surf == "jl" || surf == "kl"){
    normArea = areaVec / areaVec.Mag();
  }
  else if (surf == "iu" || surf == "ju" || surf == "ku"){
    normArea = -1.0 * (areaVec / areaVec.Mag()); //at upper surface normal should point into domain
  }


  double rho = 0;
  double normVelCellCenter = 0;

  if (bcType == "slipWall"){             //for slip wall state should be reflected across boundary face
    rho = (*this).Rho();
    vector3d<double> stateVel = (*this).Velocity();
    normVelCellCenter = stateVel.DotProd(normArea);

    vector3d<double> ghostVel (stateVel.X() - 2.0 * normArea.X() * normVelCellCenter,
			       stateVel.Y() - 2.0 * normArea.Y() * normVelCellCenter,
			       stateVel.Z() - 2.0 * normArea.Z() * normVelCellCenter);

    ghostState.SetRho(rho);
    ghostState.SetRhoU(rho*ghostVel.X());
    ghostState.SetRhoV(rho*ghostVel.Y());
    ghostState.SetRhoW(rho*ghostVel.Z());

    double pressure = (*this).Pressure(eqnState);
    double specEn = eqnState.GetSpecEnergy(pressure,rho);
    double energy = eqnState.GetEnergy(specEn, ghostVel.Mag());

    ghostState.SetRhoE(rho*energy);

    // if (surf == "iu"){
    //   cout << "norm vel cc " << normVelCellCenter  << endl;
    //   cout << "boundary vel " << stateVel << ", " << stateVel.Mag() << endl;
    //   cout << "ghost vel " << ghostVel << ", " << ghostVel.Mag() << endl;
    //   cout << "norm area " << normArea << endl;
    // }
  }
  else if (bcType == "viscousWall"){             //for viscous wall velocity at face should be 0.0
    rho = (*this).Rho();
    vector3d<double> stateVel = (*this).Velocity();
    vector3d<double> ghostVel (-1.0 * stateVel.X(),
			       -1.0 * stateVel.Y(),
			       -1.0 * stateVel.Z() );

    ghostState.SetRho(rho);
    ghostState.SetRhoU(rho*ghostVel.X());
    ghostState.SetRhoV(rho*ghostVel.Y());
    ghostState.SetRhoW(rho*ghostVel.Z());

    double pressure = (*this).Pressure(eqnState);
    double specEn = eqnState.GetSpecEnergy(pressure,rho);
    double energy = eqnState.GetEnergy(specEn, ghostVel.Mag());

    ghostState.SetRhoE(rho*energy);

  }
  else if (bcType == "subsonicInflow"){     //set velocity and density to freestream values
    rho = 1.0;
    double sos = eqnState.GetSoS( inputVars.PRef(), inputVars.RRef() );
    vector3d<double> ghostVel = inputVars.VelRef() / sos;

    ghostState.SetRho(rho);
    ghostState.SetRhoU(rho*ghostVel.X());
    ghostState.SetRhoV(rho*ghostVel.Y());
    ghostState.SetRhoW(rho*ghostVel.Z());

    double pressure = (*this).Pressure(eqnState);           //numerical bc for pressure

    double specEn = eqnState.GetSpecEnergy(pressure, rho);
    double energy = eqnState.GetEnergy(specEn, ghostVel.Mag());

    // cout << "boundary state " << (*this) << endl;
    // cout << "pressure, specEn, energy " << pressure << ", " << specEn << ", " << energy << endl;

    ghostState.SetRhoE(rho*energy);

  }
  else if (bcType == "subsonicOutflow"){     //set pressure to freestream value
    double pressure = 1.0 / eqnState.Gamma();

    //numerical bcs for density, velocity
    double sos = eqnState.GetSoS( (*this).Pressure(eqnState), (*this).Rho() );
    //rho = (*this).Rho();
    rho = eqnState.Gamma() * pressure / (sos * sos);

    double specEnergy = eqnState.GetSpecEnergy(pressure,rho);
    vector3d<double> stateVel = (*this).Velocity();

    // double normVel = stateVel.DotProd(normArea);
    // if (normVel > 0.0){ //mass flow into domain at outflow!
    //   cout << "mass flow into domain at outflow!" << endl;
    // }

    double energy = eqnState.GetEnergy(specEnergy, stateVel.Mag());

    ghostState.SetRhoE(rho*energy);

  }
  else if (bcType == "characteristic"){
      //calculate ui, u0, ci, co, mi
      vector3d<double> freeVel = inputVars.VelRef();
      double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
      freeVel = freeVel / freeSoS;
      freeSoS = freeSoS / freeSoS;

      double freeSpecEnergy = eqnState.GetSpecEnergy(1.0/eqnState.Gamma(), 1.0);
      double freeEnergy = eqnState.GetEnergy(freeSpecEnergy, freeVel.Mag());

      stateVars freeState(1.0, freeEnergy, freeVel);
      double rhoFree = freeState.Rho();
      double velFreeNorm = freeVel.DotProd(normArea); 
      double SoSFree = 1.0;

      vector3d<double> freeVelTan;                                  //freestream tangent velocity
      freeVelTan.SetX(freeVel.X() - normArea.X()*velFreeNorm);
      freeVelTan.SetY(freeVel.Y() - normArea.Y()*velFreeNorm);
      freeVelTan.SetZ(freeVel.Z() - normArea.Z()*velFreeNorm);

      vector3d<double> velInt;
      velInt.SetX((*this).Vx());
      velInt.SetY((*this).Vy());
      velInt.SetZ((*this).Vz());
      double velIntNorm;

      velIntNorm = velInt.DotProd(normArea); 
      
      double rhoInt = (*this).Rho();
      double SoSInt = eqnState.GetSoS((*this).Pressure(eqnState), rhoInt);

      vector3d<double> velIntTan;                                  //internal tangent velocity
      velIntTan.SetX(velInt.X() - normArea.X()*velIntNorm);
      velIntTan.SetY(velInt.Y() - normArea.Y()*velIntNorm);
      velIntTan.SetZ(velInt.Z() - normArea.Z()*velIntNorm);

      double machInt = fabs(velIntNorm)/SoSInt;

      double riemannInvarPlus, riemannInvarMinus;

      if ( machInt >= 1.0 && velIntNorm > 0.0 ){ //supersonic inflow

	//characteristics all go into the domain, so use freestream values for both riemann invariants
	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
	riemannInvarMinus = velFreeNorm - (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
	ghostState = freeState;
      }
      else if ( machInt >= 1.0 && velIntNorm < 0.0 ){ //supersonic outflow

	//characteristics all leave the domain, so use interior values for both riemann invariants
	riemannInvarPlus = velIntNorm + (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
	ghostState = (*this);
      }
      else if ( machInt < 1.0 && velIntNorm > 0.0 ){ //subsonic inflow

	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
	riemannInvarPlus = velFreeNorm + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus); //this is the normal velocity at the boundary
	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);

	double entropyBound = (SoSFree * SoSFree)/(eqnState.Gamma() * pow(rhoFree , eqnState.Gamma()-1.0));
	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
	vector3d<double> velBound;

	velBound.SetX( (velB * normArea.X() + freeVelTan.X()) );  
	velBound.SetY( (velB * normArea.Y() + freeVelTan.Y()) );
	velBound.SetZ( (velB * normArea.Z() + freeVelTan.Z()) );

	// cout << "inflow boundary norm vel " << velB*normArea << endl;
	// cout << "inflow norm vel " << freeVel << endl;
	// cout << "inflow tangent vel " << freeVelTan << endl;

	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

	ghostState.SetRho(rhoBound);
	ghostState.SetRhoU(rhoBound * velBound.X());
	ghostState.SetRhoV(rhoBound * velBound.Y());
	ghostState.SetRhoW(rhoBound * velBound.Z());
	ghostState.SetRhoE(rhoBound * energyBound);

      }
      else if ( machInt < 1.0 && velIntNorm < 0.0 ){ //subsonic outflow

	//characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
	riemannInvarMinus = velIntNorm - (2.0 * SoSInt)/(eqnState.Gamma() - 1.0);
	riemannInvarPlus = -1.0 * fabs(velFreeNorm) + (2.0 * SoSFree)/(eqnState.Gamma() - 1.0);
	double velB = 0.5 * (riemannInvarPlus + riemannInvarMinus);
	double SoSBound = 0.25 * (eqnState.Gamma() - 1.0) * (riemannInvarPlus - riemannInvarMinus);

	double entropyBound = (SoSInt * SoSInt)/(eqnState.Gamma() * pow(rhoInt , eqnState.Gamma()-1.0));
	double rhoBound = pow((SoSBound * SoSBound)/(eqnState.Gamma() * entropyBound) , 1.0/(eqnState.Gamma()-1.0));
	double pressBound = rhoBound * SoSBound * SoSBound / eqnState.Gamma();
	vector3d<double> velBound;

	velBound.SetX( (velB * normArea.X() + velIntTan.X()) );
	velBound.SetY( (velB * normArea.Y() + velIntTan.Y()) );
	velBound.SetZ( (velB * normArea.Z() + velIntTan.Z()) );

	double specEnBound = eqnState.GetSpecEnergy(pressBound, rhoBound);
	double energyBound = eqnState.GetEnergy(specEnBound, velBound.Mag());

	ghostState.SetRho(rhoBound);
	ghostState.SetRhoU(rhoBound * velBound.X());
	ghostState.SetRhoV(rhoBound * velBound.Y());
	ghostState.SetRhoW(rhoBound * velBound.Z());
	ghostState.SetRhoE(rhoBound * energyBound);

      }
      else {
	cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
	exit(0);
      }

  }
  else if (bcType == "supersonicInflow" || bcType == "supersonicOutflow"){
    //do nothing and return state of boundary adjacent cell
  }
  else {
    cerr << "ghost state for BC type " << bcType << " is not supported!" << endl;
    exit(0);
  }

  // cout << "surface type " << surf << " bc type " << bcType << " ghost state: " << ghostState << endl;

  return ghostState;

}

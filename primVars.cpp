#include "primVars.h"
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

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const primVars &prim){

  os << prim.rho << "   " << prim.u << "   " << prim.v << "   " << prim.w << "   " << prim.p << endl;

  return os;
}

//operator overload for addition
primVars primVars::operator + (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.rho += prim2.rho;
  prim1.u += prim2.u;
  prim1.v += prim2.v;
  prim1.w += prim2.w;
  prim1.p += prim2.p;
  return prim1;
}

//operator overload for addition with a scalar
primVars operator+ (const double &scalar, const primVars &prim2){
  primVars prim1( prim2.rho + scalar,
		  prim2.u + scalar,
		  prim2.v + scalar,
		  prim2.w + scalar,
		  prim2.p + scalar);
  return prim1;
}

//operator overload for subtraction
primVars primVars::operator - (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.rho -= prim2.rho;
  prim1.u -= prim2.u;
  prim1.v -= prim2.v;
  prim1.w -= prim2.w;
  prim1.p -= prim2.p;
  return prim1;
}

//operator overload for subtraction with a scalar
primVars operator- (const double &scalar, const primVars &prim2){
  primVars prim1( scalar - prim2.rho,
		  scalar - prim2.u,
		  scalar - prim2.v,
		  scalar - prim2.w,
		  scalar - prim2.p );
  return prim1;
}

//operator overload for elementwise multiplication
primVars primVars::operator * (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.rho *= prim2.rho;
  prim1.u *= prim2.u;
  prim1.v *= prim2.v;
  prim1.w *= prim2.w;
  prim1.p *= prim2.p;
  return prim1;
}

//member function for scalar multiplication
primVars  primVars::operator * (const double &scalar){
  primVars temp = *this;
  temp.rho *= scalar;
  temp.u *= scalar;
  temp.v *= scalar;
  temp.w *= scalar;
  temp.p *= scalar;
  return temp;
}

//member function for scalar addition
primVars  primVars::operator + (const double &scalar){
  primVars temp = *this;
  temp.rho += scalar;
  temp.u += scalar;
  temp.v += scalar;
  temp.w += scalar;
  temp.p += scalar;
  return temp;
}

//member function for scalar subtraction
primVars  primVars::operator - (const double &scalar){
  primVars temp = *this;
  temp.rho -= scalar;
  temp.u -= scalar;
  temp.v -= scalar;
  temp.w -= scalar;
  temp.p -= scalar;
  return temp;
}

//member function for scalar division
primVars  primVars::operator / (const double &scalar){
  primVars temp = *this;
  temp.rho /= scalar;
  temp.u /= scalar;
  temp.v /= scalar;
  temp.w /= scalar;
  temp.p /= scalar;
  return temp;
}

//operator overload for multiplication with a scalar
primVars operator* (const double &scalar, const primVars &prim2){
  primVars prim1( prim2.rho * scalar,
		  prim2.u * scalar,
		  prim2.v * scalar,
		  prim2.w * scalar,
		  prim2.p * scalar);
  return prim1;
}


//operator overload for elementwise division
primVars primVars::operator / (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.rho /= prim2.rho;
  prim1.u /= prim2.u;
  prim1.v /= prim2.v;
  prim1.w /= prim2.w;
  prim1.p /= prim2.p;
  return prim1;
}

//operator overload for division with a scalar
primVars operator/ (const double &scalar, const primVars &prim2){
  primVars prim1( scalar / prim2.rho,
		  scalar / prim2.u,
		  scalar / prim2.v,
		  scalar / prim2.w,
		  scalar / prim2.p);
  return prim1;
}


//member function to calculate reconstruction of primative variables from cell center to cell face
//this function uses muscle extrapolation resulting in higher order accuracy
primVars primVars::FaceReconMUSCL( const primVars &primUW2, const primVars &primDW1, const string &side, const double &kappa, const string &lim, double uw, double uw2, double dw)const{

  //primUW2 is the upwind cell furthest from the face at which the primative is being reconstructed. If the primative is being reconstructed from the left it would be at p-1
  //primUW1 is the upwind cell nearest to the face at which the primative is being reconstructed. If the primative is being reconstructed from the left it would be at p
  //         it is the object being operated on
  //primDW1 is the downwind cell. If the primative is being reconstructed from the left it would be at p+1
  //side is a string the says whether the reconstruction is from the left or right
  //kappa is the parameter that determines which scheme is implemented
  //uw is length of upwind cell
  //uw2 is length of furthest upwind cell
  //dw is length of downwind cell

  primVars facePrim;
  primVars r;
  primVars limiter;
  primVars invLimiter;
  double eps = 1.0e-25;

  primVars primUW1 = *this;

  double dPlus = (uw + dw) / (2.0 * uw);
  double dMinus = (uw + uw2) / (2.0 * uw);

  // if (dPlus < 0.6 || dPlus > 1.4){
  //   cout << "dPlus: " << dPlus << endl;
  //   cout << "uw, uw2, dw: " << uw << ", " << uw2 << ", " << dw << endl;
  // }
  // if (dMinus < 0.6 || dMinus > 1.4){
  //   cout << "dMinus: " << dMinus << endl;
  //   cout << "uw, uw2, dw: " << uw << ", " << uw2 << ", " << dw << endl;
  // }


  if (side == "left" || true){
    r =  (eps + (primDW1 - primUW1) / dPlus) / (eps + (primUW1 - primUW2) / dMinus);  //divided differences to base ratio on; eps must be listed to left of primVars
  }
  // else if (side == "right"){
  //   r = (eps + (primUW1 - primDW1)) / (eps + (primUW2 - primUW1)) ;
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
    limiter = LimiterMinmod(primUW1 - primUW2, primDW1 - primUW1, kappa);
  }
  else{
    cerr << "ERROR: Limiter " << lim << " is not recognized!" << endl;
  }

  //CHANGED last term from (primUW1 - primUW2)
  if (side == "left" || true){
    //facePrim = primUW1 + 0.25 * (primUW1 - primUW2) * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter );
    facePrim = primUW1 + ( limiter / (2.0 * (dPlus + dMinus) ) ) * ( (dPlus  - kappa * limiter) * (primUW1 - primUW2) / dMinus + 
								     (dMinus + kappa * limiter) * (primDW1 - primUW1) / dPlus );
  }
  // else if (side == "right"){
  //   facePrim = primUW1 - 0.25 * (primUW2 - primUW1) * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter );
  // }
  else {
    cerr << "ERROR: Face reconstruction from side " << side << " is not recognized!" << endl;
  }

  return facePrim;

}

//member function to calculate minmod limiter
primVars primVars::LimiterMinmod( const primVars &upwind, const primVars &downwind, const double kap)const{
  primVars limiter, sign;

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

  if (upwind.U() > 0.0 ){
    sign.SetU(1.0);
  }
  else if (upwind.U() < 0.0 ){
    sign.SetU(-1.);
  }
  else{
    sign.SetU(0.0);
  }

  if (upwind.V() > 0.0 ){
    sign.SetV(1.0);
  }
  else if (upwind.V() < 0.0 ){
    sign.SetV(-1.);
  }
  else{
    sign.SetV(0.0);
  }

  if (upwind.W() > 0.0 ){
    sign.SetW(1.0);
  }
  else if (upwind.W() < 0.0 ){
    sign.SetW(-1.);
  }
  else{
    sign.SetW(0.0);
  }

  if (upwind.P() > 0.0 ){
    sign.SetP(1.0);
  }
  else if (upwind.P() < 0.0 ){
    sign.SetP(-1.);
  }
  else{
    sign.SetP(0.0);
  }


  limiter.SetRho( sign.Rho() * max(0.0, min( fabs(upwind.Rho()), sign.Rho()*downwind.Rho()*beta ) ) );
  limiter.SetU( sign.U() * max(0.0, min( fabs(upwind.U()), sign.U()*downwind.U()*beta ) ) );
  limiter.SetV( sign.V() * max(0.0, min( fabs(upwind.V()), sign.V()*downwind.V()*beta ) ) );
  limiter.SetW( sign.W() * max(0.0, min( fabs(upwind.W()), sign.W()*downwind.W()*beta ) ) );
  limiter.SetP( sign.P() * max(0.0, min( fabs(upwind.P()), sign.P()*downwind.P()*beta ) ) );

  return limiter;
}

//member function to calculate Van Albada limiter
primVars primVars::LimiterVanAlbada( const primVars &r)const{
  primVars limiter;
  limiter = (r + r*r)/(1 + r*r);
  //if value is negative, return zero
  limiter.SetRho( max(0.0, limiter.Rho()) );
  limiter.SetU( max(0.0, limiter.U()) );
  limiter.SetV( max(0.0, limiter.V()) );
  limiter.SetW( max(0.0, limiter.W()) );
  limiter.SetP( max(0.0, limiter.P()) );
  return limiter;
}

//member function to return no limiter
primVars primVars::LimiterNone()const{
  primVars limiter;
  limiter.SetRho(1.0);
  limiter.SetU(1.0);
  limiter.SetV(1.0);
  limiter.SetW(1.0);
  limiter.SetP(1.0);
  return limiter;
}

//member function to return the state of the appropriate ghost cell
primVars primVars::GetGhostState( const string &bcType, const vector3d<double> &areaVec, const string &surf, const input &inputVars, const idealGas &eqnState, const int layer )const{

  primVars ghostState = (*this);  //set ghost state equal to boundary state to start

  //check to see that ghost layer corresponds to allowable number
  if ( !(layer == 1 || layer == 2) ){
    cerr << "ERROR: Error in primVars::GetGhostState. Requesting ghost state at a ghost layer " << layer << ". Please choose either 1 or 2" << endl;
    exit(0);
  }

  vector3d<double> normArea;

  // if (surf == "il" || surf == "jl" || surf == "kl"){
  //   normArea = -1.0 * areaVec / areaVec.Mag(); //at lower surface normal should point out of domain for ghost cell calculation
  // }
  // else if (surf == "iu" || surf == "ju" || surf == "ku"){
    normArea = areaVec / areaVec.Mag(); 
  // }

  double normVelCellCenter = 0;

  if (bcType == "slipWall"){             //for slip wall state should be reflected across boundary face, density and pressure stay equal to the boundary cell
    vector3d<double> stateVel = (*this).Velocity();
    normVelCellCenter = stateVel.DotProd(normArea);

    //for a slip wall the velocity of the boundary cell center is reflected across the boundary face to get the velocity at the ghost cell center
    vector3d<double> ghostVel (stateVel.X() - 2.0 * normArea.X() * normVelCellCenter,
			       stateVel.Y() - 2.0 * normArea.Y() * normVelCellCenter,
			       stateVel.Z() - 2.0 * normArea.Z() * normVelCellCenter);

    ghostState.SetU(ghostVel.X());
    ghostState.SetV(ghostVel.Y());
    ghostState.SetW(ghostVel.Z());

    //numerical BCs for rho and pressure, same as boundary state

  }
  else if (bcType == "viscousWall"){             //for viscous wall velocity at face should be 0.0, density and pressure stay equal to the boundary cell
    vector3d<double> stateVel = (*this).Velocity();

    //ghost cell velocity at cell center is set to opposite of velocity at boundary cell center so that velocity at face will be zero
    vector3d<double> ghostVel (-1.0 * stateVel.X(),
			       -1.0 * stateVel.Y(),
			       -1.0 * stateVel.Z() );

    ghostState.SetU(ghostVel.X());
    ghostState.SetV(ghostVel.Y());
    ghostState.SetW(ghostVel.Z());

    //numerical BCs for rho and pressure, same as boundary state

  }
  else if (bcType == "subsonicInflow"){     //set velocity and density to freestream values
    double sos = eqnState.GetSoS( inputVars.PRef(), inputVars.RRef() );
    vector3d<double> ghostVel = inputVars.VelRef() / sos;   //nondimensionalize velocity

    ghostState.SetRho(1.0);
    ghostState.SetU(ghostVel.X());
    ghostState.SetV(ghostVel.Y());
    ghostState.SetW(ghostVel.Z());

    //numerical bc for pressure, same as boundary state

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else if (bcType == "subsonicOutflow"){     //set pressure to freestream value
    ghostState.SetP( 1.0 / eqnState.Gamma() );

    //numerical bcs for density, velocity -- equal to boundary cell

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else if (bcType == "characteristic"){
    //freestream variables
    double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
    vector3d<double> freeVel = inputVars.VelRef() / freeSoS;
    primVars freeState(1.0, 1.0/eqnState.Gamma(), freeVel);
    //double velFreeNorm = freeVel.DotProd(normArea); 

    //vector3d<double> freeVelTan = freeVel - normArea * velFreeNorm;            //freestream tangent velocity

    //internal variables
    double velIntNorm = (*this).Velocity().DotProd(normArea); 
    double SoSInt = eqnState.GetSoS((*this).P(), (*this).Rho());
    //vector3d<double> velIntTan = (*this).Velocity() - normArea * velIntNorm;    //internal tangent velocity
    double machInt = fabs(velIntNorm)/SoSInt;

    if ( machInt >= 1.0 && velIntNorm < 0.0 ){ //supersonic inflow
      //characteristics all go into the domain, so use freestream values for both riemann invariants
      ghostState = freeState;

    }
    else if ( machInt >= 1.0 && velIntNorm >= 0.0 ){ //supersonic outflow
      //characteristics all leave the domain, so use interior values for both riemann invariants
      ghostState = (*this);
    }
    else if ( machInt < 1.0 && velIntNorm < 0.0 ){ //subsonic inflow

      //characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
      double rhoSoSInt = (*this).Rho() * SoSInt;
      vector3d<double> velDiff = freeState.Velocity() - (*this).Velocity();
      ghostState.SetP(0.5 * (freeState.P() + (*this).P() - rhoSoSInt * normArea.DotProd(velDiff)) );
      ghostState.SetRho(freeState.Rho() + (ghostState.P() - freeState.P()) / (SoSInt * SoSInt) );
      ghostState.SetU(freeState.U() - normArea.X() * (freeState.P() - ghostState.P()) / rhoSoSInt);
      ghostState.SetV(freeState.V() - normArea.Y() * (freeState.P() - ghostState.P()) / rhoSoSInt);
      ghostState.SetW(freeState.W() - normArea.Z() * (freeState.P() - ghostState.P()) / rhoSoSInt);

    }
    else if ( machInt < 1.0 && velIntNorm >= 0.0 ){ //subsonic outflow

      //characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
      double rhoSoSInt = (*this).Rho() * SoSInt;
      ghostState.SetP(freeState.P());
      ghostState.SetRho( (*this).Rho() + ( ghostState.P() - (*this).P() ) / (SoSInt * SoSInt)   );
      ghostState.SetU(   (*this).U() + normArea.X() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);
      ghostState.SetV(   (*this).V() + normArea.Y() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);
      ghostState.SetW(   (*this).W() + normArea.Z() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);

    }
    else {
      cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
      exit(0);
    }

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else if (bcType == "supersonicInflow"){
    //physical boundary conditions - fix everything
    double sos = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef());
    vector3d<double> vel = inputVars.VelRef() / sos;                       //nondimensional velocity

    ghostState.SetRho(1.0);
    ghostState.SetU(vel.X());
    ghostState.SetV(vel.Y());
    ghostState.SetW(vel.Z());
    ghostState.SetP(1.0 / eqnState.Gamma());

    //no matter the layer, return the same ghost state
    // if (layer == 2){ //extrapolate to get ghost state at 2nd layer
    //   ghostState = 2.0 * ghostState - (*this);
    // }

  }
  else if (bcType == "supersonicOutflow"){
    //do nothing and return boundary state -- numerical BCs for all
    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else if (bcType == "stagnationInlet"){
    double g = eqnState.Gamma() - 1.0;
    //calculate outgoing riemann invarient
    double rNeg = (*this).Velocity().DotProd(normArea) - 2.0 * (*this).SoS(eqnState) / g;

    //calculate SoS on boundary
    double cosTheta = -1.0 * (*this).Velocity().DotProd(normArea) / (*this).Velocity().Mag();
    double stagSoSsq = pow((*this).SoS(eqnState),2.0) + 0.5 * g * (*this).Velocity().MagSq();

    double sosB = -1.0 * rNeg * g / (g * cosTheta * cosTheta + 2.0) * (1.0 + cosTheta * sqrt( (g * cosTheta * cosTheta + 2.0) * stagSoSsq / 
											      (g * rNeg * rNeg) - 0.5 * g  ) );
    double tb = inputVars.StagInletT0() / inputVars.TRef() * (sosB * sosB / stagSoSsq);
    double aRef = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef());
    double pb = inputVars.StagInletP0() / (inputVars.RRef() * aRef * aRef) * pow(sosB * sosB / stagSoSsq, eqnState.Gamma()/g);
    double vbMag = sqrt(2.0 / g * (inputVars.StagInletT0() / inputVars.TRef() - tb));

    ghostState.SetRho(eqnState.GetDensityTP(tb,pb));
    ghostState.SetU(vbMag * inputVars.StagInletDx());
    ghostState.SetV(vbMag * inputVars.StagInletDy());
    ghostState.SetW(vbMag * inputVars.StagInletDz());
    ghostState.SetP(pb);

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else if (bcType == "pressureOutlet"){
    double aRef = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef());
    double pb = inputVars.PressureOutletP() / (inputVars.RRef() * aRef * aRef);

    double SoSInt = (*this).SoS(eqnState);
    double rhoSoSInt = (*this).Rho() * SoSInt;
    ghostState.SetP(pb);
    ghostState.SetRho( (*this).Rho() + ( ghostState.P() - (*this).P() ) / (SoSInt * SoSInt)   );
    ghostState.SetU(   (*this).U() + normArea.X() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);
    ghostState.SetV(   (*this).V() + normArea.Y() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);
    ghostState.SetW(   (*this).W() + normArea.Z() * ( (*this).P() - ghostState.P() ) / rhoSoSInt);

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  else {
    cerr << "ghost state for BC type " << bcType << " is not supported!" << endl;
    cerr << "surface is " << surf << endl;
    exit(0);
  }

  return ghostState;

}

//member function to take in a colMatrix of updates to the conservative variables, and update the primative variables with it.
primVars primVars::UpdateWithConsVars(const idealGas &eqnState, const colMatrix &du)const{

  //check to see that update is of proper size
  if ( du.Size() != 5 ){
    cerr << "ERROR: Error in primVars::UpdateWithConsVars. The colMatrix containing the update is not of proper size!" << endl;
    exit(0);
  }

  //convert primative to conservative and update
  colMatrix consUpdate = (*this).ConsVars(eqnState) + du;

  //convert back to primative variables
  primVars primUpdate;
  primUpdate.SetRho( consUpdate.Data(0) );
  primUpdate.SetU( consUpdate.Data(1) / consUpdate.Data(0) );
  primUpdate.SetV( consUpdate.Data(2) / consUpdate.Data(0) );
  primUpdate.SetW( consUpdate.Data(3) / consUpdate.Data(0) );
  double energy = consUpdate.Data(4) / consUpdate.Data(0);
  primUpdate.SetP( eqnState.GetPressFromEnergy( primUpdate.Rho(), energy, primUpdate.Velocity().Mag() ) );

  return primUpdate;
}

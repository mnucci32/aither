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
primVars primVars::FaceReconMUSCL( const primVars &primUW2, const primVars &primDW1, const string &side, const double &kappa, const string &lim, double upwind2face, double upwind, double central)const{

  //primUW2 is the upwind cell furthest from the face at which the primative is being reconstructed. If the primative is being reconstructed from the left it would be at p-1
  //primUW1 is the upwind cell nearest to the face at which the primative is being reconstructed. If the primative is being reconstructed from the left it would be at p
  //         it is the object being operated on
  //primDW1 is the downwind cell. If the primative is being reconstructed from the left it would be at p+1
  //side is a string the says whether the reconstruction is from the left or right
  //kappa is the parameter that determines which scheme is implemented

  primVars facePrim;
  primVars r;
  primVars limiter;
  primVars invLimiter;
  double eps = 1.0e-25;

  primVars primUW1 = *this;

  primVars num = (1.0 / central) * (eps + (primDW1 - primUW1));
  primVars denom = (1.0 / upwind) * (eps + (primUW1 - primUW2));

  if (side == "left" || "right"){
    r =  num / denom;  //divided differences to base ratio on; eps must be listed to left of primVars
  }
  // else if (side == "right"){
  //   r =  (eps + (primUW1 - primDW1)) / (eps + (primUW2 - primUW1)) ;  //divided differences to base ratio on; eps must be listed to left of primVars
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
  if (side == "left" || "right"){
    facePrim = primUW1 + 0.25 * (primUW1 - primUW2) * ( (1.0 - kappa) * limiter + (1.0 + kappa) * r * invLimiter );
  }
  // else if (side == "right"){
  //   faceState = primUW1 - 0.5 * (primUW2 - primUW1) * limiter ;
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
primVars primVars::GetGhostState( const string &bcType, const vector3d<double> &areaVec, const string &surf, const input &inputVars, const idealGas &eqnState )const{

  primVars ghostState = *this;  //set ghost state equal to boundary state to start

  vector3d<double> normArea;

  if (surf == "il" || surf == "jl" || surf == "kl"){
    normArea = areaVec / areaVec.Mag();
  }
  else if (surf == "iu" || surf == "ju" || surf == "ku"){
    normArea = -1.0 * (areaVec / areaVec.Mag()); //at upper surface normal should point into domain for ghost cell calculation
  }

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

  }
  else if (bcType == "subsonicOutflow"){     //set pressure to freestream value
    ghostState.SetP( 1.0 / eqnState.Gamma() );

    //numerical bcs for density, velocity -- equal to boundary cell

  }
  else if (bcType == "characteristic"){
      //calculate ui, u0, ci, co, mi
      vector3d<double> freeVel = inputVars.VelRef();
      double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
      freeVel = freeVel / freeSoS;
      freeSoS = freeSoS / freeSoS;

      primVars freeState(1.0, 1.0/eqnState.Gamma(), freeVel);
      double rhoFree = freeState.Rho();
      double velFreeNorm = freeVel.DotProd(normArea); 
      double SoSFree = 1.0;

      vector3d<double> freeVelTan;                                  //freestream tangent velocity
      freeVelTan.SetX(freeVel.X() - normArea.X()*velFreeNorm);
      freeVelTan.SetY(freeVel.Y() - normArea.Y()*velFreeNorm);
      freeVelTan.SetZ(freeVel.Z() - normArea.Z()*velFreeNorm);

      vector3d<double> velInt = (*this).Velocity();
      double velIntNorm;

      velIntNorm = velInt.DotProd(normArea); 
      
      double rhoInt = (*this).Rho();
      double SoSInt = eqnState.GetSoS((*this).P(), rhoInt);

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

	ghostState.SetRho(rhoBound);
	ghostState.SetU( (velB * normArea.X() + freeVelTan.X()) );
	ghostState.SetV( (velB * normArea.Y() + freeVelTan.Y()) );
	ghostState.SetW( (velB * normArea.Z() + freeVelTan.Z()) );
	ghostState.SetP(pressBound);

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

	ghostState.SetRho(rhoBound);
	ghostState.SetU( (velB * normArea.X() + velIntTan.X()) );
	ghostState.SetV( (velB * normArea.Y() + velIntTan.Y()) );
	ghostState.SetW( (velB * normArea.Z() + velIntTan.Z()) );
	ghostState.SetP(pressBound);

      }
      else {
	cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
	exit(0);
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

  }
  else if (bcType == "supersonicOutflow"){
    //do nothing and return boundary state -- numerical BCs for all
  }
  else {
    cerr << "ghost state for BC type " << bcType << " is not supported!" << endl;
    cerr << "surface is " << surf << endl;
    exit(0);
  }

  return ghostState;

}

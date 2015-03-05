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

//constructor
primVars::primVars( const genArray &a, const bool &prim, const idealGas &eqnState ){

  if (prim){ //genArray is primative variables
    for ( int ii = 0; ii < NUMVARS; ii++ ){
      data[ii] = a[ii];
    }
  }
  else{ //genArray is conserved variables
    data[0] = a[0];
    data[1] = a[1]/a[0];
    data[2] = a[2]/a[0];
    data[3] = a[3]/a[0];
    double energy = a[4]/a[0];
    data[4] = eqnState.GetPressFromEnergy(data[0], energy, (*this).Velocity().Mag() );
  }

}

//member function to initialize a state with nondimensional values
void primVars::NondimensionalInitialize(const idealGas &eos, const vector3d<double> &vel){
  data[0] = 1.0;
  data[4] = 1.0/eos.Gamma();
  data[1] = vel.X() / eos.GetSoS(data[4], data[0]);
  data[2] = vel.Y() / eos.GetSoS(data[4], data[0]);
  data[3] = vel.Z() / eos.GetSoS(data[4], data[0]);
}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const primVars &prim){
  os << prim.data[0] << ", " << prim.data[1] << ", " << prim.data[2] << ", " << prim.data[3] << ", " << prim.data[4] << endl;
  return os;
}

//operator overload for addition
primVars primVars::operator + (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.data[0] += prim2.data[0];
  prim1.data[1] += prim2.data[1];
  prim1.data[2] += prim2.data[2];
  prim1.data[3] += prim2.data[3];
  prim1.data[4] += prim2.data[4];
  return prim1;
}

//operator overload for addition with a scalar
primVars operator+ (const double &scalar, const primVars &prim2){
  primVars prim1( prim2.data[0] + scalar,
		  prim2.data[1] + scalar,
		  prim2.data[2] + scalar,
		  prim2.data[3] + scalar,
		  prim2.data[4] + scalar);
  return prim1;
}

//operator overload for subtraction
primVars primVars::operator - (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.data[0] -= prim2.data[0];
  prim1.data[1] -= prim2.data[1];
  prim1.data[2] -= prim2.data[2];
  prim1.data[3] -= prim2.data[3];
  prim1.data[4] -= prim2.data[4];
  return prim1;
}

//operator overload for subtraction with a scalar
primVars operator- (const double &scalar, const primVars &prim2){
  primVars prim1( scalar - prim2.data[0],
		  scalar - prim2.data[1],
		  scalar - prim2.data[2],
		  scalar - prim2.data[3],
		  scalar - prim2.data[4] );
  return prim1;
}

//operator overload for elementwise multiplication
primVars primVars::operator * (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.data[0] *= prim2.data[0];
  prim1.data[1] *= prim2.data[1];
  prim1.data[2] *= prim2.data[2];
  prim1.data[3] *= prim2.data[3];
  prim1.data[4] *= prim2.data[4];
  return prim1;
}

//member function for scalar multiplication
primVars  primVars::operator * (const double &scalar){
  primVars temp = *this;
  temp.data[0] *= scalar;
  temp.data[1] *= scalar;
  temp.data[2] *= scalar;
  temp.data[3] *= scalar;
  temp.data[4] *= scalar;
  return temp;
}

//member function for scalar addition
primVars  primVars::operator + (const double &scalar){
  primVars temp = *this;
  temp.data[0] += scalar;
  temp.data[1] += scalar;
  temp.data[2] += scalar;
  temp.data[3] += scalar;
  temp.data[4] += scalar;
  return temp;
}

//member function for scalar subtraction
primVars  primVars::operator - (const double &scalar){
  primVars temp = *this;
  temp.data[0] -= scalar;
  temp.data[1] -= scalar;
  temp.data[2] -= scalar;
  temp.data[3] -= scalar;
  temp.data[4] -= scalar;
  return temp;
}

//member function for scalar division
primVars  primVars::operator / (const double &scalar){
  primVars temp = *this;
  temp.data[0] /= scalar;
  temp.data[1] /= scalar;
  temp.data[2] /= scalar;
  temp.data[3] /= scalar;
  temp.data[4] /= scalar;
  return temp;
}

//operator overload for multiplication with a scalar
primVars operator* (const double &scalar, const primVars &prim2){
  primVars prim1( prim2.data[0] * scalar,
		  prim2.data[1] * scalar,
		  prim2.data[2] * scalar,
		  prim2.data[3] * scalar,
		  prim2.data[4] * scalar);
  return prim1;
}

//operator overload for elementwise division
primVars primVars::operator / (const primVars& prim2)const{
  primVars prim1 = *this;
  prim1.data[0] /= prim2.data[0];
  prim1.data[1] /= prim2.data[1];
  prim1.data[2] /= prim2.data[2];
  prim1.data[3] /= prim2.data[3];
  prim1.data[4] /= prim2.data[4];
  return prim1;
}

//operator overload for division with a scalar
primVars operator/ (const double &scalar, const primVars &prim2){
  primVars prim1( scalar / prim2.data[0],
		  scalar / prim2.data[1],
		  scalar / prim2.data[2],
		  scalar / prim2.data[3],
		  scalar / prim2.data[4]);
  return prim1;
}


//member function to calculate reconstruction of primative variables from cell center to cell face
//this function uses MUSCL extrapolation resulting in higher order accuracy
/*

____________________|_____________________
|         |       Ul|Ur        |         |
|         |         |          |         |
|   Ui-1  |   Ui    |   Ui+1   |   Ui+2  |  
|         |         |          |         |
|         |         |          |         |
|_________|_________|__________|_________|
                    |
     face to reconstruct state at
<--UW2---><--UW1---><---DW---->

The diagram above shows the stencil used to reconstruct the left and right states for a given face. For each reconstruction (right and left) only 3 points are needed, two upwind
points and one downwind point. For the left reconstruction, Ui is the first upwind point, Ui-1 is the second upwind point, and Ui+1 is the downwind point. For the right 
reconstruction Ui+1 is the first upwind point, Ui+2 is the second upwind point, and Ui is the downwind point.

For left reconstruction the MUSCL scheme goes as follows:
Ui+1/2 = Ui + 0.5 * ( (1-K) * (Ui - Ui-1) + (1+K) * (Ui+1 - Ui) )

The above equation assumes there is no limiter, and that the grid spacing is uniform. In the above equation K refers to the parameter kappa which can be varied to produce different
reconstructions. Acceptable values of K are [-1,1] (inclusive)

K = -1 -- fully upwind reconstruction - linear - results in 2nd order accuracy
K = 0 -- fromm scheme - linear - results in 2nd order accuracy
K = 0.5 -- QUICK scheme - parabolic - results in 2nd order accuracy
K = 1 -- central scheme - linear - results in 2nd order accuracy but is unstable without dissipation added
K = 1/3 -- third order - parabolic - results 2nd order accuracy with lowest error
(all order of accuracy estimates assume constant fluxes on cell faces which limits order of accuracy to 2)

With limiters the equation looks as follows:
Ui+1/2 = Ui + 0.5 * L * ( (1-K*L) * (Ui - Ui-1) + (1+K*L) * (Ui+1 - Ui) )

L represents the limiter function which ranges between 0 and 1. A value of 1 implies there is no limiter and the full accuracy of the scheme is achieved. A value of 0 implies
that the solution has been limited to first order accuracy. An in between value results in an order of accuracy between first and second.

The MUSCL scheme can be extended to nonuniform grids by adding in terms representing the difference in size between the cells. In the above diagram the values UW2, UW1, and DW
represent the length of the second upwind, first upwind, and downwind cells respectively. dP and dM represent the factors due to the change in cell size between the first
upwind to downwind and first upwind to second upwind cells.

dP = (UW1 + DW) / (2.0 * UW)
dM = (UW + UW2) / (2.0 * UW)

Ui+1/2 = Ui + 0.5 * L * (1/(dP+dM)) *  ( (dP-K*L) * (Ui - Ui-1) / dM + (dM+K*L) * (Ui+1 - Ui) / dP )

*/
primVars primVars::FaceReconMUSCL( const primVars &primUW2, const primVars &primDW1, const double &kappa, const string &lim, double uw, double uw2, double dw)const{

  //primUW2 is the upwind cell furthest from the face at which the primative is being reconstructed.
  //primUW1 is the upwind cell nearest to the face at which the primative is being reconstructed. 
  //primDW1 is the downwind cell.
  //kappa is the parameter that determines which scheme is implemented
  //uw is length of upwind cell
  //uw2 is length of furthest upwind cell
  //dw is length of downwind cell

  primVars facePrim;
  primVars r;
  primVars limiter;
  primVars invLimiter;
  double eps = 1.0e-30;

  primVars primUW1 = *this;

  double dPlus = (uw + dw) / (2.0 * uw);
  double dMinus = (uw + uw2) / (2.0 * uw);

  //divided differences to base limiter on; eps must be listed to left of primVars
  r =  (eps + (primDW1 - primUW1) / dPlus) / (eps + (primUW1 - primUW2) / dMinus);  

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

  //calculate reconstructed state at face using MUSCL method with limiter
  facePrim = primUW1 + ( limiter / (2.0 * (dPlus + dMinus) ) ) * ( (dPlus  - kappa * limiter) * (primUW1 - primUW2) / dMinus + 
								   (dMinus + kappa * limiter) * (primDW1 - primUW1) / dPlus );

  return facePrim;

}

//member function to calculate minmod limiter
primVars primVars::LimiterMinmod( const primVars &upwind, const primVars &downwind, const double kap)const{
  //upwind -- upwind state (primative)
  //downwind -- downwind state (primative)
  //kap -- MUSCL parameter kappa

  primVars limiter, sign;

  //calculate minmod parameter beta
  double beta = (3.0 - kap) / (1.0 - kap);

  //calculate minmod limiter
  if (upwind.Rho() > 0.0 ){
    sign.data[0] = 1.0;
  }
  else if (upwind.Rho() < 0.0 ){
    sign.data[0] = -1.0;
  }
  else{
    sign.data[0] = 0.0;
  }

  if (upwind.U() > 0.0 ){
    sign.data[1] = 1.0;
  }
  else if (upwind.U() < 0.0 ){
    sign.data[1] = -1.0;
  }
  else{
    sign.data[1] = 0.0;
  }

  if (upwind.V() > 0.0 ){
    sign.data[2] = 1.0;
  }
  else if (upwind.V() < 0.0 ){
    sign.data[2] = -1.0;
  }
  else{
    sign.data[2] = 0.0;
  }

  if (upwind.W() > 0.0 ){
    sign.data[3] = 1.0;
  }
  else if (upwind.W() < 0.0 ){
    sign.data[3] = -1.0;
  }
  else{
    sign.data[3] = 0.0;
  }

  if (upwind.P() > 0.0 ){
    sign.data[4] = 1.0;
  }
  else if (upwind.P() < 0.0 ){
    sign.data[4] = -1.0;
  }
  else{
    sign.data[4] = 0.0;
  }


  limiter.data[0] = sign.Rho() * max(0.0, min( fabs(upwind.Rho()), sign.Rho()*downwind.Rho()*beta ) );
  limiter.data[1] = sign.U() * max(0.0, min( fabs(upwind.U()), sign.U()*downwind.U()*beta ) );
  limiter.data[2] = sign.V() * max(0.0, min( fabs(upwind.V()), sign.V()*downwind.V()*beta ) );
  limiter.data[3] = sign.W() * max(0.0, min( fabs(upwind.W()), sign.W()*downwind.W()*beta ) );
  limiter.data[4] = sign.P() * max(0.0, min( fabs(upwind.P()), sign.P()*downwind.P()*beta ) );

  return limiter;
}

//member function to calculate Van Albada limiter
primVars primVars::LimiterVanAlbada( const primVars &r)const{
  //r -- ratio of divided differences

  primVars limiter = (r + r*r)/(1 + r*r);
  //if value is negative, return zero
  limiter.data[0] = max(0.0, limiter.Rho());
  limiter.data[1] = max(0.0, limiter.U());
  limiter.data[2] = max(0.0, limiter.V());
  limiter.data[3] = max(0.0, limiter.W());
  limiter.data[4] = max(0.0, limiter.P());
  return limiter;
}

//member function to return no limiter
primVars primVars::LimiterNone()const{
  //for no limiter return all 1s
  primVars limiter(1.0, 1.0, 1.0, 1.0, 1.0);
  return limiter;
}

//member function to return the state of the appropriate ghost cell
/*

In the diagram below Ui represents the interior state. This should be represnted by (*this). Ug1 and Ug2 are the first and second layers of ghost cells respectively. They
are the output of the function. The input "layer" determines whether Ug1 or Ug2 is output for BCs that use extrapolation. Reflected boundaries (slipWall, viscousWall) do 
not use the layer input. Instead Ui+1 must be specified as (*this) to output Ug2.

____________________|___________
|         |         |          |
|         |         |          |
|   Ug2   |   Ug1   |   Ui     |
|         |         |          |
|         |         |          |
|_________|_________|__________|
                    |
             boundary face

Currently the following boundary conditions are supported: slipWall, viscousWall, characteristic, stagnationInlet, pressureOutlet, subsonicInflow, subsonicOutflow,
supersonicInflow, supersonicOutflow

*/
primVars primVars::GetGhostState( const string &bcType, const vector3d<double> &areaVec, const string &surf, const input &inputVars, const idealGas &eqnState, const int layer )const{
  //bcType -- type of boundary condition to supply ghost cell for
  //areaVec -- area vector of boundary face
  //surf -- i, j, k surface of boundary
  //inputVar -- all input variables
  //eqnState -- equation of state
  //layer -- layer of ghost cell to return (first (closest) or second (farthest))
  //the instance of primVars being acted upon should be the interior cell bordering the boundary

  primVars ghostState = (*this);  //set ghost state equal to boundary state to start

  //check to see that ghost layer corresponds to allowable number
  if ( !(layer == 1 || layer == 2) ){
    cerr << "ERROR: Error in primVars::GetGhostState. Requesting ghost state at a ghost layer " << layer << ". Please choose either 1 or 2" << endl;
    exit(0);
  }

  //normalize area vector (should always point out of domain)
  vector3d<double> normArea;
  if (surf == "il" || surf == "jl" || surf == "kl"){
    normArea = -1.0 * areaVec / areaVec.Mag(); //at lower surface normal should point out of domain for ghost cell calculation
  }
  else if (surf == "iu" || surf == "ju" || surf == "ku"){
    normArea = areaVec / areaVec.Mag(); 
  }

  double normVelCellCenter = 0;

  //slip wall boundary condition ----------------------------------------------------------------------------------------------------------------
  //slipWall is implemented as a reflection of the interior states so that there is no mass flow across the boundary.
  if (bcType == "slipWall"){             //for slip wall state should be reflected across boundary face, density and pressure stay equal to the boundary cell
    vector3d<double> stateVel = (*this).Velocity();
    normVelCellCenter = stateVel.DotProd(normArea);

    //for a slip wall the velocity of the boundary cell center is reflected across the boundary face to get the velocity at the ghost cell center
    vector3d<double> ghostVel (stateVel.X() - 2.0 * normArea.X() * normVelCellCenter,
			       stateVel.Y() - 2.0 * normArea.Y() * normVelCellCenter,
			       stateVel.Z() - 2.0 * normArea.Z() * normVelCellCenter);

    ghostState.data[1] = ghostVel.X();
    ghostState.data[2] = ghostVel.Y();
    ghostState.data[3] = ghostVel.Z();

    //numerical BCs for rho and pressure, same as boundary state

  }
  //viscous wall boundary condition ----------------------------------------------------------------------------------------------------------------
  //viscous wall uses the interior density and pressure, but flips the sign on the velocity so that the velocity at the boundary is 0.
  else if (bcType == "viscousWall"){             //for viscous wall velocity at face should be 0.0, density and pressure stay equal to the boundary cell
    vector3d<double> stateVel = (*this).Velocity();

    //ghost cell velocity at cell center is set to opposite of velocity at boundary cell center so that velocity at face will be zero
    vector3d<double> ghostVel (-1.0 * stateVel.X(),
			       -1.0 * stateVel.Y(),
			       -1.0 * stateVel.Z() );

    ghostState.data[1] = ghostVel.X();
    ghostState.data[2] = ghostVel.Y();
    ghostState.data[3] = ghostVel.Z();

    //numerical BCs for rho and pressure, same as boundary state

  }
  //subsonic inflow boundary condition --------------------------------------------------------------------------------------------------------------
  //this boundary condition enforces density and velocity as freestream inputs (constant) and extrapolates from the interior state to get pressure
  //this is a primative implementation, stagnationInlet or characteristic are better options
  else if (bcType == "subsonicInflow"){     //set velocity and density to freestream values
    double sos = eqnState.GetSoS( inputVars.PRef(), inputVars.RRef() );
    vector3d<double> ghostVel = inputVars.VelRef() / sos;   //nondimensionalize velocity

    ghostState.data[0] = 1.0;
    ghostState.data[1] = ghostVel.X();
    ghostState.data[2] = ghostVel.Y();
    ghostState.data[3] = ghostVel.Z();

    //numerical bc for pressure, same as boundary state

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  //subsonic outflow boundary condition ------------------------------------------------------------------------------------------------------------
  //this boundary condition enforces pressure as a freestream input (constant) and extrapolates density and velocity from the interior state
  //this is a primative implementation, pressureOutlet or characteristic are better options
  else if (bcType == "subsonicOutflow"){     //set pressure to freestream value
    ghostState.data[4] = 1.0 / eqnState.Gamma();

    //numerical bcs for density, velocity -- equal to boundary cell

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  //characteristic boundary condition ---------------------------------------------------------------------------------------------------------------
  //this is a characteristic based boundary condition which is appropriate for subsonic and supersonic flow. It automatically switches between inflow
  //and outflow as needed. It also automatically determines the number of characteristics that need to be specified at the boundary (subsonic vs supersonic)
  else if (bcType == "characteristic"){
    //freestream variables
    double freeSoS = eqnState.GetSoS(inputVars.PRef(), inputVars.RRef());
    vector3d<double> freeVel = inputVars.VelRef() / freeSoS;
    primVars freeState(1.0, 1.0/eqnState.Gamma(), freeVel);

    //internal variables
    double velIntNorm = (*this).Velocity().DotProd(normArea); 
    double SoSInt = eqnState.GetSoS((*this).P(), (*this).Rho());
    double machInt = fabs(velIntNorm)/SoSInt;

    if ( machInt >= 1.0 && velIntNorm < 0.0 ){ //supersonic inflow ------------------------------------------------------------
      //characteristics all go into the domain, so use freestream values for both riemann invariants
      ghostState = freeState;

    }
    else if ( machInt >= 1.0 && velIntNorm >= 0.0 ){ //supersonic outflow -----------------------------------------------------
      //characteristics all leave the domain, so use interior values for both riemann invariants
      ghostState = (*this);
    }
    else if ( machInt < 1.0 && velIntNorm < 0.0 ){ //subsonic inflow -----------------------------------------------------------
      //characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
      double rhoSoSInt = (*this).Rho() * SoSInt;
      vector3d<double> velDiff = freeState.Velocity() - (*this).Velocity();
      ghostState.data[4] = 0.5 * (freeState.P() + (*this).P() - rhoSoSInt * normArea.DotProd(velDiff)); //plus characteristic
      //minus characteristic
      ghostState.data[0] = freeState.Rho() + (ghostState.P() - freeState.P()) / (SoSInt * SoSInt);
      ghostState.data[1] = freeState.U() - normArea.X() * (freeState.P() - ghostState.P()) / rhoSoSInt;
      ghostState.data[2] = freeState.V() - normArea.Y() * (freeState.P() - ghostState.P()) / rhoSoSInt;
      ghostState.data[3] = freeState.W() - normArea.Z() * (freeState.P() - ghostState.P()) / rhoSoSInt;
    }
    else if ( machInt < 1.0 && velIntNorm >= 0.0 ){ //subsonic outflow ----------------------------------------------------------
      //characteristics go in both directions, use interior values for plus characteristic and freestream values for minus characteristic
      double rhoSoSInt = (*this).Rho() * SoSInt;
      ghostState.data[4] = freeState.P(); //minus characteristic
      //plus characteristic
      ghostState.data[0] = (*this).Rho() + ( ghostState.P() - (*this).P() ) / (SoSInt * SoSInt);
      ghostState.data[1] = (*this).U() + normArea.X() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;
      ghostState.data[2] = (*this).V() + normArea.Y() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;
      ghostState.data[3] = (*this).W() + normArea.Z() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;
    }
    else {
      cerr << "ERROR: flow condition for characteristic BC is not recognized!" << endl;
      cerr << "Interior state: " << (*this) << endl;
      cerr << "Ghost state: " << ghostState << endl;
      exit(0);
    }

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  //supersonic inflow boundary condition ------------------------------------------------------------------------------------------------------------
  //this boundary condition enforces the entire state as the specified freestream state
  else if (bcType == "supersonicInflow"){
    //physical boundary conditions - fix everything
    double sos = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef());
    vector3d<double> vel = inputVars.VelRef() / sos;                       //nondimensional velocity
    ghostState.data[0] = 1.0; //nondimensional density
    ghostState.data[1] = vel.X();
    ghostState.data[2] = vel.Y();
    ghostState.data[3] = vel.Z();
    ghostState.data[4] = 1.0 / eqnState.Gamma(); //nondimensional pressure
  }
  //supersonic outflow boundary condition ------------------------------------------------------------------------------------------------------------
  //this boundary condition enforces the entire state as extrapolated from the interior (zeroth order extrapolation)
  else if (bcType == "supersonicOutflow"){
    //do nothing and return boundary state -- numerical BCs for all
    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }
  }
  //stagnation inlet boundary condition --------------------------------------------------------------------------------------------------------------
  //this boundary condition is appropriate for subsonic flow. It is particularly well suited for internal flows. Implementation from Blazek
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

    ghostState.data[0] = eqnState.GetDensityTP(tb,pb);
    ghostState.data[1] = vbMag * inputVars.StagInletDx();
    ghostState.data[2] = vbMag * inputVars.StagInletDy();
    ghostState.data[3] = vbMag * inputVars.StagInletDz();
    ghostState.data[4] = pb;

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  //pressure outlet boundary condition ---------------------------------------------------------------------------------------------------------------
  //this boundary condition is appropriate for subsonic flow. Implementation from Blazek
  else if (bcType == "pressureOutlet"){
    double aRef = eqnState.GetSoS(inputVars.PRef(),inputVars.RRef()); //reference speed of sound
    double pb = inputVars.PressureOutletP() / (inputVars.RRef() * aRef * aRef); //nondimensional pressure from input file

    double SoSInt = (*this).SoS(eqnState);
    double rhoSoSInt = (*this).Rho() * SoSInt;
    ghostState.data[4] = pb;
    ghostState.data[0] = (*this).Rho() + ( ghostState.P() - (*this).P() ) / (SoSInt * SoSInt);
    ghostState.data[1] = (*this).U() + normArea.X() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;
    ghostState.data[2] =  (*this).V() + normArea.Y() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;
    ghostState.data[3] =  (*this).W() + normArea.Z() * ( (*this).P() - ghostState.P() ) / rhoSoSInt;

    if (layer == 2){ //extrapolate to get ghost state at 2nd layer
      ghostState = 2.0 * ghostState - (*this);
    }

  }
  //interblock boundary condition ---------------------------------------------------------------------------------------------------------------
  //this boundary condition is appropriate for point matched interfaces between physical blocks or processor blocks
  else if (bcType == "interblock"){
    //do nothing -- assign interior state to ghost state (already done)
    //for second layer of ghost cells interior state should be 2nd interior cell

  }
  else {
    cerr << "ERROR: Error in primVars::GetGhostState ghost state for BC type " << bcType << " is not supported!" << endl;
    cerr << "surface is " << surf << endl;
    exit(0);
  }

  return ghostState;

}

//member function to take in a genArray of updates to the conservative variables, and update the primative variables with it.
//this is used in the implicit solver
primVars primVars::UpdateWithConsVars(const idealGas &eqnState, const genArray &du)const{
  // eqnState -- equation of state
  // du -- updates to conservative variables

  //convert primative to conservative and update
  genArray consUpdate = (*this).ConsVars(eqnState) + du;

  //convert back to primative variables
  primVars primUpdate(consUpdate, false, eqnState);

  return primUpdate;
}

bool primVars::IsZero()const{

  bool nonzero = false;
  for ( int ii = 0; ii < NUMVARS; ii++ ){
    if ( (*this).data[ii] != 0.0 ){
      nonzero = true;
    }
  }
  return !nonzero;
}

#include "viscBlockVars.h"
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
using std::to_string;
using std::max;

//constructors
viscBlockVars::viscBlockVars(){
  length = 6;
  numI = 2;
  numJ = 2;
  numK = 2;

  int lenCell = (numI-1)*(numJ-1)*(numK-1);

  vector<tensor<double> > vGrad (lenCell);              //dummy velocity gradient tensor for cell centers
  vector<vector3d<double> > tGrad (lenCell);              //dummy temperature gradient vector3d for cell centers

  velGrad = vGrad;
  tempGrad = tGrad;

}
//constructor -- initialize state vector with dummy variables
viscBlockVars::viscBlockVars(const plot3dBlock &blk){
  numI = blk.NumI();
  numJ = blk.NumJ();
  numK = blk.NumK();
  length = numI * numJ * numK;

  int lenCell = (numI-1)*(numJ-1)*(numK-1);

  vector<tensor<double> > vGrad (lenCell);              //dummy velocity gradient tensor for cell centers
  vector<vector3d<double> > tGrad (lenCell);              //dummy temperature gradient vector3d for cell centers

  velGrad = vGrad;
  tempGrad = tGrad;

}

//---------------------------------------------------------------------------------------------------------------//
//member functions

//member function to calculate the velocity gradient at the cell center
void viscBlockVars::CalcVelGradGG(const vector3d<double> &vil, const vector3d<double> &viu, const vector3d<double> &vjl, const vector3d<double> &vju, const vector3d<double> &vkl, const vector3d<double> &vku,
			          const vector3d<double> &ail, const vector3d<double> &aiu, const vector3d<double> &ajl, const vector3d<double> &aju, const vector3d<double> &akl, const vector3d<double> &aku,
				  const double &vol, const int &loc){

  //vil is the velocity vector at the lower i face of the cell at which the velocity gradient is being calculated
  //viu is the velocity vector at the upper i face of the cell at which the velocity gradient is being calculated
  //vjl is the velocity vector at the lower j face of the cell at which the velocity gradient is being calculated
  //vju is the velocity vector at the upper j face of the cell at which the velocity gradient is being calculated
  //vkl is the velocity vector at the lower k face of the cell at which the velocity gradient is being calculated
  //vku is the velocity vector at the upper k face of the cell at which the velocity gradient is being calculated

  //ail is the area vector at the lower i face of the cell at which the velocity gradient is being calculated
  //aiu is the area vector at the upper i face of the cell at which the velocity gradient is being calculated
  //ajl is the area vector at the lower j face of the cell at which the velocity gradient is being calculated
  //aju is the area vector at the upper j face of the cell at which the velocity gradient is being calculated
  //akl is the area vector at the lower k face of the cell at which the velocity gradient is being calculated
  //aku is the area vector at the upper k face of the cell at which the velocity gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the velocity gradient should be stored

  tensor<double> temp;
  double invVol = 1.0/vol;

  //define velocity gradient tensor
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetXX( invVol * (viu.X()*aiu.X() - vil.X()*ail.X() + vju.X()*aju.X() - vjl.X()*ajl.X() + vku.X()*aku.X() - vkl.X()*akl.X()) );
  temp.SetXY( invVol * (viu.Y()*aiu.X() - vil.Y()*ail.X() + vju.Y()*aju.X() - vjl.Y()*ajl.X() + vku.Y()*aku.X() - vkl.Y()*akl.X()) );
  temp.SetXZ( invVol * (viu.Z()*aiu.X() - vil.Z()*ail.X() + vju.Z()*aju.X() - vjl.Z()*ajl.X() + vku.Z()*aku.X() - vkl.Z()*akl.X()) );

  temp.SetYX( invVol * (viu.X()*aiu.Y() - vil.X()*ail.Y() + vju.X()*aju.Y() - vjl.X()*ajl.Y() + vku.X()*aku.Y() - vkl.X()*akl.Y()) );
  temp.SetYY( invVol * (viu.Y()*aiu.Y() - vil.Y()*ail.Y() + vju.Y()*aju.Y() - vjl.Y()*ajl.Y() + vku.Y()*aku.Y() - vkl.Y()*akl.Y()) );
  temp.SetYZ( invVol * (viu.Z()*aiu.Y() - vil.Z()*ail.Y() + vju.Z()*aju.Y() - vjl.Z()*ajl.Y() + vku.Z()*aku.Y() - vkl.Z()*akl.Y()) );

  temp.SetZX( invVol * (viu.X()*aiu.Z() - vil.X()*ail.Z() + vju.X()*aju.Z() - vjl.X()*ajl.Z() + vku.X()*aku.Z() - vkl.X()*akl.Z()) );
  temp.SetYY( invVol * (viu.Y()*aiu.Z() - vil.Y()*ail.Z() + vju.Y()*aju.Z() - vjl.Y()*ajl.Z() + vku.Y()*aku.Z() - vkl.Y()*akl.Z()) );
  temp.SetZZ( invVol * (viu.Z()*aiu.Z() - vil.Z()*ail.Z() + vju.Z()*aju.Z() - vjl.Z()*ajl.Z() + vku.Z()*aku.Z() - vkl.Z()*akl.Z()) );

  (*this).SetVelGrad(temp, loc);

}

//member function to calculate the temperature gradient at the cell center
void viscBlockVars::CalcTempGradGG(const double &til, const double &tiu, const double &tjl, const double &tju, const double &tkl, const double &tku,
				 const vector3d<double> &ail, const vector3d<double> &aiu, const vector3d<double> &ajl, const vector3d<double> &aju, const vector3d<double> &akl, const vector3d<double> &aku,
				 const double &vol, const int &loc){

  //til is the temperature at the lower i face of the cell at which the temperature gradient is being calculated
  //tiu is the temperature at the upper i face of the cell at which the temperature gradient is being calculated
  //tjl is the temperature at the lower j face of the cell at which the temperature gradient is being calculated
  //tju is the temperature at the upper j face of the cell at which the temperature gradient is being calculated
  //tkl is the temperature at the lower k face of the cell at which the temperature gradient is being calculated
  //tku is the temperature at the upper k face of the cell at which the temperature gradient is being calculated

  //ail is the area vector at the lower i face of the cell at which the temperature gradient is being calculated
  //aiu is the area vector at the upper i face of the cell at which the temperature gradient is being calculated
  //ajl is the area vector at the lower j face of the cell at which the temperature gradient is being calculated
  //aju is the area vector at the upper j face of the cell at which the temperature gradient is being calculated
  //akl is the area vector at the lower k face of the cell at which the temperature gradient is being calculated
  //aku is the area vector at the upper k face of the cell at which the temperature gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the temperature gradient should be stored

  vector3d<double> temp;
  double invVol = 1.0/vol;

  //define temperature gradient vector
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetX( invVol * (tiu*aiu.X() - til*ail.X() + tju*aju.X() - tjl*ajl.X() + tku*aku.X() - tkl*akl.X()) );
  temp.SetY( invVol * (tiu*aiu.Y() - til*ail.Y() + tju*aju.Y() - tjl*ajl.Y() + tku*aku.Y() - tkl*akl.Y()) );
  temp.SetZ( invVol * (tiu*aiu.Z() - til*ail.Z() + tju*aju.Z() - tjl*ajl.Z() + tku*aku.Z() - tkl*akl.Z()) );

  (*this).SetTempGrad(temp, loc);

}

//member function to calculate the velocity gradient on the face using central differences
template <class T>
T viscBlockVars::FaceReconCentral(const T &velU, const T &velD, const vector3d<double> &pU, const vector3d<double> &pD, const vector3d<double> &pF)const{

  //velU is the velocity at the cell center of the upwind cell
  //velD is the velocity at the cell center of the downwind cell
  //pU is the position of the cell center of the upwind cell
  //pD is the position of the cell center of the downwind cell
  //pF is the position of the face center of the face on which the reconstruction is happening

  T temp;

  double cen2cen = pU.Distance(pD);  //distance from cell center to cell center
  double up2face = pU.Distance(pF);  //distance from upwind cell center to cell face

  temp = velD * (up2face/cen2cen) + velU * (1.0 - (up2face/cen2cen));

  return temp;

}


//member function to calculate gradients at centers
void viscBlockVars::CalcCellGrads(blockVars &vars, const idealGas &eqnState, const sutherland &suth, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;

  int loc = 0;
  int iLow = 0;
  int iUp = 0;
  int jLow = 0;
  int jUp = 0;
  int kLow = 0;
  int kUp = 0;

  int ifLow = 0;
  int ifUp = 0;
  int jfLow = 0;
  int jfUp = 0;
  int kfLow = 0;
  int kfUp = 0;

  string bcName = "undefined";

  double viscConstant = 1.0;

  vector3d<double> vil, viu, vjl, vju, vkl, vku;
  double til = 0.0;
  double tiu = 0.0;
  double tjl = 0.0;
  double tju = 0.0;
  double tkl = 0.0;
  double tku = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  primVars ghostState;
  vector3d<double> ghostDistance;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	iLow = GetNeighborLowI(ii, jj, kk, imax, jmax); 
	iUp  = GetNeighborUpI(ii, jj, kk, imax, jmax);
	jLow = GetNeighborLowJ(ii, jj, kk, imax, jmax);
	jUp  = GetNeighborUpJ(ii, jj, kk, imax, jmax);
	kLow = GetNeighborLowK(ii, jj, kk, imax, jmax);
	kUp  = GetNeighborUpK(ii, jj, kk, imax, jmax);

	ifLow = GetLowerFaceI(ii, jj, kk, imax, jmax); 
	ifUp  = GetUpperFaceI(ii, jj, kk, imax, jmax);
	jfLow = GetLowerFaceJ(ii, jj, kk, imax, jmax);
	jfUp  = GetUpperFaceJ(ii, jj, kk, imax, jmax);
	kfLow = GetLowerFaceK(ii, jj, kk, imax, jmax);
	kfUp  = GetUpperFaceK(ii, jj, kk, imax, jmax);

	//test if on i lower boundary
	if (ii == 0){
	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "il");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaI(ifLow), "il", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterI(ifLow) - vars.Center(loc) ) + vars.Center(loc);

	  vil = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterI(ifLow) );
	  til = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterI(ifLow) );

	}
	else{
	  vil = (*this).FaceReconCentral( vars.State(iLow).Velocity(), vars.State(loc).Velocity(), vars.Center(iLow), vars.Center(loc), vars.FCenterI(ifLow) );
	  til = (*this).FaceReconCentral( vars.State(iLow).Temperature(eqnState), vars.State(loc).Temperature(eqnState), vars.Center(iLow), vars.Center(loc), vars.FCenterI(ifLow) );
	}

	//test if on i upper boundary
	if (ii == imax-1){
	  bcName = bound.GetBCName(ii+1, jj, kk, "iu");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaI(ifUp), "iu", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterI(ifUp) - vars.Center(loc) ) + vars.Center(loc);

	  viu = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterI(ifUp) );
	  tiu = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterI(ifUp) );

	}
	else{
	  viu = (*this).FaceReconCentral( vars.State(iUp).Velocity(),  vars.State(loc).Velocity(), vars.Center(iUp),  vars.Center(loc), vars.FCenterI(ifUp)  );
	  tiu = (*this).FaceReconCentral( vars.State(iUp).Temperature(eqnState),  vars.State(loc).Temperature(eqnState), vars.Center(iUp),  vars.Center(loc), vars.FCenterI(ifUp)  );
	}

	//test if on j lower boundary
	if (jj == 0){
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaJ(jfLow), "jl", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterJ(jfLow) - vars.Center(loc) ) + vars.Center(loc);

	  vjl = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterJ(jfLow) );
	  tjl = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterJ(jfLow) );

	}
	else{
	  vjl = (*this).FaceReconCentral( vars.State(jLow).Velocity(), vars.State(loc).Velocity(), vars.Center(jLow), vars.Center(loc), vars.FCenterJ(jfLow)  );
	  tjl = (*this).FaceReconCentral( vars.State(jLow).Temperature(eqnState), vars.State(loc).Temperature(eqnState), vars.Center(jLow), vars.Center(loc), vars.FCenterJ(jfLow)  );
	}

	//test if on j upper boundary
	if (jj == jmax-1){
	  bcName = bound.GetBCName(ii, jj+1, kk, "ju");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaJ(jfUp), "ju", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterJ(jfUp) - vars.Center(loc) ) + vars.Center(loc);

	  vju = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterJ(jfUp) );
	  tju = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterJ(jfUp) );

	}
	else{
	  vju = (*this).FaceReconCentral( vars.State(jUp).Velocity(),  vars.State(loc).Velocity(), vars.Center(jUp),  vars.Center(loc), vars.FCenterJ(jfUp)  );
	  tju = (*this).FaceReconCentral( vars.State(jUp).Temperature(eqnState),  vars.State(loc).Temperature(eqnState), vars.Center(jUp),  vars.Center(loc), vars.FCenterJ(jfUp)  );
	}

	//test if on k lower boundary
	if(kk == 0){
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaK(kfLow), "kl", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterK(kfLow) - vars.Center(loc) ) + vars.Center(loc);

	  vkl = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterK(kfLow) );
	  tkl = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterK(kfLow) );

	}
	else{
	  vkl = (*this).FaceReconCentral( vars.State(kLow).Velocity(), vars.State(loc).Velocity(), vars.Center(kLow), vars.Center(loc), vars.FCenterK(kfLow)  );
	  tkl = (*this).FaceReconCentral( vars.State(kLow).Temperature(eqnState), vars.State(loc).Temperature(eqnState), vars.Center(kLow), vars.Center(loc), vars.FCenterK(kfLow)  );
	}

	//test if on k upper boundary
	if (kk == kmax-1){
	  bcName = bound.GetBCName(ii, jj, kk+1, "ku");

	  ghostState = vars.State(loc).GetGhostState( bcName, vars.FAreaK(kfUp), "ku", inp, eqnState );
	  ghostDistance = 2.0 * ( vars.FCenterK(kfUp) - vars.Center(loc) ) + vars.Center(loc);

	  vku = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(loc).Velocity(), ghostDistance, vars.Center(loc), vars.FCenterK(kfUp) );
	  tku = (*this).FaceReconCentral( ghostState.Temperature(eqnState), vars.State(loc).Temperature(eqnState), ghostDistance, vars.Center(loc), vars.FCenterK(kfUp) );

	}
	else{
	  vku = (*this).FaceReconCentral( vars.State(kUp).Velocity(),  vars.State(loc).Velocity(), vars.Center(kUp),  vars.Center(loc), vars.FCenterK(kfUp)  );
	  tku = (*this).FaceReconCentral( vars.State(kUp).Temperature(eqnState),  vars.State(loc).Temperature(eqnState), vars.Center(kUp),  vars.Center(loc), vars.FCenterK(kfUp)  );
	}

	//calculate gradients for cell
	CalcVelGradGG(vil, viu, vjl, vju, vkl, vku, vars.FAreaI(ifLow), vars.FAreaI(ifUp), vars.FAreaJ(jfLow), vars.FAreaJ(jfUp), vars.FAreaK(kfLow), vars.FAreaK(kfUp), vars.Vol(loc), loc);
	CalcTempGradGG(til, tiu, tjl, tju, tkl, tku, vars.FAreaI(ifLow), vars.FAreaI(ifUp), vars.FAreaJ(jfLow), vars.FAreaJ(jfUp), vars.FAreaK(kfLow), vars.FAreaK(kfUp), vars.Vol(loc), loc);

	//calculate cell viscous spectral radius
	double maxViscSpeed = ViscCellSpectralRadius(vars.FAreaI(ifLow), vars.FAreaI(ifUp), vars.State(loc), eqnState, suth, vars.Vol(loc));
	maxViscSpeed += ViscCellSpectralRadius(vars.FAreaJ(jfLow), vars.FAreaJ(jfUp), vars.State(loc), eqnState, suth, vars.Vol(loc));
	maxViscSpeed += ViscCellSpectralRadius(vars.FAreaK(kfLow), vars.FAreaK(kfUp), vars.State(loc), eqnState, suth, vars.Vol(loc));
	vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(loc) + viscConstant * (mRef/Re) * maxViscSpeed, loc); 

      }
    }
  }

}

// //member function to calculate viscous fluxes on i-faces
// void viscBlockVars::CalcViscFluxI(blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

//   int imax = (*this).NumI();
//   int jmax = (*this).NumJ() - 1;
//   int kmax = (*this).NumK() - 1;

//  const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int iLow = 0;
//   int iUp = 0;
//   int iFaceUp = 0;

//   string lstr = "left";
//   string rstr = "right";

//   string bcName = "undefined";

//   primVars ghostState;
//   vector3d<double> ghostDistance;

//   double viscConstant;
//   if (inp.Kappa() == -2.0 && inp.InviscidFlux() == "roe"){  //first order upwind
//     viscConstant = 2.0;
//   }
//   else if (inp.Kappa() == -1.0 && inp.InviscidFlux() == "roe"){ //second order upwind
//     viscConstant = 1.0;
//   }
//   else{
//     viscConstant = 4.0;
//   }

//   viscousFlux tempViscFlux;

//   tensor<double> velGrad;
//   vector3d<double> tGrad;
//   vector3d<double> vel;
//   double mu = 0.0;
//   //double rho = 0.0;
//   double maxViscSpeed = 0.0;

//   double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
//   double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
//   double mRef = inp.VelRef().Mag() / aRef;

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	if (ii == 0){ //-----------------------------------------------------------------------------------------------------------------------------------------

// 	  iUp  = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);
// 	  iFaceUp  = GetNeighborUpI(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "il");

// 	  ghostState = vars.State(iUp).GetGhostState( bcName, vars.FAreaI(loc), "il", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterI(loc) - vars.Center(iUp) ) + vars.Center(iUp);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(iUp);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(iUp).Velocity(), ghostDistance, vars.Center(iUp), vars.FCenterI(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(iUp);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(iUp).Temperature(eqnState)), ghostDistance, vars.Center(iUp), vars.FCenterI(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(iUp).Rho(), ghostDistance, vars.Center(iUp), vars.FCenterI(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaI(loc) );

// 	  //at lower boundary normal points into cell, so need to subtract from residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaI(loc).Mag(), iUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaI(loc).Mag() + vars.FAreaI(iFaceUp).Mag());
// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(iUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity(vars.State(iUp).Temperature(eqnState)), avgArea, vars.Vol(iUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(iUp) + viscConstant * maxViscSpeed, iUp); 

// 	}
// 	else if (ii == imax-1){ //---------------------------------------------------------------------------------------------------------------------------------------

// 	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "iu");

// 	  ghostState = vars.State(iLow).GetGhostState( bcName, vars.FAreaI(loc), "iu", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterI(loc) - vars.Center(iLow) ) + vars.Center(iLow);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(iLow);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(iLow).Velocity(), ghostDistance, vars.Center(iLow), vars.FCenterI(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(iLow);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(iLow).Temperature(eqnState)), ghostDistance, vars.Center(iLow), vars.FCenterI(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(iLow).Rho(), ghostDistance, vars.Center(iLow), vars.FCenterI(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaI(loc) );

// 	  //at upper boundary normal points out of cell, so need to add to residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaI(loc).Mag(), iLow);

// 	  //no wave speed calculation, this is only done at the lower face of the cell
// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  // maxViscSpeed = ViscCellSpectralRadius(rho, eqnState, mu, vars.FAreaI(loc).Mag(), vars.Vol(iLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(iLow) + viscConstant * maxViscSpeed, iLow); 

// 	}
// 	else{ //------------------------------------------------------------------------------------------------------------------------------------------------

// 	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
// 	  iUp  = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);
// 	  iFaceUp  = GetNeighborUpI(ii, jj, kk, imax, jmax);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).FaceReconCentral( (*this).VelGrad(iLow), (*this).VelGrad(iUp), vars.Center(iLow), vars.Center(iUp), vars.FCenterI(loc) );
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( vars.State(iLow).Velocity(), vars.State(iUp).Velocity(), vars.Center(iLow), vars.Center(iUp), vars.FCenterI(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).FaceReconCentral( (*this).TempGrad(iLow), (*this).TempGrad(iUp), vars.Center(iLow), vars.Center(iUp), vars.FCenterI(loc) );
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity( vars.State(iLow).Temperature(eqnState) ), suth.GetViscosity( vars.State(iUp).Temperature(eqnState) ), vars.Center(iLow), vars.Center(iUp), vars.FCenterI(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( vars.State(iLow).Rho(), vars.State(iUp).Rho(), vars.Center(iLow), vars.Center(iUp), vars.FCenterI(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaI(loc) );

// 	  //area vector points from left to right, so add to left cell, subtract from right cell
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaI(loc).Mag(), iLow);
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaI(loc).Mag(), iUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaI(loc).Mag() + vars.FAreaI(iFaceUp).Mag());
// 	  // maxViscSpeed = ViscCellSpectralRadius(rho, eqnState, mu, avgArea, vars.Vol(iLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(iLow) + viscConstant * maxViscSpeed, iLow); 

// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(iUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity( vars.State(iUp).Temperature(eqnState) ), avgArea, vars.Vol(iUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(iUp) + viscConstant * maxViscSpeed, iUp); 

// 	}

//       }
//     }
//   }


// }



// //member function to calculate viscous fluxes on j-faces
// void viscBlockVars::CalcViscFluxJ(blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

//   int imax = (*this).NumI() - 1;
//   int jmax = (*this).NumJ();
//   int kmax = (*this).NumK() - 1;

//  const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int jLow = 0;
//   int jUp = 0;
//   int jFaceUp = 0;

//   string lstr = "left";
//   string rstr = "right";

//   string bcName = "undefined";
//   primVars ghostState;
//   vector3d<double> ghostDistance;

//   viscousFlux tempViscFlux;

//   double viscConstant;
//   if (inp.Kappa() == -2.0 && inp.InviscidFlux() == "roe"){  //first order upwind
//     viscConstant = 2.0;
//   }
//   else if (inp.Kappa() == -1.0 && inp.InviscidFlux() == "roe"){ //second order upwind
//     viscConstant = 1.0;
//   }
//   else{
//     viscConstant = 4.0;
//   }

//   tensor<double> velGrad;
//   vector3d<double> tGrad;
//   vector3d<double> vel;
//   double mu = 0.0;
//   //double rho = 0.0;
//   double maxViscSpeed = 0.0;

//   double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
//   double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
//   double mRef = inp.VelRef().Mag() / aRef;

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	if (jj == 0){ //-------------------------------------------------------------------------------------------------------------------------------------

// 	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);
// 	  jFaceUp  = GetNeighborUpJ(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "jl");

// 	  ghostState = vars.State(jUp).GetGhostState( bcName, vars.FAreaJ(loc), "jl", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterJ(loc) - vars.Center(jUp) ) + vars.Center(jUp);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(jUp);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(jUp).Velocity(), ghostDistance, vars.Center(jUp), vars.FCenterJ(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(jUp);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(jUp).Temperature(eqnState)), ghostDistance, vars.Center(jUp), vars.FCenterJ(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(jUp).Rho(), ghostDistance, vars.Center(jUp), vars.FCenterJ(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaJ(loc) );

// 	  //at lower boundary normal points into cell, so need to subtract from residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaJ(loc).Mag(), jUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaJ(loc).Mag() + vars.FAreaJ(jFaceUp).Mag());
// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(jUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity(vars.State(jUp).Temperature(eqnState)), avgArea, vars.Vol(jUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jUp) + viscConstant * maxViscSpeed, jUp); 

// 	  // maxViscSpeed = max(4.0/(3.0 * vars.State(jUp).Rho()), eqnState.Gamma()/vars.State(jUp).Rho()) * (mu/eqnState.GetPrandtl()) * 
// 	  //   pow(0.5 * (vars.FAreaJ(loc).Mag() + vars.FAreaJ(jFaceUp).Mag()), 2.0) / vars.Vol(jUp) ;
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jUp) + viscConstant * maxViscSpeed, jUp); 

// 	}
// 	else if (jj == jmax-1){ //--------------------------------------------------------------------------------------------------------------------------------

// 	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "ju");

// 	  ghostState = vars.State(jLow).GetGhostState( bcName, vars.FAreaJ(loc), "ju", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterJ(loc) - vars.Center(jLow) ) + vars.Center(jLow);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(jLow);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(jLow).Velocity(), ghostDistance, vars.Center(jLow), vars.FCenterJ(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(jLow);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(jLow).Temperature(eqnState)), ghostDistance, vars.Center(jLow), vars.FCenterJ(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(jLow).Rho(), ghostDistance, vars.Center(jLow), vars.FCenterJ(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaJ(loc) );

// 	  //at upper boundary normal points out of cell, so need to add to residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaJ(loc).Mag(), jLow);

// 	  //no wave speed calculation, this is only done at the lower face of the cell
// 	  // maxViscSpeed = ViscCellSpectralRadius(rho, eqnState, mu, vars.FAreaJ(loc).Mag(), vars.Vol(jLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jLow) + viscConstant * maxViscSpeed, jLow); 

// 	}
// 	else{ //-----------------------------------------------------------------------------------------------------------------------------------------------

// 	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
// 	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);
// 	  jFaceUp  = GetNeighborUpJ(ii, jj, kk, imax, jmax);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).FaceReconCentral( (*this).VelGrad(jLow), (*this).VelGrad(jUp), vars.Center(jLow), vars.Center(jUp), vars.FCenterJ(loc) );
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( vars.State(jLow).Velocity(), vars.State(jUp).Velocity(), vars.Center(jLow), vars.Center(jUp), vars.FCenterJ(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).FaceReconCentral( (*this).TempGrad(jLow), (*this).TempGrad(jUp), vars.Center(jLow), vars.Center(jUp), vars.FCenterJ(loc) );
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity( vars.State(jLow).Temperature(eqnState) ), suth.GetViscosity( vars.State(jUp).Temperature(eqnState) ), vars.Center(jLow), vars.Center(jUp), vars.FCenterJ(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( vars.State(jLow).Rho(), vars.State(jUp).Rho(), vars.Center(jLow), vars.Center(jUp), vars.FCenterJ(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaJ(loc) );

// 	  //area vector points from left to right, so add to left cell, subtract from right cell
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaJ(loc).Mag(), jLow);
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaJ(loc).Mag(), jUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaJ(loc).Mag() + vars.FAreaJ(jFaceUp).Mag());
// 	  // maxViscSpeed = ViscCellSpectralRadius(rho, eqnState, mu, avgArea, vars.Vol(jLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jLow) + viscConstant * maxViscSpeed, jLow); 

// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(jUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity( vars.State(jUp).Temperature(eqnState)), avgArea, vars.Vol(jUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jUp) + viscConstant * maxViscSpeed, jUp); 


// 	  // maxViscSpeed = max(4.0/(3.0 * vars.State(jUp).Rho()), eqnState.Gamma()/vars.State(jUp).Rho()) * (mu/eqnState.GetPrandtl()) * 
// 	  //   pow(0.5 * (vars.FAreaJ(loc).Mag() + vars.FAreaJ(jFaceUp).Mag()), 2.0) / vars.Vol(jUp) ;
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(jUp) + viscConstant * maxViscSpeed, jUp); 

// 	}

//       }
//     }
//   }


// }


// //member function to calculate viscous fluxes on j-faces
// void viscBlockVars::CalcViscFluxK(blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

//   int imax = (*this).NumI() - 1;
//   int jmax = (*this).NumJ() - 1;
//   int kmax = (*this).NumK();

//  const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int kLow = 0;
//   int kUp = 0;
//   int kFaceUp = 0;

//   string lstr = "left";
//   string rstr = "right";

//   string bcName = "undefined";
//   primVars ghostState;
//   vector3d<double> ghostDistance;

//   viscousFlux tempViscFlux;

//   double viscConstant;
//   if (inp.Kappa() == -2.0 && inp.InviscidFlux() == "roe"){  //first order upwind
//     viscConstant = 2.0;
//   }
//   else if (inp.Kappa() == -1.0 && inp.InviscidFlux() == "roe"){ //second order upwind
//     viscConstant = 1.0;
//   }
//   else{
//     viscConstant = 4.0;
//   }

//   tensor<double> velGrad;
//   vector3d<double> tGrad;
//   vector3d<double> vel;
//   double mu = 0.0;
//   //double rho = 0.0;
//   double maxViscSpeed = 0.0;

//   double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
//   double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
//   double mRef = inp.VelRef().Mag() / aRef;

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	if (kk == 0){ //-----------------------------------------------------------------------------------------------------------------------------------

// 	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);
// 	  kFaceUp  = GetNeighborUpK(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "kl");

// 	  ghostState = vars.State(kUp).GetGhostState( bcName, vars.FAreaK(loc), "kl", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterK(loc) - vars.Center(kUp) ) + vars.Center(kUp);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(kUp);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(kUp).Velocity(), ghostDistance, vars.Center(kUp), vars.FCenterK(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(kUp);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(kUp).Temperature(eqnState)), ghostDistance, vars.Center(kUp), vars.FCenterK(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(kUp).Rho(), ghostDistance, vars.Center(kUp), vars.FCenterK(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaK(loc) );

// 	  //at lower boundary normal points into cell, so need to subtract from residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaK(loc).Mag(), kUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaK(loc).Mag() + vars.FAreaK(kFaceUp).Mag());
// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(kUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity(vars.State(kUp).Temperature(eqnState)), avgArea, vars.Vol(kUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(kUp) + viscConstant * maxViscSpeed, kUp); 

// 	  // maxViscSpeed = max(4.0/(3.0 * vars.State(kUp).Rho()), eqnState.Gamma()/vars.State(kUp).Rho()) * (mu/eqnState.GetPrandtl()) * 
// 	  //   pow(0.5 * (vars.FAreaK(loc).Mag() + vars.FAreaK(kFaceUp).Mag()), 2.0) / vars.Vol(kUp) ;
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(kUp) + viscConstant * maxViscSpeed, kUp); 

// 	}
// 	else if (kk == kmax-1){ //----------------------------------------------------------------------------------------------------------------------------

// 	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

// 	  //find boundary type and get ghost state
// 	  bcName = bound.GetBCName(ii, jj, kk, "ku");

// 	  ghostState = vars.State(kLow).GetGhostState( bcName, vars.FAreaK(loc), "ku", inp, eqnState );
// 	  ghostDistance = 2.0 * ( vars.FCenterK(loc) - vars.Center(kLow) ) + vars.Center(kLow);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).VelGrad(kLow);
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( ghostState.Velocity(), vars.State(kLow).Velocity(), ghostDistance, vars.Center(kLow), vars.FCenterK(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).TempGrad(kLow);
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity(vars.State(kLow).Temperature(eqnState)), ghostDistance, vars.Center(kLow), vars.FCenterK(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( ghostState.Rho(), vars.State(kLow).Rho(), ghostDistance, vars.Center(kLow), vars.FCenterK(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaK(loc) );

// 	  //at upper boundary normal points out of cell, so need to add to residual
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaK(loc).Mag(), kLow);

// 	  //no wave speed calculation, this is only done at the lower face of the cell
// 	  // maxViscSpeed = ViscCellSpectralRadius(rho, eqnState, mu, vars.FAreaK(loc).Mag(), vars.Vol(kLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(kLow) + viscConstant * maxViscSpeed, kLow); 

// 	}
// 	else{ //-----------------------------------------------------------------------------------------------------------------------------------------------------

// 	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
// 	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);
// 	  kFaceUp  = GetNeighborUpK(ii, jj, kk, imax, jmax);

// 	  //Get velocity gradient at face
// 	  velGrad = (*this).FaceReconCentral( (*this).VelGrad(kLow), (*this).VelGrad(kUp), vars.Center(kLow), vars.Center(kUp), vars.FCenterK(loc) );
// 	  //Get velocity at face
// 	  vel = (*this).FaceReconCentral( vars.State(kLow).Velocity(), vars.State(kUp).Velocity(), vars.Center(kLow), vars.Center(kUp), vars.FCenterK(loc) );
// 	  //Get temperature gradient at face
// 	  tGrad = (*this).FaceReconCentral( (*this).TempGrad(kLow), (*this).TempGrad(kUp), vars.Center(kLow), vars.Center(kUp), vars.FCenterK(loc) );
// 	  //Get viscosity at face
// 	  mu = (*this).FaceReconCentral( suth.GetViscosity( vars.State(kLow).Temperature(eqnState) ), suth.GetViscosity( vars.State(kUp).Temperature(eqnState) ), vars.Center(kLow), vars.Center(kUp), vars.FCenterK(loc) );
// 	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
// 	  //Get density at face
// 	  //rho = (*this).FaceReconCentral( vars.State(kLow).Rho(), vars.State(kUp).Rho(), vars.Center(kLow), vars.Center(kUp), vars.FCenterK(loc) );

// 	  //calculate viscous flux
// 	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, vars.FAreaK(loc) );

// 	  //area vector points from left to right, so add to left cell, subtract from right cell
// 	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
// 	  vars.AddToResidual(-1.0 * tempViscFlux * vars.FAreaK(loc).Mag(), kLow);
// 	  vars.AddToResidual(tempViscFlux * vars.FAreaK(loc).Mag(), kUp);

// 	  //accumulate wave speed contribution
// 	  //this is done on a cell basis; when at the lower face, value for cell will be calculated
// 	  double avgArea = 0.5 * (vars.FAreaK(loc).Mag() + vars.FAreaK(kFaceUp).Mag());
// 	  // maxViscSpeed = ViscCellSpectralRadius(vars.State(kUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity( vars.State(kUp).Temperature(eqnState), avgArea, vars.Vol(kLow));
// 	  // vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(kLow) + viscConstant * maxViscSpeed, kLow); 

// 	  maxViscSpeed = ViscCellSpectralRadius(vars.State(kUp).Rho(), eqnState, (mRef/Re) * suth.GetViscosity( vars.State(kUp).Temperature(eqnState)), avgArea, vars.Vol(kUp));
// 	  vars.SetAvgWaveSpeed( vars.AvgWaveSpeed(kUp) + viscConstant * maxViscSpeed, kUp); 

// 	}

//       }
//     }
//   }


// }

//function to calculate viscous flux jacobian on i-faces
void CalcViscFluxJacI(const blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag){

  int imax = vars.NumI();
  int jmax = vars.NumJ() - 1;
  int kmax = vars.NumK() - 1;

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int iLow = 0;
  int iUp = 0;

  string bcName = "undefined";

  primVars ghostState;
  vector3d<double> centerDist;

  double mu = 0.0;
  double rho = 0.0;
  double maxViscSpeed = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (ii == 0){ //-----------------------------------------------------------------------------------------------------------------------------------------

	  iUp  = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "il");

	  ghostState = vars.State(iUp).GetGhostState( bcName, vars.FAreaI(loc), "il", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterI( loc ) - vars.Center( iUp ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(iUp).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(iUp).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaI(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(iUp, mainDiag.Data(iUp) + 0.5 * maxViscSpeed * vars.FAreaI(loc).Mag());

	}
	else if (ii == imax-1){ //---------------------------------------------------------------------------------------------------------------------------------------

	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  ghostState = vars.State(iLow).GetGhostState( bcName, vars.FAreaI(loc), "iu", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterI( loc ) - vars.Center( iLow ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(iLow).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(iLow).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaI(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(iLow, mainDiag.Data(iLow) + 0.5 * maxViscSpeed * vars.FAreaI(loc).Mag());

	}
	else{ //------------------------------------------------------------------------------------------------------------------------------------------------

	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  iUp  = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity( vars.State(iLow).Temperature(eqnState) ) + suth.GetViscosity( vars.State(iUp).Temperature(eqnState) ) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( vars.State(iLow).Rho() + vars.State(iUp).Rho() );

	  //calculate viscous flux jacobian
	  //double vol = 0.5 * (vars.Vol(iUp) + vars.Vol(iLow));
	  squareMatrix tempViscFluxJacL(5);
	  squareMatrix tempViscFluxJacR(5);
	  centerDist = vars.Center( iLow ) - vars.Center( iUp ) ;
	  //CalcTSLFluxJac( mu, eqnState, vars.FAreaI(loc), vars.State(iLow), vars.State(iUp), centerDist.Mag(), tempViscFluxJacL, tempViscFluxJacR, suth );

	  //convective flux jacobians are subtracted from lower off diagonal and added to upper off diagonal
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  //area magnitude contribution is already incorporated into flux jacobian
	  // offLowIDiag.SetData( iUp, offLowIDiag.Data(iUp) + tempViscFluxJacL * vars.FAreaI(loc).Mag());
	  // offUpIDiag.SetData( iLow, offUpIDiag.Data(iLow) - tempViscFluxJacR * vars.FAreaI(loc).Mag());

	  //accumulate wave speed contribution -- area magnitude contribution is already incorporated into maxViscSpeed
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaI(loc), centerDist);
	  mainDiag.SetData(iLow, mainDiag.Data(iLow) + 0.5 * maxViscSpeed * vars.FAreaI(loc).Mag());

	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaI(loc), centerDist);
	  mainDiag.SetData(iUp, mainDiag.Data(iUp) + 0.5 * maxViscSpeed * vars.FAreaI(loc).Mag());

	}

      }
    }
  }


}

//function to calculate viscous flux jacobian on j-faces
void CalcViscFluxJacJ(const blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag){

  int imax = vars.NumI() - 1;
  int jmax = vars.NumJ();
  int kmax = vars.NumK() - 1;

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int jLow = 0;
  int jUp = 0;

  string bcName = "undefined";

  primVars ghostState;
  vector3d<double> centerDist;

  double mu = 0.0;
  double rho = 0.0;
  double maxViscSpeed = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (jj == 0){ //-----------------------------------------------------------------------------------------------------------------------------------------

	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  ghostState = vars.State(jUp).GetGhostState( bcName, vars.FAreaJ(loc), "il", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterJ( loc ) - vars.Center( jUp ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(jUp).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(jUp).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaJ(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(jUp, mainDiag.Data(jUp) + 0.5 * maxViscSpeed * vars.FAreaJ(loc).Mag());

	}
	else if (jj == jmax-1){ //---------------------------------------------------------------------------------------------------------------------------------------

	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  ghostState = vars.State(jLow).GetGhostState( bcName, vars.FAreaJ(loc), "iu", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterJ( loc ) - vars.Center( jLow ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(jLow).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(jLow).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaJ(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(jLow, mainDiag.Data(jLow) + 0.5 * maxViscSpeed * vars.FAreaJ(loc).Mag());

	}
	else{ //------------------------------------------------------------------------------------------------------------------------------------------------

	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity( vars.State(jLow).Temperature(eqnState) ) + suth.GetViscosity( vars.State(jUp).Temperature(eqnState) ) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( vars.State(jLow).Rho() + vars.State(jUp).Rho() );

	  //calculate viscous flux jacobian
	  //double vol = 0.5 * (vars.Vol(jUp) + vars.Vol(jLow));
	  squareMatrix tempViscFluxJacL(5);
	  squareMatrix tempViscFluxJacR(5);
	  centerDist = vars.Center( jLow ) - vars.Center( jUp ) ;
	  //CalcTSLFluxJac( mu, eqnState, vars.FAreaJ(loc), vars.State(jLow), vars.State(jUp), centerDist.Mag(), tempViscFluxJacL, tempViscFluxJacR, suth );

	  //convective flux jacobians are subtracted from lower off diagonal and added to upper off diagonal
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  //area magnitude contribution is already incorporated into flux jacobian
	  // offLowJDiag.SetData( jUp, offLowJDiag.Data(jUp) + tempViscFluxJacL * vars.FAreaJ(loc).Mag());
	  // offUpJDiag.SetData( jLow, offUpJDiag.Data(jLow) - tempViscFluxJacR * vars.FAreaJ(loc).Mag());

	  //accumulate wave speed contribution -- area magnitude contribution is already incorporated into maxViscSpeed
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaJ(loc), centerDist);
	  mainDiag.SetData(jLow, mainDiag.Data(jLow) + 0.5 * maxViscSpeed * vars.FAreaJ(loc).Mag());

	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaJ(loc), centerDist);
	  mainDiag.SetData(jUp, mainDiag.Data(jUp) + 0.5 * maxViscSpeed * vars.FAreaJ(loc).Mag());

	}

      }
    }
  }


}


//function to calculate viscous flux jacobian on k-faces
void CalcViscFluxJacK(const blockVars &vars, const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag){

  int imax = vars.NumI() - 1;
  int jmax = vars.NumJ() - 1;
  int kmax = vars.NumK();

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int kLow = 0;
  int kUp = 0;

  string bcName = "undefined";

  primVars ghostState;
  vector3d<double> centerDist;

  double mu = 0.0;
  double rho = 0.0;
  double maxViscSpeed = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (kk == 0){ //-----------------------------------------------------------------------------------------------------------------------------------------

	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  ghostState = vars.State(kUp).GetGhostState( bcName, vars.FAreaK(loc), "il", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterK( loc ) - vars.Center( kUp ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(kUp).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(kUp).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaK(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(kUp, mainDiag.Data(kUp) + 0.5 * maxViscSpeed * vars.FAreaK(loc).Mag());

	}
	else if (kk == kmax-1){ //---------------------------------------------------------------------------------------------------------------------------------------

	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  ghostState = vars.State(kLow).GetGhostState( bcName, vars.FAreaK(loc), "iu", inp, eqnState );
	  centerDist = 2.0 * ( vars.FCenterK( loc ) - vars.Center( kLow ) );

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity(ghostState.Temperature(eqnState)) + suth.GetViscosity(vars.State(kLow).Temperature(eqnState)) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( ghostState.Rho() + vars.State(kLow).Rho() );

	  //accumulate wave speed contribution
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaK(loc), centerDist);

	  //boundaries only contribute to main diagonal -- area magnitude contribution is already incorporated into maxViscSpeed
	  mainDiag.SetData(kLow, mainDiag.Data(kLow) + 0.5 * maxViscSpeed * vars.FAreaK(loc).Mag());

	}
	else{ //------------------------------------------------------------------------------------------------------------------------------------------------

	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  //Get viscosity at face
	  mu = 0.5 * ( suth.GetViscosity( vars.State(kLow).Temperature(eqnState) ) + suth.GetViscosity( vars.State(kUp).Temperature(eqnState) ) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  rho = 0.5 * ( vars.State(kLow).Rho() + vars.State(kUp).Rho() );

	  //calculate viscous flux jacobian
	  //double vol = 0.5 * (vars.Vol(kUp) + vars.Vol(kLow));
	  squareMatrix tempViscFluxJacL(5);
	  squareMatrix tempViscFluxJacR(5);
	  centerDist = vars.Center( kLow ) - vars.Center( kUp ) ;
	  //CalcTSLFluxJac( mu, eqnState, vars.FAreaK(loc), vars.State(kLow), vars.State(kUp), centerDist.Mag(), tempViscFluxJacL, tempViscFluxJacR, suth );

	  //convective flux jacobians are subtracted from lower off diagonal and added to upper off diagonal
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  //area magnitude contribution is already incorporated into flux jacobian
	  // offLowKDiag.SetData( kUp, offLowKDiag.Data(kUp) + tempViscFluxJacL * vars.FAreaK(loc).Mag());
	  // offUpKDiag.SetData( kLow, offLowKDiag.Data(kLow) - tempViscFluxJacR * vars.FAreaK(loc).Mag());

	  //accumulate wave speed contribution -- area magnitude contribution is already incorporated into maxViscSpeed
	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaK(loc), centerDist);
	  mainDiag.SetData(kLow, mainDiag.Data(kLow) + 0.5 * maxViscSpeed * vars.FAreaK(loc).Mag());

	  maxViscSpeed = ViscFaceSpectralRadiusTSL(rho, eqnState, mu, vars.FAreaK(loc), centerDist);
	  mainDiag.SetData(kUp, mainDiag.Data(kUp) + 0.5 * maxViscSpeed * vars.FAreaK(loc).Mag());

	}

      }
    }
  }


}

void viscBlockVars::CalcBlockTimeStep( blockVars &vars, const input &inputVars, const double &aRef){

  int imax = vars.NumI()-1;
  int jmax = vars.NumJ()-1;
  int kmax = vars.NumK()-1;

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  for ( kk = 0; kk < kmax; kk++ ){          //loop over all cells
    for ( jj = 0; jj < jmax; jj++ ){          
      for ( ii = 0; ii < imax; ii++ ){          

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (inputVars.Dt() > 0.0){   //dt specified, use global time stepping
	  vars.SetDt(inputVars.Dt() * aRef, loc);
	}
	else if (inputVars.CFL() > 0.0){ //cfl specified, use local time stepping
	  vars.CalcCellDt(ii, jj, kk, inputVars.CFL());
	}
	else{
	  cerr << "ERROR: Neither dt or cfl was specified!" << endl;
	  exit(0);
	}

      }
    }
  }

}


//---------------------------------------------------------------------------------------------------------------//
//function declarations

//function to calculate the spectral radius on a cell face for the viscous fluxes using the thin shear layer approximation
double ViscFaceSpectralRadiusTSL(const double &rho, const idealGas &eqnState, const double &mu, const vector3d<double> &fArea, const vector3d<double> &dist){

  //rho is density at cell center
  //mu is viscoisty at cell center
  //eqnState is equation of state, used to get prandtl number
  //fMag is average face area (magnitude)
  //dist is the distance from cell center to cell center of the cells bounding the face

  vector3d<double> normArea = fArea / fArea.Mag();

  //return max(4.0/3.0, eqnState.Gamma() / eqnState.GetPrandtl()) * (mu * fMag * fMag) / (rho * dist) ;
  //return fMag * max( (4.0 * mu)/3.0, eqnState.Gamma() * mu / eqnState.GetPrandtl() ) / rho;
  //return 2.0 * fMag * fMag / ( rho * dist);
  //return max( (4.0 * mu) / (3.0 * rho * dist) , eqnState.Gamma() * mu / (eqnState.GetPrandtl() * rho * dist) );
  //return fMag * mu / (eqnState.GetPrandtl() * rho * dist) * max(4.0/3.0, eqnState.Gamma());
  //return eqnState.Gamma() * mu / (eqnState.GetPrandtl() * rho * dist);

  return 2.0 * mu / (rho * fabs(normArea.DotProd(dist)) ) ;
}

//function to calculate the spectral radius at a cell center for the viscous fluxes
double ViscCellSpectralRadius(const vector3d<double> &fAreaL, const vector3d<double> &fAreaR, const primVars &state, const idealGas &eqnState, const sutherland &suth, const double &vol){

  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  double maxTerm = max(4.0 / (3.0 * state.Rho()), eqnState.Gamma() / state.Rho()) ;
  double mu = suth.GetViscosity(state.Temperature(eqnState));
  double viscTerm = mu / eqnState.GetPrandtl();

  return maxTerm * viscTerm * fMag * fMag / vol ;
}



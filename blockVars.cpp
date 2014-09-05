#include "blockVars.h"
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
blockVars::blockVars(){
  length = 6;
  numI = 2;
  numJ = 2;
  numK = 2;

  int lenCellCen = (numI-1)*(numJ-1)*(numK-1);  

  vector<primVars> dummyState (lenCellCen);              //dummy state variable
  vector<vector3d<double> > vec1(length);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vec2(lenCellCen);             //dummy vector variable lenght of number of cells
  vector<double> scalar(lenCellCen);                      //dummy scalar variable lenght of number of cells
  colMatrix singleResid(5);
  singleResid.Zero();
  vector<colMatrix> dummyResid(lenCellCen, singleResid);

  state = dummyState;      

  vol = scalar;
  center = vec2;
  fAreaI = vec1;
  fAreaJ = vec1;
  fAreaK = vec1;
  fCenterI = vec1;
  fCenterJ = vec1;
  fCenterK = vec1;

  avgWaveSpeed = scalar;

  dt = scalar;

  residual = dummyResid;

}
//constructor -- initialize state vector with dummy variables
blockVars::blockVars(const plot3dBlock &blk){
  numI = blk.NumI();
  numJ = blk.NumJ();
  numK = blk.NumK();
  length = numI * numJ * numK;
  
  int lenCellCen = (numI-1)*(numJ-1)*(numK-1);  
  vector<primVars> dummyState (lenCellCen);              //dummy state variable
  vector<double> dummyScalar (lenCellCen);                 //dummy time variable
  colMatrix singleResid(5);
  singleResid.Zero();
  vector<colMatrix> dummyResid(lenCellCen, singleResid);

  state = dummyState;      

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  avgWaveSpeed = dummyScalar;

  dt = dummyScalar;

  residual = dummyResid;
}

//constructor -- assign passed variables to initialize state vector
blockVars::blockVars( const double density, const double pressure, const vector3d<double> vel, const plot3dBlock &blk){
  numI = blk.NumI();
  numJ = blk.NumJ();
  numK = blk.NumK();
  length = numI * numJ * numK;
  
  int lenCellCen = (numI-1)*(numJ-1)*(numK-1);  
  primVars singleState(density, pressure, vel);
  vector<primVars> initState(lenCellCen,singleState);       //initialize state vector
  vector<double> dummyScalar (lenCellCen);                  //dummy time variable

  colMatrix singleResid(5);
  singleResid.Zero();
  vector<colMatrix> dummyResid(lenCellCen, singleResid);

  state = initState;      

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  avgWaveSpeed = dummyScalar;

  dt = dummyScalar;

  residual = dummyResid;
}

//constructor -- assign passed state to initialize state vector
blockVars::blockVars( const primVars& inputState, const plot3dBlock &blk){
  numI = blk.NumI();
  numJ = blk.NumJ();
  numK = blk.NumK();
  length = numI * numJ * numK;
  
  int lenCellCen = (numI-1)*(numJ-1)*(numK-1);  
  vector<primVars> initState(lenCellCen, inputState);   //initialize state vector
  vector<double> dummyScalar (lenCellCen);                 //dummy time variable

  colMatrix singleResid(5);
  singleResid.Zero();
  vector<colMatrix> dummyResid(lenCellCen, singleResid);

  state = initState;      

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  avgWaveSpeed = dummyScalar;

  dt = dummyScalar;

  residual = dummyResid;
}

//member function to store the inviscid flux class in the place for the residual
void blockVars::AddToResidual(const inviscidFlux &flux, const int &ii){
  colMatrix temp(5);
  temp.SetData(0, flux.RhoVel());
  temp.SetData(1, flux.RhoVelU());
  temp.SetData(2, flux.RhoVelV());
  temp.SetData(3, flux.RhoVelW());
  temp.SetData(4, flux.RhoVelH());

  (*this).SetResidual( (*this).Residual(ii) + temp, ii); 
}

//member function to store the viscous flux class in the place for the residual
void blockVars::AddToResidual(const viscousFlux &flux, const int &ii){
  colMatrix temp(5);
  temp.SetData(0, 0.0);
  temp.SetData(1, flux.MomX());
  temp.SetData(2, flux.MomY());
  temp.SetData(3, flux.MomZ());
  temp.SetData(4, flux.Engy());

  (*this).SetResidual( (*this).Residual(ii) + temp, ii); 
}

//---------------------------------------------------------------------------------------------------------------//
//function declarations
//function to calculate the fluxes on the i-faces
void blockVars::CalcInvFluxI(const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

 const boundaryConditions bound = inp.BC(bb);
 const double kap = inp.Kappa();

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  string lstr = "left";
  string rstr = "right";

  double maxWS = 0.0;

  double up2faceL, upwindL, centralL, up2faceU, upwindU, centralU;

  primVars faceStateLower, faceStateUpper, ghostState;

  inviscidFlux tempFlux;

  string bcName = "undefined";

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	//find out if at a block boundary
	if ( ii == 0  ){                             //at i lower boundary ---------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "il");

	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  if (imax > 2 && kap != -2.0){ //if more than 2 faces thick, and second order, use linear extrapolation to get boundary state
	    up2faceU = (*this).Center( upperI ).Distance( (*this).FCenterI(loc) );
	    upwindU = (*this).Center( upperI ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2)) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( upperI ), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "il", maxWS, up2faceU, upwindU );
	  }
	  else{  //if not more than 2 faces thick, use cell adjacent to boundary
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( upperI ), (*this).State( upperI ), eqnState, inp, "il", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), upperI);
	}
	else if ( ii == imax-1 ){  //at i upper boundary -------------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  if (imax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( lowerI ).Distance( (*this).FCenterI(loc) );
	    upwindL = (*this).Center( lowerI ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( lowerI ), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "iu", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( lowerI ), (*this).State( lowerI ), eqnState, inp, "iu", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), lowerI);
	}
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( ii == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow
                                                                  //for slipwall need to mirror values
	                                                          //ghost state should use boundary adjacent cell and boundary normal
	  bcName = bound.GetBCName(ii-1, jj, kk, "il");           //get bc at ii=0

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( lowerI ).GetGhostState( bcName, (*this).FAreaI( GetNeighborLowI(ii, jj, kk, imax, jmax, 1) ), "il", inp, eqnState );

	  up2faceL = (*this).Center( lowerI ).Distance( (*this).FCenterI(loc) );
	  upwindL = (*this).FCenterI(loc).Distance( (*this).FCenterI( GetNeighborLowI(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( lowerI ).Distance( (*this).Center( upperI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( ghostState, (*this).State( lowerI ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperI ).Distance( (*this).FCenterI(loc) );
	  upwindU = (*this).Center( upperI ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( upperI ).Distance( (*this).Center( lowerI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ), (*this).State( lowerI ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), upperI);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), lowerI);


	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( ii == imax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii+1, jj, kk, "iu");          //get bc at ii = imax-1

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( upperI ).GetGhostState( bcName, (*this).FAreaI( GetNeighborUpI(ii, jj, kk, imax, jmax, 1) ), "iu", inp, eqnState );

	  up2faceL = (*this).Center( lowerI ).Distance( (*this).FCenterI(loc) );
	  upwindL = (*this).Center( lowerI ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( lowerI ).Distance( (*this).Center( upperI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
	  		   (*this).State( upperI ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperI ).Distance( (*this).FCenterI(loc) );
	  upwindU = (*this).FCenterI(loc).Distance( (*this).FCenterI( GetNeighborUpI(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( upperI ).Distance( (*this).Center( lowerI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( ghostState, (*this).State( lowerI ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), upperI);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), lowerI);

	}
	else{  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower i face

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  if (kap == -2.0){  //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerI ).FaceReconConst();
	    faceStateUpper = (*this).State( upperI ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( lowerI ).Distance( (*this).FCenterI(loc) );
	    upwindL = (*this).Center( lowerI ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( lowerI ).Distance( (*this).Center( upperI ) );

	    faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperI ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( upperI ).Distance( (*this).FCenterI(loc) );
	    upwindU = (*this).Center( upperI ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( upperI ).Distance( (*this).Center( lowerI ) );

	    faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerI ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), upperI);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerI) + 0.5 * maxWS * (*this).FAreaI(loc).Mag(), lowerI);

	}

      }
    }
  }


}

//function to calculate the fluxes on the j-faces
void blockVars::CalcInvFluxJ(const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);
  const double kap = inp.Kappa();

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  string lstr = "left";
  string rstr = "right";

  double maxWS = 0.0;

  double up2faceL, upwindL, centralL, up2faceU, upwindU, centralU;

  primVars faceStateLower, faceStateUpper, ghostState;

  inviscidFlux tempFlux;

  string bcName = "undefined";

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	//find out if at a block boundary
	if ( jj == 0  ){                             //at j lower boundary --------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  if (jmax > 2 && kap != -2.0){
	    up2faceU = (*this).Center( upperJ ).Distance( (*this).FCenterJ(loc) );
	    upwindU = (*this).Center( upperJ ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( upperJ ), 
				     (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "jl", maxWS, up2faceU, upwindU );
	  }
	  else {
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( upperJ ), (*this).State( upperJ ), eqnState, inp, "jl", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), upperJ);

	}
	else if ( jj == jmax-1 ){  //at j upper boundary ---------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

	  if (jmax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( lowerJ ).Distance( (*this).FCenterJ(loc) );
	    upwindL = (*this).Center( lowerJ ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( lowerJ ), 
				     (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ju", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( lowerJ ), (*this).State( lowerJ ), eqnState, inp, "ju", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), lowerJ);

	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( jj == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj-1, kk, "jl");           //get bc at jj=0

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( lowerJ ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborLowJ(ii, jj, kk, imax, jmax, 1) ), "jl", inp, eqnState );

	  up2faceL = (*this).Center( lowerJ ).Distance( (*this).FCenterJ(loc) );
	  upwindL = (*this).FCenterJ(loc).Distance( (*this).FCenterJ( GetNeighborLowJ(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( lowerJ ).Distance( (*this).Center( upperJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( ghostState, (*this).State( upperJ ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperJ ).Distance( (*this).FCenterJ(loc) );
	  upwindU = (*this).Center( upperJ ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( upperJ ).Distance( (*this).Center( lowerJ ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( lowerJ ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), upperJ);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), lowerJ);

	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( jj == jmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj+1, kk, "ju");          //get bc at jj=jmax-1

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( upperJ ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborUpJ(ii, jj, kk, imax, jmax, 1) ), "ju", inp, eqnState );

	  up2faceL = (*this).Center( lowerJ ).Distance( (*this).FCenterJ(loc) );
	  upwindL = (*this).Center( lowerJ ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( lowerJ ).Distance( (*this).Center( upperJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( upperJ ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperJ ).Distance( (*this).FCenterJ(loc) );
	  upwindU = (*this).FCenterJ(loc).Distance( (*this).FCenterJ( GetNeighborUpJ(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( upperJ ).Distance( (*this).Center( lowerJ ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( ghostState, (*this).State( lowerJ ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), upperJ);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), lowerJ);

	}
	else{ // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower j face

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerJ ).FaceReconConst();
	    faceStateUpper = (*this).State( upperJ ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( lowerJ ).Distance( (*this).FCenterJ(loc) );
	    upwindL = (*this).Center( lowerJ ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( lowerJ ).Distance( (*this).Center( upperJ ) );

	    faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperJ ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( upperJ ).Distance( (*this).FCenterJ(loc) );
	    upwindU = (*this).Center( upperJ ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( upperJ ).Distance( (*this).Center( lowerJ ) );

	    faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerJ ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );
	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), upperJ);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerJ) + 0.5 * maxWS * (*this).FAreaJ(loc).Mag(), lowerJ);

	}


      }
    }
  }


}

//function to calculate the fluxes on the k-faces
void blockVars::CalcInvFluxK(const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK();

  const boundaryConditions bound = inp.BC(bb);
  const double kap = inp.Kappa();

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  string lstr = "left";
  string rstr = "right";

  double maxWS = 0.0;

  double up2faceL, upwindL, centralL, up2faceU, upwindU, centralU;

  primVars faceStateLower, faceStateUpper, ghostState;

  inviscidFlux tempFlux;

  string bcName = "undefined";

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	//find out if at a block boundary
	if ( kk == 0  ){                             //at k lower boundary -------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  if (kmax > 2 && kap != -2.0){
	    up2faceU = (*this).Center( upperK ).Distance( (*this).FCenterK(loc) );
	    upwindU = (*this).Center( upperK ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( upperK ), 
				     (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "kl", maxWS, up2faceU, upwindU );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( upperK ), (*this).State( upperK ), eqnState, inp, "kl", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), upperK);

	}
	else if ( kk == kmax-1 ){  //at k upper boundary --------------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

	  if (kmax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( lowerK ).Distance( (*this).FCenterK(loc) );
	    upwindL = (*this).Center( lowerK ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( lowerK ), 
				     (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ku", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( lowerK ), (*this).State( lowerK ), eqnState, inp, "ku", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), lowerK);

	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( kk == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk-1, "kl");           //get bc at kk=0

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( upperK ).GetGhostState( bcName, (*this).FAreaK( GetNeighborLowK(ii, jj, kk, imax, jmax, 1) ), "kl", inp, eqnState );

	  up2faceL = (*this).Center( lowerK ).Distance( (*this).FCenterK(loc) );
	  upwindL = (*this).FCenterK(loc).Distance( (*this).FCenterK( GetNeighborLowK(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( lowerK ).Distance( (*this).Center( upperK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( ghostState, (*this).State( upperK ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperK ).Distance( (*this).FCenterK(loc) );
	  upwindU = (*this).Center( upperK ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( upperK ).Distance( (*this).Center( lowerK ) );

	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( lowerK ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), upperK);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), lowerK);

	}
	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( kk == kmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk+1, "ku");          //get bc at kk=kmax-1

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( upperK ).GetGhostState( bcName, (*this).FAreaK( GetNeighborUpK(ii, jj, kk, imax, jmax, 1) ), "ku", inp, eqnState );

	  up2faceL = (*this).Center( lowerK ).Distance( (*this).FCenterK(loc) );
	  upwindL = (*this).Center( lowerK ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( lowerK ).Distance( (*this).Center( upperK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( upperK ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( upperK ).Distance( (*this).FCenterK(loc) );
	  upwindU = (*this).FCenterK(loc).Distance( (*this).FCenterK( GetNeighborUpK(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( upperK ).Distance( (*this).Center( lowerK ) );
	  
	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( ghostState, (*this).State( lowerK ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), upperK);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), lowerK);

	}
	else{ // -------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower k face
	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerK ).FaceReconConst();
	    faceStateUpper = (*this).State( upperK ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( lowerK ).Distance( (*this).FCenterK(loc) );
	    upwindL = (*this).Center( lowerK ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( lowerK ).Distance( (*this).Center( upperK ) );

	    faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperK ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( upperK ).Distance( (*this).FCenterK(loc) );
	    upwindU = (*this).Center( upperK ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( upperK ).Distance( (*this).Center( lowerK ) );

	    faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerK ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), upperK);
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(lowerK) + 0.5 * maxWS * (*this).FAreaK(loc).Mag(), lowerK);

	}


      }
    }
  }


}

//member function to calculate the local time step. (i,j,k) are cell indices
void blockVars::CalcCellDt( const int &i, const int &j, const int &k, const double &cfl){

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;
  int loc = GetLoc1D(i, j, k, imax, jmax);

  double dt = (cfl * (*this).Vol(loc) / (*this).AvgWaveSpeed(loc)) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}

//this member function calculates the residual for an inviscid simulation
// void blockVars::CalcCellResidual(const int &ii, const int &jj, const int &kk, const int &imax, const int &jmax){

//   colMatrix resid(5);

//   int loc = GetLoc1D(ii, jj, kk, imax, jmax);
//   int iLow = GetLowerFaceI(ii, jj, kk, imax, jmax); 
//   int iUp  = GetUpperFaceI(ii, jj, kk, imax, jmax);
//   int jLow = GetLowerFaceJ(ii, jj, kk, imax, jmax);
//   int jUp  = GetUpperFaceJ(ii, jj, kk, imax, jmax);
//   int kLow = GetLowerFaceK(ii, jj, kk, imax, jmax);
//   int kUp  = GetUpperFaceK(ii, jj, kk, imax, jmax);

//   //Area vector points nominally from lower index to upper index, so the upper index fluxes must be multiplied by -1 so vector points into cell and for conservation
//   resid.SetData(0,        (*this).InvFluxI(iLow).RhoVel()  * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVel()  * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVel()  * (*this).FAreaK(kLow).Mag() 
// 		-1.0 * (*this).InvFluxI(iUp).RhoVel()   * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVel()   * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVel()   * (*this).FAreaK(kUp).Mag() );
//   resid.SetData(1,        (*this).InvFluxI(iLow).RhoVelU() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelU() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelU() * (*this).FAreaK(kLow).Mag() 
// 		-1.0 * (*this).InvFluxI(iUp).RhoVelU()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelU()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelU()  * (*this).FAreaK(kUp).Mag() );
//   resid.SetData(2,        (*this).InvFluxI(iLow).RhoVelV() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelV() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelV() * (*this).FAreaK(kLow).Mag() 
// 		-1.0 * (*this).InvFluxI(iUp).RhoVelV()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelV()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelV()  * (*this).FAreaK(kUp).Mag() );
//   resid.SetData(3,        (*this).InvFluxI(iLow).RhoVelW() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelW() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelW() * (*this).FAreaK(kLow).Mag() 
// 		-1.0 * (*this).InvFluxI(iUp).RhoVelW()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelW()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelW()  * (*this).FAreaK(kUp).Mag() );
//   resid.SetData(4,        (*this).InvFluxI(iLow).RhoVelH() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelH() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelH() * (*this).FAreaK(kLow).Mag() 
// 		-1.0 * (*this).InvFluxI(iUp).RhoVelH()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelH()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelH()  * (*this).FAreaK(kUp).Mag() );

//   (*this).SetResidual(resid, loc);

// }


void blockVars::CalcBlockTimeStep( const input &inputVars, const double &aRef){

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;
  int kmax = (*this).NumK()-1;

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  for ( kk = 0; kk < kmax; kk++ ){          //loop over all cells
    for ( jj = 0; jj < jmax; jj++ ){          
      for ( ii = 0; ii < imax; ii++ ){          

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	//(*this).CalcCellResidual(ii, jj, kk, imax, jmax);

	if (inputVars.Dt() > 0.0){   //dt specified, use global time stepping
	  (*this).SetDt(inputVars.Dt() * aRef, loc);
	}
	else if (inputVars.CFL() > 0.0){ //cfl specified, use local time stepping
	  (*this).CalcCellDt(ii, jj, kk, inputVars.CFL());
	}
	else{
	  cerr << "ERROR: Neither dt or cfl was specified!" << endl;
	  exit(0);
	}

      }
    }
  }

}

void blockVars::UpdateBlock(const input &inputVars, const int &impFlag, const idealGas &eos, const double &aRef, const int &bb, const vector<colMatrix> &du, vector<double> &l2, vector<double> &linf, int &locMaxB){

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;
  int kmax = (*this).NumK()-1;

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  if ( inputVars.TimeIntegration() != "rk4" ){
    for ( kk = 0; kk < kmax; kk++ ){          //loop over all cells
      for ( jj = 0; jj < jmax; jj++ ){          
	for ( ii = 0; ii < imax; ii++ ){          

	  loc = GetLoc1D(ii, jj, kk, imax, jmax);

	  if (inputVars.TimeIntegration() == "explicitEuler"){
	    (*this).ExplicitEulerTimeAdvance(eos, loc);
	  }
	  else if (impFlag){
	    (*this).ImplicitTimeAdvance(du[loc], eos, loc);
	  }


	  for ( unsigned int ll = 0; ll < l2.size(); ll++ ){
	    l2[ll] = l2[ll] + (*this).Residual(loc,ll) * (*this).Residual(loc,ll);

	    if ( (*this).Residual(loc,ll) > linf[4] ){
	      linf[4] = (*this).Residual(loc,ll);
	      linf[3] = (double)ll+1;
	      linf[2] = (double)kk;
	      linf[1] = (double)jj;
	      linf[0] = (double)ii;
	      locMaxB = bb;
	    }
	  }


	}
      }
    }
  }
  else if ( inputVars.TimeIntegration() == "rk4" ){
    int rr = 0;
    vector<primVars> stateN(imax*jmax*kmax);
    vector<double> dtN(imax*jmax*kmax);

    for ( rr = 0; rr < 4; rr ++ ){
      for ( kk = 0; kk < kmax; kk++ ){          //loop over all cells
	for ( jj = 0; jj < jmax; jj++ ){          
	  for ( ii = 0; ii < imax; ii++ ){          

	    loc = GetLoc1D(ii, jj, kk, imax, jmax);

	    //save state and local time step at time n
	    if (rr == 0){
	      stateN[loc] = (*this).State(loc);
	      dtN[loc] = (*this).Dt(loc);
	    }

	    (*this).RK4TimeAdvance(stateN[loc], eos, dtN[loc], loc, rr);

	    if (rr ==3){
	      for ( unsigned int ll = 0; ll < l2.size(); ll++ ){
		l2[ll] = l2[ll] + (*this).Residual(loc,ll) * (*this).Residual(loc,ll);

		if ( (*this).Residual(loc,ll) > linf[4] ){
		  linf[4] = (*this).Residual(loc,ll);
		  linf[3] = (double)ll+1;
		  linf[2] = (double)kk;
		  linf[1] = (double)jj;
		  linf[0] = (double)ii;
		  locMaxB = bb;
		}
	      }
	    }

	  }
	}
      }
      //for multistage RK4 method, calculate fluxes and residuals again
      if (rr < 3){ //no need to calculate fluxes after final RK interation
	(*this).CalcInvFluxI(eos, inputVars, bb);
	(*this).CalcInvFluxJ(eos, inputVars, bb);
	(*this).CalcInvFluxK(eos, inputVars, bb);
	(*this).CalcBlockTimeStep(inputVars, aRef);
      }


    }
  }
  else {
    cerr << "ERROR: Time integration scheme " << inputVars.TimeIntegration() << " is not recognized!" << endl;
  }



}


//member function to advance the state vector to time n+1 using explicit Euler method
void blockVars::ExplicitEulerTimeAdvance(const idealGas &eqnState, const int &loc ){

  colMatrix consVars = (*this).State(loc).ConsVars(eqnState);

  //calculate updated conserved variables
  consVars.SetData(0, consVars.Data(0) - (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc,0) );
  consVars.SetData(1, consVars.Data(1) - (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc,1) );
  consVars.SetData(2, consVars.Data(2) - (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc,2) );
  consVars.SetData(3, consVars.Data(3) - (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc,3) );
  consVars.SetData(4, consVars.Data(4) - (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc,4) );

  //calculate updated primative variables
  vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  primVars tempState (consVars.Data(0),
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  (*this).SetState(tempState, loc);

}

//member function to advance the state vector to time n+1 (for implicit methods)
void blockVars::ImplicitTimeAdvance(const colMatrix &du, const idealGas &eqnState, const int &loc ){

  colMatrix consVars = (*this).State(loc).ConsVars(eqnState);

  //calculate updated conserved variables
  for (int ii = 0; ii < du.Size(); ii++ ){
    consVars.SetData(ii, consVars.Data(ii) + du.Data(ii) );
  }

  // cout << "correction: " << endl;
  // cout << du << endl;

  //calculate updated primative variables
  vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  primVars tempState (consVars.Data(0),
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  //check for positivity
  if (tempState.Rho() < 0.0 || tempState.P() < 0.0){
    cerr << "ERROR: Density or pressure has become negative!" << endl;
    cerr << "Updated Primative variables:" << endl << tempState << endl;
    cerr << "Original Primative variables:" << endl << (*this).State(loc) << endl;
    cerr << "Correction:" << endl << du << endl;
    cerr << "Updated Conserved variables: " << consVars << endl;
    exit(1);
  }

  (*this).SetState(tempState, loc);

}


//member function to advance the state vector to time n+1 using 4th order Runge-Kutta method
void blockVars::RK4TimeAdvance( const primVars &currState, const idealGas &eqnState, const double &dt, const int &loc, const int &rk ){

  double alpha[4] = {0.25, 1.0/3.0, 0.5, 1.0};

  colMatrix consVars = currState.ConsVars(eqnState);

  //calculate updated conserved variables
  consVars.SetData(0, consVars.Data(0) - dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc,0) );
  consVars.SetData(1, consVars.Data(1) - dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc,1) );
  consVars.SetData(2, consVars.Data(2) - dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc,2) );
  consVars.SetData(3, consVars.Data(3) - dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc,3) );
  consVars.SetData(4, consVars.Data(4) - dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc,4) );

  //calculate updated primative variables
  vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  primVars tempState (consVars.Data(0),
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  (*this).SetState(tempState, loc);
}

void blockVars::ResetResidWS( ){

  colMatrix initial( (*this).Residual(0).Size() );
  initial.Zero();

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	//reset residual
	(*this).SetResidual( initial, loc ) ;
	
	//reset wave speed
	(*this).SetAvgWaveSpeed( 0.0, loc ) ;

      }
    }
  }


}


//function to calculate the flux jacobians on the i-faces
void blockVars::CalcInvFluxJacI(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag, matrixDiagonal &offLowIDiag, matrixDiagonal &offUpIDiag, const string &fluxJacType)const{

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int upperI = 0;
  int lowerI = 0;

  double maxWS = 0.0;
  double maxWSL = 0.0;
  double maxWSR = 0.0;

  primVars faceStateLower, faceStateUpper, ghostState;

  string bcName = "undefined";

  squareMatrix tempL(offLowIDiag.Data(0).Size());
  squareMatrix tempR(offLowIDiag.Data(0).Size());

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	tempL.Zero();
	tempR.Zero();

	//find out if at a block boundary
	if ( ii == 0  ){                             //at i lower boundary
	  upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "il");
	  //faceStateLower = (*this).State(upperI).GetGhostState( bcName, (*this).FAreaI(loc), "il", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateUpper = (*this).State( upperI ).FaceReconConst();

	  tempR = BoundaryFluxJacobian(bcName, (*this).FAreaI(loc), faceStateUpper, eqnState, inp, "il", fluxJacType, maxWS);

          // left flux jacobian is not needed at lower boundary (originally - maxWS)
          mainDiag.SetData(  upperI, mainDiag.Data(upperI)   + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

	  // if ( imax > 2){  //if only one cell thick, no off diagonals
	  //   lowDiag = GetMatrixDiagLowerFromMainI(upperI);
	  //   offLowIDiag.SetData(lowDiag, offLowIDiag.Data(lowDiag) + tempR * (*this).FAreaI(loc).Mag() );
	  // }

	  // if (upperI == 1){
	  //   cout << "I-flux-jacobian at cell 1, boundary " << bcName << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaI(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}
	else if ( ii == imax-1 ){  //at i upper boundary
	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  //faceStateUpper = (*this).State(lowerI).GetGhostState( bcName, (*this).FAreaI(loc), "iu", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateLower = (*this).State( lowerI ).FaceReconConst();

	  tempL = BoundaryFluxJacobian(bcName, (*this).FAreaI(loc), faceStateLower, eqnState, inp, "iu", fluxJacType, maxWS);

	  // right flux jacobian is not needed at upper boundary
          mainDiag.SetData(   lowerI, mainDiag.Data(lowerI)    + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

	  // if (imax > 2){ //if only one cell thick, no off diagonals
	  //   upDiag = GetMatrixDiagUpperFromMainI(lowerI);
	  //   offUpIDiag.SetData(upDiag, offUpIDiag.Data(upDiag) - tempL * (*this).FAreaI(loc).Mag() );
	  // }

	  // if (lowerI == 1){
	  //   cout << "I-flux-jacobian at cell 1, boundary " << bcName << ": " << endl;
	  //   cout << tempL * (*this).FAreaI(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}
	else{
	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  faceStateLower = (*this).State( lowerI ).FaceReconConst();
	  faceStateUpper = (*this).State( upperI ).FaceReconConst();

	  if ( fluxJacType == "approximateRoe" ){
	    ApproxRoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "exactRoe" ){
	    RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "laxFriedrichs" ){
	    LaxFriedrichsFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWSL, maxWSR, tempL, tempR);
	  }
	  else{
	    cerr << "ERROR: Inviscid flux jacobian type " << fluxJacType << " is not recognized!" << endl;
	    exit(1);
	  }

	  //left flux jacobian
          mainDiag.SetData(   lowerI, mainDiag.Data(lowerI)    + 0.5 * maxWSL * (*this).FAreaI(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainI(lowerI);
	  //upDiag = GetMatrixDiagUpperFromMainI(lowerI);

	  offLowIDiag.SetData(upperI, -1.0 * tempL * (*this).FAreaI(loc).Mag() );
	  // if (upDiag >= 0 && ii-1 > 0){
	  //   offUpIDiag.SetData(upDiag, offUpIDiag.Data(upDiag) - tempL * (*this).FAreaI(loc).Mag() );
	  // }

	  //right flux jacobian (originally - maxWSR)
          mainDiag.SetData(  upperI, mainDiag.Data(upperI)   + 0.5 * maxWSR * (*this).FAreaI(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainI(upperI);
	  //upDiag = GetMatrixDiagUpperFromMainI(upperI);

	  // if (lowDiag < offLowIDiag.Size() && ii+1 < imax-1 ){
	  //   offLowIDiag.SetData(lowDiag, offLowIDiag.Data(lowDiag) + tempR * (*this).FAreaI(loc).Mag() );
	  // }
	  offUpIDiag.SetData(lowerI, tempR * (*this).FAreaI(loc).Mag() );


	  // if (upperI == 1){
	  //   cout << "I-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaI(loc).Mag() << endl;
	  //   cout << "face state upperI:" << endl;
	  //   cout << faceStateUpper << endl;
	  //   cout << "face state lowerI:" << endl;
	  //   cout << faceStateLower << endl;
	  // }
	  // if (lowerI == 1){
	  //   cout << "I-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << tempL * (*this).FAreaI(loc).Mag() << endl;
	  //   cout << "face state upperI:" << endl;
	  //   cout << faceStateUpper << endl;
	  //   cout << "face state lowerI:" << endl;
	  //   cout << faceStateLower << endl;
	  // }


	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}

      }
    }
  }

}

//function to calculate the flux jacobians on the j-faces
void blockVars::CalcInvFluxJacJ(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag, matrixDiagonal &offLowJDiag, matrixDiagonal &offUpJDiag, const string &fluxJacType)const{

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int upperJ = 0;
  int lowerJ = 0;
  // int upDiag = 0;
  // int lowDiag = 0;

  double maxWS = 0.0;
  double maxWSL = 0.0;
  double maxWSR = 0.0;

  primVars faceStateLower, faceStateUpper, ghostState;

  string bcName = "undefined";

  squareMatrix tempL(offLowJDiag.Data(0).Size());
  squareMatrix tempR(offLowJDiag.Data(0).Size());

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	tempL.Zero();
	tempR.Zero();

	//find out if at a block boundary
	if ( jj == 0  ){                             //at j lower boundary
	  upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "jl");
	  faceStateLower = (*this).State(upperJ).GetGhostState( bcName, (*this).FAreaJ(loc), "jl", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();

	  tempR = BoundaryFluxJacobian(bcName, (*this).FAreaJ(loc), faceStateUpper, eqnState, inp, "jl", fluxJacType, maxWS);

          // left flux jacobian is not needed at lower boundary (originally - maxWS)
          mainDiag.SetData(  upperJ, mainDiag.Data(upperJ)   + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

	  // if (jmax > 2){  //if only one cell thick, no diagonals
	  //   lowDiag = GetMatrixDiagLowerFromMainJ(upperJ, imax);
	  //   offLowJDiag.SetData(lowDiag, offLowJDiag.Data(lowDiag) + tempR * (*this).FAreaJ(loc).Mag() );
	  // }

	  //cout << "lower diagonal storing at index (from boundary) " << lowDiag << endl;

	  // if (upperJ == 1){
	  //   cout << "J-flux-jacobian at cell 1, from face " << loc << " boundary " << bcName << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaJ(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedJ(maxWS, loc);

	}
	else if ( jj == jmax-1 ){  //at i upper boundary
	  lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  faceStateUpper = (*this).State(lowerJ).GetGhostState( bcName, (*this).FAreaJ(loc), "ju", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();

	  tempL = BoundaryFluxJacobian(bcName, (*this).FAreaJ(loc), faceStateLower, eqnState, inp, "ju", fluxJacType, maxWS);

	  // right flux jacobian is not needed at upper boundary
          mainDiag.SetData(   lowerJ, mainDiag.Data(lowerJ)    + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

	  // if (jmax > 2){ //if only one cell thick, no diagonals
	  //   upDiag = GetMatrixDiagUpperFromMainJ(lowerJ, imax);
	  //   offUpJDiag.SetData(upDiag, offUpJDiag.Data(upDiag) - tempL * (*this).FAreaJ(loc).Mag() );
	  // }

	  // if (lowerJ == 1){
	  //   cout << "J-flux-jacobian at cell 1, from face " << loc << " boundary " << bcName << ": " << endl;
	  //   cout << tempL * (*this).FAreaJ(loc).Mag() << endl;
	  // }


	  //cout << "upper diagonal storing at index (from boundary) " << upDiag << endl;

	  //(*this).SetMaxWaveSpeedJ(maxWS, loc);

	}
	else{
	  lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();
	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();

	  if ( fluxJacType == "approximateRoe" ){
	    ApproxRoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "exactRoe" ){
	    RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "laxFriedrichs" ){
	    LaxFriedrichsFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWSL, maxWSR, tempL, tempR);
	  }
	  else{
	    cerr << "ERROR: Inviscid flux jacobian type " << fluxJacType << " is not recognized!" << endl;
	    exit(1);
	  }

	  //left flux jacobian
          mainDiag.SetData(   lowerJ, mainDiag.Data(lowerJ)    + 0.5 * maxWSL * (*this).FAreaJ(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainJ(lowerJ, imax);
	  //upDiag = GetMatrixDiagUpperFromMainJ(lowerJ, imax);

	  offLowJDiag.SetData(upperJ, -1.0 * tempL * (*this).FAreaJ(loc).Mag() );
	  //cout << "lower diagonal storing at index " << lowDiag << endl;

	  // if (upDiag >= 0 && jj-1 > 0){
	  //   offUpJDiag.SetData(upDiag, offUpJDiag.Data(upDiag) - tempL * (*this).FAreaJ(loc).Mag() );
	  //   //cout << "upper diagonal storing at index " << upDiag << endl;
	  // }

	  //right flux jacobian (originally - maxWSR)
          mainDiag.SetData(  upperJ, mainDiag.Data(upperJ)   + 0.5 * maxWSR * (*this).FAreaJ(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainJ(upperJ, imax);
	  //upDiag = GetMatrixDiagUpperFromMainJ(upperJ, imax);

	  // if (lowDiag < offLowJDiag.Size() && jj+1 < jmax-1 ){
	  //   offLowJDiag.SetData(lowDiag, offLowJDiag.Data(lowDiag) + tempR * (*this).FAreaJ(loc).Mag() );
	  //   //cout << "lower diagonal storing at index " << lowDiag << endl;
	  // }

	  //cout << "upper diagonal storing at index " << upDiag << endl;
	  offUpJDiag.SetData(lowerJ, tempR * (*this).FAreaJ(loc).Mag() );


	  // if (upperJ == 1){
	  //   cout << "J-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaJ(loc).Mag() << endl;
	  // }
	  // if (lowerJ == 1){
	  //   cout << "J-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << tempL * (*this).FAreaJ(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedJ(maxWS, loc);

	}

      }
    }
  }


}

//function to calculate the flux jacobians on the k-faces
void blockVars::CalcInvFluxJacK(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag, matrixDiagonal &offLowKDiag, matrixDiagonal &offUpKDiag, const string &fluxJacType)const{

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK();

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int upperK = 0;
  int lowerK = 0;
  // int upDiag = 0;
  // int lowDiag = 0;

  double maxWS = 0.0;
  double maxWSL = 0.0;
  double maxWSR = 0.0;

  primVars faceStateLower, faceStateUpper, ghostState;

  string bcName = "undefined";

  squareMatrix tempL(offLowKDiag.Data(0).Size());
  squareMatrix tempR(offLowKDiag.Data(0).Size());

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	tempL.Zero();
	tempR.Zero();

	//find out if at a block boundary
	if ( kk == 0  ){                             //at k lower boundary
	  upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "kl");
	  faceStateLower = (*this).State(upperK).GetGhostState( bcName, (*this).FAreaK(loc), "kl", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateUpper = (*this).State( upperK ).FaceReconConst();

	  tempR = BoundaryFluxJacobian(bcName, (*this).FAreaK(loc), faceStateUpper, eqnState, inp, "kl", fluxJacType, maxWS);

          // left flux jacobian is not needed at lower boundary (originally - maxWS)
          mainDiag.SetData(  upperK, mainDiag.Data(upperK)   + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

	  // if (kmax > 2){ //if only one cell thick, no diagonals
	  //   lowDiag = GetMatrixDiagLowerFromMainK(upperK, imax, jmax);
	  //   offLowKDiag.SetData(lowDiag, offLowKDiag.Data(lowDiag) + tempR * (*this).FAreaK(loc).Mag() );
	  //   // cout << "lower diagonal storing at index (from boundary) " << lowDiag << endl;
	  // }


	  // if (upperK == 1){
	  //   cout << "K-flux-jacobian at cell 1, from face " << loc << " boundary " << bcName << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaK(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedK(maxWS, loc);

	}
	else if ( kk == kmax-1 ){  //at i upper boundary
	  lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  faceStateUpper = (*this).State(lowerK).GetGhostState( bcName, (*this).FAreaK(loc), "ku", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateLower = (*this).State( lowerK ).FaceReconConst();

	  tempL = BoundaryFluxJacobian(bcName, (*this).FAreaK(loc), faceStateLower, eqnState, inp, "ku", fluxJacType, maxWS);

	  // right flux jacobian is not needed at upper boundary
          mainDiag.SetData(   lowerK, mainDiag.Data(lowerK)    + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

	  // if (kmax > 2){ //if only one cell thick, no diagonals
	  //   upDiag = GetMatrixDiagUpperFromMainK(lowerK, imax, jmax);
	  //   offUpKDiag.SetData(upDiag, offUpKDiag.Data(upDiag) - tempL * (*this).FAreaK(loc).Mag() );
	  // }

	  // if (lowerK == 1){
	  //   cout << "K-flux-jacobian at cell 1, from face " << loc << " boundary " << bcName << ": " << endl;
	  //   cout << tempL * (*this).FAreaK(loc).Mag() << endl;
	  // }



	  // cout << "upper diagonal storing at index (from boundary) " << upDiag << endl;
	  //(*this).SetMaxWaveSpeedK(maxWS, loc);

	}
	else{
	  lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  faceStateLower = (*this).State( lowerK ).FaceReconConst();
	  faceStateUpper = (*this).State( upperK ).FaceReconConst();


	  if ( fluxJacType == "approximateRoe" ){
	    ApproxRoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "exactRoe" ){
	    RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS, tempL, tempR);
	  }
	  else if ( fluxJacType == "laxFriedrichs" ){
	    LaxFriedrichsFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWSL, maxWSR, tempL, tempR);
	  }
	  else{
	    cerr << "ERROR: Inviscid flux jacobian type " << fluxJacType << " is not recognized!" << endl;
	    exit(1);
	  }

	  //left flux jacobian
          mainDiag.SetData(   lowerK, mainDiag.Data(lowerK)    + 0.5 * maxWSL * (*this).FAreaK(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainK(lowerK, imax, jmax);
	  //upDiag = GetMatrixDiagUpperFromMainK(lowerK, imax, jmax);

	  offLowKDiag.SetData(upperK, -1.0 * tempL * (*this).FAreaK(loc).Mag() );
	  //cout << "lower diagonal storing at index " << lowDiag << endl;

	  // if (upDiag >= 0 && kk-1 > 0){
	  //   offUpKDiag.SetData(upDiag, offUpKDiag.Data(upDiag) - tempL * (*this).FAreaK(loc).Mag() );
	  //   //cout << "upper diagonal storing at index " << upDiag << endl;
	  // }

	  //right flux jacobian (originally - maxWSR)
          mainDiag.SetData(  upperK, mainDiag.Data(upperK)   + 0.5 * maxWSR * (*this).FAreaK(loc).Mag() );

	  //lowDiag = GetMatrixDiagLowerFromMainK(upperK, imax, jmax);
	  //upDiag = GetMatrixDiagUpperFromMainK(upperK, imax, jmax);

	  // if (lowDiag < offLowKDiag.Size() && kk+1 < kmax-1 ){
	  //   offLowKDiag.SetData(lowDiag, offLowKDiag.Data(lowDiag) + tempR * (*this).FAreaK(loc).Mag() );
	  //   //cout << "lower diagonal storing at index " << lowDiag << endl;
	  // }

	  //cout << "upper diagonal storing at index " << upDiag << endl;
	  offUpKDiag.SetData(lowerK, tempR * (*this).FAreaK(loc).Mag() );


	  // if (upperK == 1){
	  //   cout << "K-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << -1.0 * tempR * (*this).FAreaK(loc).Mag() << endl;
	  // }
	  // if (lowerK == 1){
	  //   cout << "K-flux-jacobian at cell 1 from face " << loc << ": " << endl;
	  //   cout << tempL * (*this).FAreaK(loc).Mag() << endl;
	  // }


	  //(*this).SetMaxWaveSpeedK(maxWS, loc);

	}

      }
    }
  }


}

//a member function to add the cell volume divided by the cell time step to the main diagonal of the implicit matrix
void blockVars::AddVolTime( colMatrix &mainDiag, const double &theta, const double &zeta) const {

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	//I = ( ( ((*this).Vol(loc) * (1.0 + zeta)) / (*this).Dt(loc) ) - theta ) * I ;
	// double I = (*this).Vol(loc)/(*this).Dt(loc) ;
	double I = ( (*this).Vol(loc) * (1.0 + zeta) ) / ( (*this).Dt(loc) * theta ) ;
	mainDiag.SetData(loc, I + mainDiag.Data(loc) );
      }
    }
  }

}

//a member function to add the cell volume divided by the cell time step to the time m - time n term
vector<colMatrix> blockVars::AddVolTime(const vector<colMatrix> &m, const vector<colMatrix> &n, const double &theta, const double &zeta) const {

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  vector<colMatrix> mMinusN = m;

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	double I = ( (*this).Vol(loc) * (1.0 + zeta) ) / ( (*this).Dt(loc) * theta ) ;
	mMinusN[loc] = I * (m[loc] - n[loc]);
      }
    }
  }
  return mMinusN;
}


//member function to calculate the delta n-1 term for the implicit bdf2 solver
void blockVars::DeltaNMinusOne(vector<colMatrix> &solDeltaNm1, const vector<colMatrix> &solTimeN, const idealGas &eqnState, const double &theta, const double &zeta){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	double coeff = ( (*this).Vol(loc) * zeta ) / ( (*this).Dt(loc) * theta ) ;
	solDeltaNm1[loc] = coeff * ( (*this).State(loc).ConsVars(eqnState) - solTimeN[loc] );
      }
    }
  }

}


//member function to return a copy of the conserved variables
vector<colMatrix> blockVars::GetCopyConsVars(const idealGas &eqnState) const {

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  vector<colMatrix> consVars(imax*jmax*kmax);

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	colMatrix temp = (*this).State(loc).ConsVars(eqnState);
	consVars[loc] = temp;
      }
    }
  }

  return consVars;

}

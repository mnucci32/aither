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
  vector<inviscidFlux> dummyFlux (length);              //dummy flux variable
  vector<vector3d<double> > vec1(length);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vec2(lenCellCen);             //dummy vector variable lenght of number of cells
  vector<double> scalar(lenCellCen);                      //dummy scalar variable lenght of number of cells
  vector<double> vec3(length);
  vector<double> singleResid(5,0.0);
  vector<vector<double> > dummyResid(lenCellCen, singleResid);

  state = dummyState;      
  invFluxI = dummyFlux;
  invFluxJ = dummyFlux;
  invFluxK = dummyFlux;

  vol = scalar;
  center = vec2;
  fAreaI = vec1;
  fAreaJ = vec1;
  fAreaK = vec1;
  fCenterI = vec1;
  fCenterJ = vec1;
  fCenterK = vec1;

  maxWaveSpeedI = vec3;
  maxWaveSpeedJ = vec3;
  maxWaveSpeedK = vec3;

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
  vector<inviscidFlux> dummyFlux (length);               //dummy flux variable
  vector<double> dummyTime (lenCellCen);                 //dummy time variable
  vector<double> singleResid(5,0.0);
  vector<vector<double> > dummyResid(lenCellCen, singleResid);


  int lengthI = numI * (numJ-1) * (numK-1);
  int lengthJ = (numI-1) * numJ * (numK-1);
  int lengthK = (numI-1) * (numJ-1) * numK;

  vector<double> vecI(lengthI,0.0);                    //dummy vectors for max wave speed
  vector<double> vecJ(lengthJ,0.0);
  vector<double> vecK(lengthK,0.0);

  state = dummyState;      
  invFluxI = dummyFlux;
  invFluxJ = dummyFlux;
  invFluxK = dummyFlux;

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  maxWaveSpeedI = vecI;
  maxWaveSpeedJ = vecJ;
  maxWaveSpeedK = vecK;

  dt = dummyTime;

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
  vector<inviscidFlux> dummyFlux (length);               //dummy flux variable
  vector<double> dummyTime (lenCellCen);                  //dummy time variable

  int lengthI = numI * (numJ-1) * (numK-1);
  int lengthJ = (numI-1) * numJ * (numK-1);
  int lengthK = (numI-1) * (numJ-1) * numK;

  vector<double> vecI(lengthI,0.0);                    //dummy vectors for max wave speed
  vector<double> vecJ(lengthJ,0.0);
  vector<double> vecK(lengthK,0.0);

  vector<double> singleResid(5,0.0);
  vector<vector<double> > dummyResid(lenCellCen, singleResid);

  state = initState;      
  invFluxI = dummyFlux;
  invFluxJ = dummyFlux;
  invFluxK = dummyFlux;

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  maxWaveSpeedI = vecI;
  maxWaveSpeedJ = vecJ;
  maxWaveSpeedK = vecK;

  dt = dummyTime;

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
  vector<inviscidFlux> dummyFlux (length);               //dummy flux variable
  vector<double> dummyTime (lenCellCen);                 //dummy time variable

  int lengthI = numI * (numJ-1) * (numK-1);
  int lengthJ = (numI-1) * numJ * (numK-1);
  int lengthK = (numI-1) * (numJ-1) * numK;

  vector<double> vecI(lengthI,0.0);                    //dummy vectors for max wave speed
  vector<double> vecJ(lengthJ,0.0);
  vector<double> vecK(lengthK,0.0);

  vector<double> singleResid(5,0.0);
  vector<vector<double> > dummyResid(lenCellCen, singleResid);

  state = initState;      
  invFluxI = dummyFlux;
  invFluxJ = dummyFlux;
  invFluxK = dummyFlux;

  vol = blk.Volume();
  center = blk.Centroid();
  fAreaI = blk.FaceAreaI();
  fAreaJ = blk.FaceAreaJ();
  fAreaK = blk.FaceAreaK();
  fCenterI = blk.FaceCenterI();
  fCenterJ = blk.FaceCenterJ();
  fCenterK = blk.FaceCenterK();

  maxWaveSpeedI = vecI;
  maxWaveSpeedJ = vecJ;
  maxWaveSpeedK = vecK;

  dt = dummyTime;

  residual = dummyResid;
}

//---------------------------------------------------------------------------------------------------------------//
//function declarations
//function to calculate the fluxes on the i-faces
void blockVars::CalcInvFluxI(const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

 const boundaryConditions bound = inp.BC()[bb];
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
	if ( ii == 0  ){                             //at i lower boundary
	  bcName = bound.GetBCName(ii, jj, kk, "il");

	  if (imax > 2 && kap != -2.0){ //if more than 2 faces thick, and second order, use linear extrapolation to get boundary state
	    up2faceU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2)) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "il", maxWS, up2faceU, upwindU );
	    //cout << "at lower boundary fAreaI is " << (*this).FAreaI()[loc]/(*this).FAreaI()[loc].Mag() << endl;
	  }
	  else{  //if not more than 2 faces thick, use cell adjacent to boundary
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "il", maxWS );
	  }

	  // cout << tempFlux << endl;

	  (*this).SetInvFluxI(tempFlux, loc);
	  (*this).SetMaxWaveSpeedI(maxWS, loc);
	}
	else if ( ii == imax-1 ){  //at i upper boundary
	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  if (imax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "iu", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "iu", maxWS );
	  }

	  (*this).SetInvFluxI(tempFlux, loc);
	  (*this).SetMaxWaveSpeedI(maxWS, loc);
	}
	else if ( ii == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow
                                                                  //for slipwall need to mirror values
	                                                          //ghost state should use boundary adjacent cell and boundary normal
	  bcName = bound.GetBCName(ii-1, jj, kk, "il");           //get bc at ii=0
	  ghostState = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaI( GetNeighborLowI(ii, jj, kk, imax, jmax, 1) ), "il", inp, eqnState );


	  up2faceL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	  upwindL = (*this).FCenterI(loc).Distance( (*this).FCenterI( GetNeighborLowI(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
	  													 (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	  upwindU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ) );

	  faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ),
	  													 (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);
	  (*this).SetInvFluxI(tempFlux, loc);
	  (*this).SetMaxWaveSpeedI(maxWS, loc);


	}
	else if ( ii == imax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii+1, jj, kk, "iu");          //get bc at ii = imax-1
	  ghostState = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaI( GetNeighborUpI(ii, jj, kk, imax, jmax, 1) ), "iu", inp, eqnState );

	  up2faceL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	  upwindL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
	  													 (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	  upwindU = (*this).FCenterI(loc).Distance( (*this).FCenterI( GetNeighborUpI(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ) );

	  faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
	  													 (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);
	  (*this).SetInvFluxI(tempFlux, loc);
	  (*this).SetMaxWaveSpeedI(maxWS, loc);


	}
	else{
	  //calculate 2 reconstructed face states for lower i face
	  if (kap == -2.0){  //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconConst();
	    faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ) );

	    faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
	    													   (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterI(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ) );

	    faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ),
	    													   (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  }
	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);
	  (*this).SetInvFluxI(tempFlux, loc);
	  (*this).SetMaxWaveSpeedI(maxWS, loc);
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

  const boundaryConditions bound = inp.BC()[bb];
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
	if ( jj == 0  ){                             //at j lower boundary
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  if (jmax > 2 && kap != -2.0){
	    up2faceU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "jl", maxWS, up2faceU, upwindU );
	  }
	  else {
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "jl", maxWS );
	  }

	  (*this).SetInvFluxJ(tempFlux, loc);
	  (*this).SetMaxWaveSpeedJ(maxWS, loc);
	}
	else if ( jj == jmax-1 ){  //at j upper boundary
	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  if (jmax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ju", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "ju", maxWS );
	  }

	  (*this).SetInvFluxJ(tempFlux, loc);
	  (*this).SetMaxWaveSpeedJ(maxWS, loc);
	}
	else if ( jj == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj-1, kk, "jl");           //get bc at jj=0
	  ghostState = (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborLowJ(ii, jj, kk, imax, jmax, 1) ), "jl", inp, eqnState );

	  up2faceL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	  upwindL = (*this).FCenterJ(loc).Distance( (*this).FCenterJ( GetNeighborLowJ(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
														 (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	  upwindU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ) );

	  faceStateUpper = (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
														 (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);
	  (*this).SetInvFluxJ(tempFlux, loc);
	  (*this).SetMaxWaveSpeedJ(maxWS, loc);

	}
	else if ( jj == jmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj+1, kk, "ju");          //get bc at jj=jmax-1
	  ghostState = (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborUpJ(ii, jj, kk, imax, jmax, 1) ), "ju", inp, eqnState );

	  up2faceL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	  upwindL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
														 (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	  upwindU = (*this).FCenterJ(loc).Distance( (*this).FCenterJ( GetNeighborUpJ(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ) );

	  faceStateUpper = (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
														 (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);
	  (*this).SetInvFluxJ(tempFlux, loc);
	  (*this).SetMaxWaveSpeedJ(maxWS, loc);
	}
	else{
	  //calculate 2 reconstructed face states for lower j face
	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).FaceReconConst();
	    faceStateUpper = (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ) );

	    faceStateLower = (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
														   (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterJ(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ) );

	    faceStateUpper = (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
														   (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );
	  }
	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);
	  (*this).SetInvFluxJ(tempFlux, loc);
	  (*this).SetMaxWaveSpeedJ(maxWS, loc);
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

  const boundaryConditions bound = inp.BC()[bb];
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
	if ( kk == 0  ){                             //at k lower boundary
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  if (kmax > 2 && kap != -2.0){
	    up2faceU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "kl", maxWS, up2faceU, upwindU );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "kl", maxWS );
	  }

	  (*this).SetInvFluxK(tempFlux, loc);
	  (*this).SetMaxWaveSpeedK(maxWS, loc);
	}
	else if ( kk == kmax-1 ){  //at k upper boundary

	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  if (kmax > 2 && kap != -2.0){
	    up2faceL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ku", maxWS, up2faceL, upwindL );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ), (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 1) ), eqnState, inp, "ku", maxWS );
	  }

	  (*this).SetInvFluxK(tempFlux, loc);
	  (*this).SetMaxWaveSpeedK(maxWS, loc);
	}
	else if ( kk == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk-1, "kl");           //get bc at kk=0
	  ghostState = (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaK( GetNeighborLowK(ii, jj, kk, imax, jmax, 1) ), "kl", inp, eqnState );

	  up2faceL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	  upwindL = (*this).FCenterK(loc).Distance( (*this).FCenterK( GetNeighborLowK(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
	  													 (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	  upwindU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );
	  centralU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ) );

	  faceStateUpper = (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
	  													 (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);
	  (*this).SetInvFluxK(tempFlux, loc);
	  (*this).SetMaxWaveSpeedK(maxWS, loc);

	}
	else if ( kk == kmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk+1, "ku");          //get bc at kk=kmax-1
	  ghostState = (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 1) ).GetGhostState( bcName, (*this).FAreaK( GetNeighborUpK(ii, jj, kk, imax, jmax, 1) ), "ku", inp, eqnState );


	  up2faceL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	  upwindL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );
	  centralL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ) );

	  faceStateLower = (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
	  													 (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	  up2faceU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	  upwindU = (*this).FCenterK(loc).Distance( (*this).FCenterK( GetNeighborUpK(ii, jj, kk, imax, jmax) ) );        //due to ghost cell set upwind distance equal to local cell length
	  centralU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ) );
	  
	  faceStateUpper = (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( ghostState,
	  													 (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);
	  (*this).SetInvFluxK(tempFlux, loc);
	  (*this).SetMaxWaveSpeedK(maxWS, loc);
	}
	else{
	  //calculate 2 reconstructed face states for lower k face
	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).FaceReconConst();
	    faceStateUpper = (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).FaceReconConst();
	  }
	  else{

	    up2faceL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	    upwindL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ) );
	    centralL = (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ) );

	    faceStateLower = (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
	    													   (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ),lstr, kap, inp.Limiter(), up2faceL, upwindL, centralL );

	    up2faceU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).FCenterK(loc) );
	    upwindU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ) );
	    centralU = (*this).Center( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).Distance( (*this).Center( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ) );


	    faceStateUpper = (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax) ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
	    													   (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax) ),rstr, kap, inp.Limiter(), up2faceU, upwindU, centralU );


	  }
	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);
	  (*this).SetInvFluxK(tempFlux, loc);
	  (*this).SetMaxWaveSpeedK(maxWS, loc);
	}


      }
    }
  }


}

//member function to calculate the local time step. (i,j,k) are cell indices
void blockVars::CalcCellDt( const int &i, const int &j, const int &k, const double &cfl){

  double wsIl, wsIu, wsJl, wsJu, wsKl, wsKu;
  double faIl, faIu, faJl, faJu, faKl, faKu;

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;

  //get cell volume
  int loc = GetLoc1D(i, j, k, imax, jmax);
  double vol = (*this).Vol(loc);

  //get cell face areas
  faIl = (*this).FAreaI( GetLowerFaceI(i, j, k, imax, jmax) ).Mag();
  faIu = (*this).FAreaI( GetUpperFaceI(i, j, k, imax, jmax) ).Mag();
  faJl = (*this).FAreaJ( GetLowerFaceJ(i, j, k, imax, jmax) ).Mag();
  faJu = (*this).FAreaJ( GetUpperFaceJ(i, j, k, imax, jmax) ).Mag();
  faKl = (*this).FAreaK( GetLowerFaceK(i, j, k, imax, jmax) ).Mag();
  faKu = (*this).FAreaK( GetUpperFaceK(i, j, k, imax, jmax) ).Mag();

  //get cell wave speeds
  wsIl = (*this).MaxWaveSpeedI( GetLowerFaceI(i, j, k, imax, jmax) );
  wsIu = (*this).MaxWaveSpeedI( GetUpperFaceI(i, j, k, imax, jmax) );
  wsJl = (*this).MaxWaveSpeedJ( GetLowerFaceJ(i, j, k, imax, jmax) );
  wsJu = (*this).MaxWaveSpeedJ( GetUpperFaceJ(i, j, k, imax, jmax) );
  wsKl = (*this).MaxWaveSpeedK( GetLowerFaceK(i, j, k, imax, jmax) );
  wsKu = (*this).MaxWaveSpeedK( GetUpperFaceK(i, j, k, imax, jmax) );


  double avgWS = 0.5 * (wsIl*faIl + wsIu*faIu + wsJl*faJl + wsJu*faJu + wsKl*faKl + wsKu*faKu);

  double dt = (cfl * vol / avgWS) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}

//this member function calculates the residual for an inviscid simulation
void blockVars::CalcCellResidual(const int &ii, const int &jj, const int &kk, const int &imax, const int &jmax){

  vector<double> resid(5);

  int loc = GetLoc1D(ii, jj, kk, imax, jmax);
  int iLow = GetLowerFaceI(ii, jj, kk, imax, jmax); 
  int iUp  = GetUpperFaceI(ii, jj, kk, imax, jmax);
  int jLow = GetLowerFaceJ(ii, jj, kk, imax, jmax);
  int jUp  = GetUpperFaceJ(ii, jj, kk, imax, jmax);
  int kLow = GetLowerFaceK(ii, jj, kk, imax, jmax);
  int kUp  = GetUpperFaceK(ii, jj, kk, imax, jmax);

  //Area vector points nominally from lower index to upper index, so the upper index fluxes must be multiplied by -1 so vector points into cell and for conservation
  resid[0] =        (*this).InvFluxI(iLow).RhoVel()  * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVel()  * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVel()  * (*this).FAreaK(kLow).Mag() 
	     -1.0 * (*this).InvFluxI(iUp).RhoVel()   * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVel()   * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVel()   * (*this).FAreaK(kUp).Mag() ;
  resid[1] =        (*this).InvFluxI(iLow).RhoVelU() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelU() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelU() * (*this).FAreaK(kLow).Mag() 
	     -1.0 * (*this).InvFluxI(iUp).RhoVelU()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelU()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelU()  * (*this).FAreaK(kUp).Mag() ;
  resid[2] =        (*this).InvFluxI(iLow).RhoVelV() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelV() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelV() * (*this).FAreaK(kLow).Mag() 
	     -1.0 * (*this).InvFluxI(iUp).RhoVelV()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelV()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelV()  * (*this).FAreaK(kUp).Mag() ;
  resid[3] =        (*this).InvFluxI(iLow).RhoVelW() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelW() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelW() * (*this).FAreaK(kLow).Mag() 
	     -1.0 * (*this).InvFluxI(iUp).RhoVelW()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelW()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelW()  * (*this).FAreaK(kUp).Mag() ;
  resid[4] =        (*this).InvFluxI(iLow).RhoVelH() * (*this).FAreaI(iLow).Mag() +     (*this).InvFluxJ(jLow).RhoVelH() * (*this).FAreaJ(jLow).Mag() +     (*this).InvFluxK(kLow).RhoVelH() * (*this).FAreaK(kLow).Mag() 
	     -1.0 * (*this).InvFluxI(iUp).RhoVelH()  * (*this).FAreaI(iUp).Mag() -1.0 * (*this).InvFluxJ(jUp).RhoVelH()  * (*this).FAreaJ(jUp).Mag() -1.0 * (*this).InvFluxK(kUp).RhoVelH()  * (*this).FAreaK(kUp).Mag() ;

  (*this).SetResidual(resid, loc);

}


void blockVars::CalcBlockResidDT( const input &inputVars, const double &aRef){

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

	(*this).CalcCellResidual(ii, jj, kk, imax, jmax);

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

void blockVars::UpdateBlock(const input &inputVars, const idealGas &eos, const double &aRef, const int &bb, vector<double> &l2, vector<double> &linf, int &locMaxB){

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;
  int kmax = (*this).NumK()-1;

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;

  if ( inputVars.TimeIntegration() == "explicitEuler" ){
    for ( kk = 0; kk < kmax; kk++ ){          //loop over all cells
      for ( jj = 0; jj < jmax; jj++ ){          
	for ( ii = 0; ii < imax; ii++ ){          

	  loc = GetLoc1D(ii, jj, kk, imax, jmax);

	  (*this).ExplicitEulerTimeAdvance(eos, loc);

	  for ( unsigned int ll = 0; ll < l2.size(); ll++ ){
	    l2[ll] = l2[ll] + (*this).Residual(loc)[ll] * (*this).Residual(loc)[ll];

	    if ( (*this).Residual(loc)[ll] > linf[4] ){
	      linf[4] = (*this).Residual(loc)[ll];
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
		l2[ll] = l2[ll] + (*this).Residual(loc)[ll] * (*this).Residual(loc)[ll];

		if ( (*this).Residual(loc)[ll] > linf[4] ){
		  linf[4] = (*this).Residual(loc)[ll];
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
	(*this).CalcBlockResidDT(inputVars, aRef);
      }


    }
  }
  else {
    cerr << "ERROR: Time integration scheme " << inputVars.TimeIntegration() << " is not recognized!" << endl;
  }



}


//member function to advance the state vector to time n+1 using explicit Euler method
void blockVars::ExplicitEulerTimeAdvance(const idealGas &eqnState, const int &loc ){

  vector<double> consVars = (*this).State(loc).ConsVars(eqnState);

  //calculate updated conserved variables
  consVars[0] = consVars[0] + (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc)[0] ;
  consVars[1] = consVars[1] + (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc)[1] ;
  consVars[2] = consVars[2] + (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc)[2] ;
  consVars[3] = consVars[3] + (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc)[3] ;
  consVars[4] = consVars[4] + (*this).Dt(loc) / (*this).Vol(loc) * (*this).Residual(loc)[4] ;

  //calculate updated primative variables
  vector3d<double> vel(consVars[1]/consVars[0], consVars[2]/consVars[0], consVars[3]/consVars[0]);

  primVars tempState (consVars[0],
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars[0], consVars[4]/consVars[0], vel.Mag() ) );

  (*this).SetState(tempState, loc);

}

//member function to advance the state vector to time n+1 using 4th order Runge-Kutta method
void blockVars::RK4TimeAdvance( const primVars &currState, const idealGas &eqnState, const double &dt, const int &loc, const int &rk ){

  double alpha[4] = {0.25, 1.0/3.0, 0.5, 1.0};

  vector<double> consVars = currState.ConsVars(eqnState);

  //calculate updated conserved variables
  consVars[0] = consVars[0] + dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc)[0] ;
  consVars[1] = consVars[1] + dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc)[1] ;
  consVars[2] = consVars[2] + dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc)[2] ;
  consVars[3] = consVars[3] + dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc)[3] ;
  consVars[4] = consVars[4] + dt / (*this).Vol(loc) * alpha[rk] * (*this).Residual(loc)[4] ;

  //calculate updated primative variables
  vector3d<double> vel(consVars[1]/consVars[0], consVars[2]/consVars[0], consVars[3]/consVars[0]);

  primVars tempState (consVars[0],
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars[0], consVars[4]/consVars[0], vel.Mag() ) );

  (*this).SetState(tempState, loc);
}

// void blockVars::TotalResidual( vector<double> &l2, vector<double> &linf, int &locMaxB, const int &bb){

//   int imax = (*this).NumI()-1;
//   int jmax = (*this).NumJ()-1;
//   int kmax = (*this).NumK()-1;

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;

//   unsigned int ll = 0;

//   for ( ii = 0; ii < imax; ii++ ){
//     for ( jj = 0; jj < jmax; jj++ ){
//       for ( kk = 0; kk < kmax; kk++ ){

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	for ( ll = 0; ll < l2.size(); ll++ ){
// 	  l2[ll] = l2[ll] + (*this).Residual(loc)[ll] * (*this).Residual(loc)[ll];

// 	  if ( (*this).Residual(loc)[ll] > linf[4] ){
// 	    linf[4] = (*this).Residual(loc)[ll];
// 	    linf[3] = (double)ll+1;
// 	    linf[2] = (double)kk;
// 	    linf[1] = (double)jj;
// 	    linf[0] = (double)ii;
// 	    locMaxB = bb;

// 	  }

// 	}
//       }
//     }
//   }

// }

//function to calculate the flux jacobians on the i-faces
void blockVars::CalcInvFluxJacI(const idealGas &eqnState, const input &inp, const int &bb, matrixDiagonal &mainDiag, matrixDiagonal &offLowIDiag, matrixDiagonal &offUpIDiag)const{

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC()[bb];

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int upperI = 0;
  int lowerI = 0;

  double maxWS = 0.0;

  primVars faceStateLower, faceStateUpper, ghostState;

  string bcName = "undefined";

  squareMatrix tempL(mainDiag.Data(0).Size());
  squareMatrix tempR(mainDiag.Data(0).Size());

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
	  faceStateLower = (*this).State(upperI).GetGhostState( bcName, (*this).FAreaI(loc), "il", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconConst();

	  RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS, tempL, tempR);

          // left flux jacobian is not needed at lower boundary
	  offLowIDiag.SetData(upperI, offLowIDiag.Data(upperI) + tempR * (*this).FAreaI(loc).Mag() );

          mainDiag.SetData(   upperI, mainDiag.Data(upperI)    - tempR * (*this).FAreaI(loc).Mag() );

	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}
	else if ( ii == imax-1 ){  //at i upper boundary
	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  faceStateUpper = (*this).State(lowerI).GetGhostState( bcName, (*this).FAreaI(loc), "iu", inp, eqnState).FaceReconConst(); //ghost state

	  faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconConst();

	  RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS, tempL, tempR);

	  // right flux jacobian is not needed at upper boundary
          mainDiag.SetData(  lowerI, mainDiag.Data(lowerI)   + tempL * (*this).FAreaI(loc).Mag() );

	  offUpIDiag.SetData(lowerI, offUpIDiag.Data(lowerI) - tempL * (*this).FAreaI(loc).Mag() );

	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}
	else{
	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  faceStateLower = (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax) ).FaceReconConst();
	  faceStateUpper = (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax) ).FaceReconConst();

	  RoeFluxJacobian(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS, tempL, tempR);

          mainDiag.SetData(   lowerI, mainDiag.Data(lowerI)    + tempL * (*this).FAreaI(loc).Mag() );
	  offLowIDiag.SetData(lowerI, offLowIDiag.Data(lowerI) + tempR * (*this).FAreaI(loc).Mag() );

          mainDiag.SetData(  upperI, mainDiag.Data(upperI)   - tempR * (*this).FAreaI(loc).Mag() );
	  offUpIDiag.SetData(lowerI, offUpIDiag.Data(lowerI) - tempL * (*this).FAreaI(loc).Mag() );

	  //(*this).SetMaxWaveSpeedI(maxWS, loc);

	}

      }
    }
  }


}


//a member function to print out the matrix structure for the block
void blockVars::PrintMatrixStructure(){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  int numCells = imax * jmax * kmax;

  for ( int rr = 0; rr < numCells; rr++ ){
    for ( int cc = 0; cc < numCells; cc++ ){

      if( cc == (rr-1) || cc == rr || cc == (rr+1) || cc == (rr-imax) || cc == (rr+imax) || cc == (rr-imax*jmax) || cc == (rr+imax*jmax) ){
	cout << "1, ";
      }
      else{
	cout << "0, ";
      }

      if (cc == numCells - 1){
	cout << endl;
      }


    }
  }

}

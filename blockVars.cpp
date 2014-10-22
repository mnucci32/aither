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
using std::min;

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

  double maxWS = 0.0;

  double upwindL, upwind2L, downwindL, upwindU, upwind2U, downwindU;

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
	  int upFaceI = GetNeighborUpI(ii, jj, kk, imax, jmax);

	  if (imax > 2 && kap != -2.0){ //if more than 2 faces thick, and second order, use linear extrapolation to get boundary state

	    int upFace2I = GetNeighborUpI(ii, jj, kk, imax, jmax, 2);

	    upwind2U =  (*this).FCenterI( upFaceI ).Distance( (*this).FCenterI( upFace2I ) );
	    upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( upperI ), (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ), 
				     eqnState, inp, "il", maxWS, upwindU, upwind2U );

	  }
	  else{  //if not more than 2 faces thick, use cell adjacent to boundary
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( upperI ), (*this).State( upperI ), eqnState, inp, "il", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + maxWS, upperI);
	}
	else if ( ii == imax-1 ){  //at i upper boundary -------------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  if (imax > 2 && kap != -2.0){

	    int lowFaceI = GetNeighborLowI(ii, jj, kk, imax, jmax);
	    int lowFace2I = GetNeighborLowI(ii, jj, kk, imax, jmax, 2);

	    upwind2L =  (*this).FCenterI( lowFaceI ).Distance( (*this).FCenterI( lowFace2I ) );
	    upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( lowerI ), (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ), 
				     eqnState, inp, "iu", maxWS, upwindL, upwind2L );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaI(loc), (*this).State( lowerI ), (*this).State( lowerI ), eqnState, inp, "iu", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);

	  //no wave speed calculation for upper faces

	}
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( ii == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow
                                                                  //for slipwall need to mirror values
	                                                          //ghost state should use boundary adjacent cell and boundary normal
	  bcName = bound.GetBCName(ii-1, jj, kk, "il");           //get bc at ii=0

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);
	  int upFaceI = GetNeighborUpI(ii, jj, kk, imax, jmax);
	  int upFace2I = GetNeighborUpI(ii, jj, kk, imax, jmax, 2);
	  int lowFaceI = GetNeighborLowI(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( lowerI ).GetGhostState( bcName, (*this).FAreaI( GetNeighborLowI(ii, jj, kk, imax, jmax, 1) ), "il", inp, eqnState);

	  upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );
	  upwind2L =  upwindL; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindL = (*this).FCenterI( loc      ).Distance( (*this).FCenterI( upFaceI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( ghostState, (*this).State( upperI ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwind2U =  (*this).FCenterI( upFaceI ).Distance( (*this).FCenterI( upFace2I ) );
	  upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );
	  downwindU = (*this).FCenterI( loc     ).Distance( (*this).FCenterI( lowFaceI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ), (*this).State( lowerI ),
								   "right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + maxWS, upperI);

	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( ii == imax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii+1, jj, kk, "iu");          //get bc at ii = imax-1

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);
	  int upFaceI = GetNeighborUpI(ii, jj, kk, imax, jmax);
	  int lowFaceI = GetNeighborLowI(ii, jj, kk, imax, jmax);
	  int lowFace2I = GetNeighborLowI(ii, jj, kk, imax, jmax, 2);

	  ghostState = (*this).State( upperI ).GetGhostState( bcName, (*this).FAreaI( GetNeighborUpI(ii, jj, kk, imax, jmax, 1) ), "iu", inp, eqnState );

	  upwind2L =  (*this).FCenterI( lowFaceI ).Distance( (*this).FCenterI( lowFace2I ) );
	  upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );
	  downwindL = (*this).FCenterI( loc      ).Distance( (*this).FCenterI( upFaceI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( upperI ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );
	  upwind2U =  upwindU; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindU = (*this).FCenterI( loc     ).Distance( (*this).FCenterI( lowFaceI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( ghostState, (*this).State( lowerI ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );
	  
	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + maxWS, upperI);

	}
	else{  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower i face

	  int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  int upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);
	  int upFaceI = GetNeighborUpI(ii, jj, kk, imax, jmax);
	  int upFace2I = GetNeighborUpI(ii, jj, kk, imax, jmax, 2);
	  int lowFaceI = GetNeighborLowI(ii, jj, kk, imax, jmax);
	  int lowFace2I = GetNeighborLowI(ii, jj, kk, imax, jmax, 2);

	  if (kap == -2.0){  //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerI ).FaceReconConst();
	    faceStateUpper = (*this).State( upperI ).FaceReconConst();
	  }
	  else{

	    upwind2L =  (*this).FCenterI( lowFaceI ).Distance( (*this).FCenterI( lowFace2I ) );
	    upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );
	    downwindL = (*this).FCenterI( loc      ).Distance( (*this).FCenterI( upFaceI ) );

	    faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerI(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperI ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	    upwind2U =  (*this).FCenterI( upFaceI ).Distance( (*this).FCenterI( upFace2I ) );
	    upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );
	    downwindU = (*this).FCenterI( loc     ).Distance( (*this).FCenterI( lowFaceI ) );

	    faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperI(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerI ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerI);
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperI);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperI) + maxWS, upperI);

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

  double maxWS = 0.0;

  double upwindL, upwind2L, downwindL, upwindU, upwind2U, downwindU;

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
	  int upFaceJ = GetNeighborUpJ(ii, jj, kk, imax, jmax);

	  if (jmax > 2 && kap != -2.0){

	    int upFace2J = GetNeighborUpJ(ii, jj, kk, imax, jmax, 2);

	    upwind2U =  (*this).FCenterJ( upFaceJ ).Distance( (*this).FCenterJ( upFace2J ) );
	    upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( upperJ ), 
				     (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "jl", maxWS, upwindU, upwind2U );

	  }
	  else {
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( upperJ ), (*this).State( upperJ ), eqnState, inp, "jl", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + maxWS, upperJ);

	}
	else if ( jj == jmax-1 ){  //at j upper boundary ---------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

	  if (jmax > 2 && kap != -2.0){

	    int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imax, jmax);
	    int lowFace2J = GetNeighborLowJ(ii, jj, kk, imax, jmax, 2);

	    upwind2L =  (*this).FCenterJ( lowFaceJ ).Distance( (*this).FCenterJ( lowFace2J ) );
	    upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( lowerJ ), 
				     (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ju", maxWS, upwindL, upwind2L );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaJ(loc), (*this).State( lowerJ ), (*this).State( lowerJ ), eqnState, inp, "ju", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);

	  //no wave speed calculation for upper faces

	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( jj == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj-1, kk, "jl");           //get bc at jj=0

	  int upFace2J = GetNeighborUpJ(ii, jj, kk, imax, jmax, 2);
	  int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imax, jmax);
	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);
	  int upFaceJ = GetNeighborUpJ(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( lowerJ ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborLowJ(ii, jj, kk, imax, jmax, 1) ), "jl", inp, eqnState );

	  upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );
	  upwind2L =  upwindL; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindL = (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( upFaceJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( ghostState, (*this).State( upperJ ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwind2U =  (*this).FCenterJ( upFaceJ ).Distance( (*this).FCenterJ( upFace2J ) );
	  upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );
	  downwindU = (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( lowFaceJ ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( lowerJ ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + maxWS, upperJ);

	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( jj == jmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj+1, kk, "ju");          //get bc at jj=jmax-1

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);
	  int upFaceJ = GetNeighborUpJ(ii, jj, kk, imax, jmax);
	  int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imax, jmax);
	  int lowFace2J = GetNeighborLowJ(ii, jj, kk, imax, jmax, 2);

	  ghostState = (*this).State( upperJ ).GetGhostState( bcName, (*this).FAreaJ( GetNeighborUpJ(ii, jj, kk, imax, jmax, 1) ), "ju", inp, eqnState );

	  upwind2L =  (*this).FCenterJ( lowFaceJ ).Distance( (*this).FCenterJ( lowFace2J ) );
	  upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );
	  downwindL = (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( upFaceJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( upperJ ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );
	  upwind2U =  upwindU; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindU = (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( lowFaceJ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( ghostState, (*this).State( lowerJ ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + maxWS, upperJ);

	}
	else{ // -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower j face

	  int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);
	  int upFaceJ = GetNeighborUpJ(ii, jj, kk, imax, jmax);
	  int upFace2J = GetNeighborUpJ(ii, jj, kk, imax, jmax, 2);
	  int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imax, jmax);
	  int lowFace2J = GetNeighborLowJ(ii, jj, kk, imax, jmax, 2);

	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerJ ).FaceReconConst();
	    faceStateUpper = (*this).State( upperJ ).FaceReconConst();
	  }
	  else{

	    upwind2L =  (*this).FCenterJ( lowFaceJ ).Distance( (*this).FCenterJ( lowFace2J ) );
	    upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );
	    downwindL = (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( upFaceJ ) );

	    faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperJ ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	    upwind2U =  (*this).FCenterJ( upFaceJ ).Distance( (*this).FCenterJ( upFace2J ) );
	    upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );
	    downwindU = (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( lowFaceJ ) );

	    faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerJ ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );
	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJ);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJ);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJ) + maxWS, upperJ);

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

  double maxWS = 0.0;

  double upwindL, upwind2L, downwindL, upwindU, upwind2U, downwindU;

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
	  int upFaceK = GetNeighborUpK(ii, jj, kk, imax, jmax);

	  if (kmax > 2 && kap != -2.0){

	    int upFace2K = GetNeighborUpK(ii, jj, kk, imax, jmax, 2);

	    upwind2U =  (*this).FCenterK( upFaceK ).Distance( (*this).FCenterK( upFace2K ) );
	    upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( upperK ), 
				     (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "kl", maxWS, upwindU, upwind2U );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( upperK ), (*this).State( upperK ), eqnState, inp, "kl", maxWS );
	  }

	  //at lower boundary normal points into cell, so need to subtract from residual
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + maxWS, upperK);

	}
	else if ( kk == kmax-1 ){  //at k upper boundary --------------------------------------------------------------------------------------------------------------------------------------
	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

	  if (kmax > 2 && kap != -2.0){
	    int lowFaceK = GetNeighborLowK(ii, jj, kk, imax, jmax);
	    int lowFace2K = GetNeighborLowK(ii, jj, kk, imax, jmax, 2);

	    upwind2L =  (*this).FCenterK( lowFaceK ).Distance( (*this).FCenterK( lowFace2K ) );
	    upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );

	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( lowerK ), 
				     (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ), eqnState, inp, "ku", maxWS, upwindL, upwind2L );
	  }
	  else{
	    tempFlux = BoundaryFlux( bcName, (*this).FAreaK(loc), (*this).State( lowerK ), (*this).State( lowerK ), eqnState, inp, "ku", maxWS );
	  }

	  //at upper boundary normal points out of cell, so need to add to residual
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);

	  //no wave speed calculation for upper faces

	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( kk == 1 && kap != -2.0){                        //lower face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk-1, "kl");           //get bc at kk=0

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);
	  int upFaceK = GetNeighborUpK(ii, jj, kk, imax, jmax);
	  int upFace2K = GetNeighborUpK(ii, jj, kk, imax, jmax, 2);
	  int lowFaceK = GetNeighborLowK(ii, jj, kk, imax, jmax);

	  ghostState = (*this).State( upperK ).GetGhostState( bcName, (*this).FAreaK( GetNeighborLowK(ii, jj, kk, imax, jmax, 1) ), "kl", inp, eqnState );

	  upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );
	  upwind2L =  upwindL; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindL = (*this).FCenterK( loc      ).Distance( (*this).FCenterK( upFaceK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( ghostState, (*this).State( upperK ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwind2U =  (*this).FCenterK( upFaceK ).Distance( (*this).FCenterK( upFace2K ) );
	  upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );
	  downwindU = (*this).FCenterK( loc     ).Distance( (*this).FCenterK( lowFaceK ) );

	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( lowerK ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + maxWS, upperK);

	}
	// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( kk == kmax-2 && kap != -2.0) {                 //upper face state reconstruction needs 1 ghost cell; set ghost cell equal to cell on boundary - works for inflow, outflow, slipwall
	  bcName = bound.GetBCName(ii, jj, kk+1, "ku");          //get bc at kk=kmax-1

	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);
	  int upFaceK = GetNeighborUpK(ii, jj, kk, imax, jmax);
	  int lowFaceK = GetNeighborLowK(ii, jj, kk, imax, jmax);
	  int lowFace2K = GetNeighborLowK(ii, jj, kk, imax, jmax, 2);

	  ghostState = (*this).State( upperK ).GetGhostState( bcName, (*this).FAreaK( GetNeighborUpK(ii, jj, kk, imax, jmax, 1) ), "ku", inp, eqnState );

	  upwind2L =  (*this).FCenterK( lowFaceK ).Distance( (*this).FCenterK( lowFace2K ) );
	  upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );
	  downwindL = (*this).FCenterK( loc      ).Distance( (*this).FCenterK( upFaceK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
								   (*this).State( upperK ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	  upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );
	  upwind2U =  upwindU; //due to ghost cell set upwind2 distance equal to upwind distance
	  downwindU = (*this).FCenterK( loc     ).Distance( (*this).FCenterK( lowFaceK) );

	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( ghostState, (*this).State( lowerK ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + maxWS, upperK);

	}
	else{ // -------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  //calculate 2 reconstructed face states for lower k face
	  int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  int upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);
	  int upFaceK = GetNeighborUpK(ii, jj, kk, imax, jmax);
	  int upFace2K = GetNeighborUpK(ii, jj, kk, imax, jmax, 2);
	  int lowFaceK = GetNeighborLowK(ii, jj, kk, imax, jmax);
	  int lowFace2K = GetNeighborLowK(ii, jj, kk, imax, jmax, 2);

	  if ( kap == -2.0 ){                         //if value is still default, use constant reconstruction
	    faceStateLower = (*this).State( lowerK ).FaceReconConst();
	    faceStateUpper = (*this).State( upperK ).FaceReconConst();
	  }
	  else{

	    upwind2L =  (*this).FCenterK( lowFaceK ).Distance( (*this).FCenterK( lowFace2K ) );
	    upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );
	    downwindL = (*this).FCenterK( loc      ).Distance( (*this).FCenterK( upFaceK ) );

	    faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( GetCellFromFaceLowerK(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( upperK ),"left", kap, inp.Limiter(), upwindL, upwind2L, downwindL );

	    upwind2U =  (*this).FCenterK( upFaceK ).Distance( (*this).FCenterK( upFace2K ) );
	    upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );
	    downwindU = (*this).FCenterK( loc     ).Distance( (*this).FCenterK( lowFaceK ) );

	    faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( GetCellFromFaceUpperK(ii, jj, kk, imax, jmax, 2) ),
								     (*this).State( lowerK ),"right", kap, inp.Limiter(), upwindU, upwind2U, downwindU );

	  }

	  tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerK);
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperK);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperK) + maxWS, upperK);

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

  double dt = cfl * ((*this).Vol(loc) / (*this).AvgWaveSpeed(loc)) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}



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

  // colMatrix consVars = (*this).State(loc).ConsVars(eqnState);

  // //calculate updated conserved variables
  // for (int ii = 0; ii < du.Size(); ii++ ){
  //   consVars.SetData(ii, consVars.Data(ii) + du.Data(ii) );
  // }

  // //calculate updated primative variables
  // vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  // primVars tempState (consVars.Data(0),
  // 		      vel.X(),
  // 		      vel.Y(),
  // 		      vel.Z(),
  // 		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  primVars tempState = (*this).State(loc).UpdateWithConsVars(eqnState, du);

  //check for positivity
  if (tempState.Rho() < 0.0 || tempState.P() < 0.0){
    cerr << "ERROR: Density or pressure has become negative!" << endl;
    cerr << "Updated Primative variables:" << endl << tempState << endl;
    cerr << "Original Primative variables:" << endl << (*this).State(loc) << endl;
    exit(0);
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


// //function to calculate the flux jacobians on the i-faces
// void blockVars::CalcInvFluxJacI(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag)const{

//   int imax = (*this).NumI();
//   int jmax = (*this).NumJ() - 1;
//   int kmax = (*this).NumK() - 1;

//   const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int upperI = 0;
//   int lowerI = 0;

//   double maxWS = 0.0;

//   primVars faceStateLower, faceStateUpper, ghostState;

//   string bcName = "undefined";

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	//find out if at a block boundary
// 	if ( ii == 0  ){                             //at i lower boundary
// 	  upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "il");
// 	  faceStateUpper = (*this).State( upperI ).FaceReconConst();

// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaI(loc), faceStateUpper, eqnState, "il", inp);

//           mainDiag.SetData(  upperI, mainDiag.Data(upperI)   + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

// 	}
// 	else if ( ii == imax-1 ){  //at i upper boundary
// 	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "iu");
// 	  faceStateLower = (*this).State( lowerI ).FaceReconConst();

// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaI(loc), faceStateLower, eqnState, "iu", inp);

//           mainDiag.SetData(   lowerI, mainDiag.Data(lowerI)    + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

// 	}
// 	else{
// 	  lowerI = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
// 	  upperI = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

// 	  faceStateLower = (*this).State( lowerI ).FaceReconConst();
// 	  faceStateUpper = (*this).State( upperI ).FaceReconConst();

// 	  maxWS = ConvSpecRad((*this).FAreaI(loc), faceStateLower, faceStateUpper, eqnState);

// 	  //left flux jacobian
//           mainDiag.SetData(   lowerI, mainDiag.Data(lowerI)    + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

// 	  //offLowIDiag.SetData(upperI, -1.0 * tempL * (*this).FAreaI(loc).Mag() );

// 	  //right flux jacobian (originally - maxWSR)
//           mainDiag.SetData(  upperI, mainDiag.Data(upperI)   + 0.5 * maxWS * (*this).FAreaI(loc).Mag() );

// 	  //offUpIDiag.SetData(lowerI, tempR * (*this).FAreaI(loc).Mag() );

// 	}

//       }
//     }
//   }

// }

// //function to calculate the flux jacobians on the j-faces
// void blockVars::CalcInvFluxJacJ(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag)const{

//   int imax = (*this).NumI() - 1;
//   int jmax = (*this).NumJ();
//   int kmax = (*this).NumK() - 1;

//   const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int upperJ = 0;
//   int lowerJ = 0;

//   double maxWS = 0.0;

//   primVars faceStateLower, faceStateUpper, ghostState;

//   string bcName = "undefined";

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	//find out if at a block boundary
// 	if ( jj == 0  ){                             //at j lower boundary
// 	  upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "jl");
// 	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();

// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaJ(loc), faceStateUpper, eqnState, "jl", inp);

//           mainDiag.SetData(  upperJ, mainDiag.Data(upperJ)   + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

// 	}
// 	else if ( jj == jmax-1 ){  //at i upper boundary
// 	  lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "ju");
// 	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();

// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaJ(loc), faceStateLower, eqnState, "ju", inp);

//           mainDiag.SetData(   lowerJ, mainDiag.Data(lowerJ)    + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

// 	}
// 	else{
// 	  lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
// 	  upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

// 	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();
// 	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();

// 	  maxWS = ConvSpecRad((*this).FAreaJ(loc), faceStateLower, faceStateUpper, eqnState);

// 	  //left flux jacobian
//           mainDiag.SetData(   lowerJ, mainDiag.Data(lowerJ)    + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

// 	  //offLowJDiag.SetData(upperJ, -1.0 * tempL * (*this).FAreaJ(loc).Mag() );

// 	  //right flux jacobian (originally - maxWSR)
//           mainDiag.SetData(  upperJ, mainDiag.Data(upperJ)   + 0.5 * maxWS * (*this).FAreaJ(loc).Mag() );

// 	  //offUpJDiag.SetData(lowerJ, tempR * (*this).FAreaJ(loc).Mag() );

// 	}

//       }
//     }
//   }


// }

// //function to calculate the flux jacobians on the k-faces
// void blockVars::CalcInvFluxJacK(const idealGas &eqnState, const input &inp, const int &bb, colMatrix &mainDiag)const{

//   int imax = (*this).NumI() - 1;
//   int jmax = (*this).NumJ() - 1;
//   int kmax = (*this).NumK();

//   const boundaryConditions bound = inp.BC(bb);

//   int ii = 0;
//   int jj = 0;
//   int kk = 0;
//   int loc = 0;
//   int upperK = 0;
//   int lowerK = 0;

//   double maxWS = 0.0;

//   primVars faceStateLower, faceStateUpper, ghostState;

//   string bcName = "undefined";

//   for ( kk = 0; kk < kmax; kk++){   
//     for ( jj = 0; jj < jmax; jj++){    
//       for ( ii = 0; ii < imax; ii++){      

// 	loc = GetLoc1D(ii, jj, kk, imax, jmax);

// 	//find out if at a block boundary
// 	if ( kk == 0  ){                             //at k lower boundary
// 	  upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "kl");
// 	  faceStateUpper = (*this).State( upperK ).FaceReconConst();

// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaK(loc), faceStateUpper, eqnState, "kl", inp);

//           mainDiag.SetData(  upperK, mainDiag.Data(upperK)   + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

// 	}
// 	else if ( kk == kmax-1 ){  //at i upper boundary
// 	  lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

// 	  bcName = bound.GetBCName(ii, jj, kk, "ku");

// 	  faceStateLower = (*this).State( lowerK ).FaceReconConst();
// 	  maxWS = BoundaryInvSpecRad(bcName, (*this).FAreaK(loc), faceStateLower, eqnState, "ku", inp);

//           mainDiag.SetData(   lowerK, mainDiag.Data(lowerK)    + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

// 	}
// 	else{
// 	  lowerK = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
// 	  upperK = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

// 	  faceStateLower = (*this).State( lowerK ).FaceReconConst();
// 	  faceStateUpper = (*this).State( upperK ).FaceReconConst();

// 	  maxWS = ConvSpecRad((*this).FAreaK(loc), faceStateLower, faceStateUpper, eqnState);

// 	  //left flux jacobian
//           mainDiag.SetData(   lowerK, mainDiag.Data(lowerK)    + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

// 	  //offLowKDiag.SetData(upperK, -1.0 * tempL * (*this).FAreaK(loc).Mag() );

// 	  //right flux jacobian (originally - maxWSR)
//           mainDiag.SetData(  upperK, mainDiag.Data(upperK)   + 0.5 * maxWS * (*this).FAreaK(loc).Mag() );

// 	  //offUpKDiag.SetData(lowerK, tempR * (*this).FAreaK(loc).Mag() );

// 	}

//       }
//     }
//   }


// }

//a member function to add the cell volume divided by the cell time step to the main diagonal of the implicit matrix
void blockVars::AddVolTime( colMatrix &mainDiag, const double &theta, const double &zeta, const double &dualCFL) const {

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;
  int loc = 0;

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){
	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	double I = ( (*this).Vol(loc) * (1.0 + zeta) ) / ( (*this).Dt(loc) * theta ) ;
	if (dualCFL > 0.0 ) { //use dual time stepping
	  double tau = (*this).AvgWaveSpeed(loc) / dualCFL; // equal to volume / tau
	  I = tau + I;
	}
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


//function to perform symmetric Gauss-Seidel relaxation to solver Ax=b
//when relax = 1.0, symmetric Gauss-Seidel is achieved. Values >1 result in symmetric successive over relaxation (SSOR)
//Values <1 result in under relaxation
double blockVars::LUSGS( const vector<vector3d<int> > &reorder, vector<colMatrix> &x, const vector<colMatrix> &solTimeMmN, const vector<colMatrix> &solDeltaNm1, const idealGas &eqnState, const input &inp, const sutherland &suth)const{

  //Aii --> block matrix of the main diagonal
  //x   --> block vector of correction
  //sweeps --> number of symmetric sweeps to perform
  //relax  --> relaxation parameter >1 is overrelaxation, <1 is underrelaxation


  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;

  //initialize x to 0
  for (unsigned int ll = 0; ll < x.size(); ll++ ){
    x[ll].Zero();
  }

  //invert main diagonal
  double AiiInv = 0.0;

  colMatrix l2Resid(x[0].Size());
  l2Resid.Zero();

  double thetaInv = 1.0 / inp.Theta();

  squareMatrix I(x[0].Size());
  I.Identity();

  colMatrix initial(x[0].Size());
  initial.Zero();

  for ( int kk = 0; kk < inp.MatrixSweeps(); kk++ ){

    vector<colMatrix> U(x.size(),initial);
    vector<colMatrix> L(x.size(),initial);

    vector<colMatrix> Uold(x.size(),initial);

    //forward sweep
    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int ilFace = GetLowerFaceI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int jlFace = GetLowerFaceJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int klFace = GetLowerFaceK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int ilFace2 = GetLowerFaceI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);
      int jlFace2 = GetLowerFaceJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);
      int klFace2 = GetLowerFaceK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);

      int il = GetNeighborLowI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int jl = GetNeighborLowJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int kl = GetNeighborLowK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      //calculate offdiagonal * update on the fly
      //only need values at offdiagonal cell 

      if ( il >=0 && il < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaI(ilFace), (*this).State(il), (*this).State(loc), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaI(ilFace2), (*this).FAreaI(ilFace), (*this).State(il).UpdateWithConsVars(eqnState, x[il]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  // double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(il), eqnState, suth, (*this).Center(il), (*this).Center(loc), (*this).FAreaI(ilFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(ilFace2), (*this).FAreaI(ilFace), (*this).State(il).UpdateWithConsVars(eqnState, x[il]), eqnState, suth, (*this).Vol(il) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(il), eqnState, (*this).FAreaI(ilFace), x[il]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaI(ilFace).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[il]) );
    }
      if ( jl >=0 && jl < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaJ(jlFace), (*this).State(jl), (*this).State(loc), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaJ(jlFace2), (*this).FAreaJ(jlFace), (*this).State(jl).UpdateWithConsVars(eqnState, x[jl]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  //double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(jl), eqnState, suth, (*this).Center(jl), (*this).Center(loc), (*this).FAreaJ(jlFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(jlFace2), (*this).FAreaJ(jlFace), (*this).State(jl).UpdateWithConsVars(eqnState, x[jl]), eqnState, suth, (*this).Vol(jl) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(jl), eqnState, (*this).FAreaJ(jlFace), x[jl]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaJ(jlFace).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[jl]) );
      }
      if ( kl >=0 && kl < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaK(klFace), (*this).State(kl), (*this).State(loc), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaK(klFace2), (*this).FAreaK(klFace), (*this).State(kl).UpdateWithConsVars(eqnState, x[kl]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  //double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(kl), eqnState, suth, (*this).Center(kl), (*this).Center(loc), (*this).FAreaK(klFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(klFace2), (*this).FAreaK(klFace), (*this).State(kl).UpdateWithConsVars(eqnState, x[kl]), eqnState, suth, (*this).Vol(kl) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(kl), eqnState, (*this).FAreaK(klFace), x[kl]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaK(klFace).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[kl]) );
      }

      double diagTimeVol = ( (*this).Vol(loc) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

      x[loc] = AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] + solTimeMmN[loc] + L[loc] - U[loc]) ; //normal at lower boundaries needs to be reversed, so add instead of subtract L


      // x[loc] = (1.0 - inp.MatrixRelaxation()) * x[loc] + inp.MatrixRelaxation() * AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] +
      // 							   solTimeMmN[loc] + L[loc]) ; //normal at lower boundaries needs to be reversed, so add i

    } //end forward sweep

    //backward sweep
    for ( int ii = (int)x.size()-1; ii >= 0; ii-- ){

      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int iuFace = GetUpperFaceI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int juFace = GetUpperFaceJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int kuFace = GetUpperFaceK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int iuFace2 = GetUpperFaceI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);
      int juFace2 = GetUpperFaceJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);
      int kuFace2 = GetUpperFaceK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax, 2);

      int iu = GetNeighborUpI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int ju = GetNeighborUpJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int ku = GetNeighborUpK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      if ( iu >=0 && iu < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaI(iuFace), (*this).State(loc), (*this).State(iu), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaI(iuFace2), (*this).FAreaI(iuFace), (*this).State(iu).UpdateWithConsVars(eqnState, x[iu]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  //double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(iu), eqnState, suth, (*this).Center(loc), (*this).Center(iu), (*this).FAreaI(iuFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(iuFace2), (*this).FAreaI(iuFace), (*this).State(iu).UpdateWithConsVars(eqnState, x[iu]), eqnState, suth, (*this).Vol(iu) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(iu), eqnState, (*this).FAreaI(iuFace), x[iu]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaI(iuFace).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[iu]) );
      }
      if ( ju >=0 && ju < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaJ(juFace), (*this).State(loc), (*this).State(ju), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaJ(juFace2), (*this).FAreaJ(juFace), (*this).State(ju).UpdateWithConsVars(eqnState, x[ju]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  //double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(ju), eqnState, suth, (*this).Center(loc), (*this).Center(ju), (*this).FAreaJ(juFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(juFace2), (*this).FAreaJ(juFace), (*this).State(ju).UpdateWithConsVars(eqnState, x[ju]), eqnState, suth, (*this).Vol(ju) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(ju), eqnState, (*this).FAreaJ(juFace), x[ju]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaJ(juFace).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[ju]) );
      }
      if ( ku >=0 && ku < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	//double specRad = ConvSpecRad( (*this).FAreaK(kuFace), (*this).State(loc), (*this).State(ku), eqnState);
	double specRad = CellSpectralRadius( (*this).FAreaK(kuFace2), (*this).FAreaK(kuFace), (*this).State(ku).UpdateWithConsVars(eqnState, x[ku]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  //double vSpecRad = (mRef/Re) * ViscFaceSpecRadTSL( (*this).State(ku), eqnState, suth, (*this).Center(loc), (*this).Center(ku), (*this).FAreaK(kuFace));
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(kuFace2), (*this).FAreaK(kuFace), (*this).State(ku).UpdateWithConsVars(eqnState, x[ku]), eqnState, suth, (*this).Vol(ku) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(ku), eqnState, (*this).FAreaK(kuFace), x[ku]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaK(kuFace).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[ku]) );
      }

      double diagTimeVol = ( (*this).Vol(loc) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

      x[loc] = x[loc] + AiiInv * ( Uold[loc] - U[loc] );

      // x[loc] = (1.0 - inp.MatrixRelaxation()) * x[loc] + inp.MatrixRelaxation() * AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] +
      //         solTimeMmN[loc] - U[loc]) ;


    } //end backward sweep

    Uold = U;

    //calculate residual
    colMatrix resid(x[0].Size());

    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      double diagTimeVol = ( (*this).Vol(loc) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      double Aii = ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() ;
      
      //normal at lower boundaries needs to be reversed, so add instead of subtract L
      resid = -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] + solTimeMmN[loc] - Aii * x[loc] + L[loc] - U[loc];
      //resid = -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] + solTimeMmN[loc] - Aii.Data(loc) * x[loc] + L[loc] - U[loc];
      l2Resid = l2Resid + resid * resid;
    }

  } //loop for sweeps

  return l2Resid.Sum();

}


//function to calculate the spectral radius on a cell face for the viscous fluxes using the thin shear layer approximation
double ViscFaceSpecRadTSL(const primVars &state, const idealGas &eqnState, const sutherland &suth, const vector3d<double> &centerL, const vector3d<double> &centerR, const vector3d<double> &fArea){

  vector3d<double> normArea = fArea / fArea.Mag();
  vector3d<double> dist = centerR - centerL;
  double mu = suth.GetViscosity(state.Temperature(eqnState));

  return 2.0 * mu / (state.Rho() * fabs(normArea.DotProd(dist)) ) ;
}

//function to return the inviscid spectral radius given a cell state, equation of state, and 2 face area vectors
double CellSpectralRadius(const vector3d<double> &fAreaL, const vector3d<double> &fAreaR, const primVars &state, const idealGas &eqnState){

  vector3d<double> normAreaL = fAreaL / fAreaL.Mag();
  vector3d<double> normAreaR = fAreaR / fAreaR.Mag();

  vector3d<double> normAvg = 0.5 * (normAreaL + normAreaR);
  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());

  return ( fabs(state.Velocity().DotProd(normAvg)) + state.SoS(eqnState) ) * fMag;

}

//function to calculate the spectral radius at a cell center for the viscous fluxes
double ViscCellSpectralRadius(const vector3d<double> &fAreaL, const vector3d<double> &fAreaR, const primVars &state, const idealGas &eqnState, const sutherland &suth, const double &vol){

  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag());
  double maxTerm = max(4.0 / (3.0 * state.Rho()), eqnState.Gamma() / state.Rho()) ;
  double mu = suth.GetViscosity(state.Temperature(eqnState));
  double viscTerm = mu / eqnState.GetPrandtl();

  return maxTerm * viscTerm * fMag * fMag / vol ;
}

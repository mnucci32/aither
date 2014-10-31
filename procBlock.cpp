#include "procBlock.h"
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
procBlock::procBlock(){
  numCells = 1;
  numVars = 5;
  numI = 1;
  numJ = 1;
  numK = 1;
  numGhosts = 0;
  parBlockStartI = 0;
  parBlockEndI = 0;
  parBlockStartJ = 0;
  parBlockEndJ = 0;
  parBlockStartK = 0;
  parBlockEndK = 0;

  int numFaces = (numI+1)*(numJ)*(numK);  

  vector<primVars> dummyState (numCells);              //dummy state variable
  vector<vector3d<double> > vec1(numFaces);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vec2(numCells);             //dummy vector variable lenght of number of cells
  vector<tensor<double> > tens(numCells);             //dummy tensor variable length of number of cells
  vector<double> scalar(numCells);                      //dummy scalar variable lenght of number of cells
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);

  state = dummyState;      

  center = vec2;
  fAreaI = vec1;
  fAreaJ = vec1;
  fAreaK = vec1;
  fCenterI = vec1;
  fCenterJ = vec1;
  fCenterK = vec1;

  residual = dummyResid;

  vol = scalar;
  avgWaveSpeed = scalar;
  dt = scalar;

  velGrad = tens;
  tempGrad = vec2;

}
//constructor -- initialize state vector with dummy variables
procBlock::procBlock(const plot3dBlock &blk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;

  if (eqnSet == "euler" || eqnSet == "navierStokes"){
    numVars = 5;
  }
  else{
    cerr << "ERROR: Error in procBlock::procBlock constructor. Equation set " << eqnSet << " is not recognized!" << endl;
    exit(0);
  }

  int viscCells = 1;
  if (eqnSet != "euler"){
    viscCells = numCells;
  }

  vector<primVars> dummyState (numCells);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);
  vector<tensor<double> > tens(viscCells);             //dummy tensor variable length of number of cells
  vector<vector3d<double> > vec(viscCells);             //dummy vector variable lenght of number of cells

  state = PadWithGhosts( dummyState, numGhosts, numI, numJ, numK );      

  vol = PadWithGhosts( blk.Volume(), numGhosts, numI, numJ, numK );
  center = PadWithGhosts( blk.Centroid(), numGhosts, numI, numJ, numK );
  fAreaI = PadWithGhosts( blk.FaceAreaI(), numGhosts, numI+1, numJ, numK );
  fAreaJ = PadWithGhosts( blk.FaceAreaJ(), numGhosts, numI, numJ+1, numK );
  fAreaK = PadWithGhosts( blk.FaceAreaK(), numGhosts, numI, numJ, numK+1 );
  fCenterI = PadWithGhosts( blk.FaceCenterI(), numGhosts, numI+1, numJ, numK );
  fCenterJ = PadWithGhosts( blk.FaceCenterJ(), numGhosts, numI, numJ+1, numK );
  fCenterK = PadWithGhosts( blk.FaceCenterK(), numGhosts, numI, numJ, numK+1 );

  avgWaveSpeed = dummyScalar;
  dt = dummyScalar;
  residual = dummyResid;

  velGrad = PadWithGhosts( tens, numGhosts, numI, numJ, numK );
  tempGrad = PadWithGhosts( vec, numGhosts, numI, numJ, numK );

}

//constructor -- assign passed variables to initialize state vector
procBlock::procBlock( const double density, const double pressure, const vector3d<double> vel, const plot3dBlock &blk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;

  if (eqnSet == "euler" || eqnSet == "navierStokes"){
    numVars = 5;
  }
  else{
    cerr << "ERROR: Error in procBlock::procBlock constructor. Equation set " << eqnSet << " is not recognized!" << endl;
    exit(0);
  }

  int viscCells = 1;
  if (eqnSet != "euler"){
    viscCells = numCells;
  }

  primVars singleState(density, pressure, vel);
  vector<primVars> dummyState (numCells, singleState);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);
  vector<tensor<double> > tens(viscCells);             //dummy tensor variable length of number of cells
  vector<vector3d<double> > vec(viscCells);             //dummy vector variable lenght of number of cells

  state = PadWithGhosts( dummyState, numGhosts, numI, numJ, numK );      

  vol = PadWithGhosts( blk.Volume(), numGhosts, numI, numJ, numK );
  center = PadWithGhosts( blk.Centroid(), numGhosts, numI, numJ, numK );
  fAreaI = PadWithGhosts( blk.FaceAreaI(), numGhosts, numI+1, numJ, numK );
  fAreaJ = PadWithGhosts( blk.FaceAreaJ(), numGhosts, numI, numJ+1, numK );
  fAreaK = PadWithGhosts( blk.FaceAreaK(), numGhosts, numI, numJ, numK+1 );
  fCenterI = PadWithGhosts( blk.FaceCenterI(), numGhosts, numI+1, numJ, numK );
  fCenterJ = PadWithGhosts( blk.FaceCenterJ(), numGhosts, numI, numJ+1, numK );
  fCenterK = PadWithGhosts( blk.FaceCenterK(), numGhosts, numI, numJ, numK+1 );

  avgWaveSpeed = dummyScalar;
  dt = dummyScalar;
  residual = dummyResid;

  velGrad = PadWithGhosts( tens, numGhosts, numI, numJ, numK );
  tempGrad = PadWithGhosts( vec, numGhosts, numI, numJ, numK );

}

//constructor -- assign passed state to initialize state vector
procBlock::procBlock( const primVars& inputState, const plot3dBlock &blk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;

  if (eqnSet == "euler" || eqnSet == "navierStokes"){
    numVars = 5;
  }
  else{
    cerr << "ERROR: Error in procBlock::procBlock constructor. Equation set " << eqnSet << " is not recognized!" << endl;
    exit(0);
  }

  int viscCells = 1;
  if (eqnSet != "euler"){
    viscCells = numCells;
  }

  vector<primVars> dummyState (numCells, inputState);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);
  vector<tensor<double> > tens(viscCells);             //dummy tensor variable length of number of cells
  vector<vector3d<double> > vec(viscCells);             //dummy vector variable lenght of number of cells

  state = PadWithGhosts( dummyState, numGhosts, numI, numJ, numK );      

  vol = PadWithGhosts( blk.Volume(), numGhosts, numI, numJ, numK );
  center = PadWithGhosts( blk.Centroid(), numGhosts, numI, numJ, numK );
  fAreaI = PadWithGhosts( blk.FaceAreaI(), numGhosts, numI+1, numJ, numK );
  fAreaJ = PadWithGhosts( blk.FaceAreaJ(), numGhosts, numI, numJ+1, numK );
  fAreaK = PadWithGhosts( blk.FaceAreaK(), numGhosts, numI, numJ, numK+1 );
  fCenterI = PadWithGhosts( blk.FaceCenterI(), numGhosts, numI+1, numJ, numK );
  fCenterJ = PadWithGhosts( blk.FaceCenterJ(), numGhosts, numI, numJ+1, numK );
  fCenterK = PadWithGhosts( blk.FaceCenterK(), numGhosts, numI, numJ, numK+1 );

  avgWaveSpeed = dummyScalar;
  dt = dummyScalar;
  residual = dummyResid;

  velGrad = PadWithGhosts( tens, numGhosts, numI, numJ, numK );
  tempGrad = PadWithGhosts( vec, numGhosts, numI, numJ, numK );

}

//member function to store the inviscid flux class in the place for the residual
void procBlock::AddToResidual(const inviscidFlux &flux, const int &ii){
  colMatrix temp(5);
  temp.SetData(0, flux.RhoVel());
  temp.SetData(1, flux.RhoVelU());
  temp.SetData(2, flux.RhoVelV());
  temp.SetData(3, flux.RhoVelW());
  temp.SetData(4, flux.RhoVelH());

  (*this).SetResidual( (*this).Residual(ii) + temp, ii); 
}

//member function to store the viscous flux class in the place for the residual
void procBlock::AddToResidual(const viscousFlux &flux, const int &ii){
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
void procBlock::CalcInvFluxI(const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI() + 1;
  int jmax = (*this).NumJ();
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

  for ( kk = 0 + (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( jj = 0 + (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( ii = 0 + (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

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

//function to calculate the fluxes on the j-faces
void procBlock::CalcInvFluxJ(const idealGas &eqnState, const input &inp, const int &bb){

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
void procBlock::CalcInvFluxK(const idealGas &eqnState, const input &inp, const int &bb){

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
void procBlock::CalcCellDt( const int &i, const int &j, const int &k, const double &cfl){

  int imax = (*this).NumI()-1;
  int jmax = (*this).NumJ()-1;
  int loc = GetLoc1D(i, j, k, imax, jmax);

  double dt = cfl * ((*this).Vol(loc) / (*this).AvgWaveSpeed(loc)) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}



void procBlock::CalcBlockTimeStep( const input &inputVars, const double &aRef){

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

void procBlock::UpdateBlock(const input &inputVars, const int &impFlag, const idealGas &eos, const double &aRef, const int &bb, const vector<colMatrix> &du, vector<double> &l2, vector<double> &linf, int &locMaxB){

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
void procBlock::ExplicitEulerTimeAdvance(const idealGas &eqnState, const int &loc ){

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
void procBlock::ImplicitTimeAdvance(const colMatrix &du, const idealGas &eqnState, const int &loc ){

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
void procBlock::RK4TimeAdvance( const primVars &currState, const idealGas &eqnState, const double &dt, const int &loc, const int &rk ){

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

void procBlock::ResetResidWS( ){

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

//a member function to add the cell volume divided by the cell time step to the time m - time n term
vector<colMatrix> procBlock::AddVolTime(const vector<colMatrix> &m, const vector<colMatrix> &n, const double &theta, const double &zeta) const {

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
void procBlock::DeltaNMinusOne(vector<colMatrix> &solDeltaNm1, const vector<colMatrix> &solTimeN, const idealGas &eqnState, const double &theta, const double &zeta){

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
vector<colMatrix> procBlock::GetCopyConsVars(const idealGas &eqnState) const {

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
double procBlock::LUSGS( const vector<vector3d<int> > &reorder, vector<colMatrix> &x, const vector<colMatrix> &solTimeMmN, const vector<colMatrix> &solDeltaNm1, const idealGas &eqnState, const input &inp, const sutherland &suth)const{

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

      x[loc] = AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) - solDeltaNm1[loc] - solTimeMmN[loc] + L[loc] ) ; //normal at lower boundaries needs to be reversed, so add instead of subtract L


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

      x[loc] = x[loc] - AiiInv * U[loc];

      // x[loc] = (1.0 - inp.MatrixRelaxation()) * x[loc] + inp.MatrixRelaxation() * AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] +
      //         solTimeMmN[loc] - U[loc]) ;


    } //end backward sweep


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

//member function to calculate the velocity gradient on the face using central differences
template <class T>
T FaceReconCentral(const T &velU, const T &velD, const vector3d<double> &pU, const vector3d<double> &pD, const vector3d<double> &pF){

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

//function to pad a vector with a specified number of ghost cells
template<class T>
vector<T> PadWithGhosts( const vector<T> &var, const int &numGhosts, const int &numI, const int &numJ, const int &numK ){

  int newI = numI + (numGhosts * 2);
  int newJ = numJ + (numGhosts * 2);
  int newK = numK + (numGhosts * 2);

  int newSize = newI * newJ * newK;

  vector<T> padBlk(newSize);

  T ghostVal;

  int newLoc = 0;
  int loc = 0;

  for ( int kk = 0; kk < newK; kk++ ){
    for ( int jj = 0; jj < newJ; jj++ ){
      for ( int ii = 0; ii < newI; ii++ ){

	newLoc = GetLoc1D(ii, jj, kk, newI, newJ);

	if ( ii < 2 || ii >= (numI + numGhosts) ) { //i ghost value
	  padBlk[newLoc] = ghostVal;
	}
	else if ( jj < 2 || jj >= (numJ + numGhosts) ) { //j ghost value
	  padBlk[newLoc] = ghostVal;
	}
	else if ( kk < 2 || kk >= (numK + numGhosts) ) { //k ghost value
	  padBlk[newLoc] = ghostVal;
	}
	else{ //physical cells

	  loc = GetLoc1D(ii-numGhosts, jj-numGhosts, kk-numGhosts, numI, numJ);
	  padBlk[newLoc] = var[loc];

	}

      }
    }
  }

  return padBlk;
}


//member function to initialize gradients
void procBlock::InitializeGrads(){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  int ii = 0;
  int jj = 0;
  int kk = 0;

  int loc = 0;

  vector3d<double> initialVector(0.0, 0.0, 0.0);
  tensor<double> initialTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	(*this).SetVelGrad( initialTensor, loc);
	(*this).SetTempGrad( initialVector, loc);

      }
    }
  }

}


//member function to calculate the velocity gradient at the cell center
void procBlock::CalcVelGradGG(const vector3d<double> &vl, const vector3d<double> &vu, const vector3d<double> &al, const vector3d<double> &au, const double &vol, const int &loc){

  //vl is the velocity vector at the lower face of the cell at which the velocity gradient is being calculated
  //vu is the velocity vector at the upper face of the cell at which the velocity gradient is being calculated

  //ail is the area vector at the lower face of the cell at which the velocity gradient is being calculated
  //aiu is the area vector at the upper face of the cell at which the velocity gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the velocity gradient should be stored

  tensor<double> temp;
  double invVol = 1.0/vol;

  //define velocity gradient tensor
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetXX( invVol * (vu.X()*au.X() - vl.X()*al.X() ) );
  temp.SetXY( invVol * (vu.Y()*au.X() - vl.Y()*al.X() ) );
  temp.SetXZ( invVol * (vu.Z()*au.X() - vl.Z()*al.X() ) );

  temp.SetYX( invVol * (vu.X()*au.Y() - vl.X()*al.Y() ) );
  temp.SetYY( invVol * (vu.Y()*au.Y() - vl.Y()*al.Y() ) );
  temp.SetYZ( invVol * (vu.Z()*au.Y() - vl.Z()*al.Y() ) );

  temp.SetZX( invVol * (vu.X()*au.Z() - vl.X()*al.Z() ) );
  temp.SetYY( invVol * (vu.Y()*au.Z() - vl.Y()*al.Z() ) );
  temp.SetZZ( invVol * (vu.Z()*au.Z() - vl.Z()*al.Z() ) );

  (*this).SetVelGrad( (*this).VelGrad(loc) + temp, loc);

}

//member function to calculate the temperature gradient at the cell center
void procBlock::CalcTempGradGG(const double &tl, const double &tu, const vector3d<double> &al, const vector3d<double> &au, const double &vol, const int &loc){

  //tl is the temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tu is the temperature at the upper face of the cell at which the temperature gradient is being calculated
  
  //al is the area vector at the lower face of the cell at which the temperature gradient is being calculated
  //au is the area vector at the upper face of the cell at which the temperature gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the temperature gradient should be stored

  vector3d<double> temp;
  double invVol = 1.0/vol;

  //define temperature gradient vector
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetX( invVol * (tu*au.X() - tl*al.X() ) );
  temp.SetY( invVol * (tu*au.Y() - tl*al.Y() ) );
  temp.SetZ( invVol * (tu*au.Z() - tl*al.Z() ) );

  (*this).SetTempGrad( (*this).TempGrad(loc) + temp, loc);

}

//member function to calculate gradients at centers
void procBlock::CalcCellGradsI(const idealGas &eqnState, const sutherland &suth, const input &inp, const int &bb){

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

  int ifLow = 0;
  int ifUp = 0;

  string bcName = "undefined";

  double viscConstant = 1.0;

  vector3d<double> vl, vu;
  double tl = 0.0;
  double tu = 0.0;

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

	ifLow = GetLowerFaceI(ii, jj, kk, imax, jmax); 
	ifUp  = GetUpperFaceI(ii, jj, kk, imax, jmax);

	//test if on i lower boundary
	if (ii == 0){
	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "il");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaI(ifLow), "il", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterI(ifLow) - (*this).Center(loc) ) + (*this).Center(loc);

	  vl = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterI(ifLow) );
	  tl = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterI(ifLow) );

	}
	else{
	  vl = FaceReconCentral( (*this).State(iLow).Velocity(), (*this).State(loc).Velocity(), (*this).Center(iLow), (*this).Center(loc), (*this).FCenterI(ifLow) );
	  tl = FaceReconCentral( (*this).State(iLow).Temperature(eqnState), (*this).State(loc).Temperature(eqnState), (*this).Center(iLow), (*this).Center(loc), (*this).FCenterI(ifLow) );
	}

	//test if on i upper boundary
	if (ii == imax-1){
	  bcName = bound.GetBCName(ii+1, jj, kk, "iu");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaI(ifUp), "iu", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterI(ifUp) - (*this).Center(loc) ) + (*this).Center(loc);

	  vu = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterI(ifUp) );
	  tu = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterI(ifUp) );

	}
	else{
	  vu = FaceReconCentral( (*this).State(iUp).Velocity(),  (*this).State(loc).Velocity(), (*this).Center(iUp),  (*this).Center(loc), (*this).FCenterI(ifUp)  );
	  tu = FaceReconCentral( (*this).State(iUp).Temperature(eqnState),  (*this).State(loc).Temperature(eqnState), (*this).Center(iUp),  (*this).Center(loc), (*this).FCenterI(ifUp)  );
	}


	//calculate gradients for cell
	CalcVelGradGG(vl, vu, (*this).FAreaI(ifLow), (*this).FAreaI(ifUp), (*this).Vol(loc), loc);
	CalcTempGradGG(tl, tu, (*this).FAreaI(ifLow), (*this).FAreaI(ifUp), (*this).Vol(loc), loc);

	//calculate cell viscous spectral radius
	double maxViscSpeed = ViscCellSpectralRadius((*this).FAreaI(ifLow), (*this).FAreaI(ifUp), (*this).State(loc), eqnState, suth, (*this).Vol(loc));
	(*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(loc) + viscConstant * (mRef/Re) * maxViscSpeed, loc); 

      }
    }
  }

}

//member function to calculate gradients at centers
void procBlock::CalcCellGradsJ(const idealGas &eqnState, const sutherland &suth, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;

  int loc = 0;
  int jLow = 0;
  int jUp = 0;

  int jfLow = 0;
  int jfUp = 0;

  string bcName = "undefined";

  double viscConstant = 1.0;

  vector3d<double> vl, vu;
  double tl = 0.0;
  double tu = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  primVars ghostState;
  vector3d<double> ghostDistance;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	jLow = GetNeighborLowJ(ii, jj, kk, imax, jmax); 
	jUp  = GetNeighborUpJ(ii, jj, kk, imax, jmax);

	jfLow = GetLowerFaceJ(ii, jj, kk, imax, jmax); 
	jfUp  = GetUpperFaceJ(ii, jj, kk, imax, jmax);

	//test if on j lower boundary
	if (jj == 0){
	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaJ(jfLow), "jl", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterJ(jfLow) - (*this).Center(loc) ) + (*this).Center(loc);

	  vl = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterJ(jfLow) );
	  tl = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterJ(jfLow) );

	}
	else{
	  vl = FaceReconCentral( (*this).State(jLow).Velocity(), (*this).State(loc).Velocity(), (*this).Center(jLow), (*this).Center(loc), (*this).FCenterJ(jfLow) );
	  tl = FaceReconCentral( (*this).State(jLow).Temperature(eqnState), (*this).State(loc).Temperature(eqnState), (*this).Center(jLow), (*this).Center(loc), (*this).FCenterJ(jfLow) );
	}

	//test if on j upper boundary
	if (jj == jmax-1){
	  bcName = bound.GetBCName(ii, jj+1, kk, "ju");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaJ(jfUp), "ju", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterJ(jfUp) - (*this).Center(loc) ) + (*this).Center(loc);

	  vu = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterJ(jfUp) );
	  tu = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterJ(jfUp) );

	}
	else{
	  vu = FaceReconCentral( (*this).State(jUp).Velocity(),  (*this).State(loc).Velocity(), (*this).Center(jUp),  (*this).Center(loc), (*this).FCenterJ(jfUp)  );
	  tu = FaceReconCentral( (*this).State(jUp).Temperature(eqnState),  (*this).State(loc).Temperature(eqnState), (*this).Center(jUp),  (*this).Center(loc), (*this).FCenterJ(jfUp)  );
	}


	//calculate gradients for cell
	CalcVelGradGG(vl, vu, (*this).FAreaJ(jfLow), (*this).FAreaJ(jfUp), (*this).Vol(loc), loc);
	CalcTempGradGG(tl, tu, (*this).FAreaJ(jfLow), (*this).FAreaJ(jfUp), (*this).Vol(loc), loc);

	//calculate cell viscous spectral radius
	double maxViscSpeed = ViscCellSpectralRadius((*this).FAreaJ(jfLow), (*this).FAreaJ(jfUp), (*this).State(loc), eqnState, suth, (*this).Vol(loc));
	(*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(loc) + viscConstant * (mRef/Re) * maxViscSpeed, loc); 

      }
    }
  }

}

//member function to calculate gradients at centers
void procBlock::CalcCellGradsK(const idealGas &eqnState, const sutherland &suth, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

  const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;

  int loc = 0;
  int kLow = 0;
  int kUp = 0;

  int kfLow = 0;
  int kfUp = 0;

  string bcName = "undefined";

  double viscConstant = 1.0;

  vector3d<double> vl, vu;
  double tl = 0.0;
  double tu = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  primVars ghostState;
  vector3d<double> ghostDistance;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);
	kLow = GetNeighborLowK(ii, jj, kk, imax, jmax); 
	kUp  = GetNeighborUpK(ii, jj, kk, imax, jmax);

	kfLow = GetLowerFaceK(ii, jj, kk, imax, jmax); 
	kfUp  = GetUpperFaceK(ii, jj, kk, imax, jmax);

	//test if on j lower boundary
	if (kk == 0){
	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaK(kfLow), "kl", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterK(kfLow) - (*this).Center(loc) ) + (*this).Center(loc);

	  vl = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterK(kfLow) );
	  tl = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterK(kfLow) );

	}
	else{
	  vl = FaceReconCentral( (*this).State(kLow).Velocity(), (*this).State(loc).Velocity(), (*this).Center(kLow), (*this).Center(loc), (*this).FCenterK(kfLow) );
	  tl = FaceReconCentral( (*this).State(kLow).Temperature(eqnState), (*this).State(loc).Temperature(eqnState), (*this).Center(kLow), (*this).Center(loc), (*this).FCenterK(kfLow) );
	}

	//test if on k upper boundary
	if (kk == kmax-1){
	  bcName = bound.GetBCName(ii, jj, kk+1, "ku");

	  ghostState = (*this).State(loc).GetGhostState( bcName, (*this).FAreaK(kfUp), "ku", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterK(kfUp) - (*this).Center(loc) ) + (*this).Center(loc);

	  vu = FaceReconCentral( ghostState.Velocity(), (*this).State(loc).Velocity(), ghostDistance, (*this).Center(loc), (*this).FCenterK(kfUp) );
	  tu = FaceReconCentral( ghostState.Temperature(eqnState), (*this).State(loc).Temperature(eqnState), ghostDistance, (*this).Center(loc), (*this).FCenterK(kfUp) );

	}
	else{
	  vu = FaceReconCentral( (*this).State(kUp).Velocity(),  (*this).State(loc).Velocity(), (*this).Center(kUp),  (*this).Center(loc), (*this).FCenterK(kfUp)  );
	  tu = FaceReconCentral( (*this).State(kUp).Temperature(eqnState),  (*this).State(loc).Temperature(eqnState), (*this).Center(kUp),  (*this).Center(loc), (*this).FCenterK(kfUp)  );
	}


	//calculate gradients for cell
	CalcVelGradGG(vl, vu, (*this).FAreaK(kfLow), (*this).FAreaK(kfUp), (*this).Vol(loc), loc);
	CalcTempGradGG(tl, tu, (*this).FAreaK(kfLow), (*this).FAreaK(kfUp), (*this).Vol(loc), loc);

	//calculate cell viscous spectral radius
	double maxViscSpeed = ViscCellSpectralRadius((*this).FAreaK(kfLow), (*this).FAreaK(kfUp), (*this).State(loc), eqnState, suth, (*this).Vol(loc));
	(*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(loc) + viscConstant * (mRef/Re) * maxViscSpeed, loc); 

      }
    }
  }

}

//member function to calculate viscous fluxes on i-faces
void procBlock::CalcViscFluxI(const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK() - 1;

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int iLow = 0;
  int iUp = 0;

  string lstr = "left";
  string rstr = "right";

  string bcName = "undefined";

  primVars ghostState;
  vector3d<double> ghostDistance;

  viscousFlux tempViscFlux;

  tensor<double> velGrad;
  vector3d<double> tGrad;
  vector3d<double> vel;
  double mu = 0.0;

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

	  ghostState = (*this).State(iUp).GetGhostState( bcName, (*this).FAreaI(loc), "il", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterI(loc) - (*this).Center(iUp) ) + (*this).Center(iUp);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(iUp);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(iUp).Velocity(), ghostDistance, (*this).Center(iUp), (*this).FCenterI(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(iUp);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(iUp).Temperature(eqnState)), ghostDistance, (*this).Center(iUp), (*this).FCenterI(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(iUp).Rho(), ghostDistance, (*this).Center(iUp), (*this).FCenterI(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaI(loc) );

	  //at lower boundary normal points into cell, so need to subtract from residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaI(loc).Mag(), iUp);

	}
	else if (ii == imax-1){ //---------------------------------------------------------------------------------------------------------------------------------------

	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "iu");

	  ghostState = (*this).State(iLow).GetGhostState( bcName, (*this).FAreaI(loc), "iu", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterI(loc) - (*this).Center(iLow) ) + (*this).Center(iLow);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(iLow);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(iLow).Velocity(), ghostDistance, (*this).Center(iLow), (*this).FCenterI(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(iLow);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(iLow).Temperature(eqnState)), ghostDistance, (*this).Center(iLow), (*this).FCenterI(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(iLow).Rho(), ghostDistance, (*this).Center(iLow), (*this).FCenterI(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaI(loc) );

	  //at upper boundary normal points out of cell, so need to add to residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaI(loc).Mag(), iLow);

	}
	else{ //------------------------------------------------------------------------------------------------------------------------------------------------

	  iLow  = GetCellFromFaceLowerI(ii, jj, kk, imax, jmax);
	  iUp  = GetCellFromFaceUpperI(ii, jj, kk, imax, jmax);

	  //Get velocity gradient at face
	  velGrad = FaceReconCentral( (*this).VelGrad(iLow), (*this).VelGrad(iUp), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );
	  //Get velocity at face
	  vel = FaceReconCentral( (*this).State(iLow).Velocity(), (*this).State(iUp).Velocity(), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );
	  //Get temperature gradient at face
	  tGrad = FaceReconCentral( (*this).TempGrad(iLow), (*this).TempGrad(iUp), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity( (*this).State(iLow).Temperature(eqnState) ), suth.GetViscosity( (*this).State(iUp).Temperature(eqnState) ), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( (*this).State(iLow).Rho(), (*this).State(iUp).Rho(), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaI(loc) );

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaI(loc).Mag(), iLow);
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaI(loc).Mag(), iUp);

	}

      }
    }
  }


}



//member function to calculate viscous fluxes on j-faces
void procBlock::CalcViscFluxJ(const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() - 1;

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int jLow = 0;
  int jUp = 0;

  string lstr = "left";
  string rstr = "right";

  string bcName = "undefined";
  primVars ghostState;
  vector3d<double> ghostDistance;

  viscousFlux tempViscFlux;

  tensor<double> velGrad;
  vector3d<double> tGrad;
  vector3d<double> vel;
  double mu = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (jj == 0){ //-------------------------------------------------------------------------------------------------------------------------------------

	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "jl");

	  ghostState = (*this).State(jUp).GetGhostState( bcName, (*this).FAreaJ(loc), "jl", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterJ(loc) - (*this).Center(jUp) ) + (*this).Center(jUp);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(jUp);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(jUp).Velocity(), ghostDistance, (*this).Center(jUp), (*this).FCenterJ(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(jUp);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(jUp).Temperature(eqnState)), ghostDistance, (*this).Center(jUp), (*this).FCenterJ(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(jUp).Rho(), ghostDistance, (*this).Center(jUp), (*this).FCenterJ(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaJ(loc) );

	  //at lower boundary normal points into cell, so need to subtract from residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaJ(loc).Mag(), jUp);

	}
	else if (jj == jmax-1){ //--------------------------------------------------------------------------------------------------------------------------------

	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "ju");

	  ghostState = (*this).State(jLow).GetGhostState( bcName, (*this).FAreaJ(loc), "ju", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterJ(loc) - (*this).Center(jLow) ) + (*this).Center(jLow);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(jLow);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(jLow).Velocity(), ghostDistance, (*this).Center(jLow), (*this).FCenterJ(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(jLow);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(jLow).Temperature(eqnState)), ghostDistance, (*this).Center(jLow), (*this).FCenterJ(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(jLow).Rho(), ghostDistance, (*this).Center(jLow), (*this).FCenterJ(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaJ(loc) );

	  //at upper boundary normal points out of cell, so need to add to residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaJ(loc).Mag(), jLow);

	}
	else{ //-----------------------------------------------------------------------------------------------------------------------------------------------

	  jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imax, jmax);
	  jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imax, jmax);

	  //Get velocity gradient at face
	  velGrad = FaceReconCentral( (*this).VelGrad(jLow), (*this).VelGrad(jUp), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );
	  //Get velocity at face
	  vel = FaceReconCentral( (*this).State(jLow).Velocity(), (*this).State(jUp).Velocity(), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );
	  //Get temperature gradient at face
	  tGrad = FaceReconCentral( (*this).TempGrad(jLow), (*this).TempGrad(jUp), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity( (*this).State(jLow).Temperature(eqnState) ), suth.GetViscosity( (*this).State(jUp).Temperature(eqnState) ), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( (*this).State(jLow).Rho(), (*this).State(jUp).Rho(), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaJ(loc) );

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaJ(loc).Mag(), jLow);
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaJ(loc).Mag(), jUp);

	}

      }
    }
  }


}


//member function to calculate viscous fluxes on j-faces
void procBlock::CalcViscFluxK(const sutherland &suth, const idealGas &eqnState, const input &inp, const int &bb){

  int imax = (*this).NumI() - 1;
  int jmax = (*this).NumJ() - 1;
  int kmax = (*this).NumK();

 const boundaryConditions bound = inp.BC(bb);

  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int kLow = 0;
  int kUp = 0;

  string lstr = "left";
  string rstr = "right";

  string bcName = "undefined";
  primVars ghostState;
  vector3d<double> ghostDistance;

  viscousFlux tempViscFlux;

  tensor<double> velGrad;
  vector3d<double> tGrad;
  vector3d<double> vel;
  double mu = 0.0;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( kk = 0; kk < kmax; kk++){   
    for ( jj = 0; jj < jmax; jj++){    
      for ( ii = 0; ii < imax; ii++){      

	loc = GetLoc1D(ii, jj, kk, imax, jmax);

	if (kk == 0){ //-----------------------------------------------------------------------------------------------------------------------------------

	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "kl");

	  ghostState = (*this).State(kUp).GetGhostState( bcName, (*this).FAreaK(loc), "kl", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterK(loc) - (*this).Center(kUp) ) + (*this).Center(kUp);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(kUp);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(kUp).Velocity(), ghostDistance, (*this).Center(kUp), (*this).FCenterK(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(kUp);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(kUp).Temperature(eqnState)), ghostDistance, (*this).Center(kUp), (*this).FCenterK(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(kUp).Rho(), ghostDistance, (*this).Center(kUp), (*this).FCenterK(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaK(loc) );

	  //at lower boundary normal points into cell, so need to subtract from residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is positive
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaK(loc).Mag(), kUp);

	}
	else if (kk == kmax-1){ //----------------------------------------------------------------------------------------------------------------------------

	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);

	  //find boundary type and get ghost state
	  bcName = bound.GetBCName(ii, jj, kk, "ku");

	  ghostState = (*this).State(kLow).GetGhostState( bcName, (*this).FAreaK(loc), "ku", inp, eqnState );
	  ghostDistance = 2.0 * ( (*this).FCenterK(loc) - (*this).Center(kLow) ) + (*this).Center(kLow);

	  //Get velocity gradient at face
	  velGrad = (*this).VelGrad(kLow);
	  //Get velocity at face
	  vel = FaceReconCentral( ghostState.Velocity(), (*this).State(kLow).Velocity(), ghostDistance, (*this).Center(kLow), (*this).FCenterK(loc) );
	  //Get temperature gradient at face
	  tGrad = (*this).TempGrad(kLow);
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity(ghostState.Temperature(eqnState)), suth.GetViscosity((*this).State(kLow).Temperature(eqnState)), ghostDistance, (*this).Center(kLow), (*this).FCenterK(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( ghostState.Rho(), (*this).State(kLow).Rho(), ghostDistance, (*this).Center(kLow), (*this).FCenterK(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaK(loc) );

	  //at upper boundary normal points out of cell, so need to add to residual
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is negative
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaK(loc).Mag(), kLow);

	}
	else{ //-----------------------------------------------------------------------------------------------------------------------------------------------------

	  kLow  = GetCellFromFaceLowerK(ii, jj, kk, imax, jmax);
	  kUp  = GetCellFromFaceUpperK(ii, jj, kk, imax, jmax);

	  //Get velocity gradient at face
	  velGrad = FaceReconCentral( (*this).VelGrad(kLow), (*this).VelGrad(kUp), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	  //Get velocity at face
	  vel = FaceReconCentral( (*this).State(kLow).Velocity(), (*this).State(kUp).Velocity(), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	  //Get temperature gradient at face
	  tGrad = FaceReconCentral( (*this).TempGrad(kLow), (*this).TempGrad(kUp), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	  //Get viscosity at face
	  mu = FaceReconCentral( suth.GetViscosity( (*this).State(kLow).Temperature(eqnState) ), suth.GetViscosity( (*this).State(kUp).Temperature(eqnState) ), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	  mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)
	  //Get density at face
	  //rho = FaceReconCentral( (*this).State(kLow).Rho(), (*this).State(kUp).Rho(), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );

	  //calculate viscous flux
	  tempViscFlux.SetFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaK(loc) );

	  //area vector points from left to right, so add to left cell, subtract from right cell
	  //but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaK(loc).Mag(), kLow);
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaK(loc).Mag(), kUp);

	}

      }
    }
  }


}

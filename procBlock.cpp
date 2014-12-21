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
  parBlock = 0;
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

}
//constructor -- initialize state vector with dummy variables
procBlock::procBlock(const plot3dBlock &blk, const int& numBlk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;
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

  vector<primVars> dummyState (numCells);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);

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

}

//constructor -- assign passed variables to initialize state vector
procBlock::procBlock( const double density, const double pressure, const vector3d<double> vel, const plot3dBlock &blk, const int &numBlk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;
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

  primVars singleState(density, pressure, vel);
  vector<primVars> dummyState (numCells,singleState);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);

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

}

//constructor -- assign passed state to initialize state vector
procBlock::procBlock( const primVars& inputState, const plot3dBlock &blk, const int &numBlk, const int &numG, const string &eqnSet){
  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;
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

  vector<double> dummyScalar (numCells);                 //dummy time variable
  vector<primVars> dummyState (numCells,inputState);              //dummy state variable
  colMatrix singleResid(numVars);
  singleResid.Zero();
  vector<colMatrix> dummyResid(numCells, singleResid);

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
void procBlock::CalcInvFluxI(const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI() + 1;
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts() + 1;
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double maxWS = 0.0;
  primVars faceStateLower,faceStateUpper;

  //loop over all physical cells
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
	int lowerING = GetCellFromFaceLowerI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int upperING = GetCellFromFaceUpperI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG);
	int upperI = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG);
	int lower2I = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG, 2);
	int upper2I = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG, 2);
	int upFaceI = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int upFace2I = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceI = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);
	int lowFace2I = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG, 2);

	if (inp.Kappa() == -2.0){  //if value is still default, use constant reconstruction
	  faceStateLower = (*this).State( lowerI ).FaceReconConst();
	  faceStateUpper = (*this).State( upperI ).FaceReconConst();
	}
	else{

	  double upwind2L =  (*this).FCenterI( lowFaceI ).Distance( (*this).FCenterI( lowFace2I ) );
	  double upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );
	  double downwindL = (*this).FCenterI( loc      ).Distance( (*this).FCenterI( upFaceI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( lower2I ), (*this).State( upperI ),
								   "left", inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  double upwind2U =  (*this).FCenterI( upFaceI ).Distance( (*this).FCenterI( upFace2I ) );
	  double upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );
	  double downwindU = (*this).FCenterI( loc     ).Distance( (*this).FCenterI( lowFaceI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( upper2I ), (*this).State( lowerI ),
								   "right", inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );

	}

	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( ii > (*this).NumGhosts() ){
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerING);
	}
	if ( ii < imax - 1 + (*this).NumGhosts() ){
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperING);
	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperING) + maxWS, upperING);
	}

      }
    }
  }

}

//function to calculate the fluxes on the j-faces
void procBlock::CalcInvFluxJ(const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() + 1;
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts() + 1;

  double maxWS = 0.0;
  primVars faceStateLower, faceStateUpper;

  //loop over all physical cells
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
	int lowerJNG = GetCellFromFaceLowerJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int upperJNG = GetCellFromFaceUpperJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG);
	int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG);
	int lower2J = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG, 2);
	int upper2J = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG, 2);
	int upFaceJ = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG);
	int upFace2J = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG);
	int lowFace2J = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG, 2);

	if ( inp.Kappa() == -2.0 ){                         //if value is still default, use constant reconstruction
	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();
	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();
	}
	else{

	  double upwind2L =  (*this).FCenterJ( lowFaceJ ).Distance( (*this).FCenterJ( lowFace2J ) );
	  double upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );
	  double downwindL = (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( upFaceJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( lower2J ),
								   (*this).State( upperJ ),"left", inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  double upwind2U =  (*this).FCenterJ( upFaceJ ).Distance( (*this).FCenterJ( upFace2J ) );
	  double upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );
	  double downwindU = (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( lowFaceJ ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( upper2J ),
								   (*this).State( lowerJ ),"right", inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );
	}

	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( jj > (*this).NumGhosts() ){
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJNG);
	}
	if ( jj < jmax - 1 + (*this).NumGhosts() ){
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJNG);
	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJNG) + maxWS, upperJNG);
	}

      }
    }
  }

}

//function to calculate the fluxes on the k-faces
void procBlock::CalcInvFluxK(const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() + 1;

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double maxWS = 0.0;
  primVars faceStateLower, faceStateUpper;

  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
	int lowerKNG = GetCellFromFaceLowerK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int upperKNG = GetCellFromFaceUpperK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG);
	int upper2K = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG, 2);
	int lower2K = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG, 2);
	int upperK = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG);
	int upFaceK = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int upFace2K = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceK = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);
	int lowFace2K = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG, 2);

	if ( inp.Kappa() == -2.0 ){                         //if value is still default, use constant reconstruction
	  faceStateLower = (*this).State( lowerK ).FaceReconConst();
	  faceStateUpper = (*this).State( upperK ).FaceReconConst();
	}
	else{

	  double upwind2L =  (*this).FCenterK( lowFaceK ).Distance( (*this).FCenterK( lowFace2K ) );
	  double upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );
	  double downwindL = (*this).FCenterK( loc      ).Distance( (*this).FCenterK( upFaceK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( lower2K ), (*this).State( upperK ),
								   "left", inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  double upwind2U =  (*this).FCenterK( upFaceK ).Distance( (*this).FCenterK( upFace2K ) );
	  double upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );
	  double downwindU = (*this).FCenterK( loc     ).Distance( (*this).FCenterK( lowFaceK ) );

	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( upper2K ), (*this).State( lowerK ),
								   "right", inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );

	}

	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( kk > (*this).NumGhosts() ){
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerKNG);
	}
	if ( kk < kmax - 1 + (*this).NumGhosts() ){
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperKNG);
	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the lower faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperKNG) + maxWS, upperKNG);
	}

      }
    }
  }

}

//member function to calculate the local time step. (i,j,k) are cell indices
void procBlock::CalcCellDt( const int &i, const int &j, const int &k, const double &cfl){

  int loc = GetLoc1D(i, j, k, (*this).NumI(), (*this).NumJ());
  int locG = GetLoc1D(i + (*this).NumGhosts(), j + (*this).NumGhosts(), k + (*this).NumGhosts()
		      , (*this).NumI() + 2 * (*this).NumGhosts(), (*this).NumJ() + 2 * (*this).NumGhosts());
  double dt = cfl * ((*this).Vol(locG) / (*this).AvgWaveSpeed(loc)) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}


void procBlock::CalcBlockTimeStep( const input &inputVars, const double &aRef){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //loop over all physical cells - no ghost cells for dt variable
  for ( int kk = 0; kk < kmax; kk++ ){          
    for ( int jj = 0; jj < jmax; jj++ ){          
      for ( int ii = 0; ii < imax; ii++ ){          

	int loc = GetLoc1D(ii, jj, kk, imax, jmax);

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

void procBlock::UpdateBlock(const input &inputVars, const int &impFlag, const idealGas &eos, const double &aRef, const int &bb, const vector<colMatrix> &du, colMatrix &l2, colMatrix &linf, int &locMaxB){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  if ( inputVars.TimeIntegration() != "rk4" ){
    for ( int kk = 0; kk < kmax; kk++ ){          //loop over all physical cells
      for ( int jj = 0; jj < jmax; jj++ ){          
	for ( int ii = 0; ii < imax; ii++ ){          

	  int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	  int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	  if (inputVars.TimeIntegration() == "explicitEuler"){
	    (*this).ExplicitEulerTimeAdvance(eos, locG, loc);
	  }
	  else if (impFlag){
	    (*this).ImplicitTimeAdvance(du[loc], eos, locG);
	  }


	  l2 = l2 + (*this).Residual(loc) * (*this).Residual(loc);
	  for ( int ll = 0; ll < l2.Size(); ll++ ){

	    if ( (*this).Residual(loc,ll) > linf.Data(4) ){
	      linf.SetData(4, (*this).Residual(loc,ll) );
	      linf.SetData(3, (double)ll+1 );
	      linf.SetData(2, (double)kk );
	      linf.SetData(1, (double)jj );
	      linf.SetData(0, (double)ii );
	      locMaxB = bb;
	    }
	  }


	}
      }
    }
  }
  else if ( inputVars.TimeIntegration() == "rk4" ){

    vector<primVars> stateN(imax*jmax*kmax);
    vector<double> dtN(imax*jmax*kmax);

    for ( int rr = 0; rr < 4; rr ++ ){
      for ( int kk = 0; kk < kmax; kk++ ){          //loop over all physical cells
	for ( int jj = 0; jj < jmax; jj++ ){          
	  for ( int ii = 0; ii < imax; ii++ ){          

	    int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	    int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	    //save state and local time step at time n
	    if (rr == 0){
	      stateN[locG] = (*this).State(locG);
	      dtN[loc] = (*this).Dt(loc);
	    }

	    (*this).RK4TimeAdvance(stateN[locG], eos, dtN[loc], locG, loc, rr);

	    if (rr ==3){
	      l2 = l2 + (*this).Residual(loc) * (*this).Residual(loc);
	      for ( int ll = 0; ll < l2.Size(); ll++ ){
	
		if ( (*this).Residual(loc,ll) > linf.Data(4) ){
		  linf.SetData(4, (*this).Residual(loc,ll) );
		  linf.SetData(3, (double)ll+1 );
		  linf.SetData(2, (double)kk );
		  linf.SetData(1, (double)jj );
		  linf.SetData(0, (double)ii );
		  locMaxB = bb;
		}
	      }
	    }

	  }
	}
      }
      //for multistage RK4 method, calculate fluxes and residuals again
      if (rr < 3){ //no need to calculate fluxes after final RK interation
	(*this).CalcInvFluxI(eos, inputVars);
	(*this).CalcInvFluxJ(eos, inputVars);
	(*this).CalcInvFluxK(eos, inputVars);
	(*this).CalcBlockTimeStep(inputVars, aRef);
      }


    }
  }
  else {
    cerr << "ERROR: Time integration scheme " << inputVars.TimeIntegration() << " is not recognized!" << endl;
  }



}


//member function to advance the state vector to time n+1 using explicit Euler method
void procBlock::ExplicitEulerTimeAdvance(const idealGas &eqnState, const int &locG, const int &loc ){

  colMatrix consVars = (*this).State(locG).ConsVars(eqnState);

  //calculate updated conserved variables
  consVars.SetData(0, consVars.Data(0) - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc,0) );
  consVars.SetData(1, consVars.Data(1) - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc,1) );
  consVars.SetData(2, consVars.Data(2) - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc,2) );
  consVars.SetData(3, consVars.Data(3) - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc,3) );
  consVars.SetData(4, consVars.Data(4) - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc,4) );

  //calculate updated primative variables
  vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  primVars tempState (consVars.Data(0),
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  (*this).SetState(tempState, locG);

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
void procBlock::RK4TimeAdvance( const primVars &currState, const idealGas &eqnState, const double &dt, const int &locG, const int &loc, const int &rk ){

  double alpha[4] = {0.25, 1.0/3.0, 0.5, 1.0};

  colMatrix consVars = currState.ConsVars(eqnState);

  //calculate updated conserved variables
  consVars.SetData(0, consVars.Data(0) - dt / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc,0) );
  consVars.SetData(1, consVars.Data(1) - dt / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc,1) );
  consVars.SetData(2, consVars.Data(2) - dt / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc,2) );
  consVars.SetData(3, consVars.Data(3) - dt / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc,3) );
  consVars.SetData(4, consVars.Data(4) - dt / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc,4) );

  //calculate updated primative variables
  vector3d<double> vel(consVars.Data(1)/consVars.Data(0), consVars.Data(2)/consVars.Data(0), consVars.Data(3)/consVars.Data(0));

  primVars tempState (consVars.Data(0),
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars.Data(0), consVars.Data(4)/consVars.Data(0), vel.Mag() ) );

  (*this).SetState(tempState, locG);
}

void procBlock::ResetResidWS( ){

  colMatrix initial( (*this).Residual(0).Size() );
  initial.Zero();

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //loop over all physical cells - no ghost cells in residual variable
  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	int loc = GetLoc1D(ii, jj, kk, imax, jmax);

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

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  vector<colMatrix> mMinusN(m.size());

  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	double I = ( (*this).Vol(locG) * (1.0 + zeta) ) / ( (*this).Dt(loc) * theta ) ;
	mMinusN[loc] = I * (m[loc] - n[loc]);
      }
    }
  }
  return mMinusN;
}


//member function to calculate the delta n-1 term for the implicit bdf2 solver
void procBlock::DeltaNMinusOne(vector<colMatrix> &solDeltaNm1, const vector<colMatrix> &solTimeN, const idealGas &eqnState, const double &theta, const double &zeta){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  //loop over physical cells
  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	double coeff = ( (*this).Vol(locG) * zeta ) / ( (*this).Dt(loc) * theta ) ;
	solDeltaNm1[loc] = coeff * ( (*this).State(locG).ConsVars(eqnState) - solTimeN[loc] );
      }
    }
  }

}


//member function to return a copy of the conserved variables
vector<colMatrix> procBlock::GetCopyConsVars(const idealGas &eqnState) const {

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();

  vector<colMatrix> consVars(imax * jmax * kmax);

  //loop over physical cells
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
      for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
	int locG = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
	int loc = GetLoc1D(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	consVars[loc] = (*this).State(locG).ConsVars(eqnState);
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


  int imax = (*this).NumI();
  int jmax = (*this).NumJ();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();


  //initialize correction (x) to 0
  for (unsigned int ll = 0; ll < x.size(); ll++ ){
    x[ll].Zero();
  }

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

    //forward sweep over all physical cells
    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      //indicies for variables without ghost cells
      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int il = GetNeighborLowI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int jl = GetNeighborLowJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int kl = GetNeighborLowK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      //indicies for variables with ghost cells
      int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      int ilFaceG = GetLowerFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int jlFaceG = GetLowerFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int klFaceG = GetLowerFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      int ilFace2G = GetLowerFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
      int jlFace2G = GetLowerFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
      int klFace2G = GetLowerFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);

      int ilG = GetNeighborLowI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int jlG = GetNeighborLowJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int klG = GetNeighborLowK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);


      if ( il >=0 && il < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaI(ilFace2G), (*this).FAreaI(ilFaceG), (*this).State(ilG).UpdateWithConsVars(eqnState, x[il]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(ilFace2G), (*this).FAreaI(ilFaceG), 
								(*this).State(ilG).UpdateWithConsVars(eqnState, x[il]), eqnState, suth, (*this).Vol(ilG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(ilG), eqnState, (*this).FAreaI(ilFaceG), x[il]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaI(ilFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[il]) );
    }
      if ( jl >=0 && jl < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaJ(jlFace2G), (*this).FAreaJ(jlFaceG), (*this).State(jlG).UpdateWithConsVars(eqnState, x[jl]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(jlFace2G), (*this).FAreaJ(jlFaceG), 
								(*this).State(jlG).UpdateWithConsVars(eqnState, x[jl]), eqnState, suth, (*this).Vol(jlG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(jlG), eqnState, (*this).FAreaJ(jlFaceG), x[jl]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaJ(jlFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[jl]) );
      }
      if ( kl >=0 && kl < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaK(klFace2G), (*this).FAreaK(klFaceG), (*this).State(klG).UpdateWithConsVars(eqnState, x[kl]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(klFace2G), (*this).FAreaK(klFaceG), 
								(*this).State(klG).UpdateWithConsVars(eqnState, x[kl]), eqnState, suth, (*this).Vol(klG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(klG), eqnState, (*this).FAreaK(klFaceG), x[kl]);

	L[loc] = L[loc] + 0.5 * ( (*this).FAreaK(klFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * I.Multiply(x[kl]) );
      }

      double diagTimeVol = ( (*this).Vol(locG) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

      x[loc] = AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) - solDeltaNm1[loc] - solTimeMmN[loc] + L[loc] ) ; //normal at lower boundaries needs to be reversed, so add instead of subtract L

    } //end forward sweep

    //backward sweep over all physical cells
    for ( int ii = (int)x.size()-1; ii >= 0; ii-- ){

      //indicies for variables without ghost cells
      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      int iu = GetNeighborUpI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int ju = GetNeighborUpJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int ku = GetNeighborUpK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

      //indicies for variables with ghost cells
      int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      int iuFaceG = GetUpperFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int juFaceG = GetUpperFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int kuFaceG = GetUpperFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      int iuFace2G = GetUpperFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
      int juFace2G = GetUpperFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
      int kuFace2G = GetUpperFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);

      int iuG = GetNeighborUpI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int juG = GetNeighborUpJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
      int kuG = GetNeighborUpK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      if ( iu >=0 && iu < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaI(iuFace2G), (*this).FAreaI(iuFaceG), (*this).State(iuG).UpdateWithConsVars(eqnState, x[iu]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(iuFace2G), (*this).FAreaI(iuFaceG), 
								(*this).State(iuG).UpdateWithConsVars(eqnState, x[iu]), eqnState, suth, (*this).Vol(iuG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(iuG), eqnState, (*this).FAreaI(iuFaceG), x[iu]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaI(iuFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[iu]) );
      }
      if ( ju >=0 && ju < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaJ(juFace2G), (*this).FAreaJ(juFaceG), (*this).State(juG).UpdateWithConsVars(eqnState, x[ju]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(juFace2G), (*this).FAreaJ(juFaceG), 
								(*this).State(juG).UpdateWithConsVars(eqnState, x[ju]), eqnState, suth, (*this).Vol(juG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(juG), eqnState, (*this).FAreaJ(juFaceG), x[ju]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaJ(juFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[ju]) );
      }
      if ( ku >=0 && ku < (int)x.size() ){
	//at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
	double specRad = CellSpectralRadius( (*this).FAreaK(kuFace2G), (*this).FAreaK(kuFaceG), (*this).State(kuG).UpdateWithConsVars(eqnState, x[ku]), eqnState);

	if (inp.EquationSet() != "euler"){ //viscous
	  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	  double mRef = inp.VelRef().Mag() / aRef;
	  double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(kuFace2G), (*this).FAreaK(kuFaceG), 
								(*this).State(kuG).UpdateWithConsVars(eqnState, x[ku]), eqnState, suth, (*this).Vol(kuG) );
	  specRad += vSpecRad;
	}

	//at given face location, call function to calculate convective flux change
	colMatrix fluxChange = ConvectiveFluxUpdate( (*this).State(kuG), eqnState, (*this).FAreaK(kuFaceG), x[ku]);

	U[loc] = U[loc] + 0.5 * ( (*this).FAreaK(kuFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * I.Multiply(x[ku]) );
      }

      double diagTimeVol = ( (*this).Vol(locG) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

      x[loc] = x[loc] - AiiInv * U[loc];

    } //end backward sweep


    //calculate residual
    colMatrix resid(x[0].Size());

    for ( int ii = 0; ii < (int)x.size(); ii++ ){

      int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
      int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

      double diagTimeVol = ( (*this).Vol(locG) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
      if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
	double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
	diagTimeVol += tau;
      }

      double Aii = ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() ;
      
      //normal at lower boundaries needs to be reversed, so add instead of subtract L
      resid = -1.0 * thetaInv * (*this).Residual(loc) + solDeltaNm1[loc] + solTimeMmN[loc] - Aii * x[loc] + L[loc] - U[loc];
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

  //loop over physical cells
  for ( int kk = numGhosts; kk < numK + numGhosts; kk++ ){
    for ( int jj = numGhosts; jj < numJ + numGhosts; jj++ ){
      for ( int ii = numGhosts; ii < numI + numGhosts; ii++ ){

	int newLoc = GetLoc1D(ii, jj, kk, newI, newJ);
	int loc = GetLoc1D(ii-numGhosts, jj-numGhosts, kk-numGhosts, numI, numJ);
	padBlk[newLoc] = var[loc];

      }
    }
  }

  return padBlk;
}

//function to pad a vector with a specified number of ghost cells
vector<primVars> PadStateWithGhosts( const primVars &var, const int &numGhosts, const int &numI, const int &numJ, const int &numK ){

  int newI = numI + (numGhosts * 2);
  int newJ = numJ + (numGhosts * 2);
  int newK = numK + (numGhosts * 2);

  int newSize = newI * newJ * newK;

  vector<primVars> padBlk(newSize,var);

  return padBlk;
}

//member function to calculate the velocity gradient at the cell center
tensor<double> procBlock::CalcVelGradGG(const vector3d<double> &vil, const vector3d<double> &viu, const vector3d<double> &vjl, const vector3d<double> &vju, 
					const vector3d<double> &vkl, const vector3d<double> &vku, const vector3d<double> &ail, const vector3d<double> &aiu,
					const vector3d<double> &ajl, const vector3d<double> &aju, const vector3d<double> &akl, const vector3d<double> &aku,
					const double &vol){

  //vil is the velocity vector at the i-lower face of the cell at which the velocity gradient is being calculated
  //viu is the velocity vector at the i-upper face of the cell at which the velocity gradient is being calculated
  //vjl is the velocity vector at the j-lower face of the cell at which the velocity gradient is being calculated
  //vju is the velocity vector at the j-upper face of the cell at which the velocity gradient is being calculated
  //vkl is the velocity vector at the k-lower face of the cell at which the velocity gradient is being calculated
  //vku is the velocity vector at the k-upper face of the cell at which the velocity gradient is being calculated

  //ail is the area vector at the lower i-face of the cell at which the velocity gradient is being calculated
  //aiu is the area vector at the upper i-face of the cell at which the velocity gradient is being calculated
  //ajl is the area vector at the lower j-face of the cell at which the velocity gradient is being calculated
  //aju is the area vector at the upper j-face of the cell at which the velocity gradient is being calculated
  //akl is the area vector at the lower k-face of the cell at which the velocity gradient is being calculated
  //aku is the area vector at the upper k-face of the cell at which the velocity gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the velocity gradient should be stored

  tensor<double> temp;
  double invVol = 1.0/vol;

  //define velocity gradient tensor
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetXX( invVol * ( viu.X()*aiu.X() - vil.X()*ail.X() + vju.X()*aju.X() - vjl.X()*ajl.X() + vku.X()*aku.X() - vkl.X()*akl.X() ) );
  temp.SetXY( invVol * ( viu.Y()*aiu.X() - vil.Y()*ail.X() + vju.Y()*aju.X() - vjl.Y()*ajl.X() + vku.Y()*aku.X() - vkl.Y()*akl.X() ) );
  temp.SetXZ( invVol * ( viu.Z()*aiu.X() - vil.Z()*ail.X() + vju.Z()*aju.X() - vjl.Z()*ajl.X() + vku.Z()*aku.X() - vkl.Z()*akl.X() ) );

  temp.SetYX( invVol * ( viu.X()*aiu.Y() - vil.X()*ail.Y() + vju.X()*aju.Y() - vjl.X()*ajl.Y() + vku.X()*aku.Y() - vkl.X()*akl.Y() ) );
  temp.SetYY( invVol * ( viu.Y()*aiu.Y() - vil.Y()*ail.Y() + vju.Y()*aju.Y() - vjl.Y()*ajl.Y() + vku.Y()*aku.Y() - vkl.Y()*akl.Y() ) );
  temp.SetYZ( invVol * ( viu.Z()*aiu.Y() - vil.Z()*ail.Y() + vju.Z()*aju.Y() - vjl.Z()*ajl.Y() + vku.Z()*aku.Y() - vkl.Z()*akl.Y() ) );

  temp.SetZX( invVol * ( viu.X()*aiu.Z() - vil.X()*ail.Z() + vju.X()*aju.Z() - vjl.X()*ajl.Z() + vku.X()*aku.Z() - vkl.X()*akl.Z() ) );
  temp.SetZY( invVol * ( viu.Y()*aiu.Z() - vil.Y()*ail.Z() + vju.Y()*aju.Z() - vjl.Y()*ajl.Z() + vku.Y()*aku.Z() - vkl.Y()*akl.Z() ) );
  temp.SetZZ( invVol * ( viu.Z()*aiu.Z() - vil.Z()*ail.Z() + vju.Z()*aju.Z() - vjl.Z()*ajl.Z() + vku.Z()*aku.Z() - vkl.Z()*akl.Z() ) );

  return temp;

}

//member function to calculate the temperature gradient at the cell center
vector3d<double> procBlock::CalcTempGradGG(const double &til, const double &tiu, const double &tjl, const double &tju, const double &tkl, const double &tku,
					   const vector3d<double> &ail, const vector3d<double> &aiu, const vector3d<double> &ajl, const vector3d<double> &aju,
					   const vector3d<double> &akl, const vector3d<double> &aku, const double &vol){

  //til is the temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tiu is the temperature at the upper face of the cell at which the temperature gradient is being calculated
  //tjl is the temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tju is the temperature at the upper face of the cell at which the temperature gradient is being calculated
  //tkl is the temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tku is the temperature at the upper face of the cell at which the temperature gradient is being calculated
  
  //ail is the area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aiu is the area vector at the upper face of the cell at which the temperature gradient is being calculated
  //ajl is the area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aju is the area vector at the upper face of the cell at which the temperature gradient is being calculated
  //akl is the area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aku is the area vector at the upper face of the cell at which the temperature gradient is being calculated

  //vol is the cell volume
  //loc is the 1D location where the temperature gradient should be stored

  vector3d<double> temp;
  double invVol = 1.0/vol;

  //define temperature gradient vector
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetX( invVol * ( tiu*aiu.X() - til*ail.X() + tju*aju.X() - tjl*ajl.X() + tku*aku.X() - tkl*akl.X() ) );
  temp.SetY( invVol * ( tiu*aiu.Y() - til*ail.Y() + tju*aju.Y() - tjl*ajl.Y() + tku*aku.Y() - tkl*akl.Y() ) );
  temp.SetZ( invVol * ( tiu*aiu.Z() - til*ail.Z() + tju*aju.Z() - tjl*ajl.Z() + tku*aku.Z() - tkl*akl.Z() ) );

  return temp;

}

//member function to calculate viscous fluxes on i-faces
void procBlock::CalcViscFluxI(const sutherland &suth, const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI() + 1;
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts() + 1;
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  //loop over all physical cells
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//face indices
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	int fUpi = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int fLowi = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);

	int fUpjUpi = GetUpperFaceJ(ii, jj, kk, imaxG - 1, jmaxG);
	int fUpjLowi = GetUpperFaceJ(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int fLowjUpi = GetLowerFaceJ(ii, jj, kk, imaxG - 1, jmaxG);
	int fLowjLowi = GetLowerFaceJ(ii - 1, jj, kk, imaxG - 1, jmaxG);

	int fUpkUpi = GetUpperFaceK(ii, jj, kk, imaxG - 1, jmaxG);
	int fUpkLowi = GetUpperFaceK(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int fLowkUpi = GetLowerFaceK(ii, jj, kk, imaxG - 1, jmaxG);
	int fLowkLowi = GetLowerFaceK(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//cell indices
	int iLow  = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG);
	int iUp  = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG);

	int jUpiUp = GetNeighborUpJ(ii, jj, kk, imaxG - 1, jmaxG);
	int jUpiLow = GetNeighborUpJ(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int jLowiUp = GetNeighborLowJ(ii, jj, kk, imaxG - 1, jmaxG);
	int jLowiLow = GetNeighborLowJ(ii - 1, jj, kk, imaxG - 1, jmaxG);

	int kUpiUp = GetNeighborUpK(ii, jj, kk, imaxG - 1, jmaxG);
	int kUpiLow = GetNeighborUpK(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int kLowiUp = GetNeighborLowK(ii, jj, kk, imaxG - 1, jmaxG);
	int kLowiLow = GetNeighborLowK(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//no ghost cell indices
	int iLowNG  = GetCellFromFaceLowerI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int iUpNG  = GetCellFromFaceUpperI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	vector3d<double> vju = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(jUpiUp).Velocity() + (*this).State(jUpiLow).Velocity() );
	vector3d<double> vjl = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(jLowiUp).Velocity() + (*this).State(jLowiLow).Velocity() );

	vector3d<double> vku = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(kUpiUp).Velocity() + (*this).State(kUpiLow).Velocity() );
	vector3d<double> vkl = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(kLowiUp).Velocity() + (*this).State(kLowiLow).Velocity() );

	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(loc) + (*this).FAreaI(fUpi) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(loc) + (*this).FAreaI(fLowi) );

	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(fUpjUpi) + (*this).FAreaJ(fUpjLowi) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(fLowjUpi) + (*this).FAreaJ(fLowjLowi) );

	vector3d<double> aku = 0.5 * ( (*this).FAreaK(fUpkUpi) + (*this).FAreaK(fUpkLowi) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(fLowkUpi) + (*this).FAreaK(fLowkLowi) );

	double vol = 0.5 * ( (*this).Vol(iLow) + (*this).Vol(iUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( (*this).State(iLow).Velocity(), (*this).State(iUp).Velocity(), vjl, vju, vkl, vku, ail, aiu, ajl, aju, akl, aku, vol);
	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(iLow).Velocity(), (*this).State(iUp).Velocity(), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );

	double tju = 0.25 * ( (*this).State(iLow).Temperature(eqnState) + (*this).State(iUp).Temperature(eqnState) + (*this).State(jUpiUp).Temperature(eqnState) +
			      (*this).State(jUpiLow).Temperature(eqnState) );
	double tjl = 0.25 * ( (*this).State(iLow).Temperature(eqnState) + (*this).State(iUp).Temperature(eqnState) + (*this).State(jLowiUp).Temperature(eqnState) +
			      (*this).State(jLowiLow).Temperature(eqnState) );

	double tku = 0.25 * ( (*this).State(iLow).Temperature(eqnState) + (*this).State(iUp).Temperature(eqnState) + (*this).State(kUpiUp).Temperature(eqnState) +
			      (*this).State(kUpiLow).Temperature(eqnState) );
	double tkl = 0.25 * ( (*this).State(iLow).Temperature(eqnState) + (*this).State(iUp).Temperature(eqnState) + (*this).State(kLowiUp).Temperature(eqnState) +
			      (*this).State(kLowiLow).Temperature(eqnState) );

	//Get temperature gradient at face
	vector3d<double> tGrad = CalcTempGradGG( (*this).State(iLow).Temperature(eqnState), (*this).State(iUp).Temperature(eqnState), tjl, tju, tkl, tku, ail, aiu, ajl, aju, akl, aku, vol);
	//Get viscosity at face
	double mu = FaceReconCentral( suth.GetViscosity( (*this).State(iLow).Temperature(eqnState) ), 
				      suth.GetViscosity( (*this).State(iUp).Temperature(eqnState) ), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );
	mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

	//calculate viscous flux
	viscousFlux tempViscFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaI(loc) );

	if (jj == (*this).NumGhosts() && ii == 32 && kk == 2){
	  cout << "vel i lower: " << (*this).State(iLow).Velocity() << endl;
	  cout << "vel i upper: " << (*this).State(iUp).Velocity() << endl;
	  cout << "vel j lower: " << vjl << endl;
	  cout << "vel j upper: " << vju << endl;
	  cout << "vel k lower: " << vkl << endl;
	  cout << "vel k upper: " << vku << endl;
	  cout << "face vel: " << vel << endl;
	  cout << "area i lower: " << ail << endl;
	  cout << "area i upper: " << aiu << endl;
	  cout << "area j lower: " << ajl << endl;
	  cout << "area j upper: " << aju << endl;
	  cout << "area k lower: " << akl << endl;
	  cout << "area k upper: " << aku << endl;
	  cout << "inverse vol: " << 1.0/vol << endl;

	  vector3d<double> viu = (*this).State(iUp).Velocity();
	  vector3d<double> vil = (*this).State(iLow).Velocity();

	  cout << 1.0/vol * ( viu * aiu.X() - vil * ail.X() )<< endl;
	  cout << 1.0/vol * ( vju * aju.X() - vjl * ajl.X() )<< endl;
	  cout << 1.0/vol * ( vku * aku.X() - vkl * akl.X() )<< endl;

	  cout << 1.0/vol * ( viu.X()*aiu.X() - vil.X()*ail.X() ) << endl;
	  cout << 1.0/vol * ( vju.X()*aju.X() - vjl.X()*ajl.X() ) << endl;

	  cout << 1.0/vol * ( viu.X()*aiu.X() ) << ", " << 1.0/vol * (-vil.X()*ail.X()) << ", " << 1.0/vol * ( vju.X()*aju.X()) << ", " 
	       << 1.0/vol *(-vjl.X()*ajl.X()) << ", " << 1.0/vol * ( vku.X()*aku.X()) << ", " << 1.0/vol * (-vkl.X()*akl.X() ) << endl;

	  cout << "vel gradient at boundary face: " << velGrad << endl;
	  cout << "temp gradient at boundary face: " << tGrad << endl;
	}

	// if ( ii == 32 && jj == 2 && kk == 2){
	//   cout << "velocity at j faces: " << vjl << ", " << vju << endl;
	//   cout << "vflux on i face: " << tempViscFlux << endl;

	// }


	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( ii > (*this).NumGhosts() ){
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaI(loc).Mag(), iLowNG);
	}
	if ( ii < imax -1 + (*this).NumGhosts() ){
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaI(loc).Mag(), iUpNG);
	}

      }
    }
  }


}



//member function to calculate viscous fluxes on j-faces
void procBlock::CalcViscFluxJ(const sutherland &suth, const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ() + 1;
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts() + 1;

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  //loop over physical cells
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//face indices
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	int fUpj = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG);
	int fLowj = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG);

	int fUpiUpj = GetUpperFaceI(ii, jj, kk, imaxG, jmaxG - 1);
	int fUpiLowj = GetUpperFaceI(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int fLowiUpj = GetLowerFaceI(ii, jj, kk, imaxG, jmaxG - 1);
	int fLowiLowj = GetLowerFaceI(ii, jj - 1, kk, imaxG, jmaxG - 1);

	int fUpkUpj = GetUpperFaceK(ii, jj, kk, imaxG, jmaxG - 1);
	int fUpkLowj = GetUpperFaceK(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int fLowkUpj = GetLowerFaceK(ii, jj, kk, imaxG, jmaxG - 1);
	int fLowkLowj = GetLowerFaceK(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//cell indices
	int jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG);
	int jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG);

	int iUpjUp = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG - 1);
	int iUpjLow = GetNeighborUpI(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int iLowjUp = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG - 1);
	int iLowjLow = GetNeighborLowI(ii, jj - 1, kk, imaxG, jmaxG - 1);

	int kUpjUp = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG - 1);
	int kUpjLow = GetNeighborUpK(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int kLowjUp = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG - 1);
	int kLowjLow = GetNeighborLowK(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//no ghost cell indices
	int jLowNG  = GetCellFromFaceLowerJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int jUpNG  = GetCellFromFaceUpperJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	vector3d<double> viu = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(iUpjUp).Velocity() + (*this).State(iUpjLow).Velocity() );
	vector3d<double> vil = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(iLowjUp).Velocity() + (*this).State(iLowjLow).Velocity() );

	vector3d<double> vku = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(kUpjUp).Velocity() + (*this).State(kUpjLow).Velocity() );
	vector3d<double> vkl = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(kLowjUp).Velocity() + (*this).State(kLowjLow).Velocity() );

	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(loc) + (*this).FAreaJ(fUpj) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(loc) + (*this).FAreaJ(fLowj) );

	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(fUpiUpj) + (*this).FAreaI(fUpiLowj) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(fLowiUpj) + (*this).FAreaI(fLowiLowj) );

	vector3d<double> aku = 0.5 * ( (*this).FAreaK(fUpkUpj) + (*this).FAreaK(fUpkLowj) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(fLowkUpj) + (*this).FAreaK(fLowkLowj) );

	double vol = 0.5 * ( (*this).Vol(jLow) + (*this).Vol(jUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( vil, viu, (*this).State(jLow).Velocity(), (*this).State(jUp).Velocity(), vkl, vku, ail, aiu, ajl, aju, akl, aku, vol);

	// if ( jj == 2 ){
	//   vector3d<double> z(0.0, 0.0, 0.0);
	//   velGrad = CalcVelGradGG( vil, viu, z, (*this).State(jUp).Velocity(), vkl, vku, ail, aiu, ajl, aju, akl, aku, vol);
	//   }

	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(jLow).Velocity(), (*this).State(jUp).Velocity(), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );

	double tiu = 0.25 * ( (*this).State(jLow).Temperature(eqnState) + (*this).State(jUp).Temperature(eqnState) + (*this).State(iUpjUp).Temperature(eqnState) +
			      (*this).State(iUpjLow).Temperature(eqnState) );
	double til = 0.25 * ( (*this).State(jLow).Temperature(eqnState) + (*this).State(jUp).Temperature(eqnState) + (*this).State(iLowjUp).Temperature(eqnState) +
			      (*this).State(iLowjLow).Temperature(eqnState) );

	double tku = 0.25 * ( (*this).State(jLow).Temperature(eqnState) + (*this).State(jUp).Temperature(eqnState) + (*this).State(kUpjUp).Temperature(eqnState) +
			      (*this).State(kUpjLow).Temperature(eqnState) );
	double tkl = 0.25 * ( (*this).State(jLow).Temperature(eqnState) + (*this).State(jUp).Temperature(eqnState) + (*this).State(kLowjUp).Temperature(eqnState) +
			      (*this).State(kLowjLow).Temperature(eqnState) );

	//Get temperature gradient at face
	vector3d<double> tGrad = CalcTempGradGG( til, tiu, (*this).State(jLow).Temperature(eqnState), (*this).State(jUp).Temperature(eqnState), tkl, tku, ail, aiu, ajl, aju, akl, aku, vol);
	//Get viscosity at face
	double mu = FaceReconCentral( suth.GetViscosity( (*this).State(jLow).Temperature(eqnState) ), 
				      suth.GetViscosity( (*this).State(jUp).Temperature(eqnState) ), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );
	mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)


	// if (jj == (*this).NumGhosts() && ii == 32 && kk == 2){
	//   cout << "vel i lower: " << vil << endl;
	//   cout << "vel i upper: " << viu << endl;
	//   cout << "vel j lower: " << (*this).State(jLow).Velocity() << endl;
	//   cout << "vel j upper: " << (*this).State(jUp).Velocity() << endl;
	//   cout << "vel k lower: " << vkl << endl;
	//   cout << "vel k upper: " << vku << endl;
	//   cout << "face vel: " << vel << endl;
	//   cout << "area i lower: " << ail << endl;
	//   cout << "area i upper: " << aiu << endl;
	//   cout << "area j lower: " << ajl << endl;
	//   cout << "area j upper: " << aju << endl;
	//   cout << "area k lower: " << akl << endl;
	//   cout << "area k upper: " << aku << endl;
	//   cout << "inverse vol: " << 1.0/vol << endl;

	//   // cout << 1.0/vol * ( viu * aiu.Y() - vil * ail.Y() )<< endl;
	//   // cout << 1.0/vol * ( (*this).State(jUp).Velocity() * aju.Y() - (*this).State(jLow).Velocity() * ajl.Y() )<< endl;
	//   // cout << 1.0/vol * ( vku * aku.Y() - vkl * akl.Y() )<< endl;

	//   cout << "vel gradient at boundary face: " << velGrad << endl;
	//   cout << "temp gradient at boundary face: " << tGrad << endl;
	// }

	//calculate viscous flux
	viscousFlux tempViscFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaJ(loc) );

	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( jj < (*this).NumGhosts() ){
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaJ(loc).Mag(), jLowNG);
	}
	if ( jj < jmax -1 + (*this).NumGhosts() ){
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaJ(loc).Mag(), jUpNG);
	}

      }
    }
  }


}


//member function to calculate viscous fluxes on j-faces
void procBlock::CalcViscFluxK(const sutherland &suth, const idealGas &eqnState, const input &inp){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() + 1;

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//face indices
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	int fUpk = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int fLowk = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);

	int fUpiUpk = GetUpperFaceI(ii, jj, kk, imaxG, jmaxG);
	int fUpiLowk = GetUpperFaceI(ii, jj, kk - 1, imaxG, jmaxG);
	int fLowiUpk = GetLowerFaceI(ii, jj, kk, imaxG, jmaxG);
	int fLowiLowk = GetLowerFaceI(ii, jj, kk - 1, imaxG, jmaxG);

	int fUpjUpk = GetUpperFaceJ(ii, jj, kk, imaxG, jmaxG);
	int fUpjLowk = GetUpperFaceJ(ii, jj, kk - 1, imaxG, jmaxG);
	int fLowjUpk = GetLowerFaceJ(ii, jj, kk, imaxG, jmaxG);
	int fLowjLowk = GetLowerFaceJ(ii, jj, kk - 1, imaxG, jmaxG);

	//cell indices
	int kLow  = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG);
	int kUp  = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG);

	int iUpkUp = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int iUpkLow = GetNeighborUpI(ii, jj, kk - 1, imaxG, jmaxG);
	int iLowkUp = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);
	int iLowkLow = GetNeighborLowI(ii, jj, kk - 1, imaxG, jmaxG);

	int jUpkUp = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int jUpkLow = GetNeighborUpK(ii, jj, kk - 1, imaxG, jmaxG);
	int jLowkUp = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);
	int jLowkLow = GetNeighborLowK(ii, jj, kk - 1, imaxG, jmaxG);

	//no ghost indices
	int kLowNG  = GetCellFromFaceLowerK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int kUpNG  = GetCellFromFaceUpperK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	vector3d<double> viu = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(iUpkUp).Velocity() + (*this).State(iUpkLow).Velocity() );
	vector3d<double> vil = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(iLowkUp).Velocity() + (*this).State(iLowkLow).Velocity() );

	vector3d<double> vju = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(jUpkUp).Velocity() + (*this).State(jUpkLow).Velocity() );
	vector3d<double> vjl = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(jLowkUp).Velocity() + (*this).State(jLowkLow).Velocity() );

	vector3d<double> aku = 0.5 * ( (*this).FAreaK(loc) + (*this).FAreaK(fUpk) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(loc) + (*this).FAreaK(fLowk) );

	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(fUpiUpk) + (*this).FAreaI(fUpiLowk) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(fLowiUpk) + (*this).FAreaI(fLowiLowk) );

	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(fUpjUpk) + (*this).FAreaJ(fUpjLowk) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(fLowjUpk) + (*this).FAreaJ(fLowjLowk) );

	double vol = 0.5 * ( (*this).Vol(kLow) + (*this).Vol(kUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( vil, viu, vjl, vju, (*this).State(kLow).Velocity(), (*this).State(kUp).Velocity(), ail, aiu, ajl, aju, akl, aku, vol);
	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(kLow).Velocity(), (*this).State(kUp).Velocity(), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	//Get temperature gradient at face

	double tiu = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(iUpkUp).Temperature(eqnState) +
			      (*this).State(iUpkLow).Temperature(eqnState) );
	double til = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(iLowkUp).Temperature(eqnState) +
			      (*this).State(iLowkLow).Temperature(eqnState) );

	double tju = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(jUpkUp).Temperature(eqnState) +
			      (*this).State(jUpkLow).Temperature(eqnState) );
	double tjl = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(jLowkUp).Temperature(eqnState) +
			      (*this).State(jLowkLow).Temperature(eqnState) );

	vector3d<double> tGrad = CalcTempGradGG( til, tiu, tjl, tju, (*this).State(kLow).Temperature(eqnState), (*this).State(kUp).Temperature(eqnState), ail, aiu, ajl, aju, akl, aku, vol);
	//Get viscosity at face
	double mu = FaceReconCentral( suth.GetViscosity( (*this).State(kLow).Temperature(eqnState) ), 
				      suth.GetViscosity( (*this).State(kUp).Temperature(eqnState) ), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

	//calculate viscous flux
	viscousFlux tempViscFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaK(loc) );

	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( kk > (*this).NumGhosts() ){
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaK(loc).Mag(), kLowNG);
	}
	if ( kk < kmax -1 + (*this).NumGhosts() ){
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaK(loc).Mag(), kUpNG);
	}

      }
    }
  }


}

//member function to assign geometric quantities such as volume, face area, etc to ghost cells
void procBlock::AssignGhostCellsGeom(){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //loop over physical I faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
      int lFaceG1_il = GetLowerFaceI(1, jj, kk, imaxG, jmaxG); 
      int lFaceG1_jl = GetLowerFaceJ(1, jj, kk, imaxG, jmaxG); 
      int lFaceG1_kl = GetLowerFaceK(1, jj, kk, imaxG, jmaxG); 

      int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);
      int lFaceG2_il = GetLowerFaceI(0, jj, kk, imaxG, jmaxG);
      int lFaceG2_jl = GetLowerFaceJ(0, jj, kk, imaxG, jmaxG);
      int lFaceG2_kl = GetLowerFaceK(0, jj, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int lFaceIn1_iu = GetUpperFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
      int lFaceIn1_jl = GetLowerFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
      int lFaceIn1_kl = GetLowerFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
      int lFaceIn2_iu = GetUpperFaceI((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
      int lFaceIn2_jl = GetLowerFaceJ((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
      int lFaceIn2_kl = GetLowerFaceK((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

      int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
      int uFaceG1_iu = GetUpperFaceI(imaxG-2, jj, kk, imaxG, jmaxG); 
      int uFaceG1_jl = GetLowerFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
      int uFaceG1_kl = GetLowerFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 

      int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);
      int uFaceG2_iu = GetUpperFaceI(imaxG-1, jj, kk, imaxG, jmaxG);
      int uFaceG2_jl = GetLowerFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);
      int uFaceG2_kl = GetLowerFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceIn1_il = GetLowerFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
      int uFaceIn1_jl = GetLowerFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
      int uFaceIn1_kl = GetLowerFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceIn2_il = GetLowerFaceI(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceIn2_jl = GetLowerFaceJ(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceIn2_kl = GetLowerFaceK(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

      int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      //mirror volume values from adjacent cells
      (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);
      (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);
      if (imax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
      }
      else{
	(*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
      }

      //mirror face area values from adjacent cells
      (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG1_jl);
      (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG1_kl);

      (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_iu);
      (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_jl);
      (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_kl);

      if (imax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG2_jl);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG2_kl);

	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_iu);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_jl);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_kl);
      }
      else{
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_jl), lFaceG2_jl);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn2_kl), lFaceG2_kl);

	(*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_iu);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_jl);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_kl);
      }

      if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG); 
	int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

	int lFaceIn1_ju = GetUpperFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int lFaceIn2_ju = GetUpperFaceJ((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	int uFaceG1_ju = GetUpperFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG2_ju = GetUpperFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);

	int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int uFaceIn2_ju = GetUpperFaceJ(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_ju);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG1_ju);

	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_ju);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG2_ju);
	}
	else{
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_ju);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_ju), uFaceG2_ju);
	}

      }

      if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG); 
	int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

	int lFaceIn1_ku = GetUpperFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int lFaceIn2_ku = GetUpperFaceK((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	int uFaceG1_ku = GetUpperFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG2_ku = GetUpperFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

	int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int uFaceIn2_ku = GetUpperFaceK(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_ku);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG1_ku);

	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_ku);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG2_ku);
	}
	else{
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_ku);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_ku), uFaceG2_ku);
	}

      }

      //cell centroid is moved interior cell width in the boundary normal direction
      vector3d<double> dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
      (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);
      dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
      (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

      if (imax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	(*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	(*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
      }
      else{
	dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
      }

      //face center is moved interior cell width in the boundary normal direction
      dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
      (*this).SetFCenterI( (*this).FCenterI(lFaceB) + dist2Move, lFaceG1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG1_jl);
      (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG1_kl);

      dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
      (*this).SetFCenterI( (*this).FCenterI(uFaceB) + dist2Move, uFaceG1_iu);
      (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG1_jl);
      (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG1_kl);

      if (imax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	(*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);
	(*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);

	dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	(*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_jl) + dist2Move, uFaceG2_jl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceG1_kl) + dist2Move, uFaceG2_kl);

      }
      else{
	dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	(*this).SetFCenterI( (*this).FCenterI(lFaceB) + dist2Move, lFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG2_jl);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG2_kl);

	dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	(*this).SetFCenterI( (*this).FCenterI(uFaceB) + dist2Move, uFaceG2_iu);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG2_jl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG2_kl);

      }

      if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG); 
	int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

	int lFaceIn1_ju = GetUpperFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	int uFaceG1_ju = GetUpperFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG2_ju = GetUpperFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);

	int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG1_ju);

	dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG1_ju);

	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_ju) + dist2Move, lFaceG2_ju);

	  dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);

	}
	else{
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG2_ju);

	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG2_ju);
	}

      }

      if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG); 
	int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

	int lFaceIn1_ku = GetUpperFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	int uFaceG1_ku = GetUpperFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG2_ku = GetUpperFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

	int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG1_ku);

	dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG1_ku);

	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	  (*this).SetFCenterK( (*this).FCenterK(lFaceG1_ku) + dist2Move, lFaceG2_ku);

	  dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	  (*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);

	}
	else{
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG2_ku);

	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG2_ku);
	}

      }

    }
  }

  //loop over physical J faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
      int lFaceG1_il = GetLowerFaceI(ii, 1, kk, imaxG, jmaxG); 
      int lFaceG1_jl = GetLowerFaceJ(ii, 1, kk, imaxG, jmaxG); 
      int lFaceG1_kl = GetLowerFaceK(ii, 1, kk, imaxG, jmaxG); 

      int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);
      int lFaceG2_il = GetLowerFaceI(ii, 0, kk, imaxG, jmaxG);
      int lFaceG2_jl = GetLowerFaceJ(ii, 0, kk, imaxG, jmaxG);
      int lFaceG2_kl = GetLowerFaceK(ii, 0, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
      int lFaceIn1_il = GetLowerFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
      int lFaceIn1_ju = GetUpperFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
      int lFaceIn1_kl = GetLowerFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

      int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
      int lFaceIn2_il = GetLowerFaceI(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
      int lFaceIn2_ju = GetUpperFaceJ(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
      int lFaceIn2_kl = GetLowerFaceK(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

      int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
      int uFaceG1_il = GetLowerFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
      int uFaceG1_ju = GetUpperFaceJ(ii, jmaxG-2, kk, imaxG, jmaxG); 
      int uFaceG1_kl = GetLowerFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 

      int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);
      int uFaceG2_il = GetLowerFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);
      int uFaceG2_ju = GetUpperFaceJ(ii, jmaxG-1, kk, imaxG, jmaxG);
      int uFaceG2_kl = GetLowerFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceIn1_il = GetLowerFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
      int uFaceIn1_jl = GetLowerFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
      int uFaceIn1_kl = GetLowerFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

      int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceIn2_il = GetLowerFaceI(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceIn2_jl = GetLowerFaceJ(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceIn2_kl = GetLowerFaceK(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

      int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

      //mirror volume values from adjacent cells
      (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);
      (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);
      if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
      }
      else{
	(*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
      }

      //mirror face area values from adjacent cells
      (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_jl);
      (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG1_il);
      (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG1_kl);

      (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_ju);
      (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_il);
      (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_kl);

      if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_jl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG2_il);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG2_kl);

	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_ju);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_il);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_kl);

      }
      else{
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_jl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn2_il), lFaceG2_il);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn2_kl), lFaceG2_kl);

	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_ju);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_il);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_kl);
      }

      if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper face areas too

	int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

	int lFaceIn1_iu = GetUpperFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int lFaceIn2_iu = GetUpperFaceI(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	int uFaceG1_iu = GetUpperFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG2_iu = GetUpperFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);

	int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int uFaceIn2_iu = GetUpperFaceI(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_iu);

	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG1_iu);

	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_iu);

	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG2_iu);
	}
	else{
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_iu);

	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_iu), uFaceG2_iu);
	}

      }

      if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

	int lFaceIn1_ku = GetUpperFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int lFaceIn2_ku = GetUpperFaceK(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	int uFaceG1_ku = GetUpperFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG2_ku = GetUpperFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

	int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int uFaceIn2_ku = GetUpperFaceK(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_ku);

	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG1_ku);

	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_ku);

	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG2_ku);
	}
	else{
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_ku);

	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_ku), uFaceG2_ku);
	}

      }

      //cell centroid is moved interior cell width in the boundary normal direction
      vector3d<double> dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
      (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);
      dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
      (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

      if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	(*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	(*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
      }
      else{
	dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
      }

      //face center is moved interior cell width in the boundary normal direction
      dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
      (*this).SetFCenterJ( (*this).FCenterJ(lFaceB) + dist2Move, lFaceG1_jl);
      (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG1_il);
      (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG1_kl);

      dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
      (*this).SetFCenterJ( (*this).FCenterJ(uFaceB) + dist2Move, uFaceG1_ju);
      (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG1_il);
      (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG1_kl);

      if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	(*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);

	dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);
	(*this).SetFCenterI( (*this).FCenterI(uFaceG1_il) + dist2Move, uFaceG2_il);
	(*this).SetFCenterK( (*this).FCenterK(uFaceG1_kl) + dist2Move, uFaceG2_kl);
      }
      else{
	dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceB) + dist2Move, lFaceG2_jl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG2_il);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG2_kl);

	dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceB) + dist2Move, uFaceG2_ju);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG2_il);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG2_kl);
      }

      if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper face areas too

	int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

	int lFaceIn1_iu = GetUpperFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	int uFaceG1_iu = GetUpperFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG2_iu = GetUpperFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);

	int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG1_iu);

	dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG1_iu);

	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	  (*this).SetFCenterI( (*this).FCenterI(lFaceG1_iu) + dist2Move, lFaceG2_iu);

	  dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	  (*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	}
	else{
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG2_iu);

	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG2_iu);
	}

      }

      if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

	int lFaceIn1_ku = GetUpperFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	int uFaceG1_ku = GetUpperFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG2_ku = GetUpperFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

	int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG1_ku);

	dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG1_ku);

	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	  (*this).SetFCenterK( (*this).FCenterK(lFaceG1_ku) + dist2Move, lFaceG2_ku);

	  dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	  (*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);
	}
	else{
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG2_ku);

	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG2_ku);
	}

      }

    }
  }

  //loop over physical K faces
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
      int lFaceG1_il = GetLowerFaceI(ii, jj, 1, imaxG, jmaxG); 
      int lFaceG1_jl = GetLowerFaceJ(ii, jj, 1, imaxG, jmaxG); 
      int lFaceG1_kl = GetLowerFaceK(ii, jj, 1, imaxG, jmaxG); 

      int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);
      int lFaceG2_il = GetLowerFaceI(ii, jj, 0, imaxG, jmaxG);
      int lFaceG2_jl = GetLowerFaceJ(ii, jj, 0, imaxG, jmaxG);
      int lFaceG2_kl = GetLowerFaceK(ii, jj, 0, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
      int lFaceIn1_il = GetLowerFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
      int lFaceIn1_jl = GetLowerFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
      int lFaceIn1_ku = GetUpperFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

      int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
      int lFaceIn2_il = GetLowerFaceI(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
      int lFaceIn2_jl = GetLowerFaceJ(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
      int lFaceIn2_ku = GetUpperFaceK(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

      int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
      int uFaceG1_il = GetLowerFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
      int uFaceG1_jl = GetLowerFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
      int uFaceG1_ku = GetUpperFaceK(ii, jj, kmaxG-2, imaxG, jmaxG); 

      int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);
      int uFaceG2_il = GetLowerFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);
      int uFaceG2_jl = GetLowerFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);
      int uFaceG2_ku = GetUpperFaceK(ii, jj, kmaxG-1, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceIn1_il = GetLowerFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
      int uFaceIn1_jl = GetLowerFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
      int uFaceIn1_kl = GetLowerFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

      int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceIn2_il = GetLowerFaceI(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceIn2_jl = GetLowerFaceJ(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceIn2_kl = GetLowerFaceK(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

      int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

      //mirror volume values from adjacent cells
      (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);
      (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);
      if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
      }
      else{
	(*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	(*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
      }

      //mirror face area values from adjacent cells
      (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_kl);
      (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG1_jl);

      (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_ku);
      (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_jl);

      if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_kl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG2_jl);

	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_ku);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_jl);
      }
      else{
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_kl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn2_il), lFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_jl), lFaceG2_jl);

	(*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_ku);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_il);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_jl);
      }

      if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper face areas too

	int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

	int lFaceIn1_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	int lFaceIn2_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);

	int uFaceIn1_iu = GetUpperFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	int uFaceIn2_iu = GetUpperFaceI(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_iu);

	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG1_iu);

	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_iu);

	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG2_iu);
	}
	else{
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_iu);

	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_iu), uFaceG2_iu);
	}

      }

      if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

	int lFaceIn1_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	int lFaceIn2_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);

	int uFaceIn1_ju = GetUpperFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	int uFaceIn2_ju = GetUpperFaceJ(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	//mirror face area values from adjacent cells
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_ju);

	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG1_ju);

	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_ju);

	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG2_ju);
	}
	else{
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_ju);

	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_ju), uFaceG2_ju);
	}

      }

      //cell centroid is moved interior cell width in the boundary normal direction
      vector3d<double> dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
      (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);
      dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
      (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

      if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	(*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	(*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
      }
      else{
	dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
      }

      //face center is moved interior cell width in the boundary normal direction
      dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
      (*this).SetFCenterK( (*this).FCenterK(lFaceB) + dist2Move, lFaceG1_kl);
      (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG1_jl);

      dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
      (*this).SetFCenterK( (*this).FCenterK(uFaceB) + dist2Move, uFaceG1_ku);
      (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG1_jl);

      if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	(*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);

	dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	(*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);
	(*this).SetFCenterI( (*this).FCenterI(uFaceG1_il) + dist2Move, uFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_jl) + dist2Move, uFaceG2_jl);
      }
      else{
	dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	(*this).SetFCenterK( (*this).FCenterK(lFaceB) + dist2Move, lFaceG2_kl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG2_jl);

	dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceB) + dist2Move, uFaceG2_ku);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG2_il);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG2_jl);
      }

      if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper face areas too

	int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

	int lFaceIn1_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);

	int uFaceIn1_iu = GetUpperFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG1_iu);

	dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG1_iu);

	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	  (*this).SetFCenterI( (*this).FCenterI(lFaceG1_iu) + dist2Move, lFaceG2_iu);

	  dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	  (*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	}
	else{
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG2_iu);

	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG2_iu);
	}

      }

      if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

	int lFaceIn1_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);

	int uFaceIn1_ju = GetUpperFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	//face center is moved interior cell width in the boundary normal direction
	dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG1_ju);

	dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG1_ju);

	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_ju) + dist2Move, lFaceG2_ju);

	  dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);
	}
	else{
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG2_ju);

	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG2_ju);
	}

      }

    }
  }

  //fill ghost cell edge lines with geometric values
  (*this).AssignGhostCellsGeomEdge();

  // cout << "geom ghost cells" << endl;
  // for( int kk = 0; kk < kmaxG; kk++ ){
  //   for( int jj = 0; jj < jmaxG; jj++ ){
  //     for( int ii = 0; ii < imaxG; ii++ ){

  // 	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
  // 	cout << ii << ", " << jj << ", " << kk << ", " << (*this).Vol(loc) << endl;

  //     }
  //     cout << endl;
  //   }
  // }


}

//member function to assign geometric quantities such as volume, face area, etc to ghost cells located on block edges
//assumes AssignGhostCellsGeom has already been run
void procBlock::AssignGhostCellsGeomEdge(){

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int j1,k1,j2,k2,je,ke;

      int gfe_j1_k1_jl, gfe_j1_k1_kl;
      int gfe_j1_k2_jl, gfe_j1_k2_kl;
      int gfe_j2_k1_jl, gfe_j2_k1_kl;
      int gfe_j2_k2_jl, gfe_j2_k2_kl;

      int gf_j1_ke_jl, gf_j1_ke_kl, gf_j1_ke_ku;
      int gf_j2_ke_jl, gf_j2_ke_kl;
      int gf_je_k1_jl, gf_je_k1_ju;
      int gf_je_k2_jl, gf_je_k2_kl;

      if ( cc == 0 ){ //at jl/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

	gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

	gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

	gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells

	gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells

	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
	gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
      }
      else if ( cc == 1 ){ //at jl/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

	gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

	gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

	gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells

	gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells

	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
	gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
      }
      else if ( cc == 2 ){ //at ju/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

	gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

	gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

	gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells

	gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells

	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
	gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

	gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

	gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

	gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells

	gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells

	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
	gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells
      }

      int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
      int gfe_j1_k1_il = GetLowerFaceI(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

      int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
      int gfe_j1_k2_il = GetLowerFaceI(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

      int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
      int gfe_j2_k1_il = GetLowerFaceI(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

      int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells
      int gfe_j2_k2_il = GetLowerFaceI(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

      int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, first layer of ghost cells
      int gf_j1_ke_il = GetLowerFaceI(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells

      int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, second layer of ghost cells
      int gf_j2_ke_il = GetLowerFaceI(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

      int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, first layer of ghost cells
      int gf_je_k1_il = GetLowerFaceI(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells

      int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, second layer of ghost cells
      int gf_je_k2_il = GetLowerFaceI(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

      //volume
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_j1_ke) + (*this).Vol(gc_je_k1) ) ,gce_j1_k1);
      (*this).SetVol( (*this).Vol(gc_j2_ke) ,gce_j2_k1);
      (*this).SetVol( (*this).Vol(gc_je_k2) ,gce_j1_k2);
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_j2_ke) + (*this).Vol(gc_je_k2) ) ,gce_j2_k2);

      //face areas
      (*this).SetFAreaI( 0.5 * ( (*this).FAreaI(gf_je_k1_il) + (*this).FAreaI(gf_j1_ke_il) ), gfe_j1_k1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_je_k1_jl), gfe_j1_k1_jl);
      (*this).SetFAreaK( (*this).FAreaK(gf_j1_ke_kl), gfe_j1_k1_kl);

      (*this).SetFAreaI( (*this).FAreaI(gf_je_k2_il), gfe_j1_k2_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_je_k2_jl), gfe_j1_k2_jl);
      (*this).SetFAreaK( (*this).FAreaK(gf_j1_ke_kl), gfe_j1_k2_kl);

      (*this).SetFAreaI( (*this).FAreaI(gf_j2_ke_il), gfe_j2_k1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_je_k1_jl), gfe_j2_k1_jl);
      (*this).SetFAreaK( (*this).FAreaK(gf_j2_ke_kl), gfe_j2_k1_kl);

      (*this).SetFAreaI( 0.5 * ( (*this).FAreaI(gf_j2_ke_il) + (*this).FAreaI(gf_je_k2_il) ), gfe_j2_k2_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_je_k2_jl), gfe_j2_k2_jl);
      (*this).SetFAreaK( (*this).FAreaK(gf_j2_ke_kl), gfe_j2_k2_kl);

      //centroids
      vector3d<double> dist2MoveK = (*this).FCenterK(gf_j1_ke_kl) - (*this).FCenterK(gf_j1_ke_ku) ;
      vector3d<double> dist2MoveJ = (*this).FCenterJ(gf_je_k1_jl) - (*this).FCenterJ(gf_je_k1_ju) ;
      (*this).SetCenter( (*this).Center(gc_j1_ke) + dist2MoveK, gce_j1_k1);
      (*this).SetCenter( (*this).Center(gc_j2_ke) + dist2MoveK, gce_j2_k1);
      (*this).SetCenter( (*this).Center(gc_je_k2) + dist2MoveJ, gce_j1_k2);
      (*this).SetCenter( (*this).Center(gc_je_k2) + 2.0 * dist2MoveJ, gce_j2_k2);

      //face centers
      (*this).SetFCenterI( (*this).FCenterI(gf_j1_ke_il) + dist2MoveK, gfe_j1_k1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_j1_ke_jl) + dist2MoveK, gfe_j1_k1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_j1_ke_kl) + dist2MoveK, gfe_j1_k1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_j2_ke_il) + dist2MoveK, gfe_j2_k1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_j2_ke_jl) + dist2MoveK, gfe_j2_k1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_j2_ke_kl) + dist2MoveK, gfe_j2_k1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_je_k2_il) + dist2MoveJ, gfe_j1_k2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_je_k2_jl) + dist2MoveJ, gfe_j1_k2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_je_k2_kl) + dist2MoveJ, gfe_j1_k2_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_je_k2_il) + 2.0 * dist2MoveJ, gfe_j2_k2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_je_k2_jl) + 2.0 * dist2MoveJ, gfe_j2_k2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_je_k2_kl) + 2.0 * dist2MoveJ, gfe_j2_k2_kl);

      //this is only done at the end of the i loop
      if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper face areas too

	int gfe_j1_k1_2il = GetLowerFaceI(ii - 1, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	vector3d<double> dist2MoveI = (*this).FCenterI(gfe_j1_k1_il) - (*this).FCenterI(gfe_j1_k1_2il);

	int gfe_j1_k1_iu = GetUpperFaceI(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	int gfe_j1_k2_iu = GetUpperFaceI(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	int gfe_j2_k1_iu = GetUpperFaceI(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	int gfe_j2_k2_iu = GetUpperFaceI(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	(*this).SetFAreaI( (*this).FAreaI(gfe_j1_k1_il), gfe_j1_k1_iu);
	(*this).SetFAreaI( (*this).FAreaI(gfe_j1_k2_il), gfe_j1_k2_iu);
	(*this).SetFAreaI( (*this).FAreaI(gfe_j2_k1_il), gfe_j2_k1_iu);
	(*this).SetFAreaI( (*this).FAreaI(gfe_j2_k2_il), gfe_j2_k2_iu);

	(*this).SetFCenterI( (*this).FCenterI(gfe_j1_k1_il) + dist2MoveI, gfe_j1_k1_iu);
	(*this).SetFCenterI( (*this).FCenterI(gfe_j1_k2_il) + dist2MoveI, gfe_j1_k2_iu);
	(*this).SetFCenterI( (*this).FCenterI(gfe_j2_k1_il) + dist2MoveI, gfe_j2_k1_iu);
	(*this).SetFCenterI( (*this).FCenterI(gfe_j2_k2_il) + dist2MoveI, gfe_j2_k2_iu);

      }

    }
  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,k1,i2,k2,ie,ke;

      int gfe_i1_k1_il, gfe_i1_k1_kl;
      int gfe_i1_k2_il, gfe_i1_k2_kl;
      int gfe_i2_k1_il, gfe_i2_k1_kl;
      int gfe_i2_k2_il, gfe_i2_k2_kl;

      int gf_i1_ke_il, gf_i1_ke_ku, gf_i1_ke_kl;
      int gf_i2_ke_il, gf_i2_ke_kl;
      int gf_ie_k1_il, gf_ie_k1_iu;
      int gf_ie_k2_il, gf_ie_k2_kl;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

	gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

	gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

	gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells

	gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

	gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

	gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

	gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells

	gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 2 ){ //at ju/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

	gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

	gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

	gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells

	gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

	gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

	gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

	gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells

	gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells
      }

      int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gfe_i1_k1_jl = GetLowerFaceJ(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

      int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gfe_i1_k2_jl = GetLowerFaceJ(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

      int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gfe_i2_k1_jl = GetLowerFaceJ(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

      int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells
      int gfe_i2_k2_jl = GetLowerFaceJ(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gf_i1_ke_jl = GetLowerFaceJ(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells

      int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gf_i2_ke_jl = GetLowerFaceJ(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

      int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gf_ie_k1_jl = GetLowerFaceJ(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells

      int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells
      int gf_ie_k2_jl = GetLowerFaceJ(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

      //volume
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_i1_ke) + (*this).Vol(gc_ie_k1) ) ,gce_i1_k1);
      (*this).SetVol( (*this).Vol(gc_i2_ke) ,gce_i2_k1);
      (*this).SetVol( (*this).Vol(gc_ie_k2) ,gce_i1_k2);
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_i2_ke) + (*this).Vol(gc_ie_k2) ) ,gce_i2_k2);

      //face areas
      (*this).SetFAreaJ( 0.5 * ( (*this).FAreaJ(gf_ie_k1_jl) + (*this).FAreaJ(gf_i1_ke_jl) ), gfe_i1_k1_jl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_k1_il), gfe_i1_k1_il);
      (*this).SetFAreaK( (*this).FAreaK(gf_i1_ke_kl), gfe_i1_k1_kl);

      (*this).SetFAreaJ( (*this).FAreaJ(gf_ie_k2_jl), gfe_i1_k2_jl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_k2_il), gfe_i1_k2_il);
      (*this).SetFAreaK( (*this).FAreaK(gf_i1_ke_kl), gfe_i1_k2_kl);

      (*this).SetFAreaJ( (*this).FAreaJ(gf_i2_ke_jl), gfe_i2_k1_jl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_k1_il), gfe_i2_k1_il);
      (*this).SetFAreaK( (*this).FAreaK(gf_i2_ke_kl), gfe_i2_k1_kl);

      (*this).SetFAreaJ( 0.5 * ( (*this).FAreaJ(gf_i2_ke_jl) + (*this).FAreaJ(gf_ie_k2_jl) ), gfe_i2_k2_jl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_k2_il), gfe_i2_k2_il);
      (*this).SetFAreaK( (*this).FAreaK(gf_i2_ke_kl), gfe_i2_k2_kl);

      //centroids
      vector3d<double> dist2MoveK = (*this).FCenterK(gf_i1_ke_kl) - (*this).FCenterK(gf_i1_ke_ku);
      vector3d<double> dist2MoveI = (*this).FCenterI(gf_ie_k1_il) - (*this).FCenterI(gf_ie_k1_iu);
      (*this).SetCenter( (*this).Center(gc_i1_ke) + dist2MoveK, gce_i1_k1);
      (*this).SetCenter( (*this).Center(gc_i2_ke) + dist2MoveK, gce_i2_k1);
      (*this).SetCenter( (*this).Center(gc_ie_k2) + dist2MoveI, gce_i1_k2);
      (*this).SetCenter( (*this).Center(gc_ie_k2) + 2.0 * dist2MoveI, gce_i2_k2);

      //face centers
      (*this).SetFCenterI( (*this).FCenterI(gf_i1_ke_il) + dist2MoveK, gfe_i1_k1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_i1_ke_jl) + dist2MoveK, gfe_i1_k1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_i1_ke_kl) + dist2MoveK, gfe_i1_k1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_i2_ke_il) + dist2MoveK, gfe_i2_k1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_i2_ke_jl) + dist2MoveK, gfe_i2_k1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_i2_ke_kl) + dist2MoveK, gfe_i2_k1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_ie_k2_il) + dist2MoveI, gfe_i1_k2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_ie_k2_jl) + dist2MoveI, gfe_i1_k2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_ie_k2_kl) + dist2MoveI, gfe_i1_k2_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_ie_k2_il) + 2.0 * dist2MoveI, gfe_i2_k2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_ie_k2_jl) + 2.0 * dist2MoveI, gfe_i2_k2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_ie_k2_kl) + 2.0 * dist2MoveI, gfe_i2_k2_kl);

      //this is only done at the end of the j loop
      if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	int gfe_i1_k1_2jl = GetLowerFaceJ(i1, jj - 1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	vector3d<double> dist2MoveJ = (*this).FCenterJ(gfe_i1_k1_jl) - (*this).FCenterJ(gfe_i1_k1_2jl);

	int gfe_i1_k1_ju = GetUpperFaceJ(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	int gfe_i1_k2_ju = GetUpperFaceJ(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	int gfe_i2_k1_ju = GetUpperFaceJ(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	int gfe_i2_k2_ju = GetUpperFaceJ(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	(*this).SetFAreaJ( (*this).FAreaJ(gfe_i1_k1_jl), gfe_i1_k1_ju);
	(*this).SetFAreaJ( (*this).FAreaJ(gfe_i1_k2_jl), gfe_i1_k2_ju);
	(*this).SetFAreaJ( (*this).FAreaJ(gfe_i2_k1_jl), gfe_i2_k1_ju);
	(*this).SetFAreaJ( (*this).FAreaJ(gfe_i2_k2_jl), gfe_i2_k2_ju);

	(*this).SetFCenterJ( (*this).FCenterJ(gfe_i1_k1_jl) + dist2MoveJ, gfe_i1_k1_ju);
	(*this).SetFCenterJ( (*this).FCenterJ(gfe_i1_k2_jl) + dist2MoveJ, gfe_i1_k2_ju);
	(*this).SetFCenterJ( (*this).FCenterJ(gfe_i2_k1_jl) + dist2MoveJ, gfe_i2_k1_ju);
	(*this).SetFCenterJ( (*this).FCenterJ(gfe_i2_k2_jl) + dist2MoveJ, gfe_i2_k2_ju);

      }

    }
  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,j1,i2,j2,ie,je;

      int gfe_i1_j1_il, gfe_i1_j1_jl;
      int gfe_i1_j2_il, gfe_i1_j2_jl;
      int gfe_i2_j1_il, gfe_i2_j1_jl;
      int gfe_i2_j2_il, gfe_i2_j2_jl;

      int gf_i1_je_il, gf_i1_je_jl, gf_i1_je_ju;
      int gf_i2_je_il, gf_i2_je_jl;
      int gf_ie_j1_il, gf_ie_j1_iu;
      int gf_ie_j2_il, gf_ie_j2_jl;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

	gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

	gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

	gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells

	gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

	gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

	gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

	gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells

	gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 2 ){ //at ju/jl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

	gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

	gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

	gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells

	gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
      }
      else if ( cc == 3 ){ //at ju/ju edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

	gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

	gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

	gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells
	gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells

	gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells

	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
	gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells
      }

      int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of j line of cells
      int gfe_i1_j1_kl = GetLowerFaceK(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

      int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of j line of cells
      int gfe_i1_j2_kl = GetLowerFaceK(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

      int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of j line of cells
      int gfe_i2_j1_kl = GetLowerFaceK(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

      int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of j line of cells
      int gfe_i2_j2_kl = GetLowerFaceK(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

      int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at j-lower end, first layer of ghost cells
      int gf_i1_je_kl = GetLowerFaceK(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, first layer of ghost cells

      int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at j-lower end, second layer of ghost cells
      int gf_i2_je_kl = GetLowerFaceK(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at j-lower end, second layer of ghost cells

      int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at i-lower end, first layer of ghost cells
      int gf_ie_j1_kl = GetLowerFaceK(ie, j1, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, first layer of ghost cells

      int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at i-lower end, second layer of ghost cells
      int gf_ie_j2_kl = GetLowerFaceK(ie, j2, kk, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at i-lower end, second layer of ghost cells

      //volume
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_i1_je) + (*this).Vol(gc_ie_j1) ) ,gce_i1_j1);
      (*this).SetVol( (*this).Vol(gc_i2_je) ,gce_i2_j1);
      (*this).SetVol( (*this).Vol(gc_ie_j2) ,gce_i1_j2);
      (*this).SetVol( 0.5 * ( (*this).Vol(gc_i2_je) + (*this).Vol(gc_ie_j2) ) ,gce_i2_j2);

      //face areas
      (*this).SetFAreaK( 0.5 * ( (*this).FAreaK(gf_ie_j1_kl) + (*this).FAreaK(gf_i1_je_kl) ), gfe_i1_j1_kl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_j1_il), gfe_i1_j1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_i1_je_jl), gfe_i1_j1_jl);

      (*this).SetFAreaK( (*this).FAreaK(gf_ie_j2_kl), gfe_i1_j2_kl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_j2_il), gfe_i1_j2_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_i1_je_jl), gfe_i1_j2_jl);

      (*this).SetFAreaK( (*this).FAreaK(gf_i2_je_kl), gfe_i2_j1_kl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_j1_il), gfe_i2_j1_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_i2_je_jl), gfe_i2_j1_jl);

      (*this).SetFAreaK( 0.5 * ( (*this).FAreaK(gf_i2_je_kl) + (*this).FAreaJ(gf_ie_j2_kl) ), gfe_i2_j2_kl);
      (*this).SetFAreaI( (*this).FAreaI(gf_ie_j2_il), gfe_i2_j2_il);
      (*this).SetFAreaJ( (*this).FAreaJ(gf_i2_je_jl), gfe_i2_j2_jl);

      //centroids
      vector3d<double> dist2MoveJ = (*this).FCenterJ(gf_i1_je_jl) - (*this).FCenterJ(gf_i1_je_ju);
      vector3d<double> dist2MoveI = (*this).FCenterI(gf_ie_j1_il) - (*this).FCenterI(gf_ie_j1_iu);
      (*this).SetCenter( (*this).Center(gc_i1_je) + dist2MoveJ, gce_i1_j1);
      (*this).SetCenter( (*this).Center(gc_i2_je) + dist2MoveJ, gce_i2_j1);
      (*this).SetCenter( (*this).Center(gc_ie_j2) + dist2MoveI, gce_i1_j2);
      (*this).SetCenter( (*this).Center(gc_ie_j2) + 2.0 * dist2MoveI, gce_i2_j2);

      //face centers
      (*this).SetFCenterI( (*this).FCenterI(gf_i1_je_il) + dist2MoveJ, gfe_i1_j1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_i1_je_jl) + dist2MoveJ, gfe_i1_j1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_i1_je_kl) + dist2MoveJ, gfe_i1_j1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_i2_je_il) + dist2MoveJ, gfe_i2_j1_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_i2_je_jl) + dist2MoveJ, gfe_i2_j1_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_i2_je_kl) + dist2MoveJ, gfe_i2_j1_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_ie_j2_il) + dist2MoveI, gfe_i1_j2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_ie_j2_jl) + dist2MoveI, gfe_i1_j2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_ie_j2_kl) + dist2MoveI, gfe_i1_j2_kl);

      (*this).SetFCenterI( (*this).FCenterI(gf_ie_j2_il) + 2.0 * dist2MoveI, gfe_i2_j2_il);
      (*this).SetFCenterJ( (*this).FCenterJ(gf_ie_j2_jl) + 2.0 * dist2MoveI, gfe_i2_j2_jl);
      (*this).SetFCenterK( (*this).FCenterK(gf_ie_j2_kl) + 2.0 * dist2MoveI, gfe_i2_j2_kl);

      //this is only done at the end of the k loop
      if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	int gfe_i1_j1_2kl = GetLowerFaceK(i1, j1, kk - 1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	vector3d<double> dist2MoveK = (*this).FCenterK(gfe_i1_j1_kl) - (*this).FCenterK(gfe_i1_j1_2kl);

	int gfe_i1_j1_ku = GetUpperFaceK(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	int gfe_i1_j2_ku = GetUpperFaceK(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	int gfe_i2_j1_ku = GetUpperFaceK(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	int gfe_i2_j2_ku = GetUpperFaceK(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	(*this).SetFAreaK( (*this).FAreaK(gfe_i1_j1_kl), gfe_i1_j1_ku);
	(*this).SetFAreaK( (*this).FAreaK(gfe_i1_j2_kl), gfe_i1_j2_ku);
	(*this).SetFAreaK( (*this).FAreaK(gfe_i2_j1_kl), gfe_i2_j1_ku);
	(*this).SetFAreaK( (*this).FAreaK(gfe_i2_j2_kl), gfe_i2_j2_ku);

	(*this).SetFCenterK( (*this).FCenterK(gfe_i1_j1_kl) + dist2MoveK, gfe_i1_j1_ku);
	(*this).SetFCenterK( (*this).FCenterK(gfe_i1_j2_kl) + dist2MoveK, gfe_i1_j2_ku);
	(*this).SetFCenterK( (*this).FCenterK(gfe_i2_j1_kl) + dist2MoveK, gfe_i2_j1_ku);
	(*this).SetFCenterK( (*this).FCenterK(gfe_i2_j2_kl) + dist2MoveK, gfe_i2_j2_ku);

      }

    }
  }

}

//member function to assign ghost cells for the inviscid flow calculation
void procBlock::AssignInviscidGhostCells(const input &inp, const idealGas &eos){

  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();
  int kmaxG = (*this).NumK() + 2 * (*this).NumGhosts();

  //loop over physical I faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
      int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(0, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "il");
      string bcNameU = bound.GetBCName(imax, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "iu");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
	bcNameU = "slipWall";
      }

      (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG1);
      (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG1);

      if (imax < 2){ //one cell thick - use one cell for both ghost cells
	(*this).SetState( (*this).State(cellLowG1), cellLowG2);
	(*this).SetState( (*this).State(cellUpG1), cellUpG2);
      }
      else{
	if (bcNameL == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 2), cellLowG2);
	}

	if (bcNameU == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 2), cellUpG2);
	}
      }

    }
  }

  //loop over physical J faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
      int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), 0, kk - (*this).NumGhosts(), "jl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jmax, kk - (*this).NumGhosts(), "ju");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
      	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
      	bcNameU = "slipWall";
      }

      (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG1);
      (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG1);

      if (jmax < 2){ //one cell thick - use once cell for both ghost cells
	(*this).SetState( (*this).State(cellLowG1), cellLowG2);
	(*this).SetState( (*this).State(cellUpG1), cellUpG2);
      }
      else{
	if (bcNameL == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 2), cellLowG2);
	}

	if (bcNameU == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 2), cellUpG2);
	}
      }

    }
  }

  //loop over physical K faces
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
      int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), 0, "kl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kmax, "ku");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
	bcNameU = "slipWall";
      }

      (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG1);
      (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG1);

      if (kmax < 2){ //one cell thick - use once cell for both ghost cells
	(*this).SetState( (*this).State(cellLowG1), cellLowG2);
	(*this).SetState( (*this).State(cellUpG1), cellUpG2);
      }
      else{
	if (bcNameL == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 2), cellLowG2);
	}

	if (bcNameU == "slipWall"){ //if slipWall, reflect second interior state over boundary face
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 2), cellUpG2);
	}
      }

    }
  }

  //assign values to edge ghost cells
  (*this).AssignInviscidGhostCellsEdge(inp, eos);

  // cout << "inviscid ghost cells" << endl;
  // for( int kk = 0; kk < kmaxG; kk++ ){
  //   for( int jj = 0; jj < jmaxG; jj++ ){
  //     for( int ii = 0; ii < imaxG; ii++ ){

  // 	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
  // 	cout << ii << ", " << jj << ", " << kk << ", " << (*this).State(loc) << endl;

  //     }
  //     cout << endl;
  //   }
  // }


}

//member function to assign inviscid ghost cells to cells on edge
//assumes AssignInviscidGhostCells has already been run
void procBlock::AssignInviscidGhostCellsEdge(const input &inp, const idealGas &eos){

  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int j1,k1,j2,k2,je,ke,je2,ke2;

      int gf_j1_ke_kl;
      int gf_j2_ke_kl;
      int gf_je_k1_jl;
      int gf_je_k2_jl;

      string bc_jl,bc_kl;
      string surfJ,surfK;

      if ( cc == 0 ){ //at jl/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if (jmax > 1){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if (kmax > 1){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "jl";
	surfK = "kl";

	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at jl/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "jl";
	surfK = "ku";

	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at ju/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "ju";
	surfK = "kl";

	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(),surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "ju";
	surfK = "ku";

	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
      int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
      int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
      int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells

      int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, first layer of ghost cells
      int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, second layer of ghost cells

      int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, first layer of ghost cells
      int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, second layer of ghost cells

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bc_jl == "viscousWall" ){
	bc_jl = "slipWall";
      }
      if ( bc_kl == "viscousWall" ){
	bc_kl = "slipWall";
      }

      if ( bc_jl == "slipWall" && !(bc_kl == "slipWall") ){  //j surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_je_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_je_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_je2_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je2_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j2_k2);
      }
      else if ( !(bc_jl == "slipWall") && bc_kl == "slipWall" ){  //k surface is a wall, but j surface is not
	(*this).SetState( (*this).State(gc_j1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_j1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_j2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k2);
      }
      else{ // both surfaces are walls or neither are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_j1_ke) + (*this).State(gc_je_k1) ) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je_k2) ,gce_j1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_j2_ke) + (*this).State(gc_je_k2) ) ,gce_j2_k2);
      }

    }
  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,k1,i2,k2,ie,ke,ie2,ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_il,bc_kl;
      string surfI,surfK;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "il";
	surfK = "kl";

	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ke;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "il";
	surfK = "ku";

	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at iu/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ke;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "iu";
	surfK = "kl";

	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at iu/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "iu";
	surfK = "ku";

	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bc_il == "viscousWall" ){
	bc_il = "slipWall";
      }
      if ( bc_kl == "viscousWall" ){
	bc_kl = "slipWall";
      }

      if ( bc_il == "slipWall" && !(bc_kl == "slipWall") ){  //i surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_ie_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_ie_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_ie2_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie2_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i2_k2);
      }
      else if ( !(bc_il == "slipWall") && bc_kl == "slipWall" ){  //k surface is a wall, but i surface is not
	(*this).SetState( (*this).State(gc_i1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_i1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_i2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k2);

	//cout << "inviscid corner " << i1 << ", " << jj << ", " << ke2 << ", " << (*this).State(gc_i1_ke2) << (*this).State(gce_i1_k2) << endl;

      }
      else{ // both surfaces are walls or neither are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_ke) + (*this).State(gc_ie_k1) ) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie_k2) ,gce_i1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_ke) + (*this).State(gc_ie_k2) ) ,gce_i2_k2);
      }

    }

  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,j1,i2,j2,ie,je,ie2,je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_il,bc_jl;
      string surfI,surfJ;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	surfI = "il";
	surfJ = "jl";

	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}
	surfI = "il";
	surfJ = "ju";

	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 2 ){ //at iu/jl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	surfI = "iu";
	surfJ = "jl";

	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 3 ){ //at iu/ju edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	surfI = "iu";
	surfJ = "ju";

	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }

      int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bc_il == "viscousWall" ){
	bc_il = "slipWall";
      }
      if ( bc_jl == "viscousWall" ){
	bc_jl = "slipWall";
      }

      if ( bc_il == "slipWall" && !(bc_jl == "slipWall") ){  //i surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_ie_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_ie_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_ie2_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie2_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i2_j2);
      }
      else if ( !(bc_il == "slipWall") && bc_jl == "slipWall" ){  //k surface is a wall, but i surface is not
	(*this).SetState( (*this).State(gc_i1_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_i1_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_i2_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j2);
      }
      else{ // both surfaces are walls or neither are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_je) + (*this).State(gc_ie_j1) ) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie_j2) ,gce_i1_j2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_je) + (*this).State(gc_ie_j2) ) ,gce_i2_j2);
      }

    }

  }

}


//member function to assign ghost cells for the viscous flow calculation
//this function assumes AssignInviscidGhostCells has been run first as it only overwrites the ghost cells associated with
//the viscousWall boundary condition
void procBlock::AssignViscousGhostCells(const input &inp, const idealGas &eos){

  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();
  int kmaxG = (*this).NumK() + 2 * (*this).NumGhosts();

  //loop over physical I faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
      int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(0, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "il");
      if (bcNameL == "viscousWall"){
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG1);
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG2);
	}
      }

      string bcNameU = bound.GetBCName(imax, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "iu");
      if (bcNameU == "viscousWall"){
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG1);
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG2);
	}
      }

    }
  }

  //loop over physical J faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
      int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), 0, kk - (*this).NumGhosts(), "jl");
      if (bcNameL == "viscousWall"){
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG1);
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG2);
	}
      }

      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jmax, kk - (*this).NumGhosts(), "ju");
      if (bcNameU == "viscousWall"){
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG1);
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG2);
	}
      }

    }
  }

  //loop over physical K faces
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

      int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
      int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

      int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);

      int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
      int uFaceB = GetUpperFaceK(ii, jj, kmax - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), 0, "kl");
      if (bcNameL == "viscousWall"){
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG1);
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG2);
	}
      }

      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kmax, "ku");
      if (bcNameU == "viscousWall"){
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG1);
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG2);
	}
      }

    }
  }

  //Assign edge ghost cells
  (*this).AssignViscousGhostCellsEdge(inp, eos);

  // cout << "viscous ghost cells" << endl;
  // for( int kk = 0; kk < kmaxG; kk++ ){
  //   for( int jj = 0; jj < jmaxG; jj++ ){
  //     for( int ii = 0; ii < imaxG; ii++ ){

  // 	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
  // 	cout << ii << ", " << jj << ", " << kk << ", " << (*this).State(loc) ;

  //     }
  //     cout << endl;
  //   }
  // }


}

//member function to assign inviscid ghost cells to cells on edge
//assumes AssignInviscidGhostCells has already been run
void procBlock::AssignViscousGhostCellsEdge(const input &inp, const idealGas &eos){

  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int j1,k1,j2,k2,je,ke,je2,ke2;

      int gf_j1_ke_kl;
      int gf_j2_ke_kl;
      int gf_je_k1_jl;
      int gf_je_k2_jl;

      string bc_jl,bc_kl;
      string surfJ,surfK;

      if ( cc == 0 ){ //at jl/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "jl";
	surfK = "kl";

	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at jl/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "jl";
	surfK = "ku";

	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at ju/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "ju";
	surfK = "kl";

	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfJ = "ju";
	surfK = "ku";

	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, first layer of ghost cells
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  //ghost face, on j-lower line of cells, at k-lower end, second layer of ghost cells

	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, first layer of ghost cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at j-lower end, second layer of ghost cells

	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
      int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
      int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
      int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells

      int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, first layer of ghost cells
      int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, second layer of ghost cells

      int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG, jmaxG);       //ghost cell, on j-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, first layer of ghost cells
      int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at j-lower end, second layer of ghost cells

      if ( bc_jl == "viscousWall" && !(bc_kl == "viscousWall") ){  //j surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_je_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_je_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_je2_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je2_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j2_k2);
      }
      else if ( !(bc_jl == "viscousWall") && bc_kl == "viscousWall" ){  //k surface is a wall, but j surface is not
	(*this).SetState( (*this).State(gc_j1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_j1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_j2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k2);
      }
      else if ( bc_jl == "viscousWall" && bc_kl == "viscousWall" ){ // both surfaces are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_j1_ke) + (*this).State(gc_je_k1) ) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je_k2) ,gce_j1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_j2_ke) + (*this).State(gc_je_k2) ) ,gce_j2_k2);
      }
      //if no boundary is a viscous wall, do nothing

    }
  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,k1,i2,k2,ie,ke,ie2,ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_il,bc_kl;
      string surfI,surfK;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "il";
	surfK = "kl";

	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "il";
	surfK = "ku";

	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at iu/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = (*this).NumGhosts() + 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "iu";
	surfK = "kl";

	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();
	if ( kmax > 1 ){
	  ke2 = kmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ke2 = ke;
	}

	surfI = "iu";
	surfK = "ku";

	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      if ( bc_il == "viscousWall" && !(bc_kl == "viscousWall") ){  //i surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_ie_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_ie_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_ie2_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie2_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i2_k2);
      }
      else if ( !(bc_il == "viscousWall") && bc_kl == "viscousWall" ){  //k surface is a wall, but i surface is not
	(*this).SetState( (*this).State(gc_i1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_i1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_i2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k2);
      }
      else if ( bc_il == "viscousWall" && bc_kl == "viscousWall" ){ // both surfaces are walls or neither are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_ke) + (*this).State(gc_ie_k1) ) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie_k2) ,gce_i1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_ke) + (*this).State(gc_ie_k2) ) ,gce_i2_k2);
      }
      //if neither surface is a wall then do nothing

    }

  }

  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){

      int i1,j1,i2,j2,ie,je,ie2,je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_il,bc_jl;
      string surfI,surfJ;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	surfI = "il";
	surfJ = "jl";

	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = (*this).NumGhosts() + 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	surfI = "il";
	surfJ = "ju";

	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je  - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 2 ){ //at iu/jl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = (*this).NumGhosts() + 1;
	}
	else{
	  je2 = je;
	}

	surfI = "iu";
	surfJ = "jl";

	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 3 ){ //at iu/ju edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();
	if ( imax > 1 ){
	  ie2 = imax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  ie2 = ie;
	}

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();
	if ( jmax > 1 ){
	  je2 = jmax - 1 + (*this).NumGhosts() - 1;
	}
	else{
	  je2 = je;
	}

	surfI = "iu";
	surfJ = "ju";

	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, first layer of ghost cells
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  //ghost face, on i-lower line of cells, at k-lower end, second layer of ghost cells

	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, first layer of ghost cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  //ghost face, on k-lower line of cells, at i-lower end, second layer of ghost cells

	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }

      int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, first layer of ghost cells
      int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG, jmaxG);       //ghost cell, on i-lower line of cells, at k-lower end, second layer of ghost cells
      int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, first layer of ghost cells
      int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG, jmaxG);       //ghost cell, on k-lower line of cells, at i-lower end, second layer of ghost cells

      if ( bc_il == "viscousWall" && !(bc_jl == "viscousWall") ){  //i surface is a wall, but k surface is not
	(*this).SetState( (*this).State(gc_ie_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_ie_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_ie2_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie2_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i2_j2);
      }
      else if ( !(bc_il == "viscousWall") && bc_jl == "viscousWall" ){  //j surface is a wall, but i surface is not
	(*this).SetState( (*this).State(gc_i1_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_i1_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_i2_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j2);
      }
      else if ( bc_il == "viscousWall" && bc_jl == "viscousWall"){ // both surfaces are walls or neither are walls
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_je) + (*this).State(gc_ie_j1) ) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie_j2) ,gce_i1_j2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_je) + (*this).State(gc_ie_j2) ) ,gce_i2_j2);
      }
      //if neither surface is a wall then do nothing

    }

  }

}

bool procBlock::AtCorner(const int& ii, const int& jj, const int& kk)const{

  bool atCorner;

  if ( ( ii < (*this).NumGhosts() || ii > (*this).NumI() - 1 + (*this).NumGhosts() ) &&
       ( jj < (*this).NumGhosts() || jj > (*this).NumJ() - 1 + (*this).NumGhosts() ) &&
       ( kk < (*this).NumGhosts() || kk > (*this).NumK() - 1 + (*this).NumGhosts() ) ){
    atCorner = true;
  }
  else{
    atCorner = false;
  }

  return atCorner;
}

bool procBlock::AtEdge(const int& ii, const int& jj, const int& kk)const{

  bool atEdge;

  if ( ( ii >= (*this).NumGhosts()     && ii <  (*this).NumI() + (*this).NumGhosts() ) &&  //at i-edge - i in physical cell range, j/k at first level of ghost cells
       ( jj == (*this).NumGhosts() - 1 || jj == (*this).NumJ() + (*this).NumGhosts() ) &&
       ( kk == (*this).NumGhosts() - 1 || kk == (*this).NumK() + (*this).NumGhosts() ) ){
    atEdge = true;
  }
  else if ( ( ii == (*this).NumGhosts() - 1 || ii == (*this).NumI() + (*this).NumGhosts() ) && //at j-edge - j in physical cell range, i/k at first level of ghost cells
	    ( jj >= (*this).NumGhosts()     && jj <  (*this).NumJ() + (*this).NumGhosts() ) &&
	    ( kk == (*this).NumGhosts() - 1 || kk == (*this).NumK() + (*this).NumGhosts() ) ){
    atEdge = true;
  }
  else if ( ( ii == (*this).NumGhosts() - 1 || ii == (*this).NumI() + (*this).NumGhosts() ) && // at k-edge - k in physical cell range, i/j at first level of ghost cells
	    ( jj == (*this).NumGhosts() - 1 || jj == (*this).NumJ() + (*this).NumGhosts() ) &&
	    ( kk >= (*this).NumGhosts()     && kk <  (*this).NumK() + (*this).NumGhosts() ) ){
    atEdge = true;
  }
  else{
    atEdge = false;
  }

  return atEdge;
}

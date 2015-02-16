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

//constructors for procBlock class
procBlock::procBlock(){
  numCells = 1;
  numVars = NUMVARS;
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
  rank = 0;
  globalPos = 0;

  int numFaces = (numI+1)*(numJ)*(numK);  

  vector<primVars> dummyState (numCells);              //dummy state variable
  vector<vector3d<double> > vec1(numFaces);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vec2(numCells);             //dummy vector variable lenght of number of cells
  vector<double> scalar(numCells);                      //dummy scalar variable lenght of number of cells
  genArray singleResid(0.0);
  vector<genArray> dummyResid(numCells, singleResid);

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

  boundaryConditions bound;
  bc = bound;

}
//constructor -- initialize state vector with dummy variables
procBlock::procBlock(const plot3dBlock &blk, const int& numBlk, const int &numG, const string &eqnSet){
  // blk -- plot3d block of which this procBlock is a subset of
  // numBlk -- the block number of blk (the parent block)
  // numG -- number of ghost cells
  // eqnSet -- which equation set is being solved

  numI = blk.NumI()-1; //i, j, k dimensions are cell-based so subtract 1
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;
  //parent block start/end are face/node based (subtract 1 from blk.Num b/c start at 0)
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;
  rank = 0;
  globalPos = 0;

  numVars = NUMVARS;

  vector<primVars> dummyState (numCells);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  genArray singleResid(0.0);
  vector<genArray> dummyResid(numCells, singleResid);

  //pad stored variable vectors with ghost cells
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

  boundaryConditions bound;
  bc = bound;

}

//constructor -- assign passed variables to initialize state vector
procBlock::procBlock( const double density, const double pressure, const vector3d<double> vel, const plot3dBlock &blk, const int &numBlk, const int &numG, const string &eqnSet,
		      const boundaryConditions &bound){
  // density -- density to initialize block with
  // pressure -- pressure to initialize block with
  // vel -- velocity to initialize block with
  // blk -- plot3d block of which this procBlock is a subset of
  // numBlk -- the block number of blk (the parent block)
  // numG -- number of ghost cells
  // eqnSet -- which equation set is being solved
  // bound -- boundary conditions for block

  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;
  //parent block start/end are face/node based (subtract 1 from blk.Num b/c start at 0)
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;
  rank = 0;
  globalPos = 0;

  bc = bound;

  numVars = NUMVARS;

  primVars singleState(density, pressure, vel);
  vector<primVars> dummyState (numCells,singleState);              //dummy state variable
  vector<double> dummyScalar (numCells);                 //dummy time variable
  genArray singleResid(0.0);
  vector<genArray> dummyResid(numCells, singleResid);

  //pad stored variable vectors with ghost cells
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
procBlock::procBlock( const primVars& inputState, const plot3dBlock &blk, const int &numBlk, const int &numG, const string &eqnSet, const boundaryConditions& bound){
  // inputState -- state to initialize block with (primative)
  // blk -- plot3d block of which this procBlock is a subset of
  // numBlk -- the block number of blk (the parent block)
  // numG -- number of ghost cells
  // eqnSet -- which equation set is being solved
  // bound -- boundary conditions for block

  numI = blk.NumI()-1;
  numJ = blk.NumJ()-1;
  numK = blk.NumK()-1;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = numBlk;

  //parent block start/end are face/node based (subtract 1 from blk.Num b/c start at 0)
  parBlockStartI = 0;
  parBlockEndI = numI;
  parBlockStartJ = 0;
  parBlockEndJ = numJ;
  parBlockStartK = 0;
  parBlockEndK = numK;
  rank = 0;
  globalPos = 0;

  bc = bound;

  numVars = NUMVARS;

  vector<double> dummyScalar (numCells);                 //dummy time variable
  vector<primVars> dummyState (numCells,inputState);              //dummy state variable
  genArray singleResid(0.0);
  vector<genArray> dummyResid(numCells, singleResid);

  //pad stored variable vectors with ghost cells
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

//constructor -- allocate space for procBlock
procBlock::procBlock( const int &ni, const int &nj, const int &nk, const int &numG ){
  // ni -- i-dimension (cell)
  // nj -- j-dimension (cell)
  // nk -- k-dimension (cell)
  // numG -- number of ghost cell layers

  numI = ni;
  numJ = nj;
  numK = nk;
  numCells = numI * numJ * numK;
  numGhosts = numG;
  parBlock = 0;

  //parent block start/end are face/node based
  parBlockStartI = 0;
  parBlockEndI = numI+1;
  parBlockStartJ = 0;
  parBlockEndJ = numJ+1;
  parBlockStartK = 0;
  parBlockEndK = numK+1;
  rank = 0;
  globalPos = 0;

  boundaryConditions bound;
  bc = bound;

  numVars = NUMVARS;
  primVars inputState;

  vector<double> dummyScalar (numCells);                                //dummy scalar variable
  vector<vector3d<double> > dummyVectorCell (numCells);                 //dummy vector3d cell variable
  vector<vector3d<double> > dummyVectorFaceI ((numI+1) * numJ * numK);  //dummy vector3d i-face variable
  vector<vector3d<double> > dummyVectorFaceJ (numI * (numJ+1) * numK);  //dummy vector3d j-face variable
  vector<vector3d<double> > dummyVectorFaceK (numI * numJ * (numK+1));  //dummy vector3d k-face variable
  vector<primVars> dummyState (numCells,inputState);                    //dummy state variable
  genArray singleResid(0.0);
  vector<genArray> dummyResid(numCells, singleResid);

  //pad stored variable vectors with ghost cells
  state = PadWithGhosts( dummyState, numGhosts, numI, numJ, numK );      

  vol = PadWithGhosts( dummyScalar, numGhosts, numI, numJ, numK );
  center = PadWithGhosts( dummyVectorCell, numGhosts, numI, numJ, numK );
  fAreaI = PadWithGhosts( dummyVectorFaceI, numGhosts, numI+1, numJ, numK );
  fAreaJ = PadWithGhosts( dummyVectorFaceJ, numGhosts, numI, numJ+1, numK );
  fAreaK = PadWithGhosts( dummyVectorFaceK, numGhosts, numI, numJ, numK+1 );
  fCenterI = PadWithGhosts( dummyVectorFaceI, numGhosts, numI+1, numJ, numK );
  fCenterJ = PadWithGhosts( dummyVectorFaceJ, numGhosts, numI, numJ+1, numK );
  fCenterK = PadWithGhosts( dummyVectorFaceK, numGhosts, numI, numJ, numK+1 );

  avgWaveSpeed = dummyScalar;
  dt = dummyScalar;
  residual = dummyResid;
}

//member function to add a member of the inviscid flux class to the residual
void procBlock::AddToResidual(const inviscidFlux &flux, const int &ii){
  // flux -- inviscid flux to add to residual
  // ii -- location of residual to add to

  genArray temp(0.0);
  temp[0] = flux.RhoVel();
  temp[1] = flux.RhoVelU();
  temp[2] = flux.RhoVelV();
  temp[3] = flux.RhoVelW();
  temp[4] = flux.RhoVelH();

  (*this).SetResidual( (*this).Residual(ii) + temp, ii); 
}

//member function to add a member of the viscous flux class to the residual
void procBlock::AddToResidual(const viscousFlux &flux, const int &ii){
  // flux -- inviscid flux to add to residual
  // ii -- location of residual to add to

  genArray temp(0.0);
  temp[0] = 0.0;
  temp[1] = flux.MomX();
  temp[2] = flux.MomY();
  temp[3] = flux.MomZ();
  temp[4] = flux.Engy();

  (*this).SetResidual( (*this).Residual(ii) + temp, ii); 
}

//---------------------------------------------------------------------------------------------------------------//
//function declarations

/* Function to calculate the inviscid fluxes on the i-faces. All phyiscal (non-ghost) i-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Ui       -->   Ui+1   |
  |            |            |
  |____________|____________|
Ui-1/2       Ui+1/2       Ui+3/2

Using the above diagram, the flux is calculated at face Ui+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Ui to Ui+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Ui, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Ui+1. 

The spectral radius in the i-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Ui+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.
*/
void procBlock::CalcInvFluxI(const idealGas &eqnState, const input &inp){
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI() + 1; //calculating fluxes on i-faces so one more face in i-direction
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts() + 1; //calculating fluxes on i-faces so one more face in i-direction
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double maxWS = 0.0;
  primVars faceStateLower,faceStateUpper;

  //loop over all physical faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG); 

	//location of lower and upper i-faces (without ghost cells included)
	int lowerING = GetCellFromFaceLowerI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax); 
	int upperING = GetCellFromFaceUpperI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//location of cells in the lower and upper i-direction from current face (with ghost cells included)
	int lowerI = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG);
	int upperI = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG);
	int lower2I = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG, 2);
	int upper2I = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG, 2);

	//location of facess in the lower and upper i-direction from current face (with ghost cells included)
	int upFaceI = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int upFace2I = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceI = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);
	int lowFace2I = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG, 2);

	if (inp.Kappa() == -2.0){  //if value is still default, use constant reconstruction (first order)
	  faceStateLower = (*this).State( lowerI ).FaceReconConst();
	  faceStateUpper = (*this).State( upperI ).FaceReconConst();
	}
	else{ //second order accuracy -- use MUSCL extrapolation

	  //length of second upwind, first upwind, and downwind cells in i-direction
	  double upwind2L =  (*this).FCenterI( lowFaceI ).Distance( (*this).FCenterI( lowFace2I ) );
	  double upwindL =   (*this).FCenterI( loc      ).Distance( (*this).FCenterI( lowFaceI ) );
	  double downwindL = (*this).FCenterI( loc      ).Distance( (*this).FCenterI( upFaceI ) );

	  faceStateLower = (*this).State( lowerI ).FaceReconMUSCL( (*this).State( lower2I ), (*this).State( upperI ),
								   inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  //length of second upwind, first upwind, and downwind cells in i-direction
	  double upwind2U =  (*this).FCenterI( upFaceI ).Distance( (*this).FCenterI( upFace2I ) );
	  double upwindU =   (*this).FCenterI( loc     ).Distance( (*this).FCenterI( upFaceI ) );
	  double downwindU = (*this).FCenterI( loc     ).Distance( (*this).FCenterI( lowFaceI ) );

	  faceStateUpper = (*this).State( upperI ).FaceReconMUSCL( (*this).State( upper2I ), (*this).State( lowerI ),
								   inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );

	}

	//calculate Roe flux at face
	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaI(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( ii > (*this).NumGhosts() ){ //at left boundary there is no left cell to add to
	  (*this).AddToResidual( tempFlux * (*this).FAreaI(loc).Mag(), lowerING);
	}
	if ( ii < imax - 1 + (*this).NumGhosts() ){ //at right boundary there is no right cell to add to
	  (*this).AddToResidual( -1.0 * tempFlux * (*this).FAreaI(loc).Mag(), upperING);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  maxWS = CellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(upFaceI), (*this).State(upperI), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperING) + maxWS, upperING);
	}

      }
    }
  }

}

/* Function to calculate the inviscid fluxes on the j-faces. All phyiscal (non-ghost) j-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Uj       -->   Uj+1   |
  |            |            |
  |____________|____________|
Uj-1/2       Uj+1/2       Uj+3/2

Using the above diagram, the flux is calculated at face Uj+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Uj to Uj+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Uj, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Uj+1. 

The spectral radius in the j-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Uj+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.
*/
void procBlock::CalcInvFluxJ(const idealGas &eqnState, const input &inp){
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ() + 1; //calculating fluxes on i-faces so one more face in j-direction
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts() + 1; //calculating fluxes on i-faces so one more face in j-direction

  double maxWS = 0.0;
  primVars faceStateLower, faceStateUpper;

  //loop over all physical faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	//location of lower and upper j-faces (without ghost cells included)
	int lowerJNG = GetCellFromFaceLowerJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int upperJNG = GetCellFromFaceUpperJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//location of cells in the lower and upper j-direction from current face (with ghost cells included)
	int lowerJ = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG);
	int upperJ = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG);
	int lower2J = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG, 2);
	int upper2J = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG, 2);

	//location of facess in the lower and upper j-direction from current face (with ghost cells included)
	int upFaceJ = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG);
	int upFace2J = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceJ = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG);
	int lowFace2J = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG, 2);

	if ( inp.Kappa() == -2.0 ){ //if value is still default, use constant reconstruction (first order)
	  faceStateLower = (*this).State( lowerJ ).FaceReconConst();
	  faceStateUpper = (*this).State( upperJ ).FaceReconConst();
	}
	else{ //second order accuracy -- use MUSCL extrapolation

	  //length of second upwind, first upwind, and downwind cells in j-direction
	  double upwind2L =  (*this).FCenterJ( lowFaceJ ).Distance( (*this).FCenterJ( lowFace2J ) );
	  double upwindL =   (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( lowFaceJ ) );
	  double downwindL = (*this).FCenterJ( loc      ).Distance( (*this).FCenterJ( upFaceJ ) );

	  faceStateLower = (*this).State( lowerJ ).FaceReconMUSCL( (*this).State( lower2J ),
								   (*this).State( upperJ ), inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  //length of second upwind, first upwind, and downwind cells in j-direction
	  double upwind2U =  (*this).FCenterJ( upFaceJ ).Distance( (*this).FCenterJ( upFace2J ) );
	  double upwindU =   (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( upFaceJ ) );
	  double downwindU = (*this).FCenterJ( loc     ).Distance( (*this).FCenterJ( lowFaceJ ) );

	  faceStateUpper = (*this).State( upperJ ).FaceReconMUSCL( (*this).State( upper2J ),
								   (*this).State( lowerJ ), inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );
	}

	//calculate Roe flux at face
	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaJ(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( jj > (*this).NumGhosts() ){ //at left boundary no left cell to add to
	  (*this).AddToResidual( tempFlux * (*this).FAreaJ(loc).Mag(), lowerJNG);
	}
	if ( jj < jmax - 1 + (*this).NumGhosts() ){ //at right boundary no right cell to add to
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaJ(loc).Mag(), upperJNG);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  maxWS = CellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(upFaceJ), (*this).State(upperJ), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperJNG) + maxWS, upperJNG);
	}

      }
    }
  }

}

/* Function to calculate the inviscid fluxes on the k-faces. All phyiscal (non-ghost) k-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Uk       -->   Uk+1   |
  |            |            |
  |____________|____________|
Uk-1/2       Uk+1/2       Uk+3/2

Using the above diagram, the flux is calculated at face Uk+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Uk to Uk+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Uk, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Uk+1. 

The spectral radius in the k-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Uk+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.
*/
void procBlock::CalcInvFluxK(const idealGas &eqnState, const input &inp){
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() + 1; //calculating fluxes on i-faces so one more face in k-direction

  //max dimensions for vectors without ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  double maxWS = 0.0;
  primVars faceStateLower, faceStateUpper;

  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	//location of lower and upper k-faces (without ghost cells included)
	int lowerKNG = GetCellFromFaceLowerK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int upperKNG = GetCellFromFaceUpperK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//location of cells in the lower and upper k-direction from current face (with ghost cells included)
	int lowerK = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG);
	int upper2K = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG, 2);
	int lower2K = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG, 2);
	int upperK = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG);

	//location of faces in the lower and upper k-direction from current face (with ghost cells included)
	int upFaceK = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int upFace2K = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG, 2);
	int lowFaceK = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);
	int lowFace2K = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG, 2);

	if ( inp.Kappa() == -2.0 ){  //if value is still default, use constant reconstruction (first order)
	  faceStateLower = (*this).State( lowerK ).FaceReconConst();
	  faceStateUpper = (*this).State( upperK ).FaceReconConst();
	}
	else{ //second order accuracy -- use MUSCL extrapolation

	  //length of second upwind, first upwind, and downwind cells in j-direction
	  double upwind2L =  (*this).FCenterK( lowFaceK ).Distance( (*this).FCenterK( lowFace2K ) );
	  double upwindL =   (*this).FCenterK( loc      ).Distance( (*this).FCenterK( lowFaceK ) );
	  double downwindL = (*this).FCenterK( loc      ).Distance( (*this).FCenterK( upFaceK ) );

	  faceStateLower = (*this).State( lowerK ).FaceReconMUSCL( (*this).State( lower2K ), (*this).State( upperK ),
								   inp.Kappa(), inp.Limiter(), upwindL, upwind2L, downwindL );

	  //length of second upwind, first upwind, and downwind cells in j-direction
	  double upwind2U =  (*this).FCenterK( upFaceK ).Distance( (*this).FCenterK( upFace2K ) );
	  double upwindU =   (*this).FCenterK( loc     ).Distance( (*this).FCenterK( upFaceK ) );
	  double downwindU = (*this).FCenterK( loc     ).Distance( (*this).FCenterK( lowFaceK ) );

	  faceStateUpper = (*this).State( upperK ).FaceReconMUSCL( (*this).State( upper2K ), (*this).State( lowerK ),
								   inp.Kappa(), inp.Limiter(), upwindU, upwind2U, downwindU );

	}

	//calculate Roe flux at face
	inviscidFlux tempFlux = RoeFlux(faceStateLower, faceStateUpper, eqnState, (*this).FAreaK(loc), maxWS);

	//area vector points from left to right, so add to left cell, subtract from right cell
	if ( kk > (*this).NumGhosts() ){ //at left boundary no left cell to add to
	  (*this).AddToResidual( tempFlux * (*this).FAreaK(loc).Mag(), lowerKNG);
	}
	if ( kk < kmax - 1 + (*this).NumGhosts() ){ //at right boundary no right cell to add to
	  (*this).AddToResidual(-1.0 * tempFlux * (*this).FAreaK(loc).Mag(), upperKNG);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  maxWS = CellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(upFaceK), (*this).State(upperK), eqnState );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(upperKNG) + maxWS, upperKNG);
	}

      }
    }
  }

}

/* Member function to calculate the local time step. (i,j,k) are cell indices. The following equation is used:

dt = CFL * V / (Lci + Lcj + Lck + C * (Lvi + Lvj + Lvk)) (Blazek 6.18)

In the above equation dt is the time step, CFL is the CFL number, V is the cell volume, Lci, Lcj, Lck are the 
convective (inviscid) spectral radii in the i, j, and k directions, C is a constant (typical value b/w 1 and 4), and 
Lvi, Lvj, Lvk are the viscous spectral radii. This function is only used when the time step isn't explicitly defined
by the user.
*/
void procBlock::CalcCellDt( const int &i, const int &j, const int &k, const double &cfl){
  // i -- i index of cell
  // j -- j index of cell
  // k -- k index of cell
  // cfl -- cfl number

  //location without ghost cells
  int loc = GetLoc1D(i, j, k, (*this).NumI(), (*this).NumJ());
  //location with ghost cells
  int locG = GetLoc1D(i + (*this).NumGhosts(), j + (*this).NumGhosts(), k + (*this).NumGhosts()
		      , (*this).NumI() + 2 * (*this).NumGhosts(), (*this).NumJ() + 2 * (*this).NumGhosts());
  double dt = cfl * ((*this).Vol(locG) / (*this).AvgWaveSpeed(loc)) ; //use nondimensional time

  (*this).SetDt(dt, loc);

}

/* Member function to calculate the time step for all cells in the procBlock. If the time step is user specified
assign that time step (after nondimensionalization) to dt variable. If time step is to be determined using CFL
number, call function to do so.
*/
void procBlock::CalcBlockTimeStep( const input &inputVars, const double &aRef){
  // inputVars -- all input variables
  // aRef -- reference speed of sound (used for time non dimensionalization)

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //loop over all physical cells - no ghost cells for dt variable
  for ( int kk = 0; kk < kmax; kk++ ){          
    for ( int jj = 0; jj < jmax; jj++ ){          
      for ( int ii = 0; ii < imax; ii++ ){          

	int loc = GetLoc1D(ii, jj, kk, imax, jmax); //current cell location

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

/* Member function to update the procBlock to advance to a new time step. For explicit methods it calls the appropriate
explicit method to update. For implicit methods it uses the correction du and calls the implicit updater.
*/
void procBlock::UpdateBlock(const input &inputVars, const int &impFlag, const idealGas &eos, const double &aRef, 
			    const vector<genArray> &du, genArray &l2, resid &linf){
  // inputVars -- all input variables
  // impFlag -- flag to determine if simulation is to be solved via explicit or implicit time stepping
  // eos -- equation of state
  // aRef -- reference speed of sound (for nondimensionalization)
  // bb
  // du -- updates to conservative variables (only used in implicit solver)
  // l2 -- l-2 norm of residual
  // linf -- l-infinity norm of residual

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  if ( inputVars.TimeIntegration() != "rk4" ){ //if not runge-kutta 4 step method for time integration
    for ( int kk = 0; kk < kmax; kk++ ){          //loop over all physical cells
      for ( int jj = 0; jj < jmax; jj++ ){          
	for ( int ii = 0; ii < imax; ii++ ){          

	  //location with and without ghost cells
	  int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	  int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	  if (inputVars.TimeIntegration() == "explicitEuler"){ //explicit euler time integration
	    (*this).ExplicitEulerTimeAdvance(eos, locG, loc);
	  }
	  else if (impFlag){ //if implicit use update (du)
	    (*this).ImplicitTimeAdvance(du[loc], eos, locG);
	  }

	  //accumulate l2 norm of residual
	  l2 = l2 + (*this).Residual(loc) * (*this).Residual(loc);

	  //if any residual is larger than previous residual, a new linf residual is found
	  for ( int ll = 0; ll < NUMVARS; ll++ ){
	    if ( (*this).Residual(loc,ll) > linf.Linf() ){
	      linf.SetLinf( (*this).Residual(loc,ll) ); //store linf residual
	      linf.SetBlock( (*this).ParentBlock() ); //store block location
	      linf.SetILoc( ii ); //store i index
	      linf.SetJLoc( jj ); //store j index
	      linf.SetKLoc( kk ); //store k index
	      linf.SetEqn( ll+1 ); //store equation number
	    }
	  }

	}
      }
    }
  }
  else if ( inputVars.TimeIntegration() == "rk4" ){ //using min storage rk4 method

    vector<primVars> stateN(imax*jmax*kmax);
    vector<double> dtN(imax*jmax*kmax);

    for ( int rr = 0; rr < 4; rr ++ ){ //loop over rk stages
      for ( int kk = 0; kk < kmax; kk++ ){          //loop over all physical cells
	for ( int jj = 0; jj < jmax; jj++ ){          
	  for ( int ii = 0; ii < imax; ii++ ){          

	    //location with and without ghost cells
	    int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	    int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	    //save state and local time step at time n
	    if (rr == 0){
	      stateN[locG] = (*this).State(locG);
	      dtN[loc] = (*this).Dt(loc);
	    }

	    //advance 1 RK stage
	    (*this).RK4TimeAdvance(stateN[locG], eos, dtN[loc], locG, loc, rr);

	    if (rr ==3){ //at last stage

	      //accumulate l2 norm of residual
	      l2 = l2 + (*this).Residual(loc) * (*this).Residual(loc);

	      for ( int ll = 0; ll < NUMVARS; ll++ ){
		if ( (*this).Residual(loc,ll) > linf.Linf() ){
		  linf.SetLinf( (*this).Residual(loc,ll) ); //store linf residual
		  linf.SetBlock( (*this).ParentBlock() ); //store block location
		  linf.SetILoc( ii ); //store i index
		  linf.SetJLoc( jj ); //store j index
		  linf.SetKLoc( kk ); //store k index
		  linf.SetEqn( ll+1 ); //store equation number
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

/* Member function to advance the state vector to time n+1 using explicit Euler method. The following equation is used:

 Un+1 = Un - dt/V * R

Un is the conserved variables at time n, Un+1 is the conserved variables at time n+1, dt is the cell's time step, V is
the cell's volume, and R is the cell's residual.
 */
void procBlock::ExplicitEulerTimeAdvance(const idealGas &eqnState, const int &locG, const int &loc ){
  // eqnState -- equation of state
  // locG -- location of cell (including ghost cells)
  // loc -- location of cell (withod ghost cells)

  //Get conserved variables for current state (time n)
  genArray consVars = (*this).State(locG).ConsVars(eqnState);
  //calculate updated conserved variables
  consVars = consVars - (*this).Dt(loc) / (*this).Vol(locG) * (*this).Residual(loc);

  //calculate updated primative variables
  vector3d<double> vel(consVars[1]/consVars[0], consVars[2]/consVars[0], consVars[3]/consVars[0]);

  primVars tempState (consVars[0],
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars[0], consVars[4]/consVars[0], vel.Mag() ) );

  //update state
  (*this).SetState(tempState, locG);
}

//member function to advance the state vector to time n+1 (for implicit methods)
void procBlock::ImplicitTimeAdvance(const genArray &du, const idealGas &eqnState, const int &loc ){
  // du -- update for a specific cell (to move from time n to n+1)
  // eqnState -- equation of state
  // loc -- location of cell

  //calculate update state (primative variables)
  primVars tempState = (*this).State(loc).UpdateWithConsVars(eqnState, du);

  //check for positivity
  if (tempState.Rho() < 0.0 || tempState.P() < 0.0){
    cerr << "ERROR: Density or pressure has become negative!" << endl;
    cerr << "Updated Primative variables:" << endl << tempState << endl;
    cerr << "Original Primative variables:" << endl << (*this).State(loc) << endl;
    exit(0);
  }

  //assign updated state
  (*this).SetState(tempState, loc);
}

/*member function to advance the state vector to time n+1 using 4th order (minimum storage) Runge-Kutta method (2nd order accurate)

 Un+1 = Un - dt/V * alpha * R

Un is the conserved variables at time n, Un+1 is the conserved variables at time n+1, dt is the cell's time step, V is
the cell's volume, alpha is the runge-kutta coefficient, and R is the cell's residual.
 */
void procBlock::RK4TimeAdvance( const primVars &currState, const idealGas &eqnState, const double &dt, const int &locG, const int &loc, const int &rk ){
  // currState -- current state (including steps within RK4) (primative)
  // eqnState -- equation of state
  // dt -- time step for cell
  // locG -- location of cell including ghost cells
  // loc -- location of cell without ghost cells
  // rk -- runge-kutta step number

  //runge-kutta step coefficients (low storage 4 step)
  double alpha[4] = {0.25, 1.0/3.0, 0.5, 1.0};

  //update conserved variables
  genArray consVars = currState.ConsVars(eqnState) - (*this).Dt(loc) / (*this).Vol(locG) * alpha[rk] * (*this).Residual(loc);

  //calculate updated primative variables
  vector3d<double> vel(consVars[1]/consVars[0], consVars[2]/consVars[0], consVars[3]/consVars[0]);

  primVars tempState (consVars[0],
		      vel.X(),
		      vel.Y(),
		      vel.Z(),
		      eqnState.GetPressFromEnergy( consVars[0], consVars[4]/consVars[0], vel.Mag() ) );

  //assign updated state
  (*this).SetState(tempState, locG);
}

//member function to reset the residual and wave speed back to zero after an iteration. This is done because the residual and wave
//speed are accumulated over many function calls.
void procBlock::ResetResidWS( ){

  //create an instance of genArray and initialize it to 0.
  genArray initial(0.0);

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //loop over all physical cells - no ghost cells in residual variable
  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	int loc = GetLoc1D(ii, jj, kk, imax, jmax); //current cell location

	//reset residual
	(*this).SetResidual( initial, loc ) ;
	
	//reset wave speed
	(*this).SetAvgWaveSpeed( 0.0, loc ) ;

      }
    }
  }

}

/* Member function to add the cell volume divided by the cell time step to the main diagonal of the time m minus time n term.

dU/dt = V/t * [ ((1 + zeta) * FD - zeta * BD) / ((1 + theta) * FD )] * Un = -Rn

The above equation shows the governing equations written in the Beam & Warming format for time integration. U is the vector of conserved variables
where n represents the time step. Theta and zeta are Beam & Warming parameters, t is the time step, V is the cell volume, and R is the residual. 
FD and BD are the forward and backward difference operators respectively. These opererators operate in the time domain. For example FD(U) = 
Un+1 - Un and BD(U) = Un - Un-1. Solving the above equation for FD(Qn) we get the following:

FD(Un) = (-t * Rn - t * theta * FD(Rn) + zeta * V * FD(Un-1)) / ((1 + zeta) * V)

FD(Rn) requires us to know the residual at time n+1, but this is unknown. To bypass this difficulty we linearize the residual using a Taylor series
expansion about time n. Rn+1 = Rn + J*FD(Un) where J is the flux jacobian dR/dU. Rearranging the above equation we get the following:

[J + (1+zeta)*V/(t*theta)] * FD(Un) = -Rn/theta + zeta*V/(t*theta) * FD(Un-1)

The above equation shows that the time m minus time n term (FD(Un)) requires a (1+zeta)V/(t*theta) term multiplied by it. That is the purpose of this
function.
*/
vector<genArray> procBlock::AddVolTime(const vector<genArray> &m, const vector<genArray> &n, const double &theta, const double &zeta) const {
  // m -- solution for block at time m
  // n -- solution for block at time n
  // theta -- Beam & Warming coefficient theta for time integration
  // zeta -- Beam & Warming coefficient zeta for time integration

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  vector<genArray> mMinusN(m.size()); //initialize a vector to hold the returned values

  //loop over all physical cells
  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	//get location of current cell with and without ghost cells
	int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	double I = ( (*this).Vol(locG) * (1.0 + zeta) ) / ( (*this).Dt(loc) * theta ) ;
	mMinusN[loc] = I * (m[loc] - n[loc]);
      }
    }
  }
  return mMinusN;
}

/* Member function to calculate the delta n-1 term for the implicit bdf2 solver.

dU/dt = V/t * [ ((1 + zeta) * FD - zeta * BD) / ((1 + theta) * FD )] * Un = -Rn

The above equation shows the governing equations written in the Beam & Warming format for time integration. U is the vector of conserved variables
where n represents the time step. Theta and zeta are Beam & Warming parameters, t is the time step, V is the cell volume, and R is the residual. 
FD and BD are the forward and backward difference operators respectively. These opererators operate in the time domain. For example FD(U) = 
Un+1 - Un and BD(U) = Un - Un-1. Solving the above equation for FD(Qn) we get the following:

FD(Un) = (-t * Rn - t * theta * FD(Rn) + zeta * V * FD(Un-1)) / ((1 + zeta) * V)

FD(Rn) requires us to know the residual at time n+1, but this is unknown. To bypass this difficulty we linearize the residual using a Taylor series
expansion about time n. Rn+1 = Rn + J*FD(Un) where J is the flux jacobian dR/dU. Rearranging the above equation we get the following:

[J + (1+zeta)*V/(t*theta)] * FD(Un) = -Rn/theta + zeta*V/(t*theta) * FD(Un-1)

The above equation shows that the time n minus time n-1 term (FD(Un-1)) requires a zeta*V/(t*theta) term multiplied by it. That is the purpose of this
function.
*/
void procBlock::DeltaNMinusOne(vector<genArray> &solDeltaNm1, const vector<genArray> &solTimeN, const idealGas &eqnState, const double &theta, const double &zeta){
  // solDeltaNm1 -- The solution at time n minus the solution at time n-1. (Un - Un-1) (output)
  // solTimeN -- The solution at time n
  // eqnState -- equation of state
  // theta -- Beam & Warming coefficient theta for time integration
  // zeta -- Beam & Warming coefficient zeta for time integration

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  //loop over physical cells
  for ( int ii = 0; ii < imax; ii++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int kk = 0; kk < kmax; kk++ ){

	//get location of current cell with and without ghost cells
	int loc =  GetLoc1D(ii, jj, kk, imax, jmax);
	int locG = GetLoc1D(ii + (*this).NumGhosts(), jj + (*this).NumGhosts(), kk + (*this).NumGhosts(), imaxG, jmaxG);

	double coeff = ( (*this).Vol(locG) * zeta ) / ( (*this).Dt(loc) * theta ) ;
	solDeltaNm1[loc] = coeff * ( (*this).State(locG).ConsVars(eqnState) - solTimeN[loc] );
      }
    }
  }

}

//Member function to return a copy of the conserved variables. This is useful for "saving" the solution at a specific time.
//It is used in the implicit solver.
vector<genArray> procBlock::GetCopyConsVars(const idealGas &eqnState) const {

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();

  vector<genArray> consVars(imax * jmax * kmax); //initialize vector to proper size

  //loop over physical cells
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
      for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

	//get location of current cell with and without ghost cells
	int locG = GetLoc1D(ii, jj, kk, imaxG, jmaxG);
	int loc = GetLoc1D(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//convert state to conservative variables
	consVars[loc] = (*this).State(locG).ConsVars(eqnState);
      }
    }
  }

  return consVars;
}

/* Member function to calculate update to solution implicitly using Lower-Upper Symmetric Gauss Seidel (LUSGS) method.

Un+1 = Un - t/V * Rn+1

The above equation shows a simple first order implicit method to calculate the solution at the next time step (n+1). The equation shows
that this method requires the residual (R) at time n+1 which is unknown. In the equation, t is the time step, and V is the volume. Since
the residual at n+1 in unknown, it must be linearized about time n as shown below.

Rn+1 = Rn + dRn/dUn * FD(Un)

In the above equation FD is the forward difference operator in time (FD(Un) = Un+1 - Un). The derivative of the residual term can be further
simplified as below.

Rn = (sum over all faces) Fni

In the above equation n refers to the time level and i refers to the face index. The summation over all faces operator will now be abbreviated
as (SF). Substituting the second and third equations into the first and rearranging we get the following.

[d(SF)Fni/dUnj + V/t] * FD(Un) = -Rn

In the above equation the index j refers to all of the cells in the stencil going into the calculation of the flux at face i. The above equation
can be simplified to A*x=b. The matrix A is an MxM block matrix where M is the number of cells in the block. Each block is an LxL block where L
is the number of equations being solved for. The sparsity of A depends on the stencil used in flux calculation. In 3D for a first order simulation
the matrix A is block pentadiagonal. For a second order approximation it would have 13 diagonals. This increases the storage requirements so in 
practice a first order approximation is used. The order of accuracy is determined by the residual calculation. The accuracy of the matrix A helps
with convergence. A poorer approximation will eventually get the correct answer with enough iteration (defect correction). Fully implicit methods
calculate and store the flux jacobians needed to populate the matrix A. The LUSGS method does not do this and instead calculates an approximate 
flux jacobian "on-the-fly" so there is no need for storage. Because an approximate flux jacobian is being used, there is no need for a highly 
accurate linear solver. Therefore the Symmetric Gauss-Seidel (SGS) method is a good candidate, and only one iteration is needed. This approximate 
flux jacobian along with the SGS linear solver form the basics of the LUSGS method (Jameson & Yoon). The LUSGS method begins by factoring the 
matrix A as shown below.

A = (D + L) * D^-1 * (D + U)

In the above equation D is the diagonal of A, L is the lower triangular portion of A, and U is the upper triangular portion of A. This allows the 
equation A*x=b to be solved in one SGS sweep as shown below.

Forward sweep:  (D + L) * FD(Un*) = -Rn
Backward sweep: (D + U) * FD(Un) = D * FD(Un*)

Another key component of the LUSGS scheme is to sweep along hyperplanes. Hyperplanes are planes of i+j+k=constant within a plot3d block. The
diagram below shows a 2D example of a block reordered to sweep along hyperplanes

           ____ ____ ____ ____ ____ ____ ____ ____
          | 20 | 26 | 32 | 37 | 41 | 44 | 46 | 47 |
          |____|____|____|____|____|____|____|____|
          | 14 | 19 | 25 | 31 | 36 | 40 | 43 | 45 |
          |____|____|____|____|____|____|____|____|
          | 9  | 13 | 18 | 24 | 30 | 35 | 39 | 42 |
   A=     |____|____|____|____|____|____|____|____|
          | 5  | 8  | 12 | 17 | 23 | 29 | 34 | 38 |
          |____|____|____|____|____|____|____|____|
          | 2  | 4  | 7  | 11 | 16 | 22 | 28 | 33 |
          |____|____|____|____|____|____|____|____|
          | 0  | 1  | 3  | 6  | 10 | 15 | 21 | 27 |
          |____|____|____|____|____|____|____|____|

This is advantageous because on the forward sweep the lower matrix L can be calculated with data from time n+1 because all of the cells contributing
to L would already have been updated. The same is true with the upper matrix U for the backward sweep. This removes the need for any storage of the 
matrix. For a given location (say A12) the matrix L would be constructed of A8 and A7, and the matrix U would be constructed of A18 and A17. This 
requires the product of the flux jacobian multiplied with the update (FD(Un)) to be calculated at these locations. For the LUSGS method the flux 
jacobians are approximated as follows:

A * S = 0.5 * (Ac * S + K * I)

A is the flux jacobian, S is the face area, Ac is the convective flux jacobian (dF/dU), K is the spectral radius multiplied by a factor, and I is
the identity matrix. The addition of the spectrial radius improves diagonal dominance which improves stability at the cost of convergence. When the 
factor multiplied by K is 1, the method is SGS, when it is < 1, it is successive overrelaxation. Reducing the factor improves convergence but hurts
stability. The product of the approximate flux jacobian with the update is calculated as shown below.

A * S * FD(Unj) = 0.5 * (Ac * FD(Unj) * S + K * I * FD(Unj)) = 0.5 * ( dFi/dUnj * FD(Unj) * S + K * I * FD(Unj)) = 0.5 * (dFi * S + K * I * FD(Unj))

The above equation shows that all that is needed to calculate the RHS of A*x=b is the update of the convective flux, and the update to the convervative
variabes (FD(Unj)) which is known due to sweeping along hyperplanes.

For viscous simulations, the viscous contribution to the spectral radius K is used, and everything else remains the same.
 */
double procBlock::LUSGS( const vector<vector3d<int> > &reorder, vector<genArray> &x, const vector<genArray> &solTimeMmN, 
			 const vector<genArray> &solDeltaNm1, const idealGas &eqnState, const input &inp, const sutherland &suth)const{
  // reorder -- order of cells to visit (this should be ordered in hyperplanes)
  // x -- correction - added to solution at time n to get to time n+1 (assumed to be zero to start)
  // solTimeMmn -- solution at time m minus n
  // solDeltaNm1 -- solution at time n minus solution at time n-1 
  // eqnState -- equation of state
  // inp -- all input variables
  // suth -- method to get temperature varying viscosity (Sutherland's law)

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  //initialize inverse to diagonal block. For LUSGS block inversion is replaced by scalar division
  double AiiInv = 0.0; 

  double thetaInv = 1.0 / inp.Theta();

  //initialize genArray to zero
  genArray initial(0.0);

  //initialize L and U matrices
  vector<genArray> U(x.size(),initial);
  vector<genArray> L(x.size(),initial);

  //------------------------------------------------------------------------------------------------------------------------------------------------------
  //forward sweep over all physical cells
  for ( int ii = 0; ii < (int)x.size(); ii++ ){

    //indicies for variables without ghost cells
    //location of cell on diagonal
    int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

    //location of cells on the lower i, j, k sides of diagonal cell 
    int il = GetNeighborLowI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
    int jl = GetNeighborLowJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
    int kl = GetNeighborLowK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

    //indicies for variables with ghost cells
    //location of cell on diagonal
    int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //location of faces of the diagonal cell that touch cells making up row of L matrix
    int ilFaceG = GetLowerFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int jlFaceG = GetLowerFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int klFaceG = GetLowerFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //location of faces of cells making up row of L matrix (used in spectral radius calculation)
    int ilFace2G = GetLowerFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
    int jlFace2G = GetLowerFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
    int klFace2G = GetLowerFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);

    //location of cells making up row of L matrix
    int ilG = GetNeighborLowI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int jlG = GetNeighborLowJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int klG = GetNeighborLowK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //if i lower diagonal cell is in physical location there is a contribution from it
    if ( il >=0 && il < (int)x.size() ){
      //at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
      double specRad = CellSpectralRadius( (*this).FAreaI(ilFace2G), (*this).FAreaI(ilFaceG), (*this).State(ilG).UpdateWithConsVars(eqnState, x[il]), eqnState);

      if (inp.EquationSet() != "euler"){ //if viscous add viscous contribution to spectral radius
	double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	double mRef = inp.VelRef().Mag() / aRef;
	double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(ilFace2G), (*this).FAreaI(ilFaceG), 
							      (*this).State(ilG).UpdateWithConsVars(eqnState, x[il]), eqnState, suth, (*this).Vol(ilG) );
	specRad += vSpecRad;
      }

      //at given face location, call function to calculate convective flux change
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(ilG), eqnState, (*this).FAreaI(ilFaceG), x[il]);

      //update L matrix
      L[loc] = L[loc] + 0.5 * ( (*this).FAreaI(ilFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * x[il] );
    }
    //if j lower diagonal cell is in physical location there is a contribution from it
    if ( jl >=0 && jl < (int)x.size() ){
      //at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
      double specRad = CellSpectralRadius( (*this).FAreaJ(jlFace2G), (*this).FAreaJ(jlFaceG), (*this).State(jlG).UpdateWithConsVars(eqnState, x[jl]), eqnState);

      if (inp.EquationSet() != "euler"){ //if viscous add viscous contribution to spectral radius
	double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	double mRef = inp.VelRef().Mag() / aRef;
	double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(jlFace2G), (*this).FAreaJ(jlFaceG), 
							      (*this).State(jlG).UpdateWithConsVars(eqnState, x[jl]), eqnState, suth, (*this).Vol(jlG) );
	specRad += vSpecRad;
      }

      //at given face location, call function to calculate convective flux change
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(jlG), eqnState, (*this).FAreaJ(jlFaceG), x[jl]);

      //update L matrix
      L[loc] = L[loc] + 0.5 * ( (*this).FAreaJ(jlFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * x[jl] );
    }
    //if k lower diagonal cell is in physical location there is a contribution from it
    if ( kl >=0 && kl < (int)x.size() ){
      //at given face location, call function to calculate spectral radius, since values are constant throughout cell, cell center values are used
      double specRad = CellSpectralRadius( (*this).FAreaK(klFace2G), (*this).FAreaK(klFaceG), (*this).State(klG).UpdateWithConsVars(eqnState, x[kl]), eqnState);

      if (inp.EquationSet() != "euler"){ //if viscous add viscous contribution to spectral radius
	double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
	double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
	double mRef = inp.VelRef().Mag() / aRef;
	double vSpecRad = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(klFace2G), (*this).FAreaK(klFaceG), 
							      (*this).State(klG).UpdateWithConsVars(eqnState, x[kl]), eqnState, suth, (*this).Vol(klG) );
	specRad += vSpecRad;
      }

      //at given face location, call function to calculate convective flux change
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(klG), eqnState, (*this).FAreaK(klFaceG), x[kl]);

      //update L matrix
      L[loc] = L[loc] + 0.5 * ( (*this).FAreaK(klFaceG).Mag() * fluxChange + inp.MatrixRelaxation() * specRad * x[kl] );
    }

    //add dual time stepping contribution to main diagonal
    double diagTimeVol = ( (*this).Vol(locG) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
    if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
      double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
      diagTimeVol += tau;
    }

    AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

    //calculate intermediate update
    x[loc] = AiiInv * ( -1.0 * thetaInv * (*this).Residual(loc) - solDeltaNm1[loc] - solTimeMmN[loc] + L[loc] ) ; //normal at lower boundaries needs to be reversed, so add instead of subtract L

  } //end forward sweep

  //------------------------------------------------------------------------------------------------------------------------------------------------------
  //backward sweep over all physical cells
  for ( int ii = (int)x.size()-1; ii >= 0; ii-- ){

    //indicies for variables without ghost cells
    //location of cell on diagonal
    int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

    //location of cells on the upper i, j, k sides of diagonal cell 
    int iu = GetNeighborUpI(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
    int ju = GetNeighborUpJ(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
    int ku = GetNeighborUpK(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);

    //indicies for variables with ghost cells
    //location of cell on diagonal
    int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //location of faces of the diagonal cell that touch cells making up row of U matrix
    int iuFaceG = GetUpperFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int juFaceG = GetUpperFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int kuFaceG = GetUpperFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //location of faces of cells making up row of U matrix (used in spectral radius calculation)
    int iuFace2G = GetUpperFaceI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
    int juFace2G = GetUpperFaceJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);
    int kuFace2G = GetUpperFaceK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG, 2);

    //location of cells making up row of U matrix
    int iuG = GetNeighborUpI(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int juG = GetNeighborUpJ(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);
    int kuG = GetNeighborUpK(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //if i upper diagonal cell is in physical location there is a contribution from it
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
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(iuG), eqnState, (*this).FAreaI(iuFaceG), x[iu]);

      //update U matrix
      U[loc] = U[loc] + 0.5 * ( (*this).FAreaI(iuFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * x[iu] );
    }
    //if j upper diagonal cell is in physical location there is a contribution from it
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
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(juG), eqnState, (*this).FAreaJ(juFaceG), x[ju]);

      //update U matrix
      U[loc] = U[loc] + 0.5 * ( (*this).FAreaJ(juFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * x[ju] );
    }
    //if k upper diagonal cell is in physical location there is a contribution from it
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
      genArray fluxChange = ConvectiveFluxUpdate( (*this).State(kuG), eqnState, (*this).FAreaK(kuFaceG), x[ku]);

      //update U matrix
      U[loc] = U[loc] + 0.5 * ( (*this).FAreaK(kuFaceG).Mag() * fluxChange - inp.MatrixRelaxation() * specRad * x[ku] );
    }

    //add dual time stepping contribution to main diagonal
    double diagTimeVol = ( (*this).Vol(locG) * (1.0 + inp.Zeta()) ) / ( (*this).Dt(loc) * inp.Theta() );
    if (inp.DualTimeCFL() > 0.0 ) { //use dual time stepping
      double tau = (*this).AvgWaveSpeed(loc) / inp.DualTimeCFL(); // equal to volume / tau
      diagTimeVol += tau;
    }

    AiiInv = 1.0 / ( ((*this).AvgWaveSpeed(loc) + diagTimeVol ) * inp.MatrixRelaxation() );

    //calculate update
    x[loc] = x[loc] - AiiInv * U[loc];

  } //end backward sweep

  //------------------------------------------------------------------------------------------------------------------------------------------------------
  //calculate residual
  //initialize LUSGS residual vector
  genArray l2Resid(0.0);

  genArray resid(0.0);

  for ( int ii = 0; ii < (int)x.size(); ii++ ){

    //location of diagonal cell with and without ghost cells
    int loc = GetLoc1D(reorder[ii].X(), reorder[ii].Y(), reorder[ii].Z(), imax, jmax);
    int locG = GetLoc1D(reorder[ii].X() + (*this).NumGhosts(), reorder[ii].Y() + (*this).NumGhosts(), reorder[ii].Z() + (*this).NumGhosts(), imaxG, jmaxG);

    //calclate dual time stepping contribution
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

  return l2Resid.Sum();

}

/*Function to return the inviscid spectral radius for one direction (i, j, or k) given a cell state, equation of state, and 2 face area vectors

L = 0.5 * (A1 + A2) * (|Vn| + SoS)

In the above equation L is the spectral radius in either the i, j, or k direction. A1 and A2 are the two face areas in that direction. Vn is the
cell velocity normal to that direction. SoS is the speed of sound at the cell
 */
double CellSpectralRadius(const vector3d<double> &fAreaL, const vector3d<double> &fAreaR, const primVars &state, const idealGas &eqnState){
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // state -- state at cell center (primative)
  // eqnState -- equation of state

  //normalize face areas
  vector3d<double> normAreaL = fAreaL / fAreaL.Mag();
  vector3d<double> normAreaR = fAreaR / fAreaR.Mag();

  vector3d<double> normAvg = 0.5 * (normAreaL + normAreaR); //average normal direction
  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag()); //average area magnitude

  //return spectral radius
  return ( fabs(state.Velocity().DotProd(normAvg)) + state.SoS(eqnState) ) * fMag;
}

/*Function to calculate the viscous spectral radius for one direction (i, j, or k).

L = max(4/(3*rho), g/rho) * mu/Pr * A^2 / V

In the above equation L is the viscous spectral radius for a given direction (i, j, or k). Rho is the density
at the cell center. G is gamma, mu is viscosity, and Pr is the Prandtl number (all at the cell center). A is the
average face area of the given direction (i, j, k), and V is the cell volume. This implementation comes from Blazek.
 */
double ViscCellSpectralRadius(const vector3d<double> &fAreaL, const vector3d<double> &fAreaR, const primVars &state, const idealGas &eqnState, const sutherland &suth, const double &vol){
  // fAreaL -- face area of lower face in either i, j, or k direction
  // fAreaR -- face area of upper face in either i, j, or k direction
  // state -- state at cell center (primative)
  // eqnState -- equation of state
  // suth -- method to the temperature varying visosity and Prandtl number (Sutherland's law)
  // vol -- cell volume

  double fMag = 0.5 * (fAreaL.Mag() + fAreaR.Mag()); //average area magnitude
  double maxTerm = max(4.0 / (3.0 * state.Rho()), eqnState.Gamma() / state.Rho()) ;
  double mu = suth.GetViscosity(state.Temperature(eqnState)); //viscosity at cell center
  double viscTerm = mu / eqnState.GetPrandtl();

  //return viscous spectral radius
  return maxTerm * viscTerm * fMag * fMag / vol ;
}

//member function to reconstruct cell variables to the face using central differences
template <class T>
T FaceReconCentral(const T &velU, const T &velD, const vector3d<double> &pU, const vector3d<double> &pD, const vector3d<double> &pF){
  // velU -- velocity at the cell center of the upwind cell
  // velD -- velocity at the cell center of the downwind cell
  // pU -- position of the cell center of the upwind cell
  // pD -- position of the cell center of the downwind cell
  // pF -- position of the face center of the face on which the reconstruction is happening

  T temp;

  double cen2cen = pU.Distance(pD);  //distance from cell center to cell center
  double up2face = pU.Distance(pF);  //distance from upwind cell center to cell face

  //reconstruct with central difference
  temp = velD * (up2face/cen2cen) + velU * (1.0 - (up2face/cen2cen));

  return temp;
}

/* Function to pad a vector with a specified number of ghost cells
           ___ ___ ___ ___ ___ ___ ___ ___
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | G | G | X | X | X | X | G | G |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|
          | E | E | G | G | G | G | E | E |
          |___|___|___|___|___|___|___|___|

In the above diagram, the cells marked with an "X" represent physical cells. The entire diagram represents the block (in 2D)
padded with 2 layers of ghost cells. The cells marked with "G" are regualar ghost cells. The cells marked with "E" are ghost
cells located along one of the 12 edges that form a plot3d block. In 3D there are also "corner" cells located at the 8 corners
that form the plot3d block. These cells are not used though. There is a place in the vector for them to make accessing the padded
vector of cells the same as for a plot3d block without ghost cells.
*/
template<class T>
vector<T> PadWithGhosts( const vector<T> &var, const int &numGhosts, const int &numI, const int &numJ, const int &numK ){
  // var -- vector of variables to pad (no ghost cells included)
  // numGhosts -- number of layers of ghost cells to pad var with
  // numI -- i dimension without ghost cells
  // numJ -- j dimension without ghost cells
  // numK -- k dimension withough ghost cells

  //max dimension for variables with ghost cells
  int newI = numI + (numGhosts * 2);
  int newJ = numJ + (numGhosts * 2);
  int newK = numK + (numGhosts * 2);

  int newSize = newI * newJ * newK; //size of vector with padded ghost cells

  vector<T> padBlk(newSize);        //initialize vector

  //loop over physical cells
  for ( int kk = numGhosts; kk < numK + numGhosts; kk++ ){
    for ( int jj = numGhosts; jj < numJ + numGhosts; jj++ ){
      for ( int ii = numGhosts; ii < numI + numGhosts; ii++ ){

	//calcualte location with and with ghost cells
	int newLoc = GetLoc1D(ii, jj, kk, newI, newJ);
	int loc = GetLoc1D(ii-numGhosts, jj-numGhosts, kk-numGhosts, numI, numJ);

	//assign the given vector of variables to the correct location within the padded vector
	padBlk[newLoc] = var[loc]; 
      }
    }
  }

  return padBlk;
}

/* Function to calculate the velocity gradient at the cell center using the Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a scalar is calculated using the Green-Gauss method. U is a scalar.
Aij is the area at face i (component j). V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out of the control volume.
 */
tensor<double> CalcVelGradGG(const vector3d<double> &vil, const vector3d<double> &viu, const vector3d<double> &vjl, const vector3d<double> &vju, 
					const vector3d<double> &vkl, const vector3d<double> &vku, const vector3d<double> &ail, const vector3d<double> &aiu,
					const vector3d<double> &ajl, const vector3d<double> &aju, const vector3d<double> &akl, const vector3d<double> &aku,
					const double &vol){

  //vil -- velocity vector at the i-lower face of the cell at which the velocity gradient is being calculated
  //viu -- velocity vector at the i-upper face of the cell at which the velocity gradient is being calculated
  //vjl -- velocity vector at the j-lower face of the cell at which the velocity gradient is being calculated
  //vju -- velocity vector at the j-upper face of the cell at which the velocity gradient is being calculated
  //vkl -- velocity vector at the k-lower face of the cell at which the velocity gradient is being calculated
  //vku -- velocity vector at the k-upper face of the cell at which the velocity gradient is being calculated

  //ail -- area vector at the lower i-face of the cell at which the velocity gradient is being calculated
  //aiu -- area vector at the upper i-face of the cell at which the velocity gradient is being calculated
  //ajl -- area vector at the lower j-face of the cell at which the velocity gradient is being calculated
  //aju -- area vector at the upper j-face of the cell at which the velocity gradient is being calculated
  //akl -- area vector at the lower k-face of the cell at which the velocity gradient is being calculated
  //aku -- area vector at the upper k-face of the cell at which the velocity gradient is being calculated

  //vol -- cell volume

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

/* Function to calculate the temperature gradient at the cell center using the Green-Gauss method

dU/dxj = (Sum)    U * Aij  / V        (j=1,2,3)
       i=1,nfaces

The above equation shows how the gradient of a scalar is calculated using the Green-Gauss method. U is a scalar.
Aij is the area at face i (component j). V is the volume of the control volume. X is the cartesian direction with
j indicating the component. The convention is for the area vectors to point out of the control volume.
 */
vector3d<double> CalcTempGradGG(const double &til, const double &tiu, const double &tjl, const double &tju, const double &tkl, const double &tku,
					   const vector3d<double> &ail, const vector3d<double> &aiu, const vector3d<double> &ajl, const vector3d<double> &aju,
					   const vector3d<double> &akl, const vector3d<double> &aku, const double &vol){

  //til -- temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tiu -- temperature at the upper face of the cell at which the temperature gradient is being calculated
  //tjl -- temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tju -- temperature at the upper face of the cell at which the temperature gradient is being calculated
  //tkl -- temperature at the lower face of the cell at which the temperature gradient is being calculated
  //tku -- temperature at the upper face of the cell at which the temperature gradient is being calculated
  
  //ail -- area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aiu -- area vector at the upper face of the cell at which the temperature gradient is being calculated
  //ajl -- area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aju -- area vector at the upper face of the cell at which the temperature gradient is being calculated
  //akl -- area vector at the lower face of the cell at which the temperature gradient is being calculated
  //aku -- area vector at the upper face of the cell at which the temperature gradient is being calculated

  //vol -- cell volume

  vector3d<double> temp;
  double invVol = 1.0/vol;

  //define temperature gradient vector
  //convention is for area vector to point out of cell, so lower values are negative, upper are positive
  temp.SetX( invVol * ( tiu*aiu.X() - til*ail.X() + tju*aju.X() - tjl*ajl.X() + tku*aku.X() - tkl*akl.X() ) );
  temp.SetY( invVol * ( tiu*aiu.Y() - til*ail.Y() + tju*aju.Y() - tjl*ajl.Y() + tku*aku.Y() - tkl*akl.Y() ) );
  temp.SetZ( invVol * ( tiu*aiu.Z() - til*ail.Z() + tju*aju.Z() - tjl*ajl.Z() + tku*aku.Z() - tkl*akl.Z() ) );

  return temp;
}

/* Function to calculate the viscous fluxes on the i-faces. All phyiscal (non-ghost) i-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Ui       -->   Ui+1   |
  |            |            |
  |____________|____________|
Ui-1/2       Ui+1/2       Ui+3/2

Using the above diagram, the flux is calculated at face Ui+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Ui to Ui+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Ui, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Ui+1. 

The spectral radius in the i-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Ui+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by constructing an alternative control volume centered
around that face as show below.

  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done. In this example the gradient is being calculated at
the face between cell Ui,j and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the alternative control volume.
The grid cells themselves cannot be used as the control volume for the gradients (averaging values at adjacent cells to get gradients at
the face) because this leads to odd/even decoupling. The face areas, volumes, and states at the center of the faces are needed for the 
alternative control volume. The left and right sides of the alternate control volume pass through the center of cells Ui,j and Ui+1,j
respectively. Therefore these values are used for the face states. The left and right face areas are calculated as the average of the face
areas of the cells that they split. For example, the left face area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j). The top and bottom 
sides of the alternative control volume pass through 4 cells each. Therefore the value of the state at the face center is the average of these
four states. For example the state at the top face is calculated as 0.25 * (Ui,j + Ui+1,j + Ui,j+1, Ui+1,j+1). The face areas of the top and
bottom sides are calculated as the average of the 2 face areas that each one passes through. For example, the top face area is calculated as
0.5 * (Ai,j+1/2 + Ai+1,j+1/2). In three dimensions each gradient calculation touches the values at 10 cells (6 shown and 4 more in/out of 
the page). The stencil for the gradients of all faces in a cell touches 15 cells. The gradient calculation with this stencil uses the "edge" 
ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxI(const sutherland &suth, const idealGas &eqnState, const input &inp){
  // suth -- method to get viscosity as a function of temperature (Sutherland's law)
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI() + 1; //calculating fluxes on i-faces so one more face in i-direction
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts() + 1; //calculating fluxes on i-faces so one more face in i-direction
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  //calculate reference Reynolds number and Mach number for nondimensionalization
  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  //coefficient for viscous spectral radii
  double vCoeff = 1.0;

  //loop over all physical faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	//location of faces in the upper and lower i-direction (with ghost cells included)
	int fUpi = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int fLowi = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);

	//location of j-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower i-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpjUpi = GetUpperFaceJ(ii, jj, kk, imaxG - 1, jmaxG);
	int fUpjLowi = GetUpperFaceJ(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int fLowjUpi = GetLowerFaceJ(ii, jj, kk, imaxG - 1, jmaxG);
	int fLowjLowi = GetLowerFaceJ(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//location of k-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower i-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpkUpi = GetUpperFaceK(ii, jj, kk, imaxG - 1, jmaxG);
	int fUpkLowi = GetUpperFaceK(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int fLowkUpi = GetLowerFaceK(ii, jj, kk, imaxG - 1, jmaxG);
	int fLowkLowi = GetLowerFaceK(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//location of cells in the upper and lower i-direction with respect to baseline face (with ghost cells included)
	int iLow  = GetCellFromFaceLowerI(ii, jj, kk, imaxG, jmaxG);
	int iUp  = GetCellFromFaceUpperI(ii, jj, kk, imaxG, jmaxG);

	//location of cells in the upper and lower j-direction and the upper 
	//and lower i-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int jUpiUp = GetNeighborUpJ(ii, jj, kk, imaxG - 1, jmaxG);
	int jUpiLow = GetNeighborUpJ(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int jLowiUp = GetNeighborLowJ(ii, jj, kk, imaxG - 1, jmaxG);
	int jLowiLow = GetNeighborLowJ(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//location of cells in the upper and lower k-direction and the upper 
	//and lower i-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int kUpiUp = GetNeighborUpK(ii, jj, kk, imaxG - 1, jmaxG);
	int kUpiLow = GetNeighborUpK(ii - 1, jj, kk, imaxG - 1, jmaxG);
	int kLowiUp = GetNeighborLowK(ii, jj, kk, imaxG - 1, jmaxG);
	int kLowiLow = GetNeighborLowK(ii - 1, jj, kk, imaxG - 1, jmaxG);

	//location of cells in the upper and lower i-direction with respect to baseline face (without ghost cells included)
	int iLowNG  = GetCellFromFaceLowerI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int iUpNG  = GetCellFromFaceUpperI(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//calculate average velocity on j and k faces of alternate control volume
	vector3d<double> vju = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(jUpiUp).Velocity() + (*this).State(jUpiLow).Velocity() );
	vector3d<double> vjl = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(jLowiUp).Velocity() + (*this).State(jLowiLow).Velocity() );

	vector3d<double> vku = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(kUpiUp).Velocity() + (*this).State(kUpiLow).Velocity() );
	vector3d<double> vkl = 0.25 * ( (*this).State(iLow).Velocity() + (*this).State(iUp).Velocity() + (*this).State(kLowiUp).Velocity() + (*this).State(kLowiLow).Velocity() );

	//calculate areas of faces in alternate control volume
	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(loc) + (*this).FAreaI(fUpi) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(loc) + (*this).FAreaI(fLowi) );

	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(fUpjUpi) + (*this).FAreaJ(fUpjLowi) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(fLowjUpi) + (*this).FAreaJ(fLowjLowi) );

	vector3d<double> aku = 0.5 * ( (*this).FAreaK(fUpkUpi) + (*this).FAreaK(fUpkLowi) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(fLowkUpi) + (*this).FAreaK(fLowkLowi) );

	//calculate volume of alternate control volume
	double vol = 0.5 * ( (*this).Vol(iLow) + (*this).Vol(iUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( (*this).State(iLow).Velocity(), (*this).State(iUp).Velocity(), vjl, vju, vkl, vku, ail, aiu, ajl, aju, akl, aku, vol);
	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(iLow).Velocity(), (*this).State(iUp).Velocity(), (*this).Center(iLow), (*this).Center(iUp), (*this).FCenterI(loc) );

	//calculate average temperature on j and k faces of alternate control volume
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

	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( ii > (*this).NumGhosts() ){ //at left boundary there is no left cell to add to
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaI(loc).Mag(), iLowNG);
	}
	if ( ii < imax -1 + (*this).NumGhosts() ){ //at right boundary there is no right cell to add to
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaI(loc).Mag(), iUpNG);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  double maxWS = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaI(loc), (*this).FAreaI(fUpi), (*this).State(iUp), eqnState, suth, (*this).Vol(iUp) );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(iUpNG) + vCoeff * maxWS, iUpNG);

	}

      }
    }
  }

}

/* Function to calculate the viscous fluxes on the j-faces. All phyiscal (non-ghost) j-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Uj       -->   Uj+1   |
  |            |            |
  |____________|____________|
Uj-1/2       Uj+1/2       Uj+3/2

Using the above diagram, the flux is calculated at face Uj+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Uj to Uj+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Uj, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Uj+1. 

The spectral radius in the j-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Uj+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by constructing an alternative control volume centered
around that face as show below.

  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done. In this example the gradient is being calculated at
the face between cell Ui,j and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the alternative control volume.
The grid cells themselves cannot be used as the control volume for the gradients (averaging values at adjacent cells to get gradients at
the face) because this leads to odd/even decoupling. The face areas, volumes, and states at the center of the faces are needed for the 
alternative control volume. The left and right sides of the alternate control volume pass through the center of cells Ui,j and Ui+1,j
respectively. Therefore these values are used for the face states. The left and right face areas are calculated as the average of the face
areas of the cells that they split. For example, the left face area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j). The top and bottom 
sides of the alternative control volume pass through 4 cells each. Therefore the value of the state at the face center is the average of these
four states. For example the state at the top face is calculated as 0.25 * (Ui,j + Ui+1,j + Ui,j+1, Ui+1,j+1). The face areas of the top and
bottom sides are calculated as the average of the 2 face areas that each one passes through. For example, the top face area is calculated as
0.5 * (Ai,j+1/2 + Ai+1,j+1/2). In three dimensions each gradient calculation touches the values at 10 cells (6 shown and 4 more in/out of 
the page). The stencil for the gradients of all faces in a cell touches 15 cells. The gradient calculation with this stencil uses the "edge" 
ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxJ(const sutherland &suth, const idealGas &eqnState, const input &inp){
  // suth -- method to get viscosity as a function of temperature (Sutherland's law)
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ() + 1; //calculating fluxes on j-faces so one more face in j-direction
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts() + 1; //calculating fluxes on j-faces so one more face in j-direction

  //calculate reference Reynolds number and Mach number for nondimensionalization
  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  //coefficient for viscous spectral radii
  double vCoeff = 1.0;

  //loop over physical faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	//location of faces in the upper and lower j-direction (with ghost cells included)
	int fUpj = GetNeighborUpJ(ii, jj, kk, imaxG, jmaxG);
	int fLowj = GetNeighborLowJ(ii, jj, kk, imaxG, jmaxG);

	//location of i-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower j-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpiUpj = GetUpperFaceI(ii, jj, kk, imaxG, jmaxG - 1);
	int fUpiLowj = GetUpperFaceI(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int fLowiUpj = GetLowerFaceI(ii, jj, kk, imaxG, jmaxG - 1);
	int fLowiLowj = GetLowerFaceI(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//location of k-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower j-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpkUpj = GetUpperFaceK(ii, jj, kk, imaxG, jmaxG - 1);
	int fUpkLowj = GetUpperFaceK(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int fLowkUpj = GetLowerFaceK(ii, jj, kk, imaxG, jmaxG - 1);
	int fLowkLowj = GetLowerFaceK(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//location of cells in the upper and lower j-direction with respect to baseline face (with ghost cells included)
	int jLow  = GetCellFromFaceLowerJ(ii, jj, kk, imaxG, jmaxG);
	int jUp  = GetCellFromFaceUpperJ(ii, jj, kk, imaxG, jmaxG);

	//location of cells in the upper and lower i-direction and the upper 
	//and lower j-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int iUpjUp = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG - 1);
	int iUpjLow = GetNeighborUpI(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int iLowjUp = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG - 1);
	int iLowjLow = GetNeighborLowI(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//location of cells in the upper and lower k-direction and the upper 
	//and lower j-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int kUpjUp = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG - 1);
	int kUpjLow = GetNeighborUpK(ii, jj - 1, kk, imaxG, jmaxG - 1);
	int kLowjUp = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG - 1);
	int kLowjLow = GetNeighborLowK(ii, jj - 1, kk, imaxG, jmaxG - 1);

	//location of cells in the upper and lower j-direction with respect to baseline face (without ghost cells included)
	int jLowNG  = GetCellFromFaceLowerJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int jUpNG  = GetCellFromFaceUpperJ(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//calculate average velocity on i and k faces of alternate control volume
	vector3d<double> viu = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(iUpjUp).Velocity() + (*this).State(iUpjLow).Velocity() );
	vector3d<double> vil = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(iLowjUp).Velocity() + (*this).State(iLowjLow).Velocity() );

	vector3d<double> vku = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(kUpjUp).Velocity() + (*this).State(kUpjLow).Velocity() );
	vector3d<double> vkl = 0.25 * ( (*this).State(jLow).Velocity() + (*this).State(jUp).Velocity() + (*this).State(kLowjUp).Velocity() + (*this).State(kLowjLow).Velocity() );

	//calculate areas of faces in alternate control volume
	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(loc) + (*this).FAreaJ(fUpj) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(loc) + (*this).FAreaJ(fLowj) );

	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(fUpiUpj) + (*this).FAreaI(fUpiLowj) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(fLowiUpj) + (*this).FAreaI(fLowiLowj) );

	vector3d<double> aku = 0.5 * ( (*this).FAreaK(fUpkUpj) + (*this).FAreaK(fUpkLowj) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(fLowkUpj) + (*this).FAreaK(fLowkLowj) );

	//calculate volume of alternate control volume
	double vol = 0.5 * ( (*this).Vol(jLow) + (*this).Vol(jUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( vil, viu, (*this).State(jLow).Velocity(), (*this).State(jUp).Velocity(), vkl, vku, ail, aiu, ajl, aju, akl, aku, vol);
	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(jLow).Velocity(), (*this).State(jUp).Velocity(), (*this).Center(jLow), (*this).Center(jUp), (*this).FCenterJ(loc) );

	//calculate average temperature on i and k faces of alternate control volume
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

	//calculate viscous flux
	viscousFlux tempViscFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaJ(loc) );

	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( jj > (*this).NumGhosts() ){ //at left boundary there is no left cell to add to
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaJ(loc).Mag(), jLowNG);
	}
	if ( jj < jmax -1 + (*this).NumGhosts() ){ //at right boundary there is no right cell to add to
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaJ(loc).Mag(), jUpNG);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  double maxWS = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaJ(loc), (*this).FAreaJ(fUpj), (*this).State(jUp), eqnState, suth, (*this).Vol(jUp) );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(jUpNG) + vCoeff * maxWS, jUpNG);

	}

      }
    }
  }


}

/* Function to calculate the viscous fluxes on the k-faces. All phyiscal (non-ghost) k-faces are looped over. The
left and right states are calculated, and then the flux at the face is calculated. The flux at the face contributes
to the residual of the cells to the left and right of the face. This contribution from the flux is added to the 
residuals and the wave speed is accumulated as well.

  ___________________________
  |            |            |
  |            |            |
  |   Uk       -->   Uk+1   |
  |            |            |
  |____________|____________|
Uk-1/2       Uk+1/2       Uk+3/2

Using the above diagram, the flux is calculated at face Uk+1/2. Since the area vector at the face always points from
lower indices to higher indices it points from Uk to Uk+1. For the residual calculation the convention is for the area
vector to point out of the cell. Therefore it is in the correct position for the residual at Uk, but the sign needs to 
be flipped when adding the contribution of the flux to the residual at Uk+1. 

The spectral radius in the k-direction is also calculated. Since this is done on a cell basis instead of a face bases, it
is only calculated for the upper cell (Uk+1 in this case). The spectral radius is added to the average wave speed variable
and is eventually used in the time step calculation if the time step isn't explicitly specified.

The velocity and temperature gradients are calculated at each cell face by constructing an alternative control volume centered
around that face as show below.

  ___________________________
  |            |            |
  |            |            |
  |   Ui,j+1   |   Ui+1,j+1 |
  |            |            |
  |_____*******|*******_____|
  |     *      |      *     |
  |     *      |      *     |
  |   Ui,j     |   Ui+1,j   |
  |     *      |      *     |
  |_____*******|*******_____|
  |            |            |
  |            |            |
  |   Ui,j-1   |   Ui+1,j-1 |
  |            |            |
  |____________|____________|

The above diagram shows a 2D schematic of how the gradient calculation is done. In this example the gradient is being calculated at
the face between cell Ui,j and Ui+1,j. The dashes represent the grid cells, and the astrisks represent the alternative control volume.
The grid cells themselves cannot be used as the control volume for the gradients (averaging values at adjacent cells to get gradients at
the face) because this leads to odd/even decoupling. The face areas, volumes, and states at the center of the faces are needed for the 
alternative control volume. The left and right sides of the alternate control volume pass through the center of cells Ui,j and Ui+1,j
respectively. Therefore these values are used for the face states. The left and right face areas are calculated as the average of the face
areas of the cells that they split. For example, the left face area would be calculated as 0.5 * (Ai+1/2,j + Ai-1/2,j). The top and bottom 
sides of the alternative control volume pass through 4 cells each. Therefore the value of the state at the face center is the average of these
four states. For example the state at the top face is calculated as 0.25 * (Ui,j + Ui+1,j + Ui,j+1, Ui+1,j+1). The face areas of the top and
bottom sides are calculated as the average of the 2 face areas that each one passes through. For example, the top face area is calculated as
0.5 * (Ai,j+1/2 + Ai+1,j+1/2). In three dimensions each gradient calculation touches the values at 10 cells (6 shown and 4 more in/out of 
the page). The stencil for the gradients of all faces in a cell touches 15 cells. The gradient calculation with this stencil uses the "edge" 
ghost cells, but not the "corner" ghost cells.
*/
void procBlock::CalcViscFluxK(const sutherland &suth, const idealGas &eqnState, const input &inp){
  // suth -- method to get viscosity as a function of temperature (Sutherland's law)
  // eqnState -- equation of state
  // inp -- all input variables

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK() + 1; //calculating fluxes on k-faces so one more face in k-direction

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();

  //calculate reference Reynolds number and Mach number for nondimensionalization
  double Re = inp.RRef() * inp.VelRef().Mag() * inp.LRef() / suth.MuRef();
  double aRef = eqnState.GetSoS( inp.PRef(), inp.RRef() );
  double mRef = inp.VelRef().Mag() / aRef;

  //coefficient for viscous spectral radii
  double vCoeff = 1.0;

  //loop over physical faces
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++){   
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++){    
      for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++){      

	//location of current face (with ghost cells included)
	int loc = GetLoc1D(ii, jj, kk, imaxG, jmaxG);

	//location of faces in the upper and lower k-direction (with ghost cells included)
	int fUpk = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int fLowk = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);

	//location of i-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower k-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpiUpk = GetUpperFaceI(ii, jj, kk, imaxG, jmaxG);
	int fUpiLowk = GetUpperFaceI(ii, jj, kk - 1, imaxG, jmaxG);
	int fLowiUpk = GetLowerFaceI(ii, jj, kk, imaxG, jmaxG);
	int fLowiLowk = GetLowerFaceI(ii, jj, kk - 1, imaxG, jmaxG);

	//location of j-faces in the upper and lower direction belonging to the cells in the upper 
	//and lower k-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int fUpjUpk = GetUpperFaceJ(ii, jj, kk, imaxG, jmaxG);
	int fUpjLowk = GetUpperFaceJ(ii, jj, kk - 1, imaxG, jmaxG);
	int fLowjUpk = GetLowerFaceJ(ii, jj, kk, imaxG, jmaxG);
	int fLowjLowk = GetLowerFaceJ(ii, jj, kk - 1, imaxG, jmaxG);

	//location of cells in the upper and lower k-direction with respect to baseline face (with ghost cells included)
	int kLow  = GetCellFromFaceLowerK(ii, jj, kk, imaxG, jmaxG);
	int kUp  = GetCellFromFaceUpperK(ii, jj, kk, imaxG, jmaxG);

	//location of cells in the upper and lower k-direction and the upper 
	//and lower i-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int iUpkUp = GetNeighborUpI(ii, jj, kk, imaxG, jmaxG);
	int iUpkLow = GetNeighborUpI(ii, jj, kk - 1, imaxG, jmaxG);
	int iLowkUp = GetNeighborLowI(ii, jj, kk, imaxG, jmaxG);
	int iLowkLow = GetNeighborLowI(ii, jj, kk - 1, imaxG, jmaxG);

	//location of cells in the upper and lower k-direction and the upper 
	//and lower j-direction of the current face (with ghost cells included) - these are used in the 
	//gradient calculation to construct the alternate control volume
	int jUpkUp = GetNeighborUpK(ii, jj, kk, imaxG, jmaxG);
	int jUpkLow = GetNeighborUpK(ii, jj, kk - 1, imaxG, jmaxG);
	int jLowkUp = GetNeighborLowK(ii, jj, kk, imaxG, jmaxG);
	int jLowkLow = GetNeighborLowK(ii, jj, kk - 1, imaxG, jmaxG);

	//location of cells in the upper and lower k-direction with respect to baseline face (without ghost cells included)
	int kLowNG  = GetCellFromFaceLowerK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);
	int kUpNG  = GetCellFromFaceUpperK(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), imax, jmax);

	//calculate average velocity on i and j faces of alternate control volume
	vector3d<double> viu = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(iUpkUp).Velocity() + (*this).State(iUpkLow).Velocity() );
	vector3d<double> vil = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(iLowkUp).Velocity() + (*this).State(iLowkLow).Velocity() );

	vector3d<double> vju = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(jUpkUp).Velocity() + (*this).State(jUpkLow).Velocity() );
	vector3d<double> vjl = 0.25 * ( (*this).State(kLow).Velocity() + (*this).State(kUp).Velocity() + (*this).State(jLowkUp).Velocity() + (*this).State(jLowkLow).Velocity() );

	//calculate areas of faces in alternate control volume
	vector3d<double> aku = 0.5 * ( (*this).FAreaK(loc) + (*this).FAreaK(fUpk) );
	vector3d<double> akl = 0.5 * ( (*this).FAreaK(loc) + (*this).FAreaK(fLowk) );

	vector3d<double> aiu = 0.5 * ( (*this).FAreaI(fUpiUpk) + (*this).FAreaI(fUpiLowk) );
	vector3d<double> ail = 0.5 * ( (*this).FAreaI(fLowiUpk) + (*this).FAreaI(fLowiLowk) );

	vector3d<double> aju = 0.5 * ( (*this).FAreaJ(fUpjUpk) + (*this).FAreaJ(fUpjLowk) );
	vector3d<double> ajl = 0.5 * ( (*this).FAreaJ(fLowjUpk) + (*this).FAreaJ(fLowjLowk) );

	//calculate volume of alternate control volume
	double vol = 0.5 * ( (*this).Vol(kLow) + (*this).Vol(kUp) );

	//Get velocity gradient at face
	tensor<double> velGrad = CalcVelGradGG( vil, viu, vjl, vju, (*this).State(kLow).Velocity(), (*this).State(kUp).Velocity(), ail, aiu, ajl, aju, akl, aku, vol);
	//Get velocity at face
	vector3d<double> vel = FaceReconCentral( (*this).State(kLow).Velocity(), (*this).State(kUp).Velocity(), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );

	//calculate average temperature on i and j faces of alternate control volume
	double tiu = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(iUpkUp).Temperature(eqnState) +
			      (*this).State(iUpkLow).Temperature(eqnState) );
	double til = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(iLowkUp).Temperature(eqnState) +
			      (*this).State(iLowkLow).Temperature(eqnState) );

	double tju = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(jUpkUp).Temperature(eqnState) +
			      (*this).State(jUpkLow).Temperature(eqnState) );
	double tjl = 0.25 * ( (*this).State(kLow).Temperature(eqnState) + (*this).State(kUp).Temperature(eqnState) + (*this).State(jLowkUp).Temperature(eqnState) +
			      (*this).State(jLowkLow).Temperature(eqnState) );

	//Get temperature gradient at face
	vector3d<double> tGrad = CalcTempGradGG( til, tiu, tjl, tju, (*this).State(kLow).Temperature(eqnState), (*this).State(kUp).Temperature(eqnState), ail, aiu, ajl, aju, akl, aku, vol);
	//Get viscosity at face
	double mu = FaceReconCentral( suth.GetViscosity( (*this).State(kLow).Temperature(eqnState) ), 
				      suth.GetViscosity( (*this).State(kUp).Temperature(eqnState) ), (*this).Center(kLow), (*this).Center(kUp), (*this).FCenterK(loc) );
	mu = mu * (mRef/Re);  //effective viscosity (due to nondimensionalization)

	//calculate viscous flux
	viscousFlux tempViscFlux( velGrad, vel, mu, suth, eqnState, tGrad, (*this).FAreaK(loc) );

	//area vector points from left to right, so add to left cell, subtract from right cell
	//but viscous fluxes are subtracted from inviscid fluxes, so sign is reversed
	if ( kk > (*this).NumGhosts() ){ //at left boundary there is no left cell to add to
	  (*this).AddToResidual(-1.0 * tempViscFlux * (*this).FAreaK(loc).Mag(), kLowNG);
	}
	if ( kk < kmax -1 + (*this).NumGhosts() ){ //at right boundary there is no right cell to add to
	  (*this).AddToResidual(tempViscFlux * (*this).FAreaK(loc).Mag(), kUpNG);

	  //calculate component of wave speed. This is done on a cell by cell basis, so only at the upper faces
	  double maxWS = (mRef/Re) * ViscCellSpectralRadius( (*this).FAreaK(loc), (*this).FAreaK(fUpk), (*this).State(kUp), eqnState, suth, (*this).Vol(kUp) );
	  (*this).SetAvgWaveSpeed( (*this).AvgWaveSpeed(kUpNG) + vCoeff * maxWS, kUpNG);
	}

      }
    }
  }

}

/* Member function to assign geometric quantities such as volume, face area, cell centroid, and face center to ghost cells. This assigns values for
regular ghost cells and "edge" ghost cells. "Corner" cells are left with no value as they are not used.
           ____ ____ ____ ____ ____ ____ ____ ____
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|

In the above diagram where X represents the physical cells, cells marked G (regular ghost cells) and E ("edge" ghost cells) are assigned geometric
values. G1 represents the first layer of ghost cells and G2 represents the second layer.
*/
void procBlock::AssignGhostCellsGeom(const input &inp){
  // inp -- all input variables

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical I faces and assign values for regular ghost cells ------------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(0, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "il");
      string bcNameU = bound.GetBCName(imax, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "iu");

      if ( bcNameL == "interblock" && bcNameU == "interblock" && imax < 2 ){
	cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot have interblock BCs surrounding a block that is 1 cell thick!" << endl;
	exit(0);
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at lower i-boundary
	int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first ghost cell at lower i-boundary
	int lFaceG1_il = GetLowerFaceI(1, jj, kk, imaxG, jmaxG); 
	int lFaceG1_jl = GetLowerFaceJ(1, jj, kk, imaxG, jmaxG); 
	int lFaceG1_kl = GetLowerFaceK(1, jj, kk, imaxG, jmaxG); 

	//location of second ghost cell at lower i-boundary
	int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first ghost cell at lower i-boundary
	int lFaceG2_il = GetLowerFaceI(0, jj, kk, imaxG, jmaxG);
	int lFaceG2_jl = GetLowerFaceJ(0, jj, kk, imaxG, jmaxG);
	int lFaceG2_kl = GetLowerFaceK(0, jj, kk, imaxG, jmaxG);

	//location of first interior cell at lower i-boundary
	int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	//location of upper i and lower j, k faces of first interior cell at lower i-boundary
	int lFaceIn1_iu = GetUpperFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int lFaceIn1_jl = GetLowerFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int lFaceIn1_kl = GetLowerFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//location of second interior cell at lower i-boundary
	int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
	//location of upper i and lower j, k faces of second interior cell at lower i-boundary
	int lFaceIn2_iu = GetUpperFaceI((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
	int lFaceIn2_jl = GetLowerFaceJ((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);
	int lFaceIn2_kl = GetLowerFaceK((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	//location of lower i-boundary face
	int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across i-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across i-boundary
	//first layer of ghost cells
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG1_jl);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG1_kl);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG2_jl);
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG2_kl);
	}
	else{
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_jl), lFaceG2_jl);
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_kl), lFaceG2_kl);
	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at lower i-boundary
	  int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG); 
	  int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

	  //location of upper j-face for interior cells at lower i-boundary
	  int lFaceIn1_ju = GetUpperFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	  int lFaceIn2_ju = GetUpperFaceJ((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_ju);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_ju);
	  }
	  else{
	    (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_ju);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at lower i-boundary
	  int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG); 
	  int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at lower i-boundary
	  int lFaceIn1_ku = GetUpperFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	  int lFaceIn2_ku = GetUpperFaceK((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_ku);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_ku);
	  }
	  else{
	    (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_ku);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	  (*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	}
	else{
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	  (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	(*this).SetFCenterI( (*this).FCenterI(lFaceB) + dist2Move, lFaceG1_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG1_jl);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG1_kl);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	  (*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);
	}
	else{
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceB) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG2_jl);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG2_kl);
	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	  //location of upper j-face for ghost cells at lower i-boundary
	  int lFaceG1_ju = GetUpperFaceJ(1, jj, kk, imaxG, jmaxG); 
	  int lFaceG2_ju = GetUpperFaceJ(0, jj, kk, imaxG, jmaxG);

	  //location of upper j-face for interior cell at lower i-boundary
	  int lFaceIn1_ju = GetUpperFaceJ((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG1_ju);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	    (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_ju) + dist2Move, lFaceG2_ju);
	  }
	  else{
	    dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	    (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG2_ju);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	  //location of upper k-face for ghost cells at lower i-boundary
	  int lFaceG1_ku = GetUpperFaceK(1, jj, kk, imaxG, jmaxG); 
	  int lFaceG2_ku = GetUpperFaceK(0, jj, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at lower i-boundary
	  int lFaceIn1_ku = GetUpperFaceK((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG1_ku);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn1_iu));
	    (*this).SetFCenterK( (*this).FCenterK(lFaceG1_ku) + dist2Move, lFaceG2_ku);
	  }
	  else{
	    dist2Move = (*this).FCenterI(lFaceB) - (*this).FCenterI(lFaceIn2_iu);
	    (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG2_ku);
	  }

	}
      }

      //upper surface ----------------------------------------------------
      if ( bcNameU != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at upper i-boundary
	int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
	//location of upper i and lower j, k faces of first ghost cell at upper i-boundary
	int uFaceG1_iu = GetUpperFaceI(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG1_jl = GetLowerFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
	int uFaceG1_kl = GetLowerFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 

	//location of second ghost cell at upper i-boundary
	int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);
	//location of upper i and lower j, k faces of second ghost cell at upper i-boundary
	int uFaceG2_iu = GetUpperFaceI(imaxG-1, jj, kk, imaxG, jmaxG);
	int uFaceG2_jl = GetLowerFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);
	int uFaceG2_kl = GetLowerFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

	//location of first interior cell at upper i-boundary
	int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first interior cell at upper i-boundary
	int uFaceIn1_il = GetLowerFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int uFaceIn1_jl = GetLowerFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	int uFaceIn1_kl = GetLowerFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//location of second interior cell at upper i-boundary
	int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of second interior cell at upper i-boundary
	int uFaceIn2_il = GetLowerFaceI(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	int uFaceIn2_jl = GetLowerFaceJ(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	int uFaceIn2_kl = GetLowerFaceK(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	//location of upper i-boundary face
	int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across i-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across i-boundary
	//first layer of ghost cells
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_iu);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_jl);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_kl);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_iu);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_jl);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_kl);
	}
	else{
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_iu);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_jl);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_kl);
	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at upper i-boundary
	  int uFaceG1_ju = GetUpperFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
	  int uFaceG2_ju = GetUpperFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);

	  //location of upper j-face for interior cells at upper i-boundary
	  int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	  int uFaceIn2_ju = GetUpperFaceJ(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG1_ju);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG2_ju);
	  }
	  else{
	    (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_ju), uFaceG2_ju);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at upper i-boundary
	  int uFaceG1_ku = GetUpperFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 
	  int uFaceG2_ku = GetUpperFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at upper i-boundary
	  int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 
	  int uFaceIn2_ku = GetUpperFaceK(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG1_ku);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG2_ku);
	  }
	  else{
	    (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_ku), uFaceG2_ku);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	  (*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
	}
	else{
	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	  (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	(*this).SetFCenterI( (*this).FCenterI(uFaceB) + dist2Move, uFaceG1_iu);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG1_jl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG1_kl);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	  (*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_jl) + dist2Move, uFaceG2_jl);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceG1_kl) + dist2Move, uFaceG2_kl);
	}
	else{
	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceB) + dist2Move, uFaceG2_iu);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG2_jl);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG2_kl);
	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper face areas too

	  //location of upper j-face for ghost cells at lower i-boundary
	  int uFaceG1_ju = GetUpperFaceJ(imaxG-2, jj, kk, imaxG, jmaxG); 
	  int uFaceG2_ju = GetUpperFaceJ(imaxG-1, jj, kk, imaxG, jmaxG);

	  //location of upper j-face for interior cell at upper i-boundary
	  int uFaceIn1_ju = GetUpperFaceJ(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG1_ju);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	    (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);
	  }
	  else{
	    dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	    (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG2_ju);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper face areas too

	  //location of upper j-face for ghost cells at upper i-boundary
	  int uFaceG1_ku = GetUpperFaceK(imaxG-2, jj, kk, imaxG, jmaxG); 
	  int uFaceG2_ku = GetUpperFaceK(imaxG-1, jj, kk, imaxG, jmaxG);

	  //location of upper j-face for interior cells at upper i-boundary
	  int uFaceIn1_ku = GetUpperFaceK(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG1_ku);

	  //second layer of ghost cells
	  if (imax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn1_il));
	    (*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);
	  }
	  else{
	    dist2Move = (*this).FCenterI(uFaceB) - (*this).FCenterI(uFaceIn2_il);
	    (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG2_ku);
	  }

	}

      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical J faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), 0, kk - (*this).NumGhosts(), "jl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jmax, kk - (*this).NumGhosts(), "ju");

      if ( bcNameL == "interblock" && bcNameU == "interblock" && jmax < 2 ){
	cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot have interblock BCs surrounding a block that is 1 cell thick!" << endl;
	exit(0);
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at lower j-boundary
	int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first ghost cell at lower j-boundary
	int lFaceG1_il = GetLowerFaceI(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG1_jl = GetLowerFaceJ(ii, 1, kk, imaxG, jmaxG); 
	int lFaceG1_kl = GetLowerFaceK(ii, 1, kk, imaxG, jmaxG); 

	//location of second ghost cell at lower i-boundary
	int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first ghost cell at lower j-boundary
	int lFaceG2_il = GetLowerFaceI(ii, 0, kk, imaxG, jmaxG);
	int lFaceG2_jl = GetLowerFaceJ(ii, 0, kk, imaxG, jmaxG);
	int lFaceG2_kl = GetLowerFaceK(ii, 0, kk, imaxG, jmaxG);

	//location of first interior cell at lower i-boundary
	int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
	//location of upper j and lower i, k faces of first interior cell at lower j-boundary
	int lFaceIn1_il = GetLowerFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int lFaceIn1_ju = GetUpperFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int lFaceIn1_kl = GetLowerFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//location of second interior cell at lower i-boundary
	int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
	//location of upper j and lower i, k faces of second interior cell at lower j-boundary
	int lFaceIn2_il = GetLowerFaceI(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
	int lFaceIn2_ju = GetUpperFaceJ(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);
	int lFaceIn2_kl = GetLowerFaceK(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	//location of lower j-boundary face
	int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across j-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across j-boundary
	//first layer of ghost cells
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_jl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG1_il);
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG1_kl);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_jl);
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG2_il);
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_kl), lFaceG2_kl);
	}
	else{
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_jl);
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_il), lFaceG2_il);
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_kl), lFaceG2_kl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at lower j-boundary
	  int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG); 
	  int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

	  //location of upper i-face for interior cells at lower j-boundary
	  int lFaceIn1_iu = GetUpperFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	  int lFaceIn2_iu = GetUpperFaceI(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_iu);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_iu);
	  }
	  else{
	    (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_iu);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at lower j-boundary
	  int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG); 
	  int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at upper j-boundary
	  int lFaceIn1_ku = GetUpperFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 
	  int lFaceIn2_ku = GetUpperFaceK(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_ku);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_ku);
	  }
	  else{
	    (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_ku);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	  (*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	}
	else{
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	  (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceB) + dist2Move, lFaceG1_jl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG1_il);
	(*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG1_kl);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);
	}
	else{
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceB) + dist2Move, lFaceG2_jl);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_kl) + dist2Move, lFaceG2_kl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at lower j-boundary
	  int lFaceG1_iu = GetUpperFaceI(ii, 1, kk, imaxG, jmaxG); 
	  int lFaceG2_iu = GetUpperFaceI(ii, 0, kk, imaxG, jmaxG);

	  //location of upper i-face for interior cells at lower j-boundary
	  int lFaceIn1_iu = GetUpperFaceI(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG1_iu);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	    (*this).SetFCenterI( (*this).FCenterI(lFaceG1_iu) + dist2Move, lFaceG2_iu);
	  }
	  else{
	    dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	    (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG2_iu);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at lower j-boundary
	  int lFaceG1_ku = GetUpperFaceK(ii, 1, kk, imaxG, jmaxG); 
	  int lFaceG2_ku = GetUpperFaceK(ii, 0, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at lower j-boundary
	  int lFaceIn1_ku = GetUpperFaceK(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG1_ku);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn1_ju));
	    (*this).SetFCenterK( (*this).FCenterK(lFaceG1_ku) + dist2Move, lFaceG2_ku);
	  }
	  else{
	    dist2Move = (*this).FCenterJ(lFaceB) - (*this).FCenterJ(lFaceIn2_ju);
	    (*this).SetFCenterK( (*this).FCenterK(lFaceIn1_ku) + dist2Move, lFaceG2_ku);
	  }

	}
      }

      //upper surface ----------------------------------------------------
      if ( bcNameU != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at upper j-boundary
	int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
	//location of upper j and lower i, k faces of first ghost cell at lower j-boundary
	int uFaceG1_il = GetLowerFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG1_ju = GetUpperFaceJ(ii, jmaxG-2, kk, imaxG, jmaxG); 
	int uFaceG1_kl = GetLowerFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 

	//location of second ghost cell at upper j-boundary
	int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);
	//location of upper j and lower i, k faces of second ghost cell at lower j-boundary
	int uFaceG2_il = GetLowerFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);
	int uFaceG2_ju = GetUpperFaceJ(ii, jmaxG-1, kk, imaxG, jmaxG);
	int uFaceG2_kl = GetLowerFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

	//location of first interior cell at upper j-boundary
	int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
	//location of lower i, j, k faces of first interior cell at lower j-boundary
	int uFaceIn1_il = GetLowerFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int uFaceIn1_jl = GetLowerFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	int uFaceIn1_kl = GetLowerFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//location of second interior cell at upper j-boundary
	int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
	//location of lower i, j, k faces of second interior cell at lower j-boundary
	int uFaceIn2_il = GetLowerFaceI(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
	int uFaceIn2_jl = GetLowerFaceJ(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);
	int uFaceIn2_kl = GetLowerFaceK(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

	//location of upper j-boundary face
	int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across j-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across j-boundary
	//first layer of ghost cells
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_ju);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_il);
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_kl);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_ju);
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_il);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_kl);
	}
	else{
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_ju);
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_il);
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_kl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at upper j-boundary
	  int uFaceG1_iu = GetUpperFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
	  int uFaceG2_iu = GetUpperFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);

	  //location of upper i-face for interior cells at upper j-boundary
	  int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	  int uFaceIn2_iu = GetUpperFaceI(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG1_iu);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG2_iu);
	  }
	  else{
	    (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_iu), uFaceG2_iu);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at upper j-boundary
	  int uFaceG1_ku = GetUpperFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 
	  int uFaceG2_ku = GetUpperFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at upper j-boundary
	  int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 
	  int uFaceIn2_ku = GetUpperFaceK(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG1_ku);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_ku), uFaceG2_ku);
	  }
	  else{
	    (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_ku), uFaceG2_ku);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	  (*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
	}
	else{
	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	  (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceB) + dist2Move, uFaceG1_ju);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG1_il);
	(*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG1_kl);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceG1_il) + dist2Move, uFaceG2_il);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceG1_kl) + dist2Move, uFaceG2_kl);
	}
	else{
	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceB) + dist2Move, uFaceG2_ju);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG2_il);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_kl) + dist2Move, uFaceG2_kl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at upper j-boundary
	  int uFaceG1_iu = GetUpperFaceI(ii, jmaxG-2, kk, imaxG, jmaxG); 
	  int uFaceG2_iu = GetUpperFaceI(ii, jmaxG-1, kk, imaxG, jmaxG);

	  //location of upper i-face for interior cells at upper j-boundary
	  int uFaceIn1_iu = GetUpperFaceI(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG1_iu);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	    (*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	  }
	  else{
	    dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	    (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG2_iu);
	  }

	}

	if ( kk == kmax - 1 + (*this).NumGhosts() ){ //at end of k-line of cells assign cell upper k-face areas too

	  //location of upper k-face for ghost cells at upper j-boundary
	  int uFaceG1_ku = GetUpperFaceK(ii, jmaxG-2, kk, imaxG, jmaxG); 
	  int uFaceG2_ku = GetUpperFaceK(ii, jmaxG-1, kk, imaxG, jmaxG);

	  //location of upper k-face for interior cells at upper j-boundary
	  int uFaceIn1_ku = GetUpperFaceK(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG1_ku);

	  //second layer of ghost cells
	  if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn1_jl));
	    (*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);
	  }
	  else{
	    dist2Move = (*this).FCenterJ(uFaceB) - (*this).FCenterJ(uFaceIn2_jl);
	    (*this).SetFCenterK( (*this).FCenterK(uFaceIn1_ku) + dist2Move, uFaceG2_ku);
	  }

	}
      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical K faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), 0, "kl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kmax, "ku");

      if ( bcNameL == "interblock" && bcNameU == "interblock" && kmax < 2 ){
	cerr << "ERROR: Error in procBlock::AssignGhostCellsGeom(). Cannot have interblock BCs surrounding a block that is 1 cell thick!" << endl;
	exit(0);
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at lower k-boundary
	int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
	//location of lower i, j, k faces of first ghost cell at lower k-boundary
	int lFaceG1_il = GetLowerFaceI(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG1_jl = GetLowerFaceJ(ii, jj, 1, imaxG, jmaxG); 
	int lFaceG1_kl = GetLowerFaceK(ii, jj, 1, imaxG, jmaxG); 

	//location of second ghost cell at lower k-boundary
	int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);
	//location of lower i, j, k faces of second ghost cell at lower k-boundary
	int lFaceG2_il = GetLowerFaceI(ii, jj, 0, imaxG, jmaxG);
	int lFaceG2_jl = GetLowerFaceJ(ii, jj, 0, imaxG, jmaxG);
	int lFaceG2_kl = GetLowerFaceK(ii, jj, 0, imaxG, jmaxG);

	//location of first interior cell at lower k-boundary
	int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
	//location of upper k and lower i, j faces of first interior cell at lower k-boundary
	int lFaceIn1_il = GetLowerFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	int lFaceIn1_jl = GetLowerFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	int lFaceIn1_ku = GetUpperFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	//location of second interior cell at lower k-boundary
	int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
	//location of upper k and lower i, j faces of second interior cell at lower k-boundary
	int lFaceIn2_il = GetLowerFaceI(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
	int lFaceIn2_jl = GetLowerFaceJ(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);
	int lFaceIn2_ku = GetUpperFaceK(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	//location of lower k-boundary face
	int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across k-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellLowIn1), cellLowG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellLowIn1), cellLowG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellLowIn2), cellLowG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across k-boundary
	//first layer of ghost cells
	(*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG1_kl);
	(*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG1_il);
	(*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG1_jl);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn1_ku), lFaceG2_kl);
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_il), lFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_jl), lFaceG2_jl);
	}
	else{
	  (*this).SetFAreaK( (*this).FAreaK(lFaceIn2_ku), lFaceG2_kl);
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_il), lFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_jl), lFaceG2_jl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at lower k-boundary
	  int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG); 
	  int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

	  //location of upper i-face for interior cells at lower k-boundary
	  int lFaceIn1_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	  int lFaceIn2_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG1_iu);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaI( (*this).FAreaI(lFaceIn1_iu), lFaceG2_iu);
	  }
	  else{
	    (*this).SetFAreaI( (*this).FAreaI(lFaceIn2_iu), lFaceG2_iu);
	  }

	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at lower k-boundary
	  int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG); 
	  int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

	  //location of upper j-face for interior cells at lower k-boundary
	  int lFaceIn1_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 
	  int lFaceIn2_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG1_ju);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn1_ju), lFaceG2_ju);
	  }
	  else{
	    (*this).SetFAreaJ( (*this).FAreaJ(lFaceIn2_ju), lFaceG2_ju);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	(*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	  (*this).SetCenter( (*this).Center(cellLowG1) + dist2Move, cellLowG2);
	}
	else{
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	  (*this).SetCenter( (*this).Center(cellLowIn1) + dist2Move, cellLowG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	(*this).SetFCenterK( (*this).FCenterK(lFaceB) + dist2Move, lFaceG1_kl);
	(*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG1_il);
	(*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG1_jl);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	  (*this).SetFCenterK( (*this).FCenterK(lFaceG1_kl) + dist2Move, lFaceG2_kl);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceG1_il) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_jl) + dist2Move, lFaceG2_jl);
	}
	else{
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	  (*this).SetFCenterK( (*this).FCenterK(lFaceB) + dist2Move, lFaceG2_kl);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_il) + dist2Move, lFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_jl) + dist2Move, lFaceG2_jl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at lower k-boundary
	  int lFaceG1_iu = GetUpperFaceI(ii, jj, 1, imaxG, jmaxG); 
	  int lFaceG2_iu = GetUpperFaceI(ii, jj, 0, imaxG, jmaxG);

	  //location of upper i-face for interior cells at lower k-boundary
	  int lFaceIn1_iu = GetUpperFaceI(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	  (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG1_iu);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	    (*this).SetFCenterI( (*this).FCenterI(lFaceG1_iu) + dist2Move, lFaceG2_iu);
	  }
	  else{
	    dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	    (*this).SetFCenterI( (*this).FCenterI(lFaceIn1_iu) + dist2Move, lFaceG2_iu);
	  }

	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at lower k-boundary
	  int lFaceG1_ju = GetUpperFaceJ(ii, jj, 1, imaxG, jmaxG); 
	  int lFaceG2_ju = GetUpperFaceJ(ii, jj, 0, imaxG, jmaxG);

	  //location of upper j-face for interior cells at lower k-boundary
	  int lFaceIn1_ju = GetUpperFaceJ(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku);
	  (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG1_ju);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn1_ku));
	    (*this).SetFCenterJ( (*this).FCenterJ(lFaceG1_ju) + dist2Move, lFaceG2_ju);
	  }
	  else{
	    dist2Move = (*this).FCenterK(lFaceB) - (*this).FCenterK(lFaceIn2_ku);
	    (*this).SetFCenterJ( (*this).FCenterJ(lFaceIn1_ju) + dist2Move, lFaceG2_ju);
	  }

	}

      }

      //upper surface -----------------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply geometry values for non interblock BCs, for interblock, do nothing

	//location of first ghost cell at upper k-boundary
	int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
	//location of upper k and lower i, j faces of first ghost cell at lower k-boundary
	int uFaceG1_il = GetLowerFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG1_jl = GetLowerFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
	int uFaceG1_ku = GetUpperFaceK(ii, jj, kmaxG-2, imaxG, jmaxG); 

	//location of second ghost cell at upper k-boundary
	int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);
	//location of upper k and lower i, j faces of second ghost cell at lower k-boundary
	int uFaceG2_il = GetLowerFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);
	int uFaceG2_jl = GetLowerFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);
	int uFaceG2_ku = GetUpperFaceK(ii, jj, kmaxG-1, imaxG, jmaxG);

	//location of first interior cell at upper k-boundary
	int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
	//location of lower i, j, k faces of first interior cell at lower k-boundary
	int uFaceIn1_il = GetLowerFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	int uFaceIn1_jl = GetLowerFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	int uFaceIn1_kl = GetLowerFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	//location of second interior cell at upper k-boundary
	int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
	//location of lower i, j, k faces of second interior cell at lower k-boundary
	int uFaceIn2_il = GetLowerFaceI(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
	int uFaceIn2_jl = GetLowerFaceJ(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);
	int uFaceIn2_kl = GetLowerFaceK(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	//location of upper k-boundary face
	int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	//Assign volume ------------------------------------------------------------------------------------------
	//mirror volume values from adjacent cells across k-boundary
	//first layer of ghost cells
	(*this).SetVol( (*this).Vol(cellUpIn1), cellUpG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetVol( (*this).Vol(cellUpIn1), cellUpG2);
	}
	else{
	  (*this).SetVol( (*this).Vol(cellUpIn2), cellUpG2);
	}

	//Assign face areas ----------------------------------------------------------------------------------------
	//mirror face area values from adjacent cells across k-boundary
	//first layer of ghost cells
	(*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG1_ku);
	(*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG1_il);
	(*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG1_jl);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn1_kl), uFaceG2_ku);
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_il), uFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_jl), uFaceG2_jl);
	}
	else{
	  (*this).SetFAreaK( (*this).FAreaK(uFaceIn2_kl), uFaceG2_ku);
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_il), uFaceG2_il);
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_jl), uFaceG2_jl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at upper k-boundary
	  int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
	  int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);

	  //location of upper i-face for interior cells at upper k-boundary
	  int uFaceIn1_iu = GetUpperFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	  int uFaceIn2_iu = GetUpperFaceI(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG1_iu);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaI( (*this).FAreaI(uFaceIn1_iu), uFaceG2_iu);
	  }
	  else{
	    (*this).SetFAreaI( (*this).FAreaI(uFaceIn2_iu), uFaceG2_iu);
	  }

	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at upper k-boundary
	  int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
	  int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);

	  //location of upper j-face for interior cells at upper k-boundary
	  int uFaceIn1_ju = GetUpperFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 
	  int uFaceIn2_ju = GetUpperFaceJ(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	  //mirror face area values from adjacent cells
	  //first layer of ghost cells
	  (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG1_ju);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn1_ju), uFaceG2_ju);
	  }
	  else{
	    (*this).SetFAreaJ( (*this).FAreaJ(uFaceIn2_ju), uFaceG2_ju);
	  }

	}

	//Assign cell centroid --------------------------------------------------------------------------------------------------------------------
	//cell centroid is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	vector3d<double> dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	(*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	  (*this).SetCenter( (*this).Center(cellUpG1) + dist2Move, cellUpG2);
	}
	else{
	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	  (*this).SetCenter( (*this).Center(cellUpIn1) + dist2Move, cellUpG2);
	}

	//Assign face centers --------------------------------------------------------------------------------------------------------------------
	//face center is moved interior cell width in the boundary normal direction
	//first layer of ghost cells
	dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	(*this).SetFCenterK( (*this).FCenterK(uFaceB) + dist2Move, uFaceG1_ku);
	(*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG1_il);
	(*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG1_jl);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	  (*this).SetFCenterK( (*this).FCenterK(uFaceG1_ku) + dist2Move, uFaceG2_ku);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceG1_il) + dist2Move, uFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_jl) + dist2Move, uFaceG2_jl);
	}
	else{
	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	  (*this).SetFCenterK( (*this).FCenterK(uFaceB) + dist2Move, uFaceG2_ku);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_il) + dist2Move, uFaceG2_il);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_jl) + dist2Move, uFaceG2_jl);
	}

	if ( ii == imax - 1 + (*this).NumGhosts() ){ //at end of i-line of cells assign cell upper i-face areas too

	  //location of upper i-face for ghost cells at upper k-boundary
	  int uFaceG1_iu = GetUpperFaceI(ii, jj, kmaxG-2, imaxG, jmaxG); 
	  int uFaceG2_iu = GetUpperFaceI(ii, jj, kmaxG-1, imaxG, jmaxG);

	  //location of upper i-face for interior cells at upper k-boundary
	  int uFaceIn1_iu = GetUpperFaceI(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	  (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG1_iu);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	    (*this).SetFCenterI( (*this).FCenterI(uFaceG1_iu) + dist2Move, uFaceG2_iu);
	  }
	  else{
	    dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	    (*this).SetFCenterI( (*this).FCenterI(uFaceIn1_iu) + dist2Move, uFaceG2_iu);
	  }

	}

	if ( jj == jmax - 1 + (*this).NumGhosts() ){ //at end of j-line of cells assign cell upper j-face areas too

	  //location of upper j-face for ghost cells at upper k-boundary
	  int uFaceG1_ju = GetUpperFaceJ(ii, jj, kmaxG-2, imaxG, jmaxG); 
	  int uFaceG2_ju = GetUpperFaceJ(ii, jj, kmaxG-1, imaxG, jmaxG);

	  //location of upper j-face for interior cells at upper k-boundary
	  int uFaceIn1_ju = GetUpperFaceJ(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	  //face center is moved interior cell width in the boundary normal direction
	  //first layer of ghost cells
	  dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl);
	  (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG1_ju);

	  //second layer of ghost cells
	  if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	    dist2Move = 2.0 * ((*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn1_kl));
	    (*this).SetFCenterJ( (*this).FCenterJ(uFaceG1_ju) + dist2Move, uFaceG2_ju);
	  }
	  else{
	    dist2Move = (*this).FCenterK(uFaceB) - (*this).FCenterK(uFaceIn2_kl);
	    (*this).SetFCenterJ( (*this).FCenterJ(uFaceIn1_ju) + dist2Move, uFaceG2_ju);
	  }

	}

      }

    }
  }

  //fill ghost cell edge lines with geometric values
  //(*this).AssignGhostCellsGeomEdge(inp);

}
 
/* Member function to assign geometric quantities such as volume, face area, cell centroid, and face center to ghost cells located on the 12 block edges. 
Assumes AssignGhostCellsGeom has already been run. 

           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J 

In the above diagram the cells marked X represent physical cells. Cells marked G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells (1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to 
the values of the ghost cells they tough at level "e". The values at edge cell 2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it touches.
*/
void procBlock::AssignGhostCellsGeomEdge( const input& inp ){

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction -------------------------------
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int j1,k1,j2,k2,je,ke;

      int gfe_j1_k1_jl, gfe_j1_k1_kl;
      int gfe_j1_k2_jl, gfe_j1_k2_kl;
      int gfe_j2_k1_jl, gfe_j2_k1_kl;
      int gfe_j2_k2_jl, gfe_j2_k2_kl;

      int gf_j1_ke_jl, gf_j1_ke_kl, gf_j1_ke_ku;
      int gf_j2_ke_jl, gf_j2_ke_kl;
      int gf_je_k1_jl, gf_je_k1_ju;
      int gf_je_k2_jl, gf_je_k2_kl;

      string bc_J, bc_K, surfJ, surfK;

      if ( cc == 0 ){ //at jl/kl edge - ghost cells are in the lower direction of both j and k, so use GetLowerFace for both
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();
	
	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, j# = j position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG); 
	gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG); 
	gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG); 
	gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG); 
	gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG);

	//ghost face, on first layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG);  
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);
	gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);

	//ghost face, on second layer of j line of cells, on non-edge layer of k line of cells
	gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG); 
	gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG); 

	//ghost face, on non-edge layer of j line of cells, on second layer of k line of cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG); 
	gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG); 

	surfJ = "jl";
	surfK = "kl";

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at jl/ku edge - ghost cells are in the lower direction of j and upper direction of k, so use GetLowerFace for J
	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, j# = j position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_jl = GetLowerFaceJ(ii, j1, k1, imaxG, jmaxG); 
	gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_jl = GetLowerFaceJ(ii, j1, k2, imaxG, jmaxG); 
	gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_jl = GetLowerFaceJ(ii, j2, k1, imaxG, jmaxG); 
	gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_jl = GetLowerFaceJ(ii, j2, k2, imaxG, jmaxG); 
	gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG); 

	//ghost face, on first layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_jl = GetLowerFaceJ(ii, j1, ke, imaxG, jmaxG); 
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG); 
	gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG); 

	//ghost face, on second layer of j line of cells, on non-edge layer of k line of cells
	gf_j2_ke_jl = GetLowerFaceJ(ii, j2, ke, imaxG, jmaxG); 
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG); 

	//ghost face, on non-edge layer of j line of cells, on first layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k1_ju = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on second layer of k line of cells
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG); 
	gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG); 

	surfJ = "jl";
	surfK = "ku";

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at ju/kl edge - ghost cells are in the lower direction of k, and upper direction of j so use GetLowerFace for k
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, j# = j position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG); 
	gfe_j1_k1_kl = GetLowerFaceK(ii, j1, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG); 
	gfe_j1_k2_kl = GetLowerFaceK(ii, j1, k2, imaxG, jmaxG);

	//ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG); 
	gfe_j2_k1_kl = GetLowerFaceK(ii, j2, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG); 
	gfe_j2_k2_kl = GetLowerFaceK(ii, j2, k2, imaxG, jmaxG);

	//ghost face, on first layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG); 
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG); 
	gf_j1_ke_ku = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG); 

	//ghost face, on second layer of j line of cells, on non-edge layer of k line of cells
	gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG); 
	gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG); 

	//ghost face, on non-edge layer of j line of cells, on second layer of k line of cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  
	gf_je_k2_kl = GetLowerFaceK(ii, je, k2, imaxG, jmaxG);  

	surfJ = "ju";
	surfK = "kl";

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(),surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at ju/ku edge - ghost cells are in the upper direction of both j and k, use GetUpperFace for both
	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, j# = j position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	gfe_j1_k1_jl = GetUpperFaceJ(ii, j1, k1, imaxG, jmaxG); 
	gfe_j1_k1_kl = GetUpperFaceK(ii, j1, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	gfe_j1_k2_jl = GetUpperFaceJ(ii, j1, k2, imaxG, jmaxG); 
	gfe_j1_k2_kl = GetUpperFaceK(ii, j1, k2, imaxG, jmaxG);

	//ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	gfe_j2_k1_jl = GetUpperFaceJ(ii, j2, k1, imaxG, jmaxG); 
	gfe_j2_k1_kl = GetUpperFaceK(ii, j2, k1, imaxG, jmaxG);

	//ghost face on edge, on second layer of j line of cells, on second layer of k line of cells
	gfe_j2_k2_jl = GetUpperFaceJ(ii, j2, k2, imaxG, jmaxG); 
	gfe_j2_k2_kl = GetUpperFaceK(ii, j2, k2, imaxG, jmaxG);

	//ghost face, on first layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_jl = GetUpperFaceJ(ii, j1, ke, imaxG, jmaxG);  
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j1_ke_ku = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  

	//ghost face, on second layer of j line of cells, on non-edge layer of k line of cells
	gf_j2_ke_jl = GetUpperFaceJ(ii, j2, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG); 
	gf_je_k1_ju = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG); 

	//ghost face, on non-edge layer of j line of cells, on second layer of k line of cells
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  
	gf_je_k2_kl = GetUpperFaceK(ii, je, k2, imaxG, jmaxG);  

	surfJ = "ju";
	surfK = "ku";

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_J != "interblock" && bc_K != "interblock" ){ 

	//cell indices and remaining face indices
	int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
	int gfe_j1_k1_il = GetLowerFaceI(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells

	int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
	int gfe_j1_k2_il = GetLowerFaceI(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells

	int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
	int gfe_j2_k1_il = GetLowerFaceI(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells

	int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells
	int gfe_j2_k2_il = GetLowerFaceI(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);          //ghost cell, on first layer of j line of cells, on non-edge layer of k line of cells
	int gf_j1_ke_il = GetLowerFaceI(ii, j1, ke, imaxG, jmaxG);  //ghost face, on first layer of j line of cells, on non-edge layer of k line of cells

	int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);          //ghost cell, on second layer of j line of cells, on non-edge layer of k line of cells
	int gf_j2_ke_il = GetLowerFaceI(ii, j2, ke, imaxG, jmaxG);  //ghost face, on second layer of j line of cells, on non-edge layer of k line of cells

	int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);          //ghost cell, on non-edge layer of j line of cells, on first layer of k line of cells
	int gf_je_k1_il = GetLowerFaceI(ii, je, k1, imaxG, jmaxG);  //ghost face, on non-edge layer of j line of cells, on first layer of k line of cells

	int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);          //ghost cell, on non-edge layer of j line of cells, on second layer of k line of cells
	int gf_je_k2_il = GetLowerFaceI(ii, je, k2, imaxG, jmaxG);  //ghost face, on non-edge layer of j line of cells, on second layer of k line of cells

	//volume -----------------------------------------------------------------------------------------------------------------------------------
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_j1_ke) + (*this).Vol(gc_je_k1) ) ,gce_j1_k1);
	(*this).SetVol( (*this).Vol(gc_j2_ke) ,gce_j2_k1);
	(*this).SetVol( (*this).Vol(gc_je_k2) ,gce_j1_k2);
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_j2_ke) + (*this).Vol(gc_je_k2) ) ,gce_j2_k2);

	//face areas --------------------------------------------------------------------------------------------------------------------------------
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

	//centroids -------------------------------------------------------------------------------------------------------------------------------------
	//edge centroid is moved distance of cell width normal to face dividing regular and edge ghost cells
	vector3d<double> dist2MoveK = (*this).FCenterK(gf_j1_ke_kl) - (*this).FCenterK(gf_j1_ke_ku) ;
	vector3d<double> dist2MoveJ = (*this).FCenterJ(gf_je_k1_jl) - (*this).FCenterJ(gf_je_k1_ju) ;
	(*this).SetCenter( (*this).Center(gc_j1_ke) + dist2MoveK, gce_j1_k1);
	(*this).SetCenter( (*this).Center(gc_j2_ke) + dist2MoveK, gce_j2_k1);
	(*this).SetCenter( (*this).Center(gc_je_k2) + dist2MoveJ, gce_j1_k2);
	(*this).SetCenter( (*this).Center(gc_je_k2) + 2.0 * dist2MoveJ, gce_j2_k2);

	//face centers -----------------------------------------------------------------------------------------------------------------------------------
	//edge face centers are moved distance of cell width normal to face dividing regular and edge ghost cells
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
	  vector3d<double> dist2MoveI = (*this).FCenterI(gfe_j1_k1_il) - (*this).FCenterI(gfe_j1_k1_2il); // i-width of cell adjacent to edge

	  int gfe_j1_k1_iu = GetUpperFaceI(ii, j1, k1, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on first layer of k line of cells
	  int gfe_j1_k2_iu = GetUpperFaceI(ii, j1, k2, imaxG, jmaxG); //ghost face on edge, on first layer of j line of cells, on second layer of k line of cells
	  int gfe_j2_k1_iu = GetUpperFaceI(ii, j2, k1, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on first layer of k line of cells
	  int gfe_j2_k2_iu = GetUpperFaceI(ii, j2, k2, imaxG, jmaxG); //ghost face on edge, on second layer of j line of cells, on second layer of k line of cells

	  //face areas
	  (*this).SetFAreaI( (*this).FAreaI(gfe_j1_k1_il), gfe_j1_k1_iu);
	  (*this).SetFAreaI( (*this).FAreaI(gfe_j1_k2_il), gfe_j1_k2_iu);
	  (*this).SetFAreaI( (*this).FAreaI(gfe_j2_k1_il), gfe_j2_k1_iu);
	  (*this).SetFAreaI( (*this).FAreaI(gfe_j2_k2_il), gfe_j2_k2_iu);

	  //face centers
	  (*this).SetFCenterI( (*this).FCenterI(gfe_j1_k1_il) + dist2MoveI, gfe_j1_k1_iu);
	  (*this).SetFCenterI( (*this).FCenterI(gfe_j1_k2_il) + dist2MoveI, gfe_j1_k2_iu);
	  (*this).SetFCenterI( (*this).FCenterI(gfe_j2_k1_il) + dist2MoveI, gfe_j2_k1_iu);
	  (*this).SetFCenterI( (*this).FCenterI(gfe_j2_k2_il) + dist2MoveI, gfe_j2_k2_iu);

	}

      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction --------------------------------
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,k1,i2,k2,ie,ke;

      int gfe_i1_k1_il, gfe_i1_k1_kl;
      int gfe_i1_k2_il, gfe_i1_k2_kl;
      int gfe_i2_k1_il, gfe_i2_k1_kl;
      int gfe_i2_k2_il, gfe_i2_k2_kl;

      int gf_i1_ke_il, gf_i1_ke_ku, gf_i1_ke_kl;
      int gf_i2_ke_il, gf_i2_ke_kl;
      int gf_ie_k1_il, gf_ie_k1_iu;
      int gf_ie_k2_il, gf_ie_k2_kl;

      string bc_I, bc_K, surfI, surfK;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both i and k, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG); 
	gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG); 
	gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG); 
	gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG);

	//ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG); 
	gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of k line of cells
	gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG); 
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG); 

	//ghost face, on non-edge layer of i line of cells, on first layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of k line of cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  
	gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);  

	surfI = "il";
	surfK = "kl";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of i and upper direction of k, so use GetLowerFace for J
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_il = GetLowerFaceI(i1, jj, k1, imaxG, jmaxG); 
	gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_il = GetLowerFaceI(i1, jj, k2, imaxG, jmaxG); 
	gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_il = GetLowerFaceI(i2, jj, k1, imaxG, jmaxG); 
	gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_il = GetLowerFaceI(i2, jj, k2, imaxG, jmaxG); 
	gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_il = GetLowerFaceI(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of k line of cells
	gf_i2_ke_il = GetLowerFaceI(i2, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k1_iu = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of k line of cells
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG); 
	gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG); 

	surfI = "il";
	surfK = "ku";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at iu/kl edge - ghost cells are in the lower direction of k, and upper direction of i so use GetLowerFace for k
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	k2 = 0;
	k1 = 1;
	ke = (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG); 
	gfe_i1_k1_kl = GetLowerFaceK(i1, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG); 
	gfe_i1_k2_kl = GetLowerFaceK(i1, jj, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG); 
	gfe_i2_k1_kl = GetLowerFaceK(i2, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG); 
	gfe_i2_k2_kl = GetLowerFaceK(i2, jj, k2, imaxG, jmaxG);

	//ghost face, on first layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG); 
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG); 
	gf_i1_ke_ku = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG); 

	//ghost face, on second layer of i line of cells, on non-edge layer of k line of cells
	gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG); 
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG); 

	//ghost face, on non-edge layer of i line of cells, on first layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of k line of cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  
	gf_ie_k2_kl = GetLowerFaceK(ie, jj, k2, imaxG, jmaxG);  

	surfI = "iu";
	surfK = "kl";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at iu/ku edge - ghost cells are in the upper direction of both i and k, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	k2 = kmaxG - 1;
	k1 = kmaxG - 2;
	ke = kmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, k# = k position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	gfe_i1_k1_il = GetUpperFaceI(i1, jj, k1, imaxG, jmaxG); 
	gfe_i1_k1_kl = GetUpperFaceK(i1, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	gfe_i1_k2_il = GetUpperFaceI(i1, jj, k2, imaxG, jmaxG); 
	gfe_i1_k2_kl = GetUpperFaceK(i1, jj, k2, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	gfe_i2_k1_il = GetUpperFaceI(i2, jj, k1, imaxG, jmaxG); 
	gfe_i2_k1_kl = GetUpperFaceK(i2, jj, k1, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of k line of cells
	gfe_i2_k2_il = GetUpperFaceI(i2, jj, k2, imaxG, jmaxG); 
	gfe_i2_k2_kl = GetUpperFaceK(i2, jj, k2, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_il = GetUpperFaceI(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i1_ke_ku = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of k line of cells
	gf_i2_ke_il = GetUpperFaceI(i2, jj, ke, imaxG, jmaxG); 
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG); 

	//ghost face, on non-edge layer of i line of cells, on first layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k1_iu = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of k line of cells
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  
	gf_ie_k2_kl = GetUpperFaceK(ie, jj, k2, imaxG, jmaxG);  

	surfI = "iu";
	surfK = "ku";

	//boundary conditioins at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_I != "interblock" && bc_K != "interblock" ){ 

	//cell indices and remaining face indices
	int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
	int gfe_i1_k1_jl = GetLowerFaceJ(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells

	int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
	int gfe_i1_k2_jl = GetLowerFaceJ(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells

	int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
	int gfe_i2_k1_jl = GetLowerFaceJ(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells

	int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells
	int gfe_i2_k2_jl = GetLowerFaceJ(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);          //ghost cell, on first layer of i line of cells, on non-edge layer of k line of cells
	int gf_i1_ke_jl = GetLowerFaceJ(i1, jj, ke, imaxG, jmaxG);  //ghost face, on first layer of i line of cells, on non-edge layer of k line of cells

	int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);          //ghost cell, on second layer of i line of cells, on non-edge layer of k line of cells
	int gf_i2_ke_jl = GetLowerFaceJ(i2, jj, ke, imaxG, jmaxG);  //ghost face, on second layer of i line of cells, on non-edge layer of k line of cells

	int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);          //ghost cell, on non-edge layer of i line of cells, on first layer of k line of cells
	int gf_ie_k1_jl = GetLowerFaceJ(ie, jj, k1, imaxG, jmaxG);  //ghost face, on non-edge layer of i line of cells, on first layer of k line of cells

	int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);          //ghost cell, on non-edge layer of i line of cells, on second layer of k line of cells
	int gf_ie_k2_jl = GetLowerFaceJ(ie, jj, k2, imaxG, jmaxG);  //ghost face, on non-edge layer of i line of cells, on second layer of k line of cells

	//volume ----------------------------------------------------------------------------------------------------------------------------------------------
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_i1_ke) + (*this).Vol(gc_ie_k1) ) ,gce_i1_k1);
	(*this).SetVol( (*this).Vol(gc_i2_ke) ,gce_i2_k1);
	(*this).SetVol( (*this).Vol(gc_ie_k2) ,gce_i1_k2);
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_i2_ke) + (*this).Vol(gc_ie_k2) ) ,gce_i2_k2);

	//face areas ------------------------------------------------------------------------------------------------------------------------------------------
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

	//centroids --------------------------------------------------------------------------------------------------------------------------------------------
	//edge centroid is moved distance of cell width normal to face dividing regular and edge ghost cells
	vector3d<double> dist2MoveK = (*this).FCenterK(gf_i1_ke_kl) - (*this).FCenterK(gf_i1_ke_ku);
	vector3d<double> dist2MoveI = (*this).FCenterI(gf_ie_k1_il) - (*this).FCenterI(gf_ie_k1_iu);
	(*this).SetCenter( (*this).Center(gc_i1_ke) + dist2MoveK, gce_i1_k1);
	(*this).SetCenter( (*this).Center(gc_i2_ke) + dist2MoveK, gce_i2_k1);
	(*this).SetCenter( (*this).Center(gc_ie_k2) + dist2MoveI, gce_i1_k2);
	(*this).SetCenter( (*this).Center(gc_ie_k2) + 2.0 * dist2MoveI, gce_i2_k2);

	//face centers ------------------------------------------------------------------------------------------------------------------------------------------
	//edge face centers are moved distance of cell width normal to face dividing regular and edge ghost cells
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
	  vector3d<double> dist2MoveJ = (*this).FCenterJ(gfe_i1_k1_jl) - (*this).FCenterJ(gfe_i1_k1_2jl); //j-width of adjacent cell

	  int gfe_i1_k1_ju = GetUpperFaceJ(i1, jj, k1, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of k line of cells
	  int gfe_i1_k2_ju = GetUpperFaceJ(i1, jj, k2, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of k line of cells
	  int gfe_i2_k1_ju = GetUpperFaceJ(i2, jj, k1, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of k line of cells
	  int gfe_i2_k2_ju = GetUpperFaceJ(i2, jj, k2, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of k line of cells

	  //face areas
	  (*this).SetFAreaJ( (*this).FAreaJ(gfe_i1_k1_jl), gfe_i1_k1_ju);
	  (*this).SetFAreaJ( (*this).FAreaJ(gfe_i1_k2_jl), gfe_i1_k2_ju);
	  (*this).SetFAreaJ( (*this).FAreaJ(gfe_i2_k1_jl), gfe_i2_k1_ju);
	  (*this).SetFAreaJ( (*this).FAreaJ(gfe_i2_k2_jl), gfe_i2_k2_ju);

	  //face centers
	  (*this).SetFCenterJ( (*this).FCenterJ(gfe_i1_k1_jl) + dist2MoveJ, gfe_i1_k1_ju);
	  (*this).SetFCenterJ( (*this).FCenterJ(gfe_i1_k2_jl) + dist2MoveJ, gfe_i1_k2_ju);
	  (*this).SetFCenterJ( (*this).FCenterJ(gfe_i2_k1_jl) + dist2MoveJ, gfe_i2_k1_ju);
	  (*this).SetFCenterJ( (*this).FCenterJ(gfe_i2_k2_jl) + dist2MoveJ, gfe_i2_k2_ju);

	}

      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the k-direction --------------------------------
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,j1,i2,j2,ie,je;

      int gfe_i1_j1_il, gfe_i1_j1_jl;
      int gfe_i1_j2_il, gfe_i1_j2_jl;
      int gfe_i2_j1_il, gfe_i2_j1_jl;
      int gfe_i2_j2_il, gfe_i2_j2_jl;

      int gf_i1_je_il, gf_i1_je_jl, gf_i1_je_ju;
      int gf_i2_je_il, gf_i2_je_jl;
      int gf_ie_j1_il, gf_ie_j1_iu;
      int gf_ie_j2_il, gf_ie_j2_jl;

      string bc_I, bc_J, surfI, surfJ;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both i and j, so use GetLowerFace for both
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, j# = j position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG); 
	gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG); 
	gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG); 
	gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG); 
	gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG); 
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG); 
	gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG); 

	//ghost face, on second layer of i line of cells, on non-edge layer of j line of cells
	gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG); 
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG); 

	//ghost face, on non-edge layer of i line of cells, on first layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG); 
	gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG); 

	//ghost face, on non-edge layer of i line of cells, on second layer of j line of cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  
	gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG);  

	surfI = "il";
	surfJ = "jl";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of i and upper direction of j, so use GetLowerFace for I
	i2 = 0;
	i1 = 1;
	ie = (*this).NumGhosts();

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, j# = j position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_il = GetLowerFaceI(i1, j1, kk, imaxG, jmaxG); 
	gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_il = GetLowerFaceI(i1, j2, kk, imaxG, jmaxG); 
	gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_il = GetLowerFaceI(i2, j1, kk, imaxG, jmaxG); 
	gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_il = GetLowerFaceI(i2, j2, kk, imaxG, jmaxG); 
	gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_il = GetLowerFaceI(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of j line of cells
	gf_i2_je_il = GetLowerFaceI(i2, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j1_iu = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of j line of cells
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  
	gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);  

	surfI = "il";
	surfJ = "ju";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 2 ){ //at iu/jl edge - ghost cells are in the upper direction of i, and lower direction of j so use GetLowerFace for J
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	j2 = 0;
	j1 = 1;
	je = (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, j# = j position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG); 
	gfe_i1_j1_jl = GetLowerFaceJ(i1, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG); 
	gfe_i1_j2_jl = GetLowerFaceJ(i1, j2, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG); 
	gfe_i2_j1_jl = GetLowerFaceJ(i2, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG); 
	gfe_i2_j2_jl = GetLowerFaceJ(i2, j2, kk, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_ju = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of j line of cells
	gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of j line of cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG); 
	gf_ie_j2_jl = GetLowerFaceJ(ie, j2, kk, imaxG, jmaxG); 

	surfI = "iu";
	surfJ = "jl";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 3 ){ //at iu/ju edge - ghost cells are in the upper direction of both i and j, use GetUpperFace for both
	i2 = imaxG - 1;
	i1 = imaxG - 2;
	ie = imax - 1 + (*this).NumGhosts();

	j2 = jmaxG - 1;
	j1 = jmaxG - 2;
	je = jmax - 1 + (*this).NumGhosts();

	//ghost edge face indices
	//naming convention - g = ghost, f = face, e = edge ghost cell, i# = i position of location, j# = j position of location, jl = j lower face
	//ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	gfe_i1_j1_il = GetUpperFaceI(i1, j1, kk, imaxG, jmaxG); 
	gfe_i1_j1_jl = GetUpperFaceJ(i1, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	gfe_i1_j2_il = GetUpperFaceI(i1, j2, kk, imaxG, jmaxG); 
	gfe_i1_j2_jl = GetUpperFaceJ(i1, j2, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	gfe_i2_j1_il = GetUpperFaceI(i2, j1, kk, imaxG, jmaxG); 
	gfe_i2_j1_jl = GetUpperFaceJ(i2, j1, kk, imaxG, jmaxG); 

	//ghost face on edge, on second layer of i line of cells, on second layer of j line of cells
	gfe_i2_j2_il = GetUpperFaceI(i2, j2, kk, imaxG, jmaxG); 
	gfe_i2_j2_jl = GetUpperFaceJ(i2, j2, kk, imaxG, jmaxG); 

	//ghost face, on first layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_il = GetUpperFaceI(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i1_je_ju = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  

	//ghost face, on second layer of i line of cells, on non-edge layer of j line of cells
	gf_i2_je_il = GetUpperFaceI(i2, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j1_iu = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on second layer of j line of cells
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  
	gf_ie_j2_jl = GetUpperFaceJ(ie, j2, kk, imaxG, jmaxG);  

	surfI = "iu";
	surfJ = "ju";

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_I != "interblock" && bc_J != "interblock" ){ 

	//ghost cell and remaining face indices
	int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of j line of cells
	int gfe_i1_j1_kl = GetLowerFaceK(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells

	int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of j line of cells
	int gfe_i1_j2_kl = GetLowerFaceK(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells

	int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of j line of cells
	int gfe_i2_j1_kl = GetLowerFaceK(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells

	int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of j line of cells
	int gfe_i2_j2_kl = GetLowerFaceK(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);          //ghost cell, on first layer of i line of cells, on non-edge layer of j line of cells
	int gf_i1_je_kl = GetLowerFaceK(i1, je, kk, imaxG, jmaxG);  //ghost face, on first layer of i line of cells, on non-edge layer of j line of cells

	int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);          //ghost cell, on second layer of i line of cells, on non-edge layer of j line of cells
	int gf_i2_je_kl = GetLowerFaceK(i2, je, kk, imaxG, jmaxG);  //ghost face, on second layer of i line of cells, on non-edge layer of j line of cells

	int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);          //ghost cell, on non-edge layer of i line of cells, on first layer of j line of cells
	int gf_ie_j1_kl = GetLowerFaceK(ie, j1, kk, imaxG, jmaxG);  //ghost face, on non-edge layer of i line of cells, on first layer of j line of cells

	int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);          //ghost cell, on non-edge layer of i line of cells, on second layer of j line of cells
	int gf_ie_j2_kl = GetLowerFaceK(ie, j2, kk, imaxG, jmaxG);  //ghost face, on non-edge layer of i line of cells, on second layer of j line of cells

	//volume -----------------------------------------------------------------------------------------------------------------------------------------------
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_i1_je) + (*this).Vol(gc_ie_j1) ) ,gce_i1_j1);
	(*this).SetVol( (*this).Vol(gc_i2_je) ,gce_i2_j1);
	(*this).SetVol( (*this).Vol(gc_ie_j2) ,gce_i1_j2);
	(*this).SetVol( 0.5 * ( (*this).Vol(gc_i2_je) + (*this).Vol(gc_ie_j2) ) ,gce_i2_j2);

	//face areas -------------------------------------------------------------------------------------------------------------------------------------------
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

	//centroids ---------------------------------------------------------------------------------------------------------------------------------------------
	//edge centroid is moved distance of cell width normal to face dividing regular and edge ghost cells
	vector3d<double> dist2MoveJ = (*this).FCenterJ(gf_i1_je_jl) - (*this).FCenterJ(gf_i1_je_ju);
	vector3d<double> dist2MoveI = (*this).FCenterI(gf_ie_j1_il) - (*this).FCenterI(gf_ie_j1_iu);
	(*this).SetCenter( (*this).Center(gc_i1_je) + dist2MoveJ, gce_i1_j1);
	(*this).SetCenter( (*this).Center(gc_i2_je) + dist2MoveJ, gce_i2_j1);
	(*this).SetCenter( (*this).Center(gc_ie_j2) + dist2MoveI, gce_i1_j2);
	(*this).SetCenter( (*this).Center(gc_ie_j2) + 2.0 * dist2MoveI, gce_i2_j2);

	//face centers ------------------------------------------------------------------------------------------------------------------------------------------
	//edge face centers are moved distance of cell width normal to face dividing regular and edge ghost cells
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
	  vector3d<double> dist2MoveK = (*this).FCenterK(gfe_i1_j1_kl) - (*this).FCenterK(gfe_i1_j1_2kl); //k-width of adjacent cell

	  int gfe_i1_j1_ku = GetUpperFaceK(i1, j1, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on first layer of j line of cells
	  int gfe_i1_j2_ku = GetUpperFaceK(i1, j2, kk, imaxG, jmaxG); //ghost face on edge, on first layer of i line of cells, on second layer of j line of cells
	  int gfe_i2_j1_ku = GetUpperFaceK(i2, j1, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on first layer of j line of cells
	  int gfe_i2_j2_ku = GetUpperFaceK(i2, j2, kk, imaxG, jmaxG); //ghost face on edge, on second layer of i line of cells, on second layer of j line of cells

	  //face areas
	  (*this).SetFAreaK( (*this).FAreaK(gfe_i1_j1_kl), gfe_i1_j1_ku);
	  (*this).SetFAreaK( (*this).FAreaK(gfe_i1_j2_kl), gfe_i1_j2_ku);
	  (*this).SetFAreaK( (*this).FAreaK(gfe_i2_j1_kl), gfe_i2_j1_ku);
	  (*this).SetFAreaK( (*this).FAreaK(gfe_i2_j2_kl), gfe_i2_j2_ku);

	  //face centers
	  (*this).SetFCenterK( (*this).FCenterK(gfe_i1_j1_kl) + dist2MoveK, gfe_i1_j1_ku);
	  (*this).SetFCenterK( (*this).FCenterK(gfe_i1_j2_kl) + dist2MoveK, gfe_i1_j2_ku);
	  (*this).SetFCenterK( (*this).FCenterK(gfe_i2_j1_kl) + dist2MoveK, gfe_i2_j1_ku);
	  (*this).SetFCenterK( (*this).FCenterK(gfe_i2_j2_kl) + dist2MoveK, gfe_i2_j2_ku);

	}
      }

    }
  }

}

/* Member function to assign values for ghost cells for the inviscid flux calculation. This function assigns values for
   regular ghost cells and "edge" ghost cells. "Corner" cells are left with no value as they are not used.
   ____ ____ ____ ____ ____ ____ ____ ____
   | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
   |____|____|____|____|____|____|____|____|
   | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | G2 | G1 | X  | X  | X  | X  | G1 | G2 |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G1 | G1 | G1 | G1 | E  | E  |
          |____|____|____|____|____|____|____|____|
          | E  | E  | G2 | G2 | G2 | G2 | E  | E  |
          |____|____|____|____|____|____|____|____|

In the above diagram where X represents the physical cells, cells marked G (regular ghost cells) and E ("edge" ghost cells) are assigned geometric
values. G1 represents the first layer of ghost cells and G2 represents the second layer.
*/
void procBlock::AssignInviscidGhostCells(const input &inp, const idealGas &eos){
  // inp -- all input variables
  // eos -- equation of state

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();
  int kmaxG = (*this).NumK() + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical I faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(0, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "il");
      string bcNameU = bound.GetBCName(imax, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "iu");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
	bcNameU = "slipWall";
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at lower i-boundary
	int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
	int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

	//location of interior cells at lower i-boundary
	int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

	//location of lower i-boundary face
	int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  if (bcNameL == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 2), cellLowG2);
	  }
	}

      }

      //upper surface -------------------------------------------------------------------------
      if ( bcNameU != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at upper i-boundary
	int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
	int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);

	//location of interior cells at upper i-boundary
	int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
	int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

	//location of upper i-boundary face
	int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  if (bcNameU == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 2), cellUpG2);
	  }
	}

      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical J faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), 0, kk - (*this).NumGhosts(), "jl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jmax, kk - (*this).NumGhosts(), "ju");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
      	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
      	bcNameU = "slipWall";
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at lower j-boundary
	int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
	int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

	//location of interior cells at lower j-boundary
	int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
	int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

	//location of lower j-boundary face
	int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use once cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  if (bcNameL == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 2), cellLowG2);
	  }
	}

      }

      //upper surface ----------------------------------------------------------
      if ( bcNameU != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at upper j-boundary
	int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
	int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);

	//location of interior cells at upper j-boundary
	int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
	int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

	//location of upper j-boundary face
	int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use once cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  if (bcNameU == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 2), cellUpG2);
	  }
	}

      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical K faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //name of boundary conditions at lower and upper boundaries
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), 0, "kl");
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kmax, "ku");

      //inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
      if ( bcNameL == "viscousWall" ){
	bcNameL = "slipWall";
      }
      if ( bcNameU == "viscousWall" ){
	bcNameU = "slipWall";
      }

      //lower surface ----------------------------------------------------------
      if ( bcNameL != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at lower k-boundary
	int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
	int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

	//location of interior cells at lower k-boundary
	int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
	int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

	//location of lower k-boundary face
	int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use once cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  if (bcNameL == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 2), cellLowG2);
	  }
	}

      }

      //upper surface ----------------------------------------------------------
      if ( bcNameU != "interblock" ) { //only supply state values for non interblock BCs, for interblock, do nothing

	//location of ghost cells at upper k-boundary
	int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
	int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);

	//location of interior cells at upper k-boundary
	int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
	int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

	//location of upper k-boundary face
	int uFaceB = GetUpperFaceK(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use once cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  if (bcNameU == "slipWall" ){ //if slipWall, reflect second interior state over boundary face instead of extrapolation
	    (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG2);
	  }
	  else{
	    (*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 2), cellUpG2);
	  }
	}

      }

    }
  }

  //assign values to edge ghost cells
  //(*this).AssignInviscidGhostCellsEdge(inp, eos);

}

/* Member function to assign values to ghost cells located on the 12 block edges for the inviscid flux calculation. Assumes AssignInviscidGhostCells
 has already been run. 

           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J 

In the above diagram the cells marked X represent physical cells. Cells marked G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells (1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to 
the values of the ghost cells they tough at level "e". The values at edge cell 2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it touches.
The exception to this rule occurs when either of the boundaries that meet at the corner are wall boundaries (slipWall, viscousWall) and the other is not. When this 
occurs the wall boundaries are "extended" into the ghost cells. This implementation is described in Blazek.
*/
void procBlock::AssignInviscidGhostCellsEdge(const input &inp, const idealGas &eos){
  // inp -- all input variables
  // eos -- equation of state

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction --------------------------------
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int j1,k1,j2,k2,je,ke,je2,ke2;

      int gf_j1_ke_kl;
      int gf_j2_ke_kl;
      int gf_je_k1_jl;
      int gf_je_k2_jl;

      string bc_J,bc_K;
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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(),surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_J = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_K = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_J != "interblock" && bc_K != "interblock" ){ 

	int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
	int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
	int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
	int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells

	int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);       //ghost cell, on first layer of j line of cells, on non-edge layer of k line of cells
	int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);       //ghost cell, on second layer of j line of cells, on non-edge layer of k line of cells
	int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);       //ghost cell, on non-edge layer of j line of cells, on first layer of k line of cells
	int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);       //ghost cell, on non-edge layer of j line of cells, on second layer of k line of cells

	int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG, jmaxG);     //ghost cell, on first layer of j line of cells, on second non-edge layer of k line of cells 
	int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG, jmaxG);     //ghost cell, on second layer of j line of cells, on second non-edge layer of k line of cells  
	int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG, jmaxG);     //ghost cell, on second non-edge layer of j line of cells, on first layer of k line of cells   
	int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG, jmaxG);     //ghost cell, on second non-edge layer of j line of cells, on second layer of k line of cells  

	//inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
	if ( bc_J == "viscousWall" ){
	  bc_J = "slipWall";
	}
	if ( bc_K == "viscousWall" ){
	  bc_K = "slipWall";
	}

	if ( bc_J == "slipWall" && !(bc_K == "slipWall") ){  //j surface is a wall, but k surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_je_k1).GetGhostState(bc_J, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j1_k1);
	  (*this).SetState( (*this).State(gc_je_k2).GetGhostState(bc_J, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j1_k2);
	  (*this).SetState( (*this).State(gc_je2_k1).GetGhostState(bc_J, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j2_k1);
	  (*this).SetState( (*this).State(gc_je2_k2).GetGhostState(bc_J, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j2_k2);
	}
	else if ( !(bc_J == "slipWall") && bc_K == "slipWall" ){  //k surface is a wall, but j surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_j1_ke).GetGhostState(bc_K, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k1);
	  (*this).SetState( (*this).State(gc_j2_ke).GetGhostState(bc_K, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k1);
	  (*this).SetState( (*this).State(gc_j1_ke2).GetGhostState(bc_K, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k2);
	  (*this).SetState( (*this).State(gc_j2_ke2).GetGhostState(bc_K, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k2);
	}
	else{ // both surfaces are walls or neither are walls - proceed as normal
	  (*this).SetState( 0.5 * ( (*this).State(gc_j1_ke) + (*this).State(gc_je_k1) ) ,gce_j1_k1);
	  (*this).SetState( (*this).State(gc_j2_ke) ,gce_j2_k1);
	  (*this).SetState( (*this).State(gc_je_k2) ,gce_j1_k2);
	  (*this).SetState( 0.5 * ( (*this).State(gc_j2_ke) + (*this).State(gc_je_k2) ) ,gce_j2_k2);
	}

      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction --------------------------------
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,k1,i2,k2,ie,ke,ie2,ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_I,bc_K;
      string surfI,surfK;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both i and k, so use GetLowerFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of i and upper direction of k, so use GetLowerFace for I
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at iu/kl edge - ghost cells are in the lower direction of k, and upper direction of i so use GetLowerFace for k
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditioins at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_K = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_I != "interblock" && bc_K != "interblock" ){ 

	//location of ghost cells
	int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
	int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
	int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
	int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

	int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);       //ghost cell, on first layer of i line of cells, on non-edge layer of k line of cells          
	int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);       //ghost cell, on second layer of i line of cells, on non-edge layer of k line of cells         
	int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on first layer of k line of cells          
	int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on second layer of k line of cells         
							                                                                                                      
	int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG, jmaxG);     //ghost cell, on first layer of i line of cells, on second non-edge layer of k line of cells   
	int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG, jmaxG);     //ghost cell, on second layer of i line of cells, on second non-edge layer of k line of cells  
	int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on first layer of k line of cells   
	int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on second layer of k line of cells  

	//inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
	if ( bc_I == "viscousWall" ){
	  bc_I = "slipWall";
	}
	if ( bc_K == "viscousWall" ){
	  bc_K = "slipWall";
	}

	if ( bc_I == "slipWall" && !(bc_K == "slipWall") ){  //i surface is a wall, but k surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_ie_k1).GetGhostState(bc_I, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i1_k1);
	  (*this).SetState( (*this).State(gc_ie_k2).GetGhostState(bc_I, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i1_k2);
	  (*this).SetState( (*this).State(gc_ie2_k1).GetGhostState(bc_I, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i2_k1);
	  (*this).SetState( (*this).State(gc_ie2_k2).GetGhostState(bc_I, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i2_k2);
	}
	else if ( !(bc_I == "slipWall") && bc_K == "slipWall" ){  //k surface is a wall, but i surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_i1_ke).GetGhostState(bc_K, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k1);
	  (*this).SetState( (*this).State(gc_i2_ke).GetGhostState(bc_K, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k1);
	  (*this).SetState( (*this).State(gc_i1_ke2).GetGhostState(bc_K, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k2);
	  (*this).SetState( (*this).State(gc_i2_ke2).GetGhostState(bc_K, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k2);

	}
	else{ // both surfaces are walls or neither are walls - proceed as normal
	  (*this).SetState( 0.5 * ( (*this).State(gc_i1_ke) + (*this).State(gc_ie_k1) ) ,gce_i1_k1);
	  (*this).SetState( (*this).State(gc_i2_ke) ,gce_i2_k1);
	  (*this).SetState( (*this).State(gc_ie_k2) ,gce_i1_k2);
	  (*this).SetState( 0.5 * ( (*this).State(gc_i2_ke) + (*this).State(gc_ie_k2) ) ,gce_i2_k2);
	}

      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the k-direction --------------------------------
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,j1,i2,j2,ie,je,ie2,je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_I,bc_J;
      string surfI,surfJ;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both i and j, so use GetLowerFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of i and upper direction of j, so use GetLowerFace for I
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 2 ){ //at iu/jl edge - ghost cells are in the lower direction of j, and upper direction of i so use GetLowerFace for j
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 3 ){ //at iu/ju edge - ghost cells are in the upper direction of both i and j, use GetUpperFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_I = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_J = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }

      //if both corner boundaries aren't interblocks, get edge cells; if bouth boundaries are interblocks, edge cells come from ghost cell exchange
      if ( bc_I != "interblock" && bc_J != "interblock" ){ 

	//location of ghost cells
	int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
	int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
	int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
	int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

	int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);       //ghost cell, on first layer of i line of cells, on non-edge layer of j line of cells          
	int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);       //ghost cell, on second layer of i line of cells, on non-edge layer of j line of cells         
	int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on first layer of j line of cells          
	int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on second layer of j line of cells         
							                                                                                                      
	int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG, jmaxG);     //ghost cell, on first layer of i line of cells, on second non-edge layer of j line of cells   
	int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG, jmaxG);     //ghost cell, on second layer of i line of cells, on second non-edge layer of j line of cells  
	int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on first layer of j line of cells   
	int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on second layer of j line of cells  

	//inviscid fluxes require different bc than viscous fluxes - treat all walls as the same
	if ( bc_I == "viscousWall" ){
	  bc_I = "slipWall";
	}
	if ( bc_J == "viscousWall" ){
	  bc_J = "slipWall";
	}

	if ( bc_I == "slipWall" && !(bc_J == "slipWall") ){  //i surface is a wall, but k surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_ie_j1).GetGhostState(bc_I, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i1_j1);
	  (*this).SetState( (*this).State(gc_ie_j2).GetGhostState(bc_I, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i1_j2);
	  (*this).SetState( (*this).State(gc_ie2_j1).GetGhostState(bc_I, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i2_j1);
	  (*this).SetState( (*this).State(gc_ie2_j2).GetGhostState(bc_I, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i2_j2);
	}
	else if ( !(bc_I == "slipWall") && bc_J == "slipWall" ){  //k surface is a wall, but i surface is not - extend wall bc
	  (*this).SetState( (*this).State(gc_i1_je).GetGhostState(bc_J, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j1);
	  (*this).SetState( (*this).State(gc_i2_je).GetGhostState(bc_J, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j1);
	  (*this).SetState( (*this).State(gc_i1_je2).GetGhostState(bc_J, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j2);
	  (*this).SetState( (*this).State(gc_i2_je2).GetGhostState(bc_J, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j2);
	}
	else{ // both surfaces are walls or neither are walls - proceed as normal
	  (*this).SetState( 0.5 * ( (*this).State(gc_i1_je) + (*this).State(gc_ie_j1) ) ,gce_i1_j1);
	  (*this).SetState( (*this).State(gc_i2_je) ,gce_i2_j1);
	  (*this).SetState( (*this).State(gc_ie_j2) ,gce_i1_j2);
	  (*this).SetState( 0.5 * ( (*this).State(gc_i2_je) + (*this).State(gc_ie_j2) ) ,gce_i2_j2);
	}

      }
    }

  }
}

/* Member function to assign ghost cells for the viscous flow calculation. This function assumes AssignInviscidGhostCells has been run first as
 it only overwrites the ghost cells associated with the viscousWall boundary condition. It overwrites both regular and edge ghost cells.
*/
void procBlock::AssignViscousGhostCells(const input &inp, const idealGas &eos){
  // inp -- all input variables
  // eos -- equation of state

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = (*this).NumI() + 2 * (*this).NumGhosts();
  int jmaxG = (*this).NumJ() + 2 * (*this).NumGhosts();
  int kmaxG = (*this).NumK() + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical I faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

      //location of ghost cells at lower i-boundary
      int cellLowG1 = GetLoc1D(1, jj, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(0, jj, kk, imaxG, jmaxG);

      //location of interior cells at lower i-boundary
      int cellLowIn1 = GetLoc1D((*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D((*this).NumGhosts() + 1, jj, kk, imaxG, jmaxG);

      //location of lower i-boundary face
      int lFaceB = GetLowerFaceI((*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      //location of ghost cells at upper i-boundary
      int cellUpG1 = GetLoc1D(imaxG-2, jj, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(imaxG-1, jj, kk, imaxG, jmaxG);

      //location of interior cells at upper i-boundary
      int cellUpIn1 = GetLoc1D(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(imaxG - 2 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG);

      //location of upper i-boundary face
      int uFaceB = GetUpperFaceI(imaxG - 1 - (*this).NumGhosts(), jj, kk, imaxG, jmaxG); 

      //boundary condition at lower boundary
      string bcNameL = bound.GetBCName(0, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "il");

      //if viscous, overwrite regular ghost cell
      if (bcNameL == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaI(lFaceB), "il", inp, eos, 1), cellLowG2);
	}
      }

      //boundary condition at upper boundary
      string bcNameU = bound.GetBCName(imax, jj - (*this).NumGhosts(), kk - (*this).NumGhosts(), "iu");

      //if viscous, overwrite regular ghost cell
      if (bcNameU == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
	if (imax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaI(uFaceB), "iu", inp, eos, 1), cellUpG2);
	}
      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical J faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //location of ghost cells at lower j-boundary
      int cellLowG1 = GetLoc1D(ii, 1, kk, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, 0, kk, imaxG, jmaxG);

      //location of interior cells at lower j-boundary
      int cellLowIn1 = GetLoc1D(ii, (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, (*this).NumGhosts() + 1, kk, imaxG, jmaxG);

      //location of lower j-boundary face
      int lFaceB = GetLowerFaceJ(ii, (*this).NumGhosts(), kk, imaxG, jmaxG); 

      //location of ghost cells at upper j-boundary
      int cellUpG1 = GetLoc1D(ii, jmaxG-2, kk, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jmaxG-1, kk, imaxG, jmaxG);

      //location of interior cells at upper j-boundary
      int cellUpIn1 = GetLoc1D(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jmaxG - 2 - (*this).NumGhosts(), kk, imaxG, jmaxG);

      //location of upper j-boundary face
      int uFaceB = GetUpperFaceJ(ii, jmaxG - 1 - (*this).NumGhosts(), kk, imaxG, jmaxG); 

      //boundary condition at lower boundary
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), 0, kk - (*this).NumGhosts(), "jl");

      //if viscous, overwrite regular ghost cell
      if (bcNameL == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaJ(lFaceB), "jl", inp, eos, 1), cellLowG2);
	}
      }

      //boundary condition at upper boundary
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jmax, kk - (*this).NumGhosts(), "ju");

      //if viscous, overwrite regular ghost cell
      if (bcNameU == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
	if (jmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellUpG1), cellUpG2);
	}
	else{
	  (*this).SetState( (*this).State(cellUpIn2).GetGhostState(bcNameU, (*this).FAreaJ(uFaceB), "ju", inp, eos, 1), cellUpG2);
	}
      }

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over physical K faces and assign values for regular ghost cells -----------------------------------------------------------------------
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){
    for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

      //location of ghost cells at lower k-boundary
      int cellLowG1 = GetLoc1D(ii, jj, 1, imaxG, jmaxG);
      int cellLowG2 = GetLoc1D(ii, jj, 0, imaxG, jmaxG);

      //location of interior cells at lower k-boundary
      int cellLowIn1 = GetLoc1D(ii, jj, (*this).NumGhosts(), imaxG, jmaxG);
      int cellLowIn2 = GetLoc1D(ii, jj, (*this).NumGhosts() + 1, imaxG, jmaxG);

      //location of lower k-boundary face
      int lFaceB = GetLowerFaceK(ii, jj, (*this).NumGhosts(), imaxG, jmaxG); 

      //location of interior cells at lower k-boundary
      int cellUpG1 = GetLoc1D(ii, jj, kmaxG-2, imaxG, jmaxG);
      int cellUpG2 = GetLoc1D(ii, jj, kmaxG-1, imaxG, jmaxG);

      //location of interior cells at upper k-boundary
      int cellUpIn1 = GetLoc1D(ii, jj, kmaxG - 1 - (*this).NumGhosts(), imaxG, jmaxG);
      int cellUpIn2 = GetLoc1D(ii, jj, kmaxG - 2 - (*this).NumGhosts(), imaxG, jmaxG);

      //location of upper k-boundary face
      int uFaceB = GetUpperFaceK(ii, jj, kmax - 1 - (*this).NumGhosts(), imaxG, jmaxG); 

      //name of boundary condition at lower boundary
      string bcNameL = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), 0, "kl");

      //if viscous, overwrite regular ghost cells
      if (bcNameL == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellLowIn1).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG1);

	//second layer of ghost cells
	if (kmax < 2){ //one cell thick - use one cell for both ghost cells
	  (*this).SetState( (*this).State(cellLowG1), cellLowG2);
	}
	else{
	  (*this).SetState( (*this).State(cellLowIn2).GetGhostState(bcNameL, (*this).FAreaK(lFaceB), "kl", inp, eos, 1), cellLowG2);
	}
      }

      //boundary condition at upper boundary
      string bcNameU = bound.GetBCName(ii - (*this).NumGhosts(), jj - (*this).NumGhosts(), kmax, "ku");

      //if viscous, overwrite regular ghost cells
      if (bcNameU == "viscousWall"){
	//first layer of ghost cells
	(*this).SetState( (*this).State(cellUpIn1).GetGhostState(bcNameU, (*this).FAreaK(uFaceB), "ku", inp, eos, 1), cellUpG1);

	//second layer of ghost cells
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

}

/* Member function to assign values to ghost cells located on the 12 block edges for the viscous flux calculation. Assumes AssignViscousGhostCells
 has already been run. Only overwrites edge ghost cells if one of the boundaries at the corner is a viscousWall boundary condition.

           ____ ____ ____ ____ ____ ____ ____ ____
          | G2 | G1 | X  | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         K| G2 | G1 | X  | X  | X  | X  | X  | X  |
         ^|____|____|____|____|____|____|____|____|
         || G2 | G1 | X  | X  | X  | X  | X  | X  |
         ||____|____|____|____|____|____|____|____|
         e| G2 | G1 | X* | X  | X  | X  | X  | X  |
          |____|____|____|____|____|____|____|____|
         1| E  | E  | G1 | G1 | G1 | G1 | G1 | G1 |
          |____|____|____|____|____|____|____|____|
         2| E  | E  | G2 | G2 | G2 | G2 | G2 | G2 |
          |____|____|____|____|____|____|____|____|
            2    1     e ----> J 

In the above diagram the cells marked X represent physical cells. Cells marked G1 and G2 represent the first and second layer of ghost cells respectively. The
cells marked E are the edge ghost cells that need to be assigned values. At each corner location (X*) there are 4 edge ghost cells that need to be filled. The
axes on the side of the diagram indicate the coordinates of the edge ghost cells (1, 2) as well as the coordinates of the adjacent regualar ghost cells (e).

The values at edge cell 1,1 are the average of the values at the two ghost cells it touches at level "e". The values at edge cells 1,2 and 2,1 are identical to 
the values of the ghost cells they tough at level "e". The values at edge cell 2,2 are the average of the values at the two (1,2 & 2,1) edge ghost cells it touches.
The exception to this rule occurs when either of the boundaries that meet at the corner are viscousWall boundaries and the other is not. When this 
occurs the viscousWall boundaries are "extended" into the ghost cells. This implementation is described in Blazek.
*/
void procBlock::AssignViscousGhostCellsEdge(const input &inp, const idealGas &eos){
  // inp -- all input variables
  // eos -- equation of state

  //get boundary conditions for block
  const boundaryConditions bound = inp.BC( (*this).ParentBlock() );

  //max dimensions for vectors without ghost cells
  int imax = (*this).NumI();
  int jmax = (*this).NumJ();
  int kmax = (*this).NumK();

  //max dimensions for vectors with ghost cells
  int imaxG = imax + 2 * (*this).NumGhosts();
  int jmaxG = jmax + 2 * (*this).NumGhosts();
  int kmaxG = kmax + 2 * (*this).NumGhosts();

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper j sides of block - this will include 4 edges that run in the i-direction --------------------------------
  //edges at the jl/kl, jl/ku, ju/kl, ju/ku sides will be accounted for in this loop
  for ( int ii = (*this).NumGhosts(); ii < imax + (*this).NumGhosts(); ii++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetLowerFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetLowerFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetLowerFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetLowerFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
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

	//ghost face, on first/second layer of j line of cells, on non-edge layer of k line of cells
	gf_j1_ke_kl = GetUpperFaceK(ii, j1, ke, imaxG, jmaxG);  
	gf_j2_ke_kl = GetUpperFaceK(ii, j2, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of j line of cells, on first/second layer of k line of cells
	gf_je_k1_jl = GetUpperFaceJ(ii, je, k1, imaxG, jmaxG);  
	gf_je_k2_jl = GetUpperFaceJ(ii, je, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_jl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, ke - (*this).NumGhosts(), surfJ);
	bc_kl = bound.GetBCName(ii - (*this).NumGhosts(), je - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //location of ghost cells
      int gce_j1_k1 = GetLoc1D(ii, j1, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of j line of cells, on first layer of k line of cells
      int gce_j1_k2 = GetLoc1D(ii, j1, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of j line of cells, on second layer of k line of cells
      int gce_j2_k1 = GetLoc1D(ii, j2, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on first layer of k line of cells
      int gce_j2_k2 = GetLoc1D(ii, j2, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of j line of cells, on second layer of k line of cells

      int gc_j1_ke = GetLoc1D(ii, j1, ke, imaxG, jmaxG);      //ghost cell, on first layer of j line of cells, on non-edge layer of k line of cells          
      int gc_j2_ke = GetLoc1D(ii, j2, ke, imaxG, jmaxG);      //ghost cell, on second layer of j line of cells, on non-edge layer of k line of cells         
      int gc_je_k1 = GetLoc1D(ii, je, k1, imaxG, jmaxG);      //ghost cell, on non-edge layer of j line of cells, on first layer of k line of cells          
      int gc_je_k2 = GetLoc1D(ii, je, k2, imaxG, jmaxG);      //ghost cell, on non-edge layer of j line of cells, on second layer of k line of cells         
							                                                                                                     
      int gc_j1_ke2 = GetLoc1D(ii, j1, ke2, imaxG, jmaxG);    //ghost cell, on first layer of j line of cells, on second non-edge layer of k line of cells   
      int gc_j2_ke2 = GetLoc1D(ii, j2, ke2, imaxG, jmaxG);    //ghost cell, on second layer of j line of cells, on second non-edge layer of k line of cells  
      int gc_je2_k1 = GetLoc1D(ii, je2, k1, imaxG, jmaxG);    //ghost cell, on second non-edge layer of j line of cells, on first layer of k line of cells   
      int gc_je2_k2 = GetLoc1D(ii, je2, k2, imaxG, jmaxG);    //ghost cell, on second non-edge layer of j line of cells, on second layer of k line of cells  

      if ( bc_jl == "viscousWall" && !(bc_kl == "viscousWall") ){  //j surface is a viscous wall, but k surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_je_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_je_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_je2_k1).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k1_jl), surfJ, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je2_k2).GetGhostState(bc_jl, (*this).FAreaJ(gf_je_k2_jl), surfJ, inp, eos, 1) ,gce_j2_k2);
      }
      else if ( !(bc_jl == "viscousWall") && bc_kl == "viscousWall" ){  //k surface is a viscous wall, but j surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_j1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_j1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j1_ke_kl), surfK, inp, eos, 1) ,gce_j1_k2);
	(*this).SetState( (*this).State(gc_j2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_j2_ke_kl), surfK, inp, eos, 1) ,gce_j2_k2);
      }
      else if ( bc_jl == "viscousWall" && bc_kl == "viscousWall" ){ // both surfaces are viscous walls - proceed as normal
	(*this).SetState( 0.5 * ( (*this).State(gc_j1_ke) + (*this).State(gc_je_k1) ) ,gce_j1_k1);
	(*this).SetState( (*this).State(gc_j2_ke) ,gce_j2_k1);
	(*this).SetState( (*this).State(gc_je_k2) ,gce_j1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_j2_ke) + (*this).State(gc_je_k2) ) ,gce_j2_k2);
      }
      //if no boundary is a viscous wall, do nothing

    }
  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the j-direction --------------------------------
  //edges at the il/kl, il/ku, iu/kl, iu/ku sides will be accounted for in this loop
  for ( int jj = (*this).NumGhosts(); jj < jmax + (*this).NumGhosts(); jj++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,k1,i2,k2,ie,ke,ie2,ke2;

      int gf_i1_ke_kl;
      int gf_i2_ke_kl;
      int gf_ie_k1_il;
      int gf_ie_k2_il;

      string bc_il,bc_kl;
      string surfI,surfK;

      if ( cc == 0 ){ //at il/kl edge - ghost cells are in the lower direction of both i and k, so use GetLowerFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 1 ){ //at il/ku edge - ghost cells are in the lower direction of i and upper direction of k, so use GetLowerFace for I
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetLowerFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetLowerFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }
      else if ( cc == 2 ){ //at iu/kl edge - ghost cells are in the lower direction of k, and upper direction of i so use GetLowerFace for k
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetLowerFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetLowerFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfK);

      }
      else if ( cc == 3 ){ //at iu/ku edge - ghost cells are in the upper direction of both i and k, use GetUpperFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of k line of cells
	gf_i1_ke_kl = GetUpperFaceK(i1, jj, ke, imaxG, jmaxG);  
	gf_i2_ke_kl = GetUpperFaceK(i2, jj, ke, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of k line of cells
	gf_ie_k1_il = GetUpperFaceI(ie, jj, k1, imaxG, jmaxG);  
	gf_ie_k2_il = GetUpperFaceI(ie, jj, k2, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, jj - (*this).NumGhosts(), ke - (*this).NumGhosts(), surfI);
	bc_kl = bound.GetBCName(ie - (*this).NumGhosts(), jj - (*this).NumGhosts(), ke - (*this).NumGhosts() + 1, surfK);

      }

      //location of ghost cells
      int gce_i1_k1 = GetLoc1D(i1, jj, k1, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_k2 = GetLoc1D(i1, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_k1 = GetLoc1D(i2, jj, k1, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_k2 = GetLoc1D(i2, jj, k2, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_ke = GetLoc1D(i1, jj, ke, imaxG, jmaxG);       //ghost cell, on first layer of i line of cells, on non-edge layer of k line of cells          
      int gc_i2_ke = GetLoc1D(i2, jj, ke, imaxG, jmaxG);       //ghost cell, on second layer of i line of cells, on non-edge layer of k line of cells         
      int gc_ie_k1 = GetLoc1D(ie, jj, k1, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on first layer of k line of cells          
      int gc_ie_k2 = GetLoc1D(ie, jj, k2, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on second layer of k line of cells         
							                                                                                                      
      int gc_i1_ke2 = GetLoc1D(i1, jj, ke2, imaxG, jmaxG);     //ghost cell, on first layer of i line of cells, on second non-edge layer of k line of cells   
      int gc_i2_ke2 = GetLoc1D(i2, jj, ke2, imaxG, jmaxG);     //ghost cell, on second layer of i line of cells, on second non-edge layer of k line of cells  
      int gc_ie2_k1 = GetLoc1D(ie2, jj, k1, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on first layer of k line of cells   
      int gc_ie2_k2 = GetLoc1D(ie2, jj, k2, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on second layer of k line of cells  

      if ( bc_il == "viscousWall" && !(bc_kl == "viscousWall") ){  //i surface is a viscous wall, but k surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_ie_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_ie_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_ie2_k1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k1_il), surfI, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie2_k2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_k2_il), surfI, inp, eos, 1) ,gce_i2_k2);
      }
      else if ( !(bc_il == "viscousWall") && bc_kl == "viscousWall" ){  //k surface is a viscous wall, but i surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_i1_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_i1_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i1_ke_kl), surfK, inp, eos, 1) ,gce_i1_k2);
	(*this).SetState( (*this).State(gc_i2_ke2).GetGhostState(bc_kl, (*this).FAreaK(gf_i2_ke_kl), surfK, inp, eos, 1) ,gce_i2_k2);
      }
      else if ( bc_il == "viscousWall" && bc_kl == "viscousWall" ){ // both surfaces are viscous walls - proceed as normal
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_ke) + (*this).State(gc_ie_k1) ) ,gce_i1_k1);
	(*this).SetState( (*this).State(gc_i2_ke) ,gce_i2_k1);
	(*this).SetState( (*this).State(gc_ie_k2) ,gce_i1_k2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_ke) + (*this).State(gc_ie_k2) ) ,gce_i2_k2);
      }
      //if neither surface is a wall then do nothing

    }

  }

  //--------------------------------------------------------------------------------------------------------------------------------------------
  //loop over edges at lower and upper i sides of block - this will include 4 edges that run in the k-direction --------------------------------
  //edges at the il/jl, il/ju, iu/jl, iu/ju sides will be accounted for in this loop
  for ( int kk = (*this).NumGhosts(); kk < kmax + (*this).NumGhosts(); kk++ ){

    for ( int cc = 0; cc < 4; cc++ ){ //loop over 4 edges

      int i1,j1,i2,j2,ie,je,ie2,je2;

      int gf_i1_je_jl;
      int gf_i2_je_jl;
      int gf_ie_j1_il;
      int gf_ie_j2_il;

      string bc_il,bc_jl;
      string surfI,surfJ;

      if ( cc == 0 ){ //at il/jl edge - ghost cells are in the lower direction of both i and j, so use GetLowerFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 1 ){ //at il/ju edge - ghost cells are in the lower direction of i and upper direction of j, so use GetLowerFace for I
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetLowerFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetLowerFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je  - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 2 ){ //at iu/jl edge - ghost cells are in the lower direction of j, and upper direction of i so use GetLowerFace for J
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetLowerFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetLowerFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfJ);

      }
      else if ( cc == 3 ){ //at iu/ju edge - ghost cells are in the upper direction of both i and j, use GetUpperFace for both
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

	//ghost face, on first/second layer of i line of cells, on non-edge layer of j line of cells
	gf_i1_je_jl = GetUpperFaceJ(i1, je, kk, imaxG, jmaxG);  
	gf_i2_je_jl = GetUpperFaceJ(i2, je, kk, imaxG, jmaxG);  

	//ghost face, on non-edge layer of i line of cells, on first/second layer of j line of cells
	gf_ie_j1_il = GetUpperFaceI(ie, j1, kk, imaxG, jmaxG);  
	gf_ie_j2_il = GetUpperFaceI(ie, j2, kk, imaxG, jmaxG);  

	//boundary conditions at corner
	bc_il = bound.GetBCName(ie - (*this).NumGhosts() + 1, je - (*this).NumGhosts(), kk - (*this).NumGhosts(), surfI);
	bc_jl = bound.GetBCName(ie - (*this).NumGhosts(), je - (*this).NumGhosts() + 1, kk - (*this).NumGhosts(), surfJ);

      }

      //location of ghost cells
      int gce_i1_j1 = GetLoc1D(i1, j1, kk, imaxG, jmaxG);      //ghost cell on edge, on first layer of i line of cells, on first layer of k line of cells
      int gce_i1_j2 = GetLoc1D(i1, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on first layer of i line of cells, on second layer of k line of cells
      int gce_i2_j1 = GetLoc1D(i2, j1, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on first layer of k line of cells
      int gce_i2_j2 = GetLoc1D(i2, j2, kk, imaxG, jmaxG);         //ghost cell on edge, on second layer of i line of cells, on second layer of k line of cells

      int gc_i1_je = GetLoc1D(i1, je, kk, imaxG, jmaxG);       //ghost cell, on first layer of i line of cells, on non-edge layer of j line of cells          
      int gc_i2_je = GetLoc1D(i2, je, kk, imaxG, jmaxG);       //ghost cell, on second layer of i line of cells, on non-edge layer of j line of cells         
      int gc_ie_j1 = GetLoc1D(ie, j1, kk, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on first layer of j line of cells          
      int gc_ie_j2 = GetLoc1D(ie, j2, kk, imaxG, jmaxG);       //ghost cell, on non-edge layer of i line of cells, on second layer of j line of cells         
							                                                                                                      
      int gc_i1_je2 = GetLoc1D(i1, je2, kk, imaxG, jmaxG);     //ghost cell, on first layer of i line of cells, on second non-edge layer of j line of cells   
      int gc_i2_je2 = GetLoc1D(i2, je2, kk, imaxG, jmaxG);     //ghost cell, on second layer of i line of cells, on second non-edge layer of j line of cells  
      int gc_ie2_j1 = GetLoc1D(ie2, j1, kk, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on first layer of j line of cells   
      int gc_ie2_j2 = GetLoc1D(ie2, j2, kk, imaxG, jmaxG);     //ghost cell, on second non-edge layer of i line of cells, on second layer of j line of cells  

      if ( bc_il == "viscousWall" && !(bc_jl == "viscousWall") ){  //i surface is a viscous wall, but k surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_ie_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_ie_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_ie2_j1).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j1_il), surfI, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie2_j2).GetGhostState(bc_il, (*this).FAreaI(gf_ie_j2_il), surfI, inp, eos, 1) ,gce_i2_j2);
      }
      else if ( !(bc_il == "viscousWall") && bc_jl == "viscousWall" ){  //j surface is a viscous wall, but i surface is not - extend wall bc
	(*this).SetState( (*this).State(gc_i1_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_i1_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i1_je_jl), surfJ, inp, eos, 1) ,gce_i1_j2);
	(*this).SetState( (*this).State(gc_i2_je2).GetGhostState(bc_jl, (*this).FAreaJ(gf_i2_je_jl), surfJ, inp, eos, 1) ,gce_i2_j2);
      }
      else if ( bc_il == "viscousWall" && bc_jl == "viscousWall"){ // both surfaces are viscous walls - proceed as normal
	(*this).SetState( 0.5 * ( (*this).State(gc_i1_je) + (*this).State(gc_ie_j1) ) ,gce_i1_j1);
	(*this).SetState( (*this).State(gc_i2_je) ,gce_i2_j1);
	(*this).SetState( (*this).State(gc_ie_j2) ,gce_i1_j2);
	(*this).SetState( 0.5 * ( (*this).State(gc_i2_je) + (*this).State(gc_ie_j2) ) ,gce_i2_j2);
      }
      //if neither surface is a viscous wall then do nothing

    }

  }

}

/* Member function to determine where in padded plot3dBlock an index is located. It takes in an i, j, k cell location and returns a boolean indicating
if the given i, j, k location corresponds to a corner location. Corner locations are not used at all. This function is NOT USED in the code but is useful
for debugging purposes.
 */
bool procBlock::AtCorner(const int& ii, const int& jj, const int& kk)const{
  // ii -- i index of location to test
  // jj -- j index of location to test
  // kk -- k index of location to test

  bool atCorner;

  //if all (i, j, & k) are outside of the limits of physical cells, location is a corner location
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

/* Member function to determine where in padded plot3dBlock an index is located. It takes in an i, j, k cell location and returns a boolean indicating
if the given i, j, k location corresponds to a edge location. Edge locations are used in the gradient calculations. This function is NOT USED in the 
code but is useful for debugging purposes.
 */
bool procBlock::AtEdge(const int& ii, const int& jj, const int& kk)const{
  // ii -- i index of location to test
  // jj -- j index of location to test
  // kk -- k index of location to test

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


/* Function to swap ghost cell geometry between two blocks at an interblock boundary. Slices are removed from the physical cells (extending into ghost cells at the edges)
of one block and inserted into the ghost cells of its partner block. The reverse is also true. The slices are taken in the coordinate system orientation of their parent
block. 

         Interior Cells       Ghost Cells                Ghost Cells       Interior Cells
         ________ ________|_________ _________        _________ _________|_________ _________
     Ui-3/2    Ui-1/2     |      Uj+1/2    Uj-1/2  Ui-3/2    Ui-1/2      |      Uj+1/2     Uj+3/2
        |        |        |         |         |       |        |         |         |         |
        | Ui-1   |  Ui    |   Uj+1  |  Uj     |       |  Ui-1  |   Ui    |  Uj     |  Uj+1   |
        |        |        |         |         |       |        |         |         |         |
        |________|________|_________|_________|       |________|_________|_________|_________|
                          |                                              |

The above diagram shows the resulting values after the ghost cell swap. The logic ensures that the ghost cells at the interblock boundary exactly match their partner block 
as if there were no separation in the grid.

Only 3 faces at each ghost cell need to be swapped (i.e. the lower face for Ui is the upper face for Ui-1). At the end of a line (i-line, j-line or k-line), both the upper
and lower faces need to be swapped.
*/
void SwapSlice( const interblock &inter, procBlock &blk1, procBlock &blk2, const bool& geom ){
  // inter -- interblock boundary information
  // blk1 -- first block involved in interblock boundary
  // blk2 -- second block involved in interblock boundary
  // geom -- boolean to determine whether to swap geometry or states

  //Get indices for slice coming from first block to swap
  int is1, ie1, js1, je1, ks1, ke1;

  //if at upper boundary no need to adjust for ghost cells as constant surface is already at the interior cells when acounting for ghost cells
  //if at the lower boundary adjust the constant surface by the number of ghost cells to get to the first interior cell 
  int upLowFac = ( inter.BoundaryFirst() % 2 == 0 ) ? 0 : blk1.NumGhosts();

  if ( inter.BoundaryFirst() == 1 || inter.BoundaryFirst() == 2 ){ //direction 3 is i
    //extend min/maxes to cover ghost cells
    is1 = inter.ConstSurfaceFirst() + upLowFac;
    ie1 = is1 + blk1.NumGhosts() - 1;
    //direction 1 is j
    js1 = inter.Dir1StartFirst();
    je1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();
    //direction 2 is k
    ks1 = inter.Dir2StartFirst();
    ke1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  }
  else if ( inter.BoundaryFirst() == 3 || inter.BoundaryFirst() == 4 ){ //direction 3 is j
    //extend min/maxes to cover ghost cells
    js1 = inter.ConstSurfaceFirst() + upLowFac;
    je1 = js1 + blk1.NumGhosts() - 1;
    //direction 1 is k
    ks1 = inter.Dir1StartFirst();
    ke1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();
    //direction 2 is i
    is1 = inter.Dir2StartFirst();
    ie1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  }
  else if ( inter.BoundaryFirst() == 5 || inter.BoundaryFirst() == 6 ){ //direction 3 is k
    //extend min/maxes to cover ghost cells
    ks1 = inter.ConstSurfaceFirst() + upLowFac;
    ke1 = ks1 + blk1.NumGhosts() - 1;
    //direction 1 is i
    is1 = inter.Dir1StartFirst();
    ie1 = inter.Dir1EndFirst() - 1 + 2 * blk1.NumGhosts();
    //direction 2 is j
    js1 = inter.Dir2StartFirst();
    je1 = inter.Dir2EndFirst() - 1 + 2 * blk1.NumGhosts();
  }
  else{
    cerr << "ERROR: Error in procBlock::SwapSlice(). Surface boundary " << inter.BoundaryFirst() << " is not recognized!" << endl;
  }

  //Get indices for slice coming from second block to swap
  int is2, ie2, js2, je2, ks2, ke2;

  //if at upper boundary no need to adjust for ghost cells as constant surface is already at the interior cells when acounting for ghost cells
  //if at the lower boundary adjust the constant surface by the number of ghost cells to get to the first interior cell 
  upLowFac = ( inter.BoundarySecond() % 2 == 0 ) ? 0 : blk2.NumGhosts();

  if ( inter.BoundarySecond() == 1 || inter.BoundarySecond() == 2 ){ //direction 3 is i
    //extend min/maxes to cover ghost cells
    is2 = inter.ConstSurfaceSecond() + upLowFac;
    ie2 = is2 + blk2.NumGhosts() - 1;
    //direction 1 is j
    js2 = inter.Dir1StartSecond();
    je2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();
    //direction 2 is k
    ks2 = inter.Dir2StartSecond();
    ke2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  }
  else if ( inter.BoundarySecond() == 3 || inter.BoundarySecond() == 4 ){ //direction 3 is j
    //extend min/maxes to cover ghost cells
    js2 = inter.ConstSurfaceSecond() + upLowFac;
    je2 = js2 + blk2.NumGhosts() - 1;
    //direction 1 is k
    ks2 = inter.Dir1StartSecond();
    ke2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();
    //direction 2 is i
    is2 = inter.Dir2StartSecond();
    ie2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  }
  else if ( inter.BoundarySecond() == 5 || inter.BoundarySecond() == 6 ){ //direction 3 is k
    //extend min/maxes to cover ghost cells
    ks2 = inter.ConstSurfaceSecond() + upLowFac;
    ke2 = ks2 + blk2.NumGhosts() - 1;
    //direction 1 is i
    is2 = inter.Dir1StartSecond();
    ie2 = inter.Dir1EndSecond() - 1 + 2 * blk2.NumGhosts();
    //direction 2 is j
    js2 = inter.Dir2StartSecond();
    je2 = inter.Dir2EndSecond() - 1 + 2 * blk2.NumGhosts();
  }
  else{
    cerr << "ERROR: Error in procBlock::SwapSlice(). Surface boundary " << inter.BoundarySecond() << " is not recognized!" << endl;
  }

  geomSlice geom1, geom2;
  stateSlice state1, state2;
  if (geom){ //get geomSlices to swap
    geom1 = blk1.GetGeomSlice(is1, ie1, js1, je1, ks1, ke1);
    geom2 = blk2.GetGeomSlice(is2, ie2, js2, je2, ks2, ke2);
  }
  else{   //get stateSlices to swap
    state1 = blk1.GetStateSlice(is1, ie1, js1, je1, ks1, ke1);
    state2 = blk2.GetStateSlice(is2, ie2, js2, je2, ks2, ke2);
  }

  //change interblocks to work with slice and ghosts
  interblock inter1 = inter;
  interblock inter2 = inter;

  //if at an upper surface, start block at upper boundary (after including ghosts), if at lower surface, start block at 0
  int blkStart = (inter1.BoundarySecond() % 2 == 0) ? inter1.ConstSurfaceSecond() + blk1.NumGhosts() : 0;
  inter1.SetConstSurfaceFirst(0); //slice always starts at 0
  inter1.SetConstSurfaceSecond( blkStart );
  //adjust direction 1 start and end for ghost cells
  inter1.SetDir1EndFirst( inter1.Dir1EndFirst() - inter1.Dir1StartFirst() + 2 * blk1.NumGhosts());
  inter1.SetDir1EndSecond( inter1.Dir1EndSecond() + 2 * blk1.NumGhosts());
  inter1.SetDir1StartFirst(0); //slice always starts at 0
  //adjust direction 2 start and end for ghost cells
  inter1.SetDir2EndFirst( inter1.Dir2EndFirst() - inter1.Dir2StartFirst() + 2 * blk1.NumGhosts());
  inter1.SetDir2EndSecond( inter1.Dir2EndSecond() + 2 * blk1.NumGhosts());
  inter1.SetDir2StartFirst(0); //slice always starts at 0
  inter1.SwapOrder(); //have block be first entry, slice second

  //if at an upper surface, start block at upper boundary (after including ghosts), if at lower surface, start block at 0
  blkStart = (inter2.BoundaryFirst() % 2 == 0) ? inter2.ConstSurfaceFirst() + blk2.NumGhosts() : 0;
  inter2.SetConstSurfaceSecond(0); //slice always starts at 0
  inter2.SetConstSurfaceFirst( blkStart );
  //adjust direction 1 start and end for ghost cells
  inter2.SetDir1EndSecond( inter2.Dir1EndSecond() - inter2.Dir1StartSecond() + 2 * blk2.NumGhosts());
  inter2.SetDir1EndFirst( inter2.Dir1EndFirst() + 2 * blk2.NumGhosts());
  inter2.SetDir1StartSecond(0); //slice always starts at 0
  //adjust direction 2 start and end for ghost cells
  inter2.SetDir2EndSecond( inter2.Dir2EndSecond() - inter2.Dir2StartSecond() + 2 * blk2.NumGhosts());
  inter2.SetDir2EndFirst( inter2.Dir2EndFirst() + 2 * blk2.NumGhosts());
  inter2.SetDir2StartSecond(0); //slice always starts at 0

  //put slices in proper blocks
  if (geom){ //put geomSlices in procBlock
    blk1.PutGeomSlice(geom2, inter2, blk2.NumGhosts());
    blk2.PutGeomSlice(geom1, inter1, blk1.NumGhosts());
  }
  else{ //put stateSlices in procBlock
    blk1.PutStateSlice(state2, inter2, blk2.NumGhosts());
    blk2.PutStateSlice(state1, inter1, blk1.NumGhosts());
  }

}

/* Function to swap slice using MPI. This is similar to the SwapSlice member function, but is called when the neighboring procBlocks are on different
processors.
*/
void procBlock::SwapSliceMPI( const interblock &inter, const int &rank, const MPI_Datatype &MPI_cellData ){
  // inter -- interblock boundary information
  // rank -- processor rank
  // MPI_cellData -- MPI datatype for passing primVars, genArray

  //Get indices for slice coming from block to swap
  int is, ie, js, je, ks, ke;

  if ( rank == inter.RankFirst() ){ //local block is first in interblock
    //if at upper boundary no need to adjust for ghost cells as constant surface is already at the interior cells when acounting for ghost cells
    //if at the lower boundary adjust the constant surface by the number of ghost cells to get to the first interior cell 
    int upLowFac = ( inter.BoundaryFirst() % 2 == 0 ) ? 0 : (*this).NumGhosts();

    if ( inter.BoundaryFirst() == 1 || inter.BoundaryFirst() == 2 ){ //direction 3 is i
      //extend min/maxes to cover ghost cells
      is = inter.ConstSurfaceFirst() + upLowFac;
      ie = is + (*this).NumGhosts() - 1;
      //direction 1 is j
      js = inter.Dir1StartFirst();
      je = inter.Dir1EndFirst() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is k
      ks = inter.Dir2StartFirst();
      ke = inter.Dir2EndFirst() - 1 + 2 * (*this).NumGhosts();
    }
    else if ( inter.BoundaryFirst() == 3 || inter.BoundaryFirst() == 4 ){ //direction 3 is j
      //extend min/maxes to cover ghost cells
      js = inter.ConstSurfaceFirst() + upLowFac;
      je = js + (*this).NumGhosts() - 1;
      //direction 1 is k
      ks = inter.Dir1StartFirst();
      ke = inter.Dir1EndFirst() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is i
      is = inter.Dir2StartFirst();
      ie = inter.Dir2EndFirst() - 1 + 2 * (*this).NumGhosts();
    }
    else if ( inter.BoundaryFirst() == 5 || inter.BoundaryFirst() == 6 ){ //direction 3 is k
      //extend min/maxes to cover ghost cells
      ks = inter.ConstSurfaceFirst() + upLowFac;
      ke = ks + (*this).NumGhosts() - 1;
      //direction 1 is i
      is = inter.Dir1StartFirst();
      ie = inter.Dir1EndFirst() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is j
      js = inter.Dir2StartFirst();
      je = inter.Dir2EndFirst() - 1 + 2 * (*this).NumGhosts();
    }
    else{
      cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Surface boundary " << inter.BoundaryFirst() << " is not recognized!" << endl;
    }
  }
  else if( rank == inter.RankSecond() ){ //local block is second in interblock
    //if at upper boundary no need to adjust for ghost cells as constant surface is already at the interior cells when acounting for ghost cells
    //if at the lower boundary adjust the constant surface by the number of ghost cells to get to the first interior cell 
    int upLowFac = ( inter.BoundarySecond() % 2 == 0 ) ? 0 : (*this).NumGhosts();

    if ( inter.BoundarySecond() == 1 || inter.BoundarySecond() == 2 ){ //direction 3 is i
      //extend min/maxes to cover ghost cells
      is = inter.ConstSurfaceSecond() + upLowFac;
      ie = is + (*this).NumGhosts() - 1;
      //direction 1 is j
      js = inter.Dir1StartSecond();
      je = inter.Dir1EndSecond() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is k
      ks = inter.Dir2StartSecond();
      ke = inter.Dir2EndSecond() - 1 + 2 * (*this).NumGhosts();
    }
    else if ( inter.BoundarySecond() == 3 || inter.BoundarySecond() == 4 ){ //direction 3 is j
      //extend min/maxes to cover ghost cells
      js = inter.ConstSurfaceSecond() + upLowFac;
      je = js + (*this).NumGhosts() - 1;
      //direction 1 is k
      ks = inter.Dir1StartSecond();
      ke = inter.Dir1EndSecond() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is i
      is = inter.Dir2StartSecond();
      ie = inter.Dir2EndSecond() - 1 + 2 * (*this).NumGhosts();
    }
    else if ( inter.BoundarySecond() == 5 || inter.BoundarySecond() == 6 ){ //direction 3 is k
      //extend min/maxes to cover ghost cells
      ks = inter.ConstSurfaceSecond() + upLowFac;
      ke = ks + (*this).NumGhosts() - 1;
      //direction 1 is i
      is = inter.Dir1StartSecond();
      ie = inter.Dir1EndSecond() - 1 + 2 * (*this).NumGhosts();
      //direction 2 is j
      js = inter.Dir2StartSecond();
      je = inter.Dir2EndSecond() - 1 + 2 * (*this).NumGhosts();
    }
    else{
      cerr << "ERROR: Error in procBlock::SwapSliceMPI(). Surface boundary " << inter.BoundarySecond() << " is not recognized!" << endl;
    }
  }
  else{
    cerr << "ERROR: Error in procBlock::SwapSliceMPIMPI(). Processor rank does not match either of interblock ranks!" << endl;
    exit(0);
  }

  //get local state slice to swap
  stateSlice state = (*this).GetStateSlice(is, ie, js, je, ks, ke);

  //swap with mpi_send_recv_replace
  //pack data into buffer, but first get size
  int bufSize = 0;
  int tempSize = 0;
  MPI_Pack_size(state.NumCells(), MPI_cellData, MPI_COMM_WORLD, &tempSize); //add size for states
  bufSize += tempSize;
  MPI_Pack_size(11, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for ints in class stateSlice
  bufSize += tempSize;

  char *buffer = new char[bufSize]; //allocate buffer to pack data into

  //pack data into buffer
  int position = 0;
  MPI_Pack(&state.state[0], state.NumCells(), MPI_cellData, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.numCells, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.numI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.numJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.numK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlock, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockStartI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockEndI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockStartJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockEndJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockStartK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&state.parBlockEndK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);

  MPI_Status status;
  if ( rank == inter.RankFirst() ){ //send/recv with second entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankSecond(), 1, inter.RankSecond(), 1, MPI_COMM_WORLD, &status);
  }
  else{ //send/recv with first entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankFirst(), 1, inter.RankFirst(), 1, MPI_COMM_WORLD, &status);
  }

  //put slice back into procBlock
  position = 0;
  MPI_Unpack(buffer, bufSize, &position, &state.state[0], state.NumCells(), MPI_cellData, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.numCells, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.numI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.numJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.numK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlock, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockStartI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockEndI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockStartJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockEndJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockStartK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &state.parBlockEndK, 1, MPI_INT, MPI_COMM_WORLD);

  delete [] buffer;

  //change interblocks to work with slice and ghosts
  interblock interAdj = inter;

  if ( rank == inter.RankSecond() ) { //block to insert into is first in interblock
    //if at an upper surface, start block at upper boundary (after including ghosts), if at lower surface, start block at 0
    int blkStart = (interAdj.BoundarySecond() % 2 == 0) ? interAdj.ConstSurfaceSecond() + (*this).NumGhosts() : 0;
    interAdj.SetConstSurfaceFirst(0); //slice always starts at 0
    interAdj.SetConstSurfaceSecond( blkStart );
    //adjust direction 1 start and end for ghost cells
    interAdj.SetDir1EndFirst( interAdj.Dir1EndFirst() - interAdj.Dir1StartFirst() + 2 * (*this).NumGhosts());
    interAdj.SetDir1EndSecond( interAdj.Dir1EndSecond() + 2 * (*this).NumGhosts());
    interAdj.SetDir1StartFirst(0); //slice always starts at 0
    //adjust direction 2 start and end for ghost cells
    interAdj.SetDir2EndFirst( interAdj.Dir2EndFirst() - interAdj.Dir2StartFirst() + 2 * (*this).NumGhosts());
    interAdj.SetDir2EndSecond( interAdj.Dir2EndSecond() + 2 * (*this).NumGhosts());
    interAdj.SetDir2StartFirst(0); //slice always starts at 0
    interAdj.SwapOrder(); //have block be first entry, slice second
  }
  else{ //block to insert into is second in interblock, so pass swapped version
    //if at an upper surface, start block at upper boundary (after including ghosts), if at lower surface, start block at 0
    int blkStart = (interAdj.BoundaryFirst() % 2 == 0) ? interAdj.ConstSurfaceFirst() + (*this).NumGhosts() : 0;
    interAdj.SetConstSurfaceSecond(0); //slice always starts at 0
    interAdj.SetConstSurfaceFirst( blkStart );
    //adjust direction 1 start and end for ghost cells
    interAdj.SetDir1EndSecond( interAdj.Dir1EndSecond() - interAdj.Dir1StartSecond() + 2 * (*this).NumGhosts());
    interAdj.SetDir1EndFirst( interAdj.Dir1EndFirst() + 2 * (*this).NumGhosts());
    interAdj.SetDir1StartSecond(0); //slice always starts at 0
    //adjust direction 2 start and end for ghost cells
    interAdj.SetDir2EndSecond( interAdj.Dir2EndSecond() - interAdj.Dir2StartSecond() + 2 * (*this).NumGhosts());
    interAdj.SetDir2EndFirst( interAdj.Dir2EndFirst() + 2 * (*this).NumGhosts());
    interAdj.SetDir2StartSecond(0); //slice always starts at 0
  }

  //interert stateSlice into procBlock
  (*this).PutStateSlice(state, interAdj, (*this).NumGhosts());

}

/* Function to return a vector of location indicies for ghost cells at an interblock boundary. The vector is formatted as shown below:

  vector = [i j k]

The vector will contain 3 entries corresponding to the i, j, and k locations of either the first or second pair in the interblock, depending on
what is specified in the pairID variable. The indices returned will correspond to cell locations and will take into account the orientation of the
patches that comprise the interblock with relation to each other.
*/
vector<int> GetSwapLoc( const int &l1, const int &l2, const int &l3, const interblock &inter, const bool &pairID ){
  // l1 -- index of direction 1 within slice to insert
  // l2 -- index of direction 2 within slice to insert
  // l3 -- index of direction 3 within slice to insert
  // inter -- interblock boundary condition
  // pairID -- returning index for first or second block in interblock match

  //preallocate vector to return
  vector<int> loc(3);

  if (pairID) { //working on first in pair -----------------------------------------------------------------------------------------------
    //first patch in pair is calculated using orientation 1
    if ( inter.BoundaryFirst() == 1 || inter.BoundaryFirst() == 2 ){ //i-patch
      //get direction 1 length
      loc[1] = inter.Dir1StartFirst() + l1; //direction 1 is j
      loc[2] = inter.Dir2StartFirst() + l2; //direction 2 is k
      loc[0] = inter.ConstSurfaceFirst() + l3 ; //add l3 to get to ghost cells (cell index instead of face)
    }
    else if ( inter.BoundaryFirst() == 3 || inter.BoundaryFirst() == 4 ){ //j-patch
      //get direction 1 length
      loc[2] = inter.Dir1StartFirst() + l1; //direction 1 is k
      loc[0] = inter.Dir2StartFirst() + l2; //direction 2 is i
      loc[1] = inter.ConstSurfaceFirst() + l3 ; //add l3 to get to ghost cells (cell index instead of face)
    }
    else if ( inter.BoundaryFirst() == 5 || inter.BoundaryFirst() == 6 ){ //k-patch
      //get direction 1 length
      loc[0] = inter.Dir1StartFirst() + l1; //direction 1 is i
      loc[1] = inter.Dir2StartFirst() + l2; //direction 2 is j
      loc[2] = inter.ConstSurfaceFirst() + l3 ; //add l3 to get to ghost cells (cell index instead of face)
    }
    else{ 
      cerr << "ERROR: Error in procBlock:GetSwapLoc(). Boundary surface " << inter.BoundaryFirst() << " is not recognized!" << endl;
      cerr << "Please choose a number between 1-6." << endl;
    }
  }
  //--------------------------------------------------------------------------------------------------------------------------------------------
  //need to use orientation for second in pair
  else{ //working on second in pair

    //-------------------------------------------------------------------------------------------------------
    if ( inter.BoundarySecond() == 1 || inter.BoundarySecond() == 2){ //i-patch

      if ( inter.Orientation() == 2 || inter.Orientation() == 4 || inter.Orientation() == 5 || inter.Orientation() == 7 ){ //swap dir 1 and 2
	//direction 1 is j (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[2] = ( inter.Orientation() == 5 || inter.Orientation() == 7 ) ? inter.Dir2EndSecond() - 1 - l1 : inter.Dir2StartSecond() + l1 ;

	//direction 2 is k (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[1] = ( inter.Orientation() == 4 || inter.Orientation() == 7 ) ? inter.Dir1EndSecond() - 1 - l2 : inter.Dir1StartSecond() + l2 ;
      }
      else{ //no direction swap
	//direction 1 is j -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[1] = ( inter.Orientation() == 6 || inter.Orientation() == 8 ) ? inter.Dir1EndSecond() - 1 - l1 : inter.Dir1StartSecond() + l1 ;

	//direction 1 is k -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[2] = ( inter.Orientation() == 3 || inter.Orientation() == 8 ) ? inter.Dir2EndSecond() - 1 - l2 : inter.Dir2StartSecond() + l2 ;
      }

      //calculate index for all ghost layers
      loc[0] = inter.ConstSurfaceSecond() + l3 ; //add l3 to get to ghost cells
    }
    //-------------------------------------------------------------------------------------------------------
    else if ( inter.BoundarySecond() == 3 || inter.BoundarySecond() == 4 ){ //j-patch

      if ( inter.Orientation() == 2 || inter.Orientation() == 4 || inter.Orientation() == 5 || inter.Orientation() == 7 ){ //swap dir 1 and 2
	//direction 1 is k (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[0] = ( inter.Orientation() == 5 || inter.Orientation() == 7 ) ? inter.Dir2EndSecond() - 1 - l1 : inter.Dir2StartSecond() + l1 ;

	//direction 2 is i (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[2] = ( inter.Orientation() == 4 || inter.Orientation() == 7 ) ? inter.Dir1EndSecond() - 1 - l2 : inter.Dir1StartSecond() + l2 ;
      }
      else{ //no direction swap
	//direction 1 is k -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[2] = ( inter.Orientation() == 3 || inter.Orientation() == 8 ) ? inter.Dir1EndSecond() - 1 - l1 : inter.Dir1StartSecond() + l1 ;

	//direction 2 is i -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[0] = ( inter.Orientation() == 6 || inter.Orientation() == 8 ) ? inter.Dir2EndSecond() - 1 - l2 : inter.Dir2StartSecond() + l2 ;
      }

      //calculate index for all ghost layers
      loc[1] = inter.ConstSurfaceSecond() + l3 ; //add l3 to get to ghost cells
    }
    //-------------------------------------------------------------------------------------------------------
    else if ( inter.BoundarySecond() == 5 || inter.BoundarySecond() == 6 ){ //k-patch

      if ( inter.Orientation() == 2 || inter.Orientation() == 4 || inter.Orientation() == 5 || inter.Orientation() == 7 ){ //swap dir 1 and 2
	//direction 1 is i (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[1] = ( inter.Orientation() == 5 || inter.Orientation() == 7 ) ? inter.Dir2EndSecond() - 1 - l1 : inter.Dir2StartSecond() + l1 ;

	//direction 2 is j (swapped) -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[0] = ( inter.Orientation() == 4 || inter.Orientation() == 7 ) ? inter.Dir1EndSecond() - 1 - l2 : inter.Dir1StartSecond() + l2 ;
      }
      else{ //no direction swap
	//direction 1 is i -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[0] = ( inter.Orientation() == 3 || inter.Orientation() == 8 ) ? inter.Dir1EndSecond() - 1 - l1 : inter.Dir1StartSecond() + l1 ;

	//direction 2 is j -- if true direction reversed -- subtract 1 from End to get to cell index
	loc[1] = ( inter.Orientation() == 6 || inter.Orientation() == 8 ) ? inter.Dir2EndSecond() - 1 - l2 : inter.Dir2StartSecond() + l2 ;
      }

      //calculate index for all ghost layers
      loc[2] = inter.ConstSurfaceSecond() + l3 ; //add l3 to get to ghost cells
    }
    //-------------------------------------------------------------------------------------------------------
    else{
      cerr << "ERROR: Error in procBlock.cpp:GetSwapLoc(). Boundary surface of " << inter.BoundarySecond() << " is not recognized!" << endl;
      cerr << "Valid numbers are between 1 and 6." << endl;
      exit(0);
    }
  }

  return loc;
}

/* Function to populate ghost cells with proper cell states for inviscid flow calculation. This function operates on the entire grid and uses interblock
boundaries to pass the correct data between grid blocks.
*/
void GetBoundaryConditions(vector<procBlock> &states, const input &inp, const idealGas &eos, const vector<interblock> &connections, const int &rank, const MPI_Datatype &MPI_cellData){
  // states -- vector of all procBlocks in the solution domain
  // inp -- all input variables
  // eos -- equation of state
  // connections -- vector of interblock connections

  //loop over all blocks and assign inviscid ghost cells
  for ( unsigned int ii = 0; ii < states.size(); ii++ ){
    states[ii].AssignInviscidGhostCells(inp, eos);
  }

  //loop over connections and swap ghost cells where needed
  for ( unsigned int ii = 0; ii < connections.size(); ii++ ){
    if ( connections[ii].RankFirst() == rank && connections[ii].RankSecond() == rank ) { //both sides of interblock are on same processor, swap w/o mpi
      SwapSlice( connections[ii], states[connections[ii].LocalBlockFirst()], states[connections[ii].LocalBlockSecond()], false);
    }
    else if ( connections[ii].RankFirst() == rank ){ //rank matches rank of one side of interblock, swap over mpi
      states[connections[ii].LocalBlockFirst()].SwapSliceMPI( connections[ii], rank, MPI_cellData ); 
    }
    else if ( connections[ii].RankSecond() == rank ) { //rank matches rank of one side of interblock, swap over mpi
      states[connections[ii].LocalBlockSecond()].SwapSliceMPI( connections[ii], rank, MPI_cellData );
    }
    //if rank doesn't match either side of interblock, then do nothing and move on to the next interblock
  }

  //loop over all blocks and get ghost cell edge data
  for ( unsigned int ii = 0; ii < states.size(); ii++) {
    states[ii].AssignInviscidGhostCellsEdge(inp, eos);
  }

}

/* Member function to get a slice (portion) of the geometry of a procBlock. The geom slice remains in the orientation of its parent block unless
the revI, revJ, or revK flags are activated to reverse either the i, j, or k directions respectively. This function essentially cuts out a portion
of a procBlock and returns a geomSlice with that portion of the block geometry.
*/
geomSlice procBlock::GetGeomSlice(const int &is, const int &ie, const int &js, const int &je, const int &ks, const int &ke, 
				  const bool revI, const bool revJ, const bool revK ) const {
  // is -- starting i-cell index for slice
  // ie -- ending i-cell index for slice
  // js -- starting j-cell index for slice
  // je -- ending j-cell index for slice
  // ks -- starting k-cell index for slice
  // ke -- ending k-cell index for slice
  // revI -- flag to reverse i direction of indices (default false)
  // revJ -- flag to reverse j direction of indices (default false)
  // revK -- flag to reverse k direction of indices (default false)

  int sizeI = ie - is + 1;
  int sizeJ = je - js + 1;
  int sizeK = ke - ks + 1;

  //initialize slice
  geomSlice slice(sizeI, sizeJ, sizeK, (*this).ParentBlock(), is, ie+1, js, je+1, ks, ke+1 );

  //get parent block maxes
  int imaxPar = (*this).NumI() + 2.0 * (*this).NumGhosts();
  int jmaxPar = (*this).NumJ() + 2.0 * (*this).NumGhosts();

  //loop over all cells in slice and populate
  for ( int kk = 0; kk < sizeK; kk++ ){
    for ( int jj = 0; jj < sizeJ; jj++ ){
      for ( int ii = 0; ii < sizeI; ii++ ){

	//determine if direction needs to be reversed
	int k = revK ? sizeK - 1 - kk : kk;
	double kFac = revK ? -1.0 : 1.0;

	int j = revJ ? sizeJ - 1 - jj : jj;
	double jFac = revJ ? -1.0 : 1.0;

	int i = revI ? sizeI - 1 - ii : ii;
	double iFac = revI ? -1.0 : 1.0;

	//cell locations
	int locPar = GetLoc1D(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int loc = GetLoc1D(ii, jj, kk, sizeI, sizeJ);

	//lower i-face locations
	int lowIPar = GetLowerFaceI(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int lowI = GetLowerFaceI(ii, jj, kk, sizeI, sizeJ);

	//upper i-face locations
	int upIPar = GetUpperFaceI(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int upI = GetUpperFaceI(ii, jj, kk, sizeI, sizeJ);

	//lower j-face locations
	int lowJPar = GetLowerFaceJ(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int lowJ = GetLowerFaceJ(ii, jj, kk, sizeI, sizeJ);

	//upper j-face locations
	int upJPar = GetUpperFaceJ(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int upJ = GetUpperFaceJ(ii, jj, kk, sizeI, sizeJ);

	//lower k-face locations
	int lowKPar = GetLowerFaceK(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int lowK = GetLowerFaceK(ii, jj, kk, sizeI, sizeJ);

	//upper k-face locations
	int upKPar = GetUpperFaceK(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int upK = GetUpperFaceK(ii, jj, kk, sizeI, sizeJ);

	//assign cell variables
	slice.SetVol( (*this).Vol(locPar), loc );
	slice.SetCenter( (*this).Center(locPar), loc );

	//assign i-face variables
	if ( revI ){ //if direction is reversed, upper/lower faces need to be swapped
	  slice.SetFAreaI( iFac * (*this).FAreaI(upIPar), lowI );
	  slice.SetFCenterI( (*this).FCenterI(upIPar), lowI );

	  if ( ii == sizeI - 1 ){ //at end of i-line assign upper face values too
	    slice.SetFAreaI( iFac * (*this).FAreaI(lowIPar), upI );
	    slice.SetFCenterI( (*this).FCenterI(locPar), upI );
	  }
	}
	else{
	  slice.SetFAreaI( iFac * (*this).FAreaI(lowIPar), lowI );
	  slice.SetFCenterI( (*this).FCenterI(lowIPar), lowI );

	  if ( ii == sizeI - 1 ){ //at end of i-line assign upper face values too
	    slice.SetFAreaI( iFac * (*this).FAreaI(upIPar), upI );
	    slice.SetFCenterI( (*this).FCenterI(upIPar), upI );
	  }
	}

	//assign j-face variables
	if ( revJ ){ //if direction is reversed, upper/lower faces need to be swapped
	  slice.SetFAreaJ( jFac * (*this).FAreaJ(upJPar), lowJ );
	  slice.SetFCenterJ( (*this).FCenterJ(upJPar), lowJ );

	  if ( jj == sizeJ - 1 ){ //at end of j-line assign upper face values too
	    slice.SetFAreaJ( jFac * (*this).FAreaJ(lowJPar), upJ );
	    slice.SetFCenterJ( (*this).FCenterJ(lowJPar), upJ );
	  }
	}
	else{
	  slice.SetFAreaJ( jFac * (*this).FAreaJ(lowJPar), lowJ );
	  slice.SetFCenterJ( (*this).FCenterJ(lowJPar), lowJ );

	  if ( jj == sizeJ - 1 ){ //at end of j-line assign upper face values too
	    slice.SetFAreaJ( jFac * (*this).FAreaJ(upJPar), upJ );
	    slice.SetFCenterJ( (*this).FCenterJ(upJPar), upJ );
	  }
	}

	//assign k-face variables
	if ( revK ){ //if direction is reversed, upper/lower faces need to be swapped
	  slice.SetFAreaK( kFac * (*this).FAreaK(upKPar), lowK ); 
	  slice.SetFCenterK( (*this).FCenterK(upKPar), lowK );

	  if ( kk == sizeK - 1 ){ //at end of k-line assign upper face values too
	    slice.SetFAreaK( kFac * (*this).FAreaK(lowKPar), upK );
	    slice.SetFCenterK( (*this).FCenterK(lowKPar), upK );
	  }
	}
	else{
	  slice.SetFAreaK( kFac * (*this).FAreaK(lowKPar), lowK );
	  slice.SetFCenterK( (*this).FCenterK(lowKPar), lowK );

	  if ( kk == sizeK - 1 ){ //at end of k-line assign upper face values too
	    slice.SetFAreaK( kFac * (*this).FAreaK(upKPar), upK );
	    slice.SetFCenterK( (*this).FCenterK(upKPar), upK );
	  }
	}

      }
    }
  }

  return slice;
}

/* Member function to get a slice (portion) of the states of a procBlock. The state slice remains in the orientation of its parent block unless
the revI, revJ, or revK flags are activated to reverse either the i, j, or k directions respectively. This function essentially cuts out a portion
of a procBlock and returns a stateSlice with that portion of the block states.
*/
stateSlice procBlock::GetStateSlice(const int &is, const int &ie, const int &js, const int &je, const int &ks, const int &ke,
				    const bool revI, const bool revJ, const bool revK ) const {
  // is -- starting i-cell index for slice
  // ie -- ending i-cell index for slice
  // js -- starting j-cell index for slice
  // je -- ending j-cell index for slice
  // ks -- starting k-cell index for slice
  // ke -- ending k-cell index for slice
  // revI -- flag to reverse i direction of indices (default false)
  // revJ -- flag to reverse j direction of indices (default false)
  // revK -- flag to reverse k direction of indices (default false)

  int sizeI = ie - is + 1;
  int sizeJ = je - js + 1;
  int sizeK = ke - ks + 1;

  //initialize state slice
  stateSlice states(sizeI, sizeJ, sizeK, (*this).ParentBlock(), is, ie+1, js, je+1, ks, ke+1 );

  //get parent block maxes
  int imaxPar = (*this).NumI() + 2.0 * (*this).NumGhosts();
  int jmaxPar = (*this).NumJ() + 2.0 * (*this).NumGhosts();

  //loop over all cells in slice and populate
  for ( int kk = 0; kk < sizeK; kk++ ){
    for ( int jj = 0; jj < sizeJ; jj++ ){
      for ( int ii = 0; ii < sizeI; ii++ ){

	//determine if direction needs to be reversed
	int k = revK ? sizeK - 1 - kk : kk;
	int j = revJ ? sizeJ - 1 - jj : jj;
	int i = revI ? sizeI - 1 - ii : ii;

	//cell locations
	int locPar = GetLoc1D(is+i, js+j, ks+k, imaxPar, jmaxPar);
	int loc = GetLoc1D(ii, jj, kk, sizeI, sizeJ);

	//assign cell variables
	states.SetState( (*this).State(locPar), loc);

      }
    }
  }

  return states;
}

/* Member function to overwrite a section of a procBlock's geometry with a geomSlice. The function uses the orientation supplied in the interblock to orient the 
geomSlice relative to the procBlock. It assumes that the procBlock is listed first, and the geomSlice second in the interblock data structure.
*/
void procBlock::PutGeomSlice( const geomSlice &slice, const interblock& inter, const int &d3){
  // slice -- geomSlice to insert int procBlock
  // inter -- interblock data structure describing the patches and their orientation
  // d3 -- distance of direction normal to patch to insert

  //check that number of cells to insert matches
  int blkCell = (inter.Dir1EndFirst() - inter.Dir1StartFirst()) * (inter.Dir2EndFirst() - inter.Dir2StartFirst()) * d3;
  if ( blkCell != slice.NumCells() ){
    cerr << "ERROR: Error in procBlock::PutGeomSlice(). Number of cells being inserted does not match designated space to insert to." << endl;
    cerr << "Direction 1, 2, 3 of procBlock: " << inter.Dir1EndFirst() - inter.Dir1StartFirst() << ", " << inter.Dir2EndFirst() - inter.Dir2StartFirst()
	 << ", " << d3 << endl;
    cerr << "Direction I, J, K of geomSlice: " << slice.NumI() << ", " << slice.NumJ() << ", " << slice.NumK() << endl;
    exit(0);
  }

  //get procBlock maxes
  int imaxB = (*this).NumI() + 2.0 * (*this).NumGhosts();
  int jmaxB = (*this).NumJ() + 2.0 * (*this).NumGhosts();

  //get slice maxes
  int imaxS = slice.NumI();
  int jmaxS = slice.NumJ();

  //determine if area direction needs to be reversed
  double aFac3 = ( (inter.BoundaryFirst() + inter.BoundarySecond()) % 2 == 0 ) ? -1.0 : 1.0 ;
  double aFac1 = ( inter.Orientation() == 3 || inter.Orientation() == 4 || inter.Orientation() == 7 || inter.Orientation() == 8 ) ? -1.0 : 1.0 ;
  double aFac2 = ( inter.Orientation() == 5 || inter.Orientation() == 6 || inter.Orientation() == 7 || inter.Orientation() == 8 ) ? -1.0 : 1.0 ;

  //loop over cells to insert
  for ( int l3 = 0; l3 < d3; l3++ ){
    for ( int l2 = 0; l2 < (inter.Dir2EndFirst() - inter.Dir2StartFirst()); l2++ ){
      for ( int l1 = 0; l1 < (inter.Dir1EndFirst() - inter.Dir1StartFirst()); l1++ ){

	//get block and slice indices
	vector<int> indB = GetSwapLoc(l1, l2, l3, inter, true);
	vector<int> indS = GetSwapLoc(l1, l2, l3, inter, false);

	//get cell locations
	int locB = GetLoc1D(indB[0], indB[1], indB[2], imaxB, jmaxB);
	int locS = GetLoc1D(indS[0], indS[1], indS[2], imaxS, jmaxS);

	//get i-face locations
	int IlowB = GetLowerFaceI(indB[0], indB[1], indB[2], imaxB, jmaxB);
	int IupB  = GetUpperFaceI(indB[0], indB[1], indB[2], imaxB, jmaxB);

	int IlowS = GetLowerFaceI(indS[0], indS[1], indS[2], imaxS, jmaxS);
	int IupS  = GetUpperFaceI(indS[0], indS[1], indS[2], imaxS, jmaxS);

	//get j-face locations
	int JlowB = GetLowerFaceJ(indB[0], indB[1], indB[2], imaxB, jmaxB);
	int JupB  = GetUpperFaceJ(indB[0], indB[1], indB[2], imaxB, jmaxB);

	int JlowS = GetLowerFaceJ(indS[0], indS[1], indS[2], imaxS, jmaxS);
	int JupS  = GetUpperFaceJ(indS[0], indS[1], indS[2], imaxS, jmaxS);

	//get k-face locations
	int KlowB = GetLowerFaceK(indB[0], indB[1], indB[2], imaxB, jmaxB);
	int KupB  = GetUpperFaceK(indB[0], indB[1], indB[2], imaxB, jmaxB);

	int KlowS = GetLowerFaceK(indS[0], indS[1], indS[2], imaxS, jmaxS);
	int KupS  = GetUpperFaceK(indS[0], indS[1], indS[2], imaxS, jmaxS);

	//swap cell data
	(*this).SetVol( slice.Vol(locS), locB);
	(*this).SetCenter( slice.Center(locS), locB);

	//----------------------------------------------------------------------------------------------------------------------------------------
	//swap face data
	if ( inter.BoundaryFirst() <= 2 && inter.BoundarySecond() <= 2 ){ //both patches i, i to i, j to j, k to k

	  //swap face data for direction 3
	  (*this).SetFCenterI( slice.FCenterI(IlowS) , IlowB);
	  (*this).SetFAreaI( aFac3 * slice.FAreaI(IlowS) , IlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterI( slice.FCenterI(IupS) , IupB);
	    (*this).SetFAreaI( aFac3 * slice.FAreaI(IupS) , IupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterJ(JlowS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaJ(JlowS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterJ(JupS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaJ(JupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterJ(JupS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaJ(JupS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterJ(JlowS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaJ(JlowS), JupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterK(KlowS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaK(KlowS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterK(KupS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaK(KupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterK(KupS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaK(KupS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterK(KlowS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaK(KlowS), KupB );
	    }
	  }

	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 2 && inter.BoundaryFirst() <= 4 && inter.BoundarySecond() > 2 && inter.BoundarySecond() <= 4 ){ //both patches j, j to j, k to k, i to i

	  //swap face data for direction 3
	  (*this).SetFCenterJ( slice.FCenterJ(JlowS) , JlowB);
	  (*this).SetFAreaJ( aFac3 * slice.FAreaJ(JlowS) , JlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterJ( slice.FCenterJ(JupS) , JupB);
	    (*this).SetFAreaJ( aFac3 * slice.FAreaJ(JupS) , JupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterK(KlowS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaK(KlowS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterK(KupS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaK(KupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterK(KupS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaK(KupS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterK(KlowS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaK(KlowS), KupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterI(IlowS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaI(IlowS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterI(IupS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaI(IupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterI(IupS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaI(IupS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterI(IlowS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaI(IlowS), IupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 4 && inter.BoundaryFirst() <= 6 && inter.BoundarySecond() > 4 && inter.BoundarySecond() <= 6 ){ //both patches k, k to k, i to i, j to j

	  //swap face data for direction 3
	  (*this).SetFCenterK( slice.FCenterK(KlowS) , KlowB);
	  (*this).SetFAreaK( aFac3 * slice.FAreaK(KlowS) , KlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterK( slice.FCenterK(KupS) , KupB);
	    (*this).SetFAreaK( aFac3 * slice.FAreaK(KupS) , KupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterI(IlowS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaI(IlowS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterI(IupS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaI(IupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterI(IupS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaI(IupS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterI(IlowS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaI(IlowS), IupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterJ(JlowS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaJ(JlowS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterJ(JupS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaJ(JupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterJ(JupS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaJ(JupS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterJ(JlowS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaJ(JlowS), JupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() <= 2 && inter.BoundarySecond() > 2 && inter.BoundarySecond() <= 4){ //patches are i/j  - i to j, j to k, k to i

	  //swap face data for direction 3
	  (*this).SetFCenterI( slice.FCenterJ(JlowS) , IlowB);
	  (*this).SetFAreaI( aFac3 * slice.FAreaJ(JlowS) , IlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterI( slice.FCenterJ(JupS) , IupB);
	    (*this).SetFAreaI( aFac3 * slice.FAreaJ(JupS) , IupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterK(KlowS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaK(KlowS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterK(KupS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaK(KupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterK(KupS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaK(KupS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterK(KlowS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaK(KlowS), JupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterI(IlowS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaI(IlowS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterI(IupS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaI(IupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterI(IupS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaI(IupS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterI(IlowS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaI(IlowS), KupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() <= 2 && inter.BoundarySecond() > 4 && inter.BoundarySecond() <= 6){ //patches are i/k  - i to k, j to i, k to j

	  //swap face data for direction 3
	  (*this).SetFCenterI( slice.FCenterK(KlowS) , IlowB);
	  (*this).SetFAreaI( aFac3 * slice.FAreaK(KlowS) , IlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterI( slice.FCenterK(KupS) , IupB);
	    (*this).SetFAreaI( aFac3 * slice.FAreaK(KupS) , IupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterI(IlowS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaI(IlowS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterI(IupS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaI(IupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterI(IupS), JlowB );
	    (*this).SetFAreaJ( aFac1 * slice.FAreaI(IupS), JlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterJ( slice.FCenterI(IlowS), JupB );
	      (*this).SetFAreaJ( aFac1 * slice.FAreaI(IlowS), JupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterJ(JlowS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaJ(JlowS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterJ(JupS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaJ(JupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterJ(JupS), KlowB );
	    (*this).SetFAreaK( aFac2 * slice.FAreaJ(JupS), KlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterK( slice.FCenterJ(JlowS), KupB );
	      (*this).SetFAreaK( aFac2 * slice.FAreaJ(JlowS), KupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 2 && inter.BoundaryFirst() <= 4 && inter.BoundarySecond() <= 2 ){ //patches are j/i, j to i, k to j, i to k

	  //swap face data for direction 3
	  (*this).SetFCenterJ( slice.FCenterI(IlowS) , JlowB);
	  (*this).SetFAreaJ( aFac3 * slice.FAreaI(IlowS) , JlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterJ( slice.FCenterI(IupS) , JupB);
	    (*this).SetFAreaJ( aFac3 * slice.FAreaI(IupS) , JupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterJ(JlowS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaJ(JlowS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterJ(JupS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaJ(JupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterJ(JupS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaJ(JupS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterJ(JlowS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaJ(JlowS), KupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterK(KlowS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaK(KlowS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterK(KupS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaK(KupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterK(KupS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaK(KupS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterK(KlowS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaK(KlowS), IupB );
	    }
	  }

	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 2 && inter.BoundaryFirst() <= 4 && inter.BoundarySecond() > 4 && inter.BoundarySecond() <= 6 ){ //patches are j/k, j to k, k to i, i to j

	  //swap face data for direction 3
	  (*this).SetFCenterJ( slice.FCenterK(KlowS) , JlowB);
	  (*this).SetFAreaJ( aFac3 * slice.FAreaK(KlowS) , JlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterJ( slice.FCenterJ(JupS) , JupB);
	    (*this).SetFAreaJ( aFac3 * slice.FAreaJ(JupS) , JupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterK( slice.FCenterI(IlowS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaI(IlowS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterI(IupS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaI(IupS), KupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterK( slice.FCenterI(IupS), KlowB );
	    (*this).SetFAreaK( aFac1 * slice.FAreaI(IupS), KlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterK( slice.FCenterI(IlowS), KupB );
	      (*this).SetFAreaK( aFac1 * slice.FAreaI(IlowS), KupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterJ(JlowS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaJ(JlowS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterJ(JupS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaJ(JupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterJ(JupS), IlowB );
	    (*this).SetFAreaI( aFac2 * slice.FAreaJ(JupS), IlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterI( slice.FCenterJ(JlowS), IupB );
	      (*this).SetFAreaI( aFac2 * slice.FAreaJ(JlowS), IupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 4 && inter.BoundaryFirst() <= 6 && inter.BoundarySecond() <= 2 ){ //patches are k/i, k to i, i to j, j to k

	  //swap face data for direction 3
	  (*this).SetFCenterK( slice.FCenterI(IlowS) , KlowB);
	  (*this).SetFAreaK( aFac3 * slice.FAreaI(IlowS) , KlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterK( slice.FCenterI(IupS) , KupB);
	    (*this).SetFAreaK( aFac3 * slice.FAreaI(IupS) , KupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterJ(JlowS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaJ(JlowS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterJ(JupS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaJ(JupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterJ(JupS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaJ(JupS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterJ(JlowS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaJ(JlowS), IupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterK(KlowS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaK(KlowS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterK(KupS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaK(KupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterK(KupS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaK(KupS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterK(KlowS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaK(KlowS), JupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else if ( inter.BoundaryFirst() > 4 && inter.BoundaryFirst() <= 6 && inter.BoundarySecond() <= 2 ){ //patches are k/j, k to j, i to k, j to i

	  //swap face data for direction 3
	  (*this).SetFCenterK( slice.FCenterJ(JlowS) , KlowB);
	  (*this).SetFAreaK( aFac3 * slice.FAreaJ(JlowS) , KlowB);

	  if ( l3 == (d3 - 1) ){ //at end of direction 3 line
	    (*this).SetFCenterK( slice.FCenterJ(JupS) , KupB);
	    (*this).SetFAreaK( aFac3 * slice.FAreaJ(JupS) , KupB);
	  }

	  //swap face data for direction 1
	  if ( aFac1 == 1.0 ){
	    (*this).SetFCenterI( slice.FCenterK(KlowS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaK(KlowS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterK(KupS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaK(KupS), IupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterI( slice.FCenterK(KupS), IlowB );
	    (*this).SetFAreaI( aFac1 * slice.FAreaK(KupS), IlowB );

	    if ( l1 == (inter.Dir1EndFirst() - inter.Dir1StartFirst() - 1) ){//at end of direction 1 line
	      (*this).SetFCenterI( slice.FCenterK(KlowS), IupB );
	      (*this).SetFAreaI( aFac1 * slice.FAreaK(KlowS), IupB );
	    }
	  }

	  //swap face data for direction 2
	  if ( aFac2 == 1.0 ){
	    (*this).SetFCenterJ( slice.FCenterI(IlowS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaI(IlowS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterI(IupS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaI(IupS), JupB );
	    }
	  }
	  else{ //if direction is reversed, upper/lower faces need to be swapped
	    (*this).SetFCenterJ( slice.FCenterI(IupS), JlowB );
	    (*this).SetFAreaJ( aFac2 * slice.FAreaI(IupS), JlowB );

	    if ( l2 == (inter.Dir2EndFirst() - inter.Dir2StartFirst() - 1) ){//at end of direction 2 line
	      (*this).SetFCenterJ( slice.FCenterI(IlowS), JupB );
	      (*this).SetFAreaJ( aFac2 * slice.FAreaI(IlowS), JupB );
	    }
	  }
	}
	//----------------------------------------------------------------------------------------------------------------------------------------
	else{
	  cerr << "ERROR: Error in procBlock::PutGeomSlice(). Unable to swap face quantities because behavior for interface with boundary pair "
	       << inter.BoundaryFirst() << ", " << inter.BoundarySecond() << " is not defined." << endl;
	  exit(0);
	}

      }
    }
  }

}

/* Member function to overwrite a section of a procBlock's states with a stateSlice. The function uses the orientation supplied in the interblock to orient the 
stateSlice relative to the procBlock. It assumes that the procBlock is listed first, and the stateSlice second in the interblock data structure.
*/
void procBlock::PutStateSlice( const stateSlice &slice, const interblock &inter, const int &d3){
  // slice -- geomSlice to insert int procBlock
  // inter -- interblock data structure defining the patches and their orientation
  // d3 -- distance of direction normal to patch to insert

  //check that number of cells to insert matches
  int blkCell = (inter.Dir1EndFirst() - inter.Dir1StartFirst()) * (inter.Dir2EndFirst() - inter.Dir2StartFirst()) * d3;
  if ( blkCell != slice.NumCells() ){
    cerr << "ERROR: Error in procBlock::PutStateSlice(). Number of cells being inserted does not match designated space to insert to." << endl;
    cerr << "Direction 1, 2, 3 of procBlock: " << inter.Dir1EndFirst() - inter.Dir1StartFirst() << ", " << inter.Dir2EndFirst() - inter.Dir2StartFirst()
	 << ", " << d3 << endl;
    cerr << "Direction I, J, K of geomSlice: " << slice.NumI() << ", " << slice.NumJ() << ", " << slice.NumK() << endl;
    exit(0);
  }

  //get block maxes
  int imaxB = (*this).NumI() + 2.0 * (*this).NumGhosts();
  int jmaxB = (*this).NumJ() + 2.0 * (*this).NumGhosts();

  //get slice maxes
  int imaxS = slice.NumI();
  int jmaxS = slice.NumJ();

  //loop over cells to insert
  for ( int l3 = 0; l3 < d3; l3++ ){
    for ( int l2 = 0; l2 < (inter.Dir2EndFirst() - inter.Dir2StartFirst()); l2++ ){
      for ( int l1 = 0; l1 < (inter.Dir1EndFirst() - inter.Dir1StartFirst()); l1++ ){

	//get block and slice indices
	vector<int> indB = GetSwapLoc(l1, l2, l3, inter, true);
	vector<int> indS = GetSwapLoc(l1, l2, l3, inter, false);

	//get cell locations
	int locB = GetLoc1D(indB[0], indB[1], indB[2], imaxB, jmaxB);
	int locS = GetLoc1D(indS[0], indS[1], indS[2], imaxS, jmaxS);

	//swap cell data
	(*this).SetState( slice.State(locS), locB);

      }
    }
  }

}

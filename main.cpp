#include <iostream>      //cout, cerr, endl
#include "plot3d.h"
#include "vector3d.h"
#include "tensor.h"
#include "input.h"
#include <vector>        //stl vector
#include <algorithm>     //max_element
#include <numeric>       //accumulate
#include "procBlock.h"
#include "inviscidFlux.h"
#include "viscousFlux.h"
#include "primVars.h"
#include "eos.h"
#include "boundaryConditions.h"
#include "output.h"
#include "matrix.h"
#include "parallel.h"
#include <fenv.h>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::clock_t;
using std::clock;

int main( int argc, char *argv[] ) {

  //initialize MPI and make calls to get number of processors and rank of each processor
  int numProcs, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //get MPI version
  int version, subversion;
  MPI_Get_version(&version, &subversion);
  if (rank == ROOT ){
    cout << "Using MPI Version " << version << "." << subversion << endl;
  }

  cout << "Hello from processor " << rank << " of " << numProcs << "!" << endl;

  //start clock to time simulation
  clock_t start;
  double duration;
  start = clock();

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); //enable exceptions so code won't run with NANs

  double totalCells = 0.0;
  input inputVars;
  vector<int> loadBal;
  int numProcBlock = 0;

  string inputFile = argv[1];  //name of input file is the second argument (the executable being the first)

  //broadcast input file name to all names for portability
  BroadcastString(inputFile);

  //Parse input file
  inputVars.ReadInput(inputFile, rank);

  //Determine number of equations
  int numEqns = 0;
  if ( (inputVars.EquationSet() == "euler") || (inputVars.EquationSet() == "navierStokes") ){
    numEqns = 5;
  }
  else{
    cerr << "ERROR: Equations set is not recognized. Cannot determine number of equations!" << endl;
  }

  //determine number of ghost cells
  int numGhost = 2;

  //get equation of state
  idealGas eos(inputVars.Gamma(),inputVars.R());

  //Initialize state vector with nondimensional variables
  //get reference speed of sound
  double aRef = eos.GetSoS(inputVars.PRef(),inputVars.RRef());

  //get a reference velocity
  vector3d<double> velRef = inputVars.VelRef();

  //initialize a single state
  primVars state(1.0, 1.0/eos.Gamma(), velRef/aRef);

  //initialize sutherland's law for viscosity
  sutherland suth(inputVars.TRef());

  vector<plot3dBlock> mesh;
  vector<interblock> connections;
  vector<procBlock> stateBlocks;

  if (rank == ROOT ){

    cout << "Number of equations: " << numEqns << endl << endl;

    //Read grid
    mesh = ReadP3dGrid(inputVars.GridName(), totalCells);

    //split blocks and BCs 
    vector<boundaryConditions> bcs = inputVars.AllBC();
    plot3dBlock blk = mesh[0].Split("j",30);
    cout << "block split" << endl;
    procBlock lowSplit(state, mesh[0], 0, numGhost, inputVars.EquationSet(), bcs[0]);
    vector<procBlock> ls(1,lowSplit);
    WriteCellCenter("lowSplit",ls);
    procBlock upSplit(state, blk, 0, numGhost, inputVars.EquationSet(), bcs[0]);
    vector<procBlock> us(1,upSplit);
    WriteCellCenter("upSplit",us);

    mesh[0].Join(blk, "j");
    cout << "block joined" << endl;
    procBlock joined(state, mesh[0], 0, numGhost, inputVars.EquationSet(), bcs[0]);
    vector<procBlock> splitJoin(1,joined);
    WriteCellCenter("splitJoin",splitJoin);
    exit(0);

    //Get interblock BCs
    connections = GetInterblockBCs( bcs, mesh );

    //Could send proc3dblocks to processors here, or initialize all on ROOT processor


    //initialize the whole mesh with one state and assign ghost cells geometry ------------------
    stateBlocks.resize( mesh.size() );
    for ( int ll = 0; ll < (int)mesh.size(); ll++) {
      stateBlocks[ll] = procBlock(state, mesh[ll], ll, numGhost, inputVars.EquationSet(), bcs[ll]);
      stateBlocks[ll].AssignGhostCellsGeom(inputVars);
    }

    //swap geometry for interblock BCs
    for ( unsigned int ii = 0; ii < connections.size(); ii++ ){
      SwapSlice( connections[ii], stateBlocks[connections[ii].BlockFirst()], stateBlocks[connections[ii].BlockSecond()], true);
    }
    //Get ghost cell edge data
    for ( int ll = 0; ll < (int)mesh.size(); ll++) {
      stateBlocks[ll].AssignGhostCellsGeomEdge(inputVars);
    }

    cout << endl << "Solution Initialized" << endl;
    //----------------------------------------------------------------------------------------------

    //decompose grid
    loadBal = CubicDecomposition(stateBlocks, numProcs, connections, totalCells);

  }

  //set MPI datatypes
  MPI_Datatype MPI_vec3d, MPI_cellData, MPI_procBlockInts, MPI_interblock, MPI_DOUBLE_5INT;
  SetDataTypesMPI(numEqns, MPI_vec3d, MPI_cellData, MPI_procBlockInts, MPI_interblock, MPI_DOUBLE_5INT );

  //send number of procBlocks to all processors
  SendNumProcBlocks( loadBal, rank, numProcBlock);

  //send procBlocks to appropriate processor
  vector<procBlock> localStateBlocks = SendProcBlocks(stateBlocks, rank, numProcBlock, MPI_cellData, MPI_vec3d);

  //send connections to all processors
  SendConnections( connections, MPI_interblock );

  //create operation
  MPI_Op MPI_MAX_LINF;
  MPI_Op_create( (MPI_User_function *) MaxLinf, true, &MPI_MAX_LINF);

  //----------------------------------------------------------------------------------------------

  //determine if implict or explicit
  bool implicitFlag = false;
  if ( inputVars.TimeIntegration() == "implicitEuler" || inputVars.TimeIntegration() == "crankNicholson" || inputVars.TimeIntegration() == "bdf2" ){
    implicitFlag = true;
  }

  //column matrix of zeros
  genArray initial(0.0);

  //preallocate vectors for old solution
  //outermost vector for blocks, inner vector for cell is blocks, genArray for variables in cell
  vector<vector<genArray> > solTimeN(numProcBlock);
  vector<vector<genArray> > solDeltaNm1(numProcBlock);

  //initialize residual variables
  genArray residL2(0.0); //l2 norm residuals
  genArray residL2First(0.0); //l2 norm residuals to normalize by
  resid residLinf; //linf residuals
  double matrixResid = 0.0; //matrix inversion residual (only for implicit runs)

  if ( rank == ROOT) {
    //Write out cell centers grid file
    WriteCellCenter(inputVars.GridName(),stateBlocks);
    //WriteCellCenterGhost(inputVars.GridName(),stateBlocks);
    //Write out initial results
    WriteFun(inputVars.GridName(),stateBlocks, eos, 0, inputVars.RRef(), aRef, inputVars.TRef());
    WriteRes(inputVars.GridName(), 0, inputVars.OutputFrequency());
  }

  for ( int nn = 0; nn < inputVars.Iterations(); nn++ ){            //loop over time

    //calculate cfl number
    inputVars.CalcCFL(nn);

    for ( int mm = 0; mm < inputVars.NonlinearIterations(); mm++ ){    //loop over nonlinear iterations

      //Get boundary conditions for all blocks
      GetBoundaryConditions(localStateBlocks, inputVars, eos, connections, rank, MPI_cellData);

      for ( int bb = 0; bb < (int)localStateBlocks.size(); bb++ ){             //loop over number of blocks

	//calculate inviscid fluxes
	localStateBlocks[bb].CalcInvFluxI(eos, inputVars);
	localStateBlocks[bb].CalcInvFluxJ(eos, inputVars);
	localStateBlocks[bb].CalcInvFluxK(eos, inputVars);

	//if viscous change ghost cells and calculate viscous fluxes
	if (inputVars.EquationSet() == "navierStokes"){

	  //determine ghost cell values for viscous fluxes
	  localStateBlocks[bb].AssignViscousGhostCells(inputVars, eos);

	  //calculate viscous fluxes
	  localStateBlocks[bb].CalcViscFluxI(suth, eos, inputVars);
	  localStateBlocks[bb].CalcViscFluxJ(suth, eos, inputVars);
	  localStateBlocks[bb].CalcViscFluxK(suth, eos, inputVars);
	}

	//calculate the time step to use in the simulation (either user specified or derived from CFL)
	localStateBlocks[bb].CalcBlockTimeStep(inputVars, aRef);

	//if implicit get old solution, reorder block, and use linear solver
	if (implicitFlag){

	  //store time-n solution
	  if (mm == 0){ //first nonlinear iteration, save solution
	    solTimeN[bb] = localStateBlocks[bb].GetCopyConsVars(eos);
	    if (nn == 0){ //at first iteration, resize vector for old solution and calculate solution at time n=-1
	      solDeltaNm1[bb].resize(solTimeN[bb].size());
	      localStateBlocks[bb].DeltaNMinusOne(solDeltaNm1[bb], solTimeN[bb], eos, inputVars.Theta(), inputVars.Zeta());
	    }
	  }

	  //add volume divided by time step term to time m minus time n values
	  vector<genArray> solTimeMmN = localStateBlocks[bb].AddVolTime(localStateBlocks[bb].GetCopyConsVars(eos), solTimeN[bb], 
								    inputVars.Theta(), inputVars.Zeta());

	  //reorder block (by hyperplanes) for lusgs
	  vector<vector3d<int> > reorder = HyperplaneReorder(localStateBlocks[bb].NumI(), localStateBlocks[bb].NumJ(), localStateBlocks[bb].NumK());

	  //reserve space for correction du
	  vector<genArray> du(localStateBlocks[bb].NumCells(),initial);

	  //calculate correction (du)
	  matrixResid += localStateBlocks[bb].LUSGS(reorder, du, solTimeMmN, solDeltaNm1[bb], eos, inputVars, suth );

	  //update solution
	  localStateBlocks[bb].UpdateBlock(inputVars, implicitFlag, eos, aRef, du, residL2, residLinf );

	  //assign time n to time n-1 at end of nonlinear iterations
	  if (inputVars.TimeIntegration() == "bdf2" && mm == inputVars.NonlinearIterations()-1 ){
	    localStateBlocks[bb].DeltaNMinusOne(solDeltaNm1[bb], solTimeN[bb], eos, inputVars.Theta(), inputVars.Zeta());
	  }

	} //conditional for implicit solver

	else{ //explicit
	  //update solution
	  vector<genArray> dummyCorrection(1); //not used in explicit update
	  localStateBlocks[bb].UpdateBlock(inputVars, implicitFlag, eos, aRef, dummyCorrection, residL2, residLinf );
	}

	//zero residuals and wave speed
	localStateBlocks[bb].ResetResidWS();

      } //loop for blocks --------------------------------------------------------------------------------------------------

      //Get residuals from all processors
      residL2.GlobalReduceMPI(rank, numEqns);
      residLinf.GlobalReduceMPI(rank, MPI_DOUBLE_5INT, MPI_MAX_LINF);

      //Get matrix residuals from all processors
      if ( rank == ROOT ){
      	MPI_Reduce(MPI_IN_PLACE, &matrixResid, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
      }
      else{
	MPI_Reduce(&matrixResid, &matrixResid, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
      }


      if ( rank == ROOT ){
	//finish calculation of L2 norm of residual
	for ( int cc = 0; cc < numEqns; cc++ ){
	  residL2[cc] = sqrt(residL2[cc]);
	}
	//finish calculation of matrix residual
	matrixResid = sqrt(matrixResid/(totalCells * numEqns));

	//print out run information
	WriteResiduals(inputVars, residL2First, residL2, residLinf, matrixResid, nn, mm);
      }

      //reset residuals
      residL2.Zero();
      residLinf.Zero();
      matrixResid = 0.0;

    } //loop for nonlinear iterations --------------------------------------------------------------------------------------


    if ( (nn+1)  % inputVars.OutputFrequency() == 0 ){ //write out function file

      //send/recv solutions
      GetProcBlocks( stateBlocks, localStateBlocks, rank, MPI_cellData );

      if ( rank == ROOT ){
	cout << "writing out function file at iteration " << nn << endl;
	//Write out function file
	WriteFun(inputVars.GridName(),stateBlocks, eos, (double) (nn+1), inputVars.RRef(), aRef, inputVars.TRef());
	WriteRes(inputVars.GridName(), (nn+1), inputVars.OutputFrequency());
      }
    }

  } //loop for time ---------------------------------------------------------------------------------------------------------


  if ( rank == ROOT ){
    cout << endl;
    cout << "Program Complete" << endl;
    PrintTime();

    duration = (clock() - start)/(double) CLOCKS_PER_SEC;
    cout << "Total Time: " << duration << " seconds" << endl;
  }

  MPI_Finalize();

  return 0;
}

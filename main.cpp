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
#include <fenv.h>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::clock_t;
using std::clock;

int main( int argc, char *argv[] ) {

  //start clock to time simulation
  clock_t start;
  double duration;
  start = clock();

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); //enable exceptions so code won't run with NANs

  const string inputFile = argv[1];  //name of input file is the second argument (the executable being the first)

  //Parse input file
  double totalCells = 0.0;
  input inputVars = ReadInput(inputFile);

  //Determine number of equations
  int numEqns = 0;
  if ( (inputVars.EquationSet() == "euler") || (inputVars.EquationSet() == "navierStokes") ){
    numEqns = 5;
  }
  else{
    cerr << "ERROR: Equations set is not recognized. Cannot determine number of equations!" << endl;
  }

  cout << "Number of equations: " << numEqns << endl << endl;

  //determine number of ghost cells
  int numGhost = 2;

  //Read grid
  vector<plot3dBlock> mesh = ReadP3dGrid(inputVars.GridName(), totalCells);

  //Get interblock BCs
  vector<interblock> connections = GetInterblockBCs( inputVars.AllBC(), mesh );

  //Initialize state vector with nondimensional variables
  //get reference speed of sound
  idealGas eos(inputVars.Gamma(),inputVars.R());                          //create an equation of state
  double aRef = eos.GetSoS(inputVars.PRef(),inputVars.RRef());

  //get a reference velocity
  vector3d<double> velRef = inputVars.VelRef();

  //initialize a single state
  primVars state(1.0, 1.0/eos.Gamma(), velRef/aRef);

  //initialize sutherland's law for viscosity
  sutherland suth(inputVars.TRef());

  //initialize the whole mesh with one state and assign ghost cells geometry
  vector<procBlock> stateBlocks( mesh.size() );
  for ( int ll = 0; ll < (int)mesh.size(); ll++) {
    stateBlocks[ll] = procBlock(state, mesh[ll], ll, numGhost, inputVars.EquationSet());
    stateBlocks[ll].AssignGhostCellsGeom();
  }

  //Populate interblock boundaries with correct geometry
  for ( unsigned int ii = 0; ii < connections.size(); ii++ ){
    SwapGhostGeom( connections[ii], stateBlocks[connections[ii].BlockFirst()], stateBlocks[connections[ii].BlockSecond()]);
  }
  //Get ghost cell edge data
  for ( int ll = 0; ll < (int)mesh.size(); ll++) {
    stateBlocks[ll].AssignGhostCellsGeomEdge();
  }

  cout << endl << "Solution Initialized" << endl;

  //determine if implict or explicit
  bool implicitFlag = false;
  if ( inputVars.TimeIntegration() == "implicitEuler" || inputVars.TimeIntegration() == "crankNicholson" || inputVars.TimeIntegration() == "bdf2" ){
    implicitFlag = true;
  }

  //column matrix of zeros
  colMatrix initial(numEqns);
  initial.Zero();

  //preallocate vectors for old solution and solution correction
  //outermost vector for blocks, inner vector for cell is blocks, colMatrix for variables in cell
  vector<vector<colMatrix> > solTimeN(mesh.size());
  vector<vector<colMatrix> > solDeltaNm1(mesh.size());
  vector<vector<colMatrix> > du(mesh.size());

  //initialize residual variables
  int locMaxB = 0; //block with max residual
  colMatrix residL2 = initial; //l2 norm residuals
  colMatrix residL2First = initial; //l2 norm residuals to normalize by
  colMatrix residLinf = initial; //linf residuals
  double matrixResid = 0.0; //matrix inversion residual (only for implicit runs)

  //Write out cell centers grid file
  WriteCellCenter(inputVars.GridName(),stateBlocks);
  //WriteCellCenterGhost(inputVars.GridName(),stateBlocks);
  //Write out initial results
  WriteFun(inputVars.GridName(),stateBlocks, eos, 0, inputVars.RRef(), aRef, inputVars.TRef());
  WriteRes(inputVars.GridName(), 0, inputVars.OutputFrequency());

  for ( int nn = 0; nn < inputVars.Iterations(); nn++ ){            //loop over time

    //calculate cfl number
    inputVars.CalcCFL(nn);

    for ( int mm = 0; mm < inputVars.NonlinearIterations(); mm++ ){    //loop over nonlinear iterations

      //Get boundary conditions
      GetBoundaryConditions(stateBlocks, inputVars, eos, connections);

      for ( int bb = 0; bb < (int)mesh.size(); bb++ ){             //loop over number of blocks

	//initialize implicit matrix
	if (implicitFlag){

	  //reserve space for correction du
	  //only need to do this if it is first iteration
	  if (nn == 0){
	    du[bb].resize(stateBlocks[bb].NumCells(), initial);
	  }

	} //end of implicit conditional

	//calculate inviscid fluxes
	stateBlocks[bb].CalcInvFluxI(eos, inputVars);
	stateBlocks[bb].CalcInvFluxJ(eos, inputVars);
	stateBlocks[bb].CalcInvFluxK(eos, inputVars);

	//if viscous change ghost cells and calculate viscous fluxes
	if (inputVars.EquationSet() == "navierStokes"){

	  //determine ghost cell values for viscous fluxes
	  stateBlocks[bb].AssignViscousGhostCells(inputVars, eos);

	  //calculate viscous fluxes
	  stateBlocks[bb].CalcViscFluxI(suth, eos, inputVars);
	  stateBlocks[bb].CalcViscFluxJ(suth, eos, inputVars);
	  stateBlocks[bb].CalcViscFluxK(suth, eos, inputVars);
	}

	//calculate the time step to use in the simulation (either user specified or derived from CFL)
	stateBlocks[bb].CalcBlockTimeStep(inputVars, aRef);

	//if implicit get old solution, reorder block, and use linear solver
	if (implicitFlag){

	  //store time-n solution
	  if (mm == 0){ //first nonlinear iteration, save solution
	    solTimeN[bb] = stateBlocks[bb].GetCopyConsVars(eos);
	    if (nn == 0){ //at first iteration, resize vector for old solution and calculate solution at time n=-1
	      solDeltaNm1[bb].resize(solTimeN[bb].size());
	      stateBlocks[bb].DeltaNMinusOne(solDeltaNm1[bb], solTimeN[bb], eos, inputVars.Theta(), inputVars.Zeta());
	    }
	  }

	  //add volume divided by time step term to time m minus time n values
	  vector<colMatrix> solTimeMmN = stateBlocks[bb].AddVolTime(stateBlocks[bb].GetCopyConsVars(eos), solTimeN[bb], inputVars.Theta(), inputVars.Zeta());

	  //reorder block (by hyperplanes) for lusgs
	  vector<vector3d<int> > reorder = HyperplaneReorder(stateBlocks[bb].NumI(), stateBlocks[bb].NumJ(), stateBlocks[bb].NumK());

	  //calculate correction (du)
	  matrixResid += stateBlocks[bb].LUSGS(reorder, du[bb], solTimeMmN, solDeltaNm1[bb], eos, inputVars, suth );

	} //conditional for implicit solver

      } //loop for blocks

      //after update for all blocks has been calculated and stored, update all blocks
      //blocks cannot be updated within block loop because ghost cells for connection boundaries are determined from adjacent blocks
      //this would result in the ghost cell being at time n+1 when it should be at time n
      for ( int dd = 0; dd < (int)mesh.size(); dd++ ){             //loop over number of blocks

	//update solution
	stateBlocks[dd].UpdateBlock(inputVars, implicitFlag, eos, aRef, du[dd], residL2, residLinf, locMaxB);

	//if implicit, assign time n to time n-1 at end of nonlinear iterations
	if (implicitFlag && inputVars.TimeIntegration() == "bdf2" && mm == inputVars.NonlinearIterations()-1 ){
	  stateBlocks[dd].DeltaNMinusOne(solDeltaNm1[dd], solTimeN[dd], eos, inputVars.Theta(), inputVars.Zeta());
	}

	//zero residuals and wave speed
	stateBlocks[dd].ResetResidWS();

      } //loop for blocks


      //finish calculation of L2 norm of residual
      for ( int cc = 0; cc < residL2.Size(); cc++ ){
	residL2.SetData(cc, sqrt(residL2.Data(cc)) );
      }
      //finish calculation of matrix residual
      matrixResid = sqrt(matrixResid/(totalCells * numEqns));

      //print out run information
      WriteResiduals(inputVars, residL2First, residL2, residLinf, matrixResid, locMaxB, nn, mm);

      //reset residuals
      residL2 = initial;
      residLinf = initial;
      locMaxB = 0;
      matrixResid = 0.0;

    } //loop for nonlinear iterations

    if ( (nn+1)  % inputVars.OutputFrequency() == 0 ){ //write out function file
      cout << "writing out function file at iteration " << nn << endl;
      //Write out function file
      WriteFun(inputVars.GridName(),stateBlocks, eos, (double) (nn+1), inputVars.RRef(), aRef, inputVars.TRef());
      WriteRes(inputVars.GridName(), (nn+1), inputVars.OutputFrequency());
    }


  } //loop for time


  cout << endl;
  cout << "Program Complete" << endl;
  PrintTime();

  duration = (clock() - start)/(double) CLOCKS_PER_SEC;
  cout << "Total Time: " << duration << " seconds" << endl;

  return 0;
}

#include <iostream>      //cout, cerr, endl
#include "plot3d.h"
#include "vector3d.h"
#include "tensor.h"
#include "input.h"
#include <vector>        //stl vector
#include <algorithm>     //max_element
#include <numeric>       //accumulate
#include "blockVars.h"
#include "inviscidFlux.h"
#include "viscousFlux.h"
#include "viscBlockVars.h"
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

  clock_t start;
  double duration;
  start = clock();

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  const string inputFile = argv[1];  //name of input file is the second argument (the executable being the first)

  const double eps = 1.0e-30;

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

  //Read grid
  plot3dMesh mesh = ReadP3dGrid(inputVars.GridName(), totalCells);

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

  //initialize the whole mesh
  vector<blockVars> stateBlocks( mesh.NumBlocks() );
  vector<viscBlockVars> viscBlocks( mesh.NumBlocks() );
  int ll = 0;
  for ( ll = 0; ll < mesh.NumBlocks(); ll++) {
    stateBlocks[ll] = blockVars(state, mesh.Blocks(ll));
    if (inputVars.EquationSet() == "navierStokes"){
      viscBlocks[ll] = viscBlockVars(mesh.Blocks(ll));
    }
    else if (inputVars.EquationSet() == "euler"){
      //do nothing extra
    }
    else{
      cerr << "ERROR: Equation set " << inputVars.EquationSet() << " is not recognized!" << endl;
      exit(0);
    }

  }

  cout << endl << "Solution Initialized" << endl;

  //determine if implict or explicit
  bool implicitFlag = false;
  if ( inputVars.TimeIntegration() == "implicitEuler" || inputVars.TimeIntegration() == "crankNicholson" || inputVars.TimeIntegration() == "bdf2" ){
    implicitFlag = true;
  }

  colMatrix initial(numEqns);
  initial.Zero();

  vector<vector<colMatrix> > solTimeN(mesh.NumBlocks());
  vector<vector<colMatrix> > solDeltaNm1(mesh.NumBlocks());
  vector<vector<colMatrix> > du(mesh.NumBlocks());

  int bb = 0;
  unsigned int cc = 0;
  int nn = 0;
  int mm = 0;

  int locMaxB = 0;

  vector<double> residL2(numEqns, 0.0);
  vector<double> residL2First(numEqns, 0.0);
  vector<double> residLinf(numEqns, 0.0);
  double matrixResid = 0.0;

  //Write out cell centers grid file
  WriteCellCenter(inputVars.GridName(),stateBlocks);

  for ( nn = 0; nn < inputVars.Iterations(); nn++ ){            //loop over time

    //calculate cfl number
    inputVars.CalcCFL(nn);

    for ( mm = 0; mm < inputVars.NonlinearIterations(); mm++ ){    //loop over nonlinear iterations

      for ( bb = 0; bb < mesh.NumBlocks(); bb++ ){             //loop over number of blocks

	int numElems = (mesh.Blocks(bb).NumI() - 1) * (mesh.Blocks(bb).NumJ() - 1) * (mesh.Blocks(bb).NumK() - 1);

	//initialize implicit matrix
	if (implicitFlag){

	  //reserve space for correction du
	  //only need to do this if it is first iteration
	  if (nn == 0){
	    du[bb].resize(numElems,initial);
	  }

	} //end of implicit conditional

	//calculate inviscid fluxes
	stateBlocks[bb].CalcInvFluxI(eos, inputVars, bb);
	stateBlocks[bb].CalcInvFluxJ(eos, inputVars, bb);
	stateBlocks[bb].CalcInvFluxK(eos, inputVars, bb);

	//if viscous, calculate gradients and viscous fluxes
	if (inputVars.EquationSet() == "navierStokes"){
	  viscBlocks[bb].InitializeGrads(stateBlocks[bb]);

	  viscBlocks[bb].CalcCellGradsI(stateBlocks[bb], eos, suth, inputVars, bb);
	  viscBlocks[bb].CalcCellGradsJ(stateBlocks[bb], eos, suth, inputVars, bb);
	  viscBlocks[bb].CalcCellGradsK(stateBlocks[bb], eos, suth, inputVars, bb);

	  viscBlocks[bb].CalcViscFluxI(stateBlocks[bb], suth, eos, inputVars, bb);
	  viscBlocks[bb].CalcViscFluxJ(stateBlocks[bb], suth, eos, inputVars, bb);
	  viscBlocks[bb].CalcViscFluxK(stateBlocks[bb], suth, eos, inputVars, bb);
	}

	//calculate cell time step
	// if (inputVars.EquationSet() == "navierStokes"){
	//   viscBlocks[bb].CalcBlockTimeStep(stateBlocks[bb], inputVars, aRef);
	// }
	// else{
	  stateBlocks[bb].CalcBlockTimeStep(inputVars, aRef);
	// }


	//if implicit calculate flux jacobians and assembly matrix
	if (implicitFlag){

	  //store time-n solution
	  if (mm == 0){
	    solTimeN[bb] = stateBlocks[bb].GetCopyConsVars(eos);
	    if (nn == 0){
	      solDeltaNm1[bb].resize(solTimeN[bb].size());
	      stateBlocks[bb].DeltaNMinusOne(solDeltaNm1[bb], solTimeN[bb], eos, inputVars.Theta(), inputVars.Zeta());
	    }
	  }

	  // stateBlocks[bb].CalcInvFluxJacI( eos, inputVars, bb, mainDiag);
	  // stateBlocks[bb].CalcInvFluxJacJ( eos, inputVars, bb, mainDiag);
	  // stateBlocks[bb].CalcInvFluxJacK( eos, inputVars, bb, mainDiag);


	  // if (inputVars.EquationSet() == "navierStokes" ){
	  //   CalcViscFluxJacI(stateBlocks[bb], suth, eos, inputVars, bb, mainDiag);
	  //   CalcViscFluxJacJ(stateBlocks[bb], suth, eos, inputVars, bb, mainDiag);
	  //   CalcViscFluxJacK(stateBlocks[bb], suth, eos, inputVars, bb, mainDiag);
	  // }

	  //add volume divided by time step term to main diagonal
	  //stateBlocks[bb].AddVolTime(mainDiag, inputVars.Theta(), inputVars.Zeta(), inputVars.DualTimeCFL());

	  //add volume divided by time step term time m - time n term
	  vector<colMatrix> solTimeMmN = stateBlocks[bb].AddVolTime(stateBlocks[bb].GetCopyConsVars(eos), solTimeN[bb], inputVars.Theta(), inputVars.Zeta());

	  //reorder block for lusgs
	  vector<vector3d<int> > reorder = HyperplaneReorder(stateBlocks[bb].NumI()-1, stateBlocks[bb].NumJ()-1, stateBlocks[bb].NumK()-1);

	  //calculate correction (du)
	  matrixResid += stateBlocks[bb].LUSGS(reorder, du[bb], solTimeMmN, solDeltaNm1[bb], eos, inputVars, suth );


	} //conditional for implicit solver

      } //loop for blocks

      //after update for all blocks has been calculated and stored, update all blocks
      //blocks cannot be updated within block loop because ghost cells for connection boundaries are determined from adjacent blocks
      //this would result in the ghost cell being at time n+1 when it should be at time n
      for ( int dd = 0; dd < mesh.NumBlocks(); dd++ ){             //loop over number of blocks

	//update solution
	stateBlocks[dd].UpdateBlock(inputVars, implicitFlag, eos, aRef, dd, du[dd], residL2, residLinf, locMaxB);

	//if implicit, assign time n to time n-1 at end of nonlinear iterations
	if (implicitFlag && inputVars.TimeIntegration() == "bdf2" && mm == inputVars.NonlinearIterations()-1 ){
	  stateBlocks[dd].DeltaNMinusOne(solDeltaNm1[dd], solTimeN[dd], eos, inputVars.Theta(), inputVars.Zeta());
	}

	//zero residuals and wave speed
	stateBlocks[dd].ResetResidWS();

      } //loop for blocks


      //finish calculation of L2 norm of residual
      for ( cc = 0; cc < residL2.size(); cc++ ){
	residL2[cc] = sqrt(residL2[cc]);

	if (nn == 0 && mm == 0){
	  residL2First[cc] = residL2[cc];
	  cout << residL2[cc] << endl;
	}

	//normalize residuals
	residL2[cc] = (residL2[cc]+eps) / (residL2First[cc]+eps) ;
      }

      matrixResid = sqrt(matrixResid/(totalCells * numEqns));


      //print out run information
      if (nn%100 == 0 && mm == 0){  
	if (inputVars.Dt() > 0.0){
	  cout << "STEP    NONLINEAR     DT     RES-Mass     Res-Mom-X     Res-Mom-Y     Res-Mom-Z     Res-Energy    Max Res Eqn    Max Res Blk    Max Res I    Max Res J    Max Res K    Max Res    Res-Matrix" << endl;
	}
	else if (inputVars.CFL() > 0.0){
	  cout << "STEP    NONLINEAR     CFL     RES-Mass     Res-Mom-X     Res-Mom-Y     Res-Mom-Z     Res-Energy   Max Res Eqn    Max Res Blk    Max Res I    Max Res J    Max Res K    Max Res    Res-Matrix" << endl;
	}

      }
      if (inputVars.Dt() > 0.0){
	cout << nn << "     " << mm << "     " << inputVars.Dt() << "     " << residL2[0] <<  "     " << residL2[1] << "     " << residL2[2] << "     " << residL2[3] << "     " << residL2[4] << "     " 
	     << residLinf[3] << "     " << locMaxB << "     " << residLinf[0] <<"     " << residLinf[1] << "     " << residLinf[2] << "     " << residLinf[4] << "     " << matrixResid << endl;
      }
      else if (inputVars.CFL() > 0.0){
	cout << nn << "     " << mm << "     " << inputVars.CFL() << "     " << residL2[0] <<  "     " << residL2[1] << "     " << residL2[2] << "     " << residL2[3] << "     " << residL2[4] << "     " 
	     << residLinf[3] << "     " << locMaxB << "     " << residLinf[0] <<"     " << residLinf[1] << "     " << residLinf[2] << "     " << residLinf[4] << "     " << matrixResid << endl;
      }

      //reset residuals
      for ( cc = 0; cc < residL2.size(); cc++ ){
	residL2[cc] = 0.0;
	residLinf[cc] = 0.0;
      }
      locMaxB = 0;
      matrixResid = 0.0;

    } //loop for nonlinear iterations

    if ( (nn+1)  % inputVars.OutputFrequency() == 0 ){ //write out function file
      cout << "write out function file at iteration " << nn << endl;
      //Write out function file
      WriteFun(inputVars.GridName(),stateBlocks, viscBlocks, eos, (double) (nn+1), inputVars.RRef(), aRef, inputVars.TRef());
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

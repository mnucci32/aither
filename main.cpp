/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <iostream>      // cout, cerr, endl
#include <vector>        // stl vector
#include <ctime>         // clock
#include <string>        // stl string
#include <memory>        // unique_ptr

#ifdef __linux__
#include <cfenv>         // exceptions
#elif __APPLE__
#include <xmmintrin.h>
#endif

#include "plot3d.hpp"
#include "vector3d.hpp"
#include "input.hpp"
#include "procBlock.hpp"
#include "primVars.hpp"
#include "eos.hpp"
#include "boundaryConditions.hpp"
#include "output.hpp"
#include "genArray.hpp"
#include "parallel.hpp"
#include "turbulence.hpp"
#include "gradients.hpp"
#include "resid.hpp"
#include "multiArray3d.hpp"
#include "kdtree.hpp"
#include "fluxJacobian.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::clock_t;
using std::clock;

int main(int argc, char *argv[]) {
  // Initialize MPI and make calls to get number
  // of processors and rank of each processor
  auto numProcs = 1;
  auto rank = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get MPI version
  auto version = 3;
  auto subversion = 0;
  MPI_Get_version(&version, &subversion);
  if ( rank == ROOTP ) {
    cout << "Aither version " << MAJORVERSION << "." << MINORVERSION << "."
         << PATCHNUMBER << endl;
    cout << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
    cout << "Using MPI Version " << version << "." << subversion << endl;
    cout << "Using " << numProcs << " processors" << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Start clock to time simulation
  auto start = clock();

  // Enable exceptions so code won't run with NANs
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#elif __APPLE__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);    
#endif
  
  // Name of input file is the second argument (the executable being the first)
  string inputFile = argv[1];

  // Broadcast input file name to all names for portability
  BroadcastString(inputFile);

  auto totalCells = 0.0;
  input inputVars(inputFile);
  decomposition decomp;
  auto numProcBlock = 0;

  // Parse input file
  inputVars.ReadInput(rank);

  // Determine number of ghost cells
  auto numGhost = 2;

  // Get equation of state
  idealGas eos(inputVars.Gamma(), inputVars.R());

  // Initialize sutherland's law for viscosity
  sutherland suth(inputVars.TRef(), inputVars.RRef(), inputVars.LRef(),
                  inputVars.PRef(), inputVars.VelRef(), eos);

  // Initialize state vector with nondimensional variables
  // Get reference speed of sound
  auto aRef = eos.SoS(inputVars.PRef(), inputVars.RRef());
  primVars state;
  state.NondimensionalInitialize(eos, aRef, inputVars, suth);

  // Get turbulence model
  auto turb = inputVars.AssignTurbulenceModel();

  vector<plot3dBlock> mesh;
  vector<interblock> connections;
  vector<procBlock> stateBlocks;
  vector<vector3d<double>> viscFaces;

  if (rank == ROOTP) {
    cout << "Number of equations: " << inputVars.NumEquations() << endl << endl;

    // Read grid
    mesh = ReadP3dGrid(inputVars.GridName(), inputVars.LRef(), totalCells);

    // Get BCs for blocks
    auto bcs = inputVars.AllBC();

    // Decompose grid
    if (inputVars.DecompMethod() == "manual") {
      decomp = ManualDecomposition(mesh, bcs, numProcs);
    } else if (inputVars.DecompMethod() == "cubic") {
      decomp = CubicDecomposition(mesh, bcs, numProcs);
    } else {
      cerr << "ERROR: Domain decomposition method " << inputVars.DecompMethod()
           << " is not recognized!" << endl;
      exit(0);
    }

    // Get interblock BCs
    connections = GetInterblockBCs(bcs, mesh, decomp);

    // Initialize the whole mesh with one state and assign ghost cells geometry
    stateBlocks.resize(mesh.size());
    for (auto ll = 0; ll < static_cast<int>(mesh.size()); ll++) {
      stateBlocks[ll] = procBlock(state, mesh[ll], decomp.ParentBlock(ll),
                                  numGhost, bcs[ll], ll, decomp.Rank(ll),
                                  decomp.LocalPosition(ll) );
      stateBlocks[ll].AssignGhostCellsGeom();
    }

    // Swap geometry for interblock BCs
    for (auto ii = 0; ii < static_cast<int>(connections.size()); ii++) {
      SwapSlice(connections[ii], stateBlocks[connections[ii].BlockFirst()],
                stateBlocks[connections[ii].BlockSecond()], true);
    }
    // Get ghost cell edge data
    for (auto ll = 0; ll < static_cast<int>(mesh.size()); ll++) {
      stateBlocks[ll].AssignGhostCellsGeomEdge();
    }

    // Get face centers of faces with viscous wall BC
    viscFaces = GetViscousFaceCenters(stateBlocks);

    cout << "Solution Initialized" << endl << endl;
    //---------------------------------------------------------------------
  }

  // Set MPI datatypes
  MPI_Datatype MPI_vec3d, MPI_cellData, MPI_procBlockInts,
      MPI_interblock, MPI_DOUBLE_5INT, MPI_vec3dMag;
  SetDataTypesMPI(MPI_vec3d, MPI_cellData, MPI_procBlockInts,
                  MPI_interblock, MPI_DOUBLE_5INT, MPI_vec3dMag);

  // Send number of procBlocks to all processors
  SendNumProcBlocks(decomp.NumBlocksOnAllProc(), numProcBlock);

  // Send procBlocks to appropriate processor
  auto localStateBlocks = SendProcBlocks(stateBlocks, rank, numProcBlock,
                                         MPI_cellData, MPI_vec3d, MPI_vec3dMag);

  // Send connections to all processors
  SendConnections(connections, MPI_interblock);

  // Broadcast viscous face centers to all processors
  BroadcastViscFaces(MPI_vec3d, viscFaces);

  // Create operation
  MPI_Op MPI_MAX_LINF;
  MPI_Op_create(reinterpret_cast<MPI_User_function *> (MaxLinf), true,
                &MPI_MAX_LINF);

  //-----------------------------------------------------------------------
  // wall distance calculation

  auto wallStart = clock();

  if (rank == ROOTP) {
    cout << "Starting wall distance calculation..." << endl;
    cout << "Building k-d tree..." << endl;
  }

  // Construct k-d tree for wall distance calculation
  kdtree tree(viscFaces);

  if (rank == ROOTP) {
    auto kdDuration = (clock() - wallStart) /
        static_cast<double> (CLOCKS_PER_SEC);
    cout << "K-d tree complete after " << kdDuration << " seconds" << endl;
  }

  if (tree.Size() > 0) {
    CalcWallDistance(localStateBlocks, tree);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == ROOTP) {
    auto wallDuration = (clock() - wallStart) /
        static_cast<double> (CLOCKS_PER_SEC);
    cout << "Wall distance calculation finished after " << wallDuration
         << " seconds" << endl << endl;
  }

  //-----------------------------------------------------------------------
  // Preallocate vectors for old solution
  // Outermost vector for blocks, genArray for variables in cell
  vector<multiArray3d<genArray>> solTimeN(numProcBlock);
  vector<multiArray3d<genArray>> solDeltaNm1(numProcBlock);
  vector<multiArray3d<genArray>> solDeltaMmN(numProcBlock);
  vector<multiArray3d<genArray>> du(numProcBlock);
  // Allocate array for flux jacobian
  vector<multiArray3d<fluxJacobian>> mainDiagonal(numProcBlock);

  genArray residL2First(0.0);  // l2 norm residuals to normalize by

  // Send/recv solutions - necessary to get wall distances
  GetProcBlocks(stateBlocks, localStateBlocks, rank, MPI_cellData);

  if (rank == ROOTP) {
    // Write out cell centers grid file
    WriteCellCenter(inputVars.GridName(), stateBlocks, decomp,
                    inputVars.LRef());

    // Write out initial results
    WriteFun(stateBlocks, eos, suth, 0, decomp, inputVars, turb);
    WriteRes(inputVars.SimNameRoot(), 0, inputVars.OutputFrequency());
  }

  // ----------------------------------------------------------------------
  // ----------------------- Start Main Loop ------------------------------
  // ----------------------------------------------------------------------
  for (auto nn = 0; nn < inputVars.Iterations(); nn++) {   // loop over time
    // Calculate cfl number
    inputVars.CalcCFL(nn);

    // Store time-n solution, for time integration methods that require it
    if (inputVars.IsImplicit() || inputVars.TimeIntegration() == "rk4") {
      solTimeN = GetCopyConsVars(localStateBlocks, eos);
    }

    // loop over nonlinear iterations
    for (auto mm = 0; mm < inputVars.NonlinearIterations(); mm++) {
      // Get boundary conditions for all blocks
      GetBoundaryConditions(localStateBlocks, inputVars, eos, suth, turb,
                            connections, rank, MPI_cellData);

      // Loop over number of blocks
      for (auto bb = 0; bb < numProcBlock; bb++) {
        // Get solution at time M - N if implicit
        if (inputVars.IsImplicit()) {
          solDeltaMmN[bb] = localStateBlocks[bb].SolTimeMMinusN(
              solTimeN[bb], eos, inputVars, mm);
          if (nn == 0 && mm == 0) {
            // At first iteration, resize array for old solution
            solDeltaNm1[bb].ClearResize(localStateBlocks[bb].NumI(),
                                        localStateBlocks[bb].NumJ(),
                                        localStateBlocks[bb].NumK());
            mainDiagonal[bb].ClearResize(localStateBlocks[bb].NumI(),
                                         localStateBlocks[bb].NumJ(),
                                         localStateBlocks[bb].NumK());
            du[bb].ClearResize(localStateBlocks[bb].NumIG(),
                               localStateBlocks[bb].NumJG(),
                               localStateBlocks[bb].NumKG(), genArray(0.0));
          }
        }

        // Calculate residual (RHS)
        localStateBlocks[bb].CalcResidual(suth, eos, inputVars, turb,
                                          mainDiagonal[bb]);

        // Calculate the time step to use in the simulation
        // (either user specified or derived from CFL)
        localStateBlocks[bb].CalcBlockTimeStep(inputVars, aRef);
      }

      // Initialize residual variables
      genArray residL2(0.0);  // l2 norm residuals
      resid residLinf;  // linf residuals
      auto matrixResid = 0.0;
      if (inputVars.IsImplicit()) {
        matrixResid = ImplicitUpdate(localStateBlocks, mainDiagonal, du,
                                     inputVars, eos, aRef, suth, solTimeN,
                                     solDeltaMmN, solDeltaNm1, turb, mm,
                                     residL2, residLinf);
      } else {  // explicit time integration
        ExplicitUpdate(localStateBlocks, inputVars, eos, aRef, suth, solTimeN,
                       turb, mm, residL2, residLinf);
      }

      // ----------------------------------------------------------------------
      // Get residuals from all processors
      residL2.GlobalReduceMPI(rank, inputVars.NumEquations());
      residLinf.GlobalReduceMPI(rank, MPI_DOUBLE_5INT, MPI_MAX_LINF);

      // Get matrix residuals from all processors
      if (rank == ROOTP) {
        MPI_Reduce(MPI_IN_PLACE, &matrixResid, 1, MPI_DOUBLE, MPI_SUM,
                   ROOTP, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(&matrixResid, &matrixResid, 1, MPI_DOUBLE, MPI_SUM,
                   ROOTP, MPI_COMM_WORLD);
      }

      if (rank == ROOTP) {
        // Finish calculation of L2 norm of residual
	residL2.SquareRoot();

        // Finish calculation of matrix residual
        matrixResid = sqrt(matrixResid/(totalCells * inputVars.NumEquations()));

        // Print out run information
        WriteResiduals(inputVars, residL2First, residL2, residLinf, matrixResid,
                       nn, mm);
      }
    }  // loop for nonlinear iterations ---------------------------------------

    if ((nn+1) % inputVars.OutputFrequency() == 0) {  // write out function file
      // Send/recv solutions
      GetProcBlocks(stateBlocks, localStateBlocks, rank, MPI_cellData);

      if (rank == ROOTP) {
        cout << "writing out function file at iteration " << nn << endl;
        // Write out function file
        WriteFun(stateBlocks, eos, suth, (nn+1), decomp, inputVars, turb);
        WriteRes(inputVars.SimNameRoot(), (nn+1), inputVars.OutputFrequency());
      }
    }
  }  // loop for time step -----------------------------------------------------

  if (rank == ROOTP) {
    cout << endl << "Program Complete" << endl;
    PrintTime();

    auto duration = (clock() - start) / static_cast<double> (CLOCKS_PER_SEC);
    cout << "Total Time: " << duration << " seconds" << endl;
  }

  // Free datatypes previously created
  FreeDataTypesMPI(MPI_vec3d, MPI_cellData, MPI_procBlockInts,
                   MPI_interblock, MPI_DOUBLE_5INT, MPI_vec3dMag);

  MPI_Finalize();

  return 0;
}

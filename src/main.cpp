/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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
#include <chrono>        // clock
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
#include "gridLevel.hpp"
#include "primitive.hpp"
#include "physicsModels.hpp"
#include "thermodynamic.hpp"
#include "boundaryConditions.hpp"
#include "output.hpp"
#include "varArray.hpp"
#include "parallel.hpp"
#include "resid.hpp"
#include "multiArray3d.hpp"
#include "kdtree.hpp"
#include "fluxJacobian.hpp"
#include "utility.hpp"
#include "matMultiArray3d.hpp"
#include "mgSolution.hpp"
#include "logFileManager.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

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

  // Enable exceptions so code won't run with NANs
#ifdef __linux__
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
#elif __APPLE__
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#endif

  // Check command line inputs
  // Name of input file is the second argument (the executable being the first)
  if (rank == ROOTP && !(argc == 2 || argc == 3)) {
    cerr << "USAGE: <mpirun -np n> aither inputFile.inp <restartFile.rst>" << endl;
    cerr << "       Arguments in <> are optional." << endl;
    cerr << "       If not invoked with mpirun, 1 processor will be used." << endl;
    cerr << "       If no restart file specified, none will be used." << endl;
    exit(EXIT_FAILURE);
  }
  string inputFile = argv[1];
  string restartFile = (argc == 3) ? argv[2] : "none";

  // Broadcast input/restart file names to all processors for portability
  BroadcastString(inputFile);
  BroadcastString(restartFile);

  auto totalCells = 0.0;
  input inp(inputFile, restartFile);
  decomposition decomp;
  auto numProcBlock = 0;

  // Parse input file
  inp.ReadInput(rank);
  logFileManager logs(inp, rank);

  // nondimensionalize fluid data
  inp.NondimensionalizeFluid();

  // Get physics models
  const auto phys = inp.AssignPhysicsModels();

  // Nondimensionalize BC & IC data
  inp.NondimensionalizeStateData(phys.EoS());

  mgSolution solution;  // only keep finest grid level globally
  vector<vector3d<double>> viscFaces;

  if (rank == ROOTP) {
    cout << "Number of equations: " << inp.NumEquations() << endl << endl;

    // Read grid
    auto mesh = ReadP3dGrid(inp.GridName(), inp.LRef(), totalCells);
    // Get BCs for blocks
    auto bcs = inp.AllBC();

    // Decompose grid
    if (inp.DecompMethod() == "manual") {
      decomp = ManualDecomposition(mesh, bcs, numProcs);
    } else if (inp.DecompMethod() == "cubic") {
      decomp = CubicDecomposition(mesh, bcs, numProcs);
    } else {
      cerr << "ERROR: Domain decomposition method " << inp.DecompMethod()
           << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }

    solution.ConstructFinestLevel(mesh, bcs, decomp, phys, restartFile, inp,
                                  logs.L2First());

    // Get face centers of faces with viscous wall BC
    viscFaces = GetViscousFaceCenters(solution.Finest().Blocks());

    cout << "Solution Initialized" << endl << endl;
    //---------------------------------------------------------------------
  }

  // Set MPI datatypes
  MPI_Datatype MPI_vec3d, MPI_procBlockInts, MPI_connection, MPI_DOUBLE_5INT,
      MPI_vec3dMag, MPI_uncoupledScalar, MPI_tensorDouble;
  SetDataTypesMPI(MPI_vec3d, MPI_procBlockInts, MPI_connection, MPI_DOUBLE_5INT,
                  MPI_vec3dMag, MPI_uncoupledScalar, MPI_tensorDouble);

  // Broadcast decomposition to all processors
  decomp.Broadcast();

  // Send number of procBlocks to all processors
  SendNumProcBlocks(decomp.NumBlocksOnAllProc(), numProcBlock);

  // Send finest gridLevel to appropriate processor
  auto localSolution = solution.SendFinestGridLevel(
      rank, numProcBlock, MPI_vec3d, MPI_vec3dMag, MPI_connection, inp);
  localSolution.ConstructMultigrids(decomp, inp, phys, rank, MPI_connection,
                                    MPI_vec3d, MPI_vec3dMag);

  // Update auxillary variables (temperature, viscosity, etc), cell widths
  localSolution.AuxillaryAndWidths(phys);

  // Broadcast viscous face centers to all processors
  BroadcastViscFaces(MPI_vec3d, viscFaces);

  // Create operation
  MPI_Op MPI_MAX_LINF;
  MPI_Op_create(reinterpret_cast<MPI_User_function *> (MaxLinf), true,
                &MPI_MAX_LINF);

  //-----------------------------------------------------------------------
  // wall distance calculation

  const auto wallStart = std::chrono::high_resolution_clock::now();

  if (rank == ROOTP) {
    cout << "Starting wall distance calculation..." << endl;
    cout << "Building k-d tree..." << endl;
  }

  // Construct k-d tree for wall distance calculation
  // Using finest grid level faces for all levels
  kdtree tree(viscFaces);

  if (rank == ROOTP) {
    const auto kdEnd = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> kdDuration = kdEnd - wallStart;
    cout << "K-d tree complete after " << kdDuration.count() << " seconds"
         << endl;
  }

  if (tree.Size() > 0) {
    localSolution.CalcWallDistance(tree);
    localSolution.SwapWallDist(rank, inp.NumberGhostLayers());
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == ROOTP) {
    const auto wallEnd = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> wallDuration = wallEnd - wallStart;
    cout << "Wall distance calculation finished after " << wallDuration.count()
         << " seconds" << endl << endl;
  }

  //-----------------------------------------------------------------------
  // Send/recv solutions - necessary to get wall distances
  solution.GetFinestGridLevel(localSolution, rank, MPI_uncoupledScalar,
                              MPI_vec3d, MPI_tensorDouble, inp);

  if (rank == ROOTP) {
    // Write out cell centers grid file
    WriteCellCenter(inp.GridName(), solution.Finest().Blocks(), decomp, inp);

    // Write out initial results
    WriteOutput(solution.Finest().Blocks(), phys, inp.IterationStart(), decomp,
                inp);
  }

  // ----------------------------------------------------------------------
  // ----------------------- Start Main Loop ------------------------------
  // ----------------------------------------------------------------------
  // loop over time
  for (auto nn = 0; nn < inp.Iterations(); ++nn) {
    MPI_Barrier(MPI_COMM_WORLD);
    logs.GetIterStart();

    // Calculate cfl number
    inp.CalcCFL(nn);

    // Store time-n solution, for time integration methods that require it
    localSolution.StoreOldSolution(inp, phys, nn);

    // loop over nonlinear iterations
    for (auto mm = 0; mm < inp.NonlinearIterations(); ++mm) {
      // Initialize residual variables
      // l2 norm residuals
      residual residL2(inp.NumEquations(), inp.NumSpecies());
      resid residLinf;  // linf residuals

      // advance an iteration
      auto matrixResid = localSolution.Iterate(
          inp, phys, MPI_tensorDouble, MPI_vec3d, mm, rank, residL2, residLinf);

      // ----------------------------------------------------------------------
      // Get residuals from all processors
      residL2.GlobalReduceMPI(rank);
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
        matrixResid = sqrt(matrixResid/(totalCells * inp.NumEquations()));

        // Print out run information
        logs.WriteResiduals(inp, residL2, residLinf, matrixResid,
                            nn + inp.IterationStart(), mm);
      }
    }  // loop for nonlinear iterations ---------------------------------------

    // write out function file
    if (inp.WriteOutput(nn) || inp.WriteRestart(nn)) {
      // Send/recv solutions
      solution.GetFinestGridLevel(localSolution, rank, MPI_uncoupledScalar,
                                  MPI_vec3d, MPI_tensorDouble, inp);

      if (rank == ROOTP && inp.WriteOutput(nn)) {
        cout << "writing out function file at iteration "
             << nn + inp.IterationStart()<< endl;
        // Write out function file
        WriteOutput(solution.Finest().Blocks(), phys,
                    (nn + inp.IterationStart() + 1), decomp, inp);
      }
      if (rank == ROOTP && inp.WriteRestart(nn)) {
        cout << "writing out restart file at iteration "
             << nn + inp.IterationStart()<< endl;
        // Write out restart file
        WriteRestart(solution.Finest().Blocks(), phys,
                     (nn + inp.IterationStart() + 1), decomp, inp,
                     logs.L2First());
      }
    }
    logs.WriteTime(nn);
  }  // loop for time step -----------------------------------------------------

  if (rank == ROOTP) {
    cout << endl << "Program Complete" << endl;
    PrintTime();
  }

  // Free datatypes previously created
  FreeDataTypesMPI(MPI_vec3d, MPI_procBlockInts, MPI_connection,
                   MPI_DOUBLE_5INT, MPI_vec3dMag, MPI_uncoupledScalar,
                   MPI_tensorDouble);
  // Free operation previously created
  MPI_Op_free(&MPI_MAX_LINF);

  MPI_Finalize();

  return EXIT_SUCCESS;
}

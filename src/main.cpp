/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

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

  // Start clock to time simulation
  const auto start = std::chrono::high_resolution_clock::now();

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

  // nondimensionalize fluid data
  inp.NondimensionalizeFluid();

  // Get physics models
  const auto phys = inp.AssignPhysicsModels();

  // Nondimensionalize BC & IC data
  inp.NondimensionalizeStateData(phys.EoS());

  vector<plot3dBlock> mesh;
  vector<connection> connections;
  vector<procBlock> stateBlocks;
  vector<vector3d<double>> viscFaces;
  // l2 norm residuals to normalize by
  residual residL2First(inp.NumEquations(), inp.NumSpecies());

  if (rank == ROOTP) {
    cout << "Number of equations: " << inp.NumEquations() << endl << endl;

    // Read grid
    mesh = ReadP3dGrid(inp.GridName(), inp.LRef(), totalCells);
    vector<vector3d<int>> gridSizes(mesh.size());
    for (auto ii = 0U; ii < mesh.size(); ++ii) {
      gridSizes[ii] = {mesh[ii].NumCellsI(), mesh[ii].NumCellsJ(),
                       mesh[ii].NumCellsK()};
    }

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

    // Get connection BCs
    connections = GetConnectionBCs(bcs, mesh, decomp, inp);

    // Initialize the whole mesh with ICs and assign ghost cells geometry
    stateBlocks.resize(mesh.size());
    for (auto ll = 0U; ll < mesh.size(); ll++) {
      stateBlocks[ll] =
          procBlock(mesh[ll], decomp.ParentBlock(ll), bcs[ll], ll,
                    decomp.Rank(ll), decomp.LocalPosition(ll), inp);
      stateBlocks[ll].InitializeStates(inp, phys);
      stateBlocks[ll].AssignGhostCellsGeom();
    }
    // if restart, get data from restart file
    if (inp.IsRestart()) {
      ReadRestart(stateBlocks, restartFile, decomp, inp, phys, residL2First,
                  gridSizes);
    }

    // Swap geometry for connection BCs
    for (auto &conn : connections) {
      SwapGeomSlice(conn, stateBlocks[conn.BlockFirst()],
                    stateBlocks[conn.BlockSecond()]);
    }
    // Get ghost cell edge data
    for (auto &block : stateBlocks) {
      block.AssignGhostCellsGeomEdge();
    }

    // Get face centers of faces with viscous wall BC
    viscFaces = GetViscousFaceCenters(stateBlocks);

    cout << "Solution Initialized" << endl << endl;
    //---------------------------------------------------------------------
  }

  // Set MPI datatypes
  MPI_Datatype MPI_vec3d, MPI_procBlockInts, MPI_connection, MPI_DOUBLE_5INT,
      MPI_vec3dMag, MPI_uncoupledScalar, MPI_tensorDouble;
  SetDataTypesMPI(MPI_vec3d, MPI_procBlockInts, MPI_connection, MPI_DOUBLE_5INT,
                  MPI_vec3dMag, MPI_uncoupledScalar, MPI_tensorDouble);

  // Send number of procBlocks to all processors
  SendNumProcBlocks(decomp.NumBlocksOnAllProc(), numProcBlock);

  // Send procBlocks to appropriate processor
  auto localStateBlocks =
      SendProcBlocks(stateBlocks, rank, numProcBlock, MPI_vec3d, MPI_vec3dMag,
                     inp);

  // Update auxillary variables (temperature, viscosity, etc), cell widths
  for (auto &block : localStateBlocks) {
    block.UpdateAuxillaryVariables(phys, false);
    block.CalcCellWidths();
  }

  // Send connections to all processors
  SendConnections(connections, MPI_connection);

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
  kdtree tree(viscFaces);

  if (rank == ROOTP) {
    const auto kdEnd = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> kdDuration = kdEnd - wallStart;
    cout << "K-d tree complete after " << kdDuration.count() << " seconds"
         << endl;
  }

  if (tree.Size() > 0) {
    CalcWallDistance(localStateBlocks, tree);
    SwapWallDist(localStateBlocks, connections, rank, inp.NumberGhostLayers());
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == ROOTP) {
    const auto wallEnd = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> wallDuration = wallEnd - wallStart;
    cout << "Wall distance calculation finished after " << wallDuration.count()
         << " seconds" << endl << endl;
  }

  //-----------------------------------------------------------------------
  // Allocate array for flux jacobian
  vector<matMultiArray3d> mainDiagonal(numProcBlock);
  if (inp.IsImplicit()) {
    ResizeArrays(localStateBlocks, inp, mainDiagonal);
  }

  // Send/recv solutions - necessary to get wall distances
  GetProcBlocks(stateBlocks, localStateBlocks, rank, MPI_uncoupledScalar,
                MPI_vec3d, MPI_tensorDouble, inp);

  ofstream resFile;
  if (rank == ROOTP) {
    // Open residual file
    if (inp.IsRestart()) {
      resFile.open(inp.SimNameRoot() + ".resid", ios::app);
    } else {
      resFile.open(inp.SimNameRoot() + ".resid", ios::out);
    }
    if (resFile.fail()) {
      cerr << "ERROR: Could not open residual file!" << endl;
      exit(EXIT_FAILURE);
    }

    // Write out cell centers grid file
    WriteCellCenter(inp.GridName(), stateBlocks, decomp, inp);

    // Write out initial results
    WriteFun(stateBlocks, phys, inp.IterationStart(), decomp, inp);
    WriteMeta(inp, inp.IterationStart());
  }

  // ----------------------------------------------------------------------
  // ----------------------- Start Main Loop ------------------------------
  // ----------------------------------------------------------------------
  // loop over time
  for (auto nn = 0; nn < inp.Iterations(); nn++) {
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate cfl number
    inp.CalcCFL(nn);

    // Store time-n solution, for time integration methods that require it
    if (inp.NeedToStoreTimeN()) {
      AssignSolToTimeN(localStateBlocks, phys);
      if (!inp.IsRestart() && inp.IsMultilevelInTime() && nn == 0) {
        AssignSolToTimeNm1(localStateBlocks);
      }
    }

    // loop over nonlinear iterations
    for (auto mm = 0; mm < inp.NonlinearIterations(); mm++) {
      // Get boundary conditions for all blocks
      GetBoundaryConditions(localStateBlocks, inp, phys, connections, rank);

      // Calculate residual (RHS)
      CalcResidual(localStateBlocks, mainDiagonal, phys, inp, connections, rank,
                   MPI_tensorDouble, MPI_vec3d);

      // Calculate time step
      CalcTimeStep(localStateBlocks, inp);

      // Initialize residual variables
      // l2 norm residuals
      residual residL2(inp.NumEquations(), inp.NumSpecies());
      resid residLinf;  // linf residuals
      auto matrixResid = 0.0;
      if (inp.IsImplicit()) {
        matrixResid = ImplicitUpdate(localStateBlocks, mainDiagonal, inp, phys,
                                     mm, residL2, residLinf, connections, rank);
      } else {  // explicit time integration
        ExplicitUpdate(localStateBlocks, inp, phys, mm, residL2, residLinf);
      }

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
        WriteResiduals(inp, residL2First, residL2, residLinf, matrixResid,
                       nn + inp.IterationStart(), mm, resFile);
      }
    }  // loop for nonlinear iterations ---------------------------------------

    // write out function file
    if (inp.WriteOutput(nn) || inp.WriteRestart(nn)) {
      // Send/recv solutions
      GetProcBlocks(stateBlocks, localStateBlocks, rank, MPI_uncoupledScalar,
                    MPI_vec3d, MPI_tensorDouble, inp);

      if (rank == ROOTP && inp.WriteOutput(nn)) {
        cout << "writing out function file at iteration "
             << nn + inp.IterationStart()<< endl;
        // Write out function file
        WriteFun(stateBlocks, phys, (nn + inp.IterationStart() + 1), decomp,
                 inp);
        WriteMeta(inp, (nn + inp.IterationStart() + 1));
      }
      if (rank == ROOTP && inp.WriteRestart(nn)) {
        cout << "writing out restart file at iteration "
             << nn + inp.IterationStart()<< endl;
        // Write out restart file
        WriteRestart(stateBlocks, phys, (nn + inp.IterationStart() + 1), decomp,
                     inp, residL2First);
      }
    }
  }  // loop for time step -----------------------------------------------------

  if (rank == ROOTP) {
    // close residual file
    resFile.close();

    cout << endl << "Program Complete" << endl;
    PrintTime();

    const auto simEnd = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> duration = simEnd - start;
    cout << "Total Time: " << duration.count() << " seconds" << endl;
  }

  // Free datatypes previously created
  FreeDataTypesMPI(MPI_vec3d, MPI_procBlockInts, MPI_connection,
                   MPI_DOUBLE_5INT, MPI_vec3dMag, MPI_uncoupledScalar,
                   MPI_tensorDouble);

  MPI_Finalize();

  return EXIT_SUCCESS;
}

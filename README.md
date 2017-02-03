# Accurate Implicit Three-dimensional Efficient RANS (AITHER)

[![Build Status](https://travis-ci.org/mnucci32/aither.svg?branch=master)](https://travis-ci.org/mnucci32/aither)

### About The code
This code is for a 3D Navier-Stokes computational fluid dynamics solver. It is a cell centered, structured solver, using mulit-block structured grids in Plot3D format. It uses explicit and implicit time integration methods. It uses MUSCL extrapolation to reconstruct the primative variables from the cell centers to the cell faces. The code uses the Roe flux difference splitting scheme for the inviscid fluxes, and a central scheme for the viscous fluxes. It is second order accurate in both space and time.

### Current Status
The code is 2nd order accurate in space and time. Available explicit time integration methods are forward euler (1st order) and a minimum storage four stage Runge-Kutta method (2nd order). The implicit solver (LU-SGS, BLU-SGS, DPLUR, BDPLUR) is implemented for implicit time integration. Dual time stepping is implemented for time accuracy in the implicit solver. Available implicit time integrations methods come from the Beam and Warming family of methods and are the implicit euler (1st order), Crank-Nicholson (2nd order), and BDF2 (2nd order) methods. The code has been thoroughly commented. It has been made parallel using MPI. Currently the Wilcox K-Omega 2006 and SST 2003 turbulence models are available.

### To Do List
* Add SST-DES turbulence model
* Add WALE and Smagorinsky subgrid scale models for LES
* Add Couette flow regression test for isothermal wall, moving wall, periodic boundary conditions
* Add restart capability
* Add wall functions for turbulence models
* Add multigrid scheme for improved convergence

### Dependencies
* MPI - OpenMPI and MPICH have both been used in the past. Aither is currently developed with OpenMPI
* C++ compiler with C++14 support
* Cmake - Cmake only depends on a C++ compiler

### How To compile
Aither is compiled and installed with the standard cmake process.

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation /path/to/source
make
make install
```

Cmake will automatically look for an MPI package. To specify a specific installation, set *-DMPI_DIR* to the MPI installation directory.

### How To Run
```bash
mpirun -np 1 aither inputFile.inp <restartFile.rst> >outputFile.out 2>errorFile.err &
```
The restart file argument is optional.

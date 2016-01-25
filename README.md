# Accurate Implicit Three-dimensional Efficient RANS (AITHER)

### About The code
This code is for a 3D Navier-Stokes computational fluid dynamics solver. It is a cell centered, structured solver, using mulit-block structured grids in Plot3D format. It uses explicit and implicit time integration methods. It uses MUSCL extrapolation to reconstruct the primative variables from the cell centers to the cell faces. The code uses the Roe flux difference splitting scheme for the inviscid fluxes, and a central scheme for the viscous fluxes. It is second order accurate in both space and time.

### Current Status
The code is 2nd order accurate in space and time. Available explicit time integration methods are forward euler (1st order) and a minimum storage four stage Runge-Kutta method (2nd order). The implicit solver (LU-SGS) is implemented for implicit time integration. Dual time stepping is implemented for time accuracy in the implicit solver. Available implicit time integrations methods come from the Beam and Warming family of methods and are the implicit euler (1st order), Crank-Nicholson (2nd order), and BDF2 (2nd order) methods. The code has been thoroughly commented. It has been made parallel using MPI. Currently the Wilcox K-Omega 2006 and SST 2003 turbulence models are available.

### To Do List
* Implement block LUSGS implicit method
* Implement higher order state reconstruction for DES / Hybrid RANS/LES
* Implement tubulence model for DES / Hybrid RANS/LES

### How To compile
Assuming you have g++ (or equivalent) with c++14 support and some version of MPI (I use MPICH), just type 'make'.

### How To Run
```bash
mpirun -np 1 aither inputFile.inp >outputFile.out 2>errorFile.err &
```
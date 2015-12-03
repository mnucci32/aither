---
layout: post
title: "Compiling Aither & Running A Simulation"
date: 2015-11-24 22:30
tags: [CFD, Aither, MPI, C++, Compiler]
comments: true
---
# Compiling Aither
To compile Aither you must have a C++ compiler with C++14 support. [G++](https://gcc.gnu.org/) version 5 and up should work, as well as [Clang](http://clang.llvm.org/) 3.4 and up. Before compiling, you must have an MPI implementation installed. Aither should work with any standard implementation, and has been used with both [OpenMPI](http://www.open-mpi.org/) and [MPICH](https://www.mpich.org/). Once these dependencies are installed, all that is necessary to compile the code is to navigate to the source code directory and type *make*.

# Running A Simulation
The Aither repository comes with a suite of test cases to run. For example to run the transonic flow over a bump in a channel case, navigate to the *testCases/transonicBump* directory and run the simulation with the following command.

{% highlight bash %}
mpirun -np 1 aither transonicBump.inp >transonicBump.out 2>transonicBump.err
{% endhighlight %}

The results will be **_center.fun* files written to the run directory. The residuals will be output to standard out which the above command has placed in the *transonicBump.out* file. Any warnings or errors will be written to standard error which is placed in the *transonicBump.err* file. Aither is a cell-centered code, so the solution is calculated at the centroid of each cell. This is why *center* is attached to the results files. The centroids of the grid file will also be written out as **_center.xyz* for visualization purposes. The PLOT3D file of the cell centroids and simulation results can be read into [Paraview](http://www.paraview.org/) for visualization.
#ifndef PARALLELHEADERDEF             //only if the macro PARALLELHEADERDEF is not defined execute these lines of code
#define PARALLELHEADERDEF             //define the macro

/* This header contains the function declarations for file output.

It contains function headers to write out the grid at the cell centers in Plot3D format, as well as the Plot3D function files.
It also writes out a master file in Ensight format to name the Plot3D functions.  */

#include "mpi.h" //parallelism
#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "plot3d.h" //plot3d
#include "primVars.h" //primVars
#include "procBlock.h" //procBlock
#include "geomSlice.h" //geomSlice, stateSlice
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::cout;
using std::endl;
using std::cerr;

#define ROOT 0

//function definitions
void ManualDecomposition(vector<procBlock>&, const int&);

void SetDataTypesMPI(const int&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&);

vector<procBlock> SendProcBlocks( const vector<procBlock>&, const int&, const MPI_Datatype&, const MPI_Datatype&, const MPI_Datatype& );

#endif

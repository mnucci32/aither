#ifndef PARALLELHEADERDEF             //only if the macro PARALLELHEADERDEF is not defined execute these lines of code
#define PARALLELHEADERDEF             //define the macro

/* This header contains the function declarations for many of the parallel functions in the code */

#include "mpi.h" //parallelism
#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "plot3d.h" //plot3d
#include "primVars.h" //primVars
#include "procBlock.h" //procBlock
#include "boundaryConditions.h" //interblock
#include <iostream>
#include "macros.h"

using std::vector;
using std::string;
using std::ios;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;


//function definitions
decomposition ManualDecomposition(vector<plot3dBlock>&, vector<boundaryConditions>&, const int& );
decomposition CubicDecomposition(vector<plot3dBlock>&, vector<boundaryConditions>&, const int& );

void SendNumProcBlocks(const vector<int>&, const int&, int&);
void SendConnections(vector<interblock>&, const MPI_Datatype& );

void SetDataTypesMPI(const int&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype& );

vector<procBlock> SendProcBlocks( const vector<procBlock>&, const int&, const int&, const MPI_Datatype&, const MPI_Datatype& );
void GetProcBlocks( vector<procBlock>&, const vector<procBlock>&, const int&, const MPI_Datatype& );

void MaxLinf( resid*, resid*, int*, MPI_Datatype* );

void BroadcastString( string &str );
#endif

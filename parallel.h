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

class decomposition {
  vector<int> rank;                 //rank of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> parBlock;             //parent block of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> localPos;             //local position of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> splitHistBlkLow;      //lower block of split (vector size equals number of splits)
  vector<int> splitHistBlkUp;       //upper block of split (vector size equals number of splits)
  vector<string> splitHistDir;      //direction of split (vector size eqals number of splits)
  int numProcs;                     //number of processors

 public:
  //constructor
  decomposition();
  decomposition( const int&, const int& );

  //member functions
  double IdealLoad( const vector<plot3dBlock>&, const int& );   //determines ideal load given number of processors
  double MaxLoad( const vector<plot3dBlock>&);                  //determines maximum load
  int NumBlocksOnProc( const int&);                             //determines number of blocks on a given processor
  vector<int> NumBlocksOnAllProc();                  //determines number of blocks on all processors
  void SendToProc( const int&, const int&);                     //sends an entire block to a given processor
  void Split( const int&, const int&, const string& );          //adds data for a split

  //destructor
  ~decomposition() {}

};

//function definitions
vector<int> ManualDecomposition(vector<plot3dBlock>&, vector<vector3d<int> >&, const int&, const int& );
vector<int> CubicDecomposition(vector<plot3dBlock>&, vector<vector3d<int> >&, vector<boundaryConditions>&, const int&, const int& );

void SendNumProcBlocks(const vector<int>&, const int&, int&);
void SendConnections(vector<interblock>&, const MPI_Datatype& );

void SetDataTypesMPI(const int&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype& );

vector<procBlock> SendProcBlocks( const vector<procBlock>&, const int&, const int&, const MPI_Datatype&, const MPI_Datatype& );
void GetProcBlocks( vector<procBlock>&, const vector<procBlock>&, const int&, const MPI_Datatype& );

void MaxLinf( resid*, resid*, int*, MPI_Datatype* );

void BroadcastString( string &str );
#endif

#include "procBlock.h"
#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;


//constructors for geomSlice class
geomSlice::geomSlice(){
  numCells = 1;
  numI = 1;
  numJ = 1;
  numK = 1;
  parBlock = 0;
  parBlockStartI = 0;
  parBlockEndI = 0;
  parBlockStartJ = 0;
  parBlockEndJ = 0;
  parBlockStartK = 0;
  parBlockEndK = 0;

  int numFaces = (numI+1)*(numJ)*(numK);  

  vector<vector3d<double> > vec1(numFaces);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vec2(numCells);                 //dummy vector variable length of number of cells
  vector<double> scalar(numCells);                          //dummy scalar variable length of number of cells

  center = vec2;
  fAreaI = vec1;
  fAreaJ = vec1;
  fAreaK = vec1;
  fCenterI = vec1;
  fCenterJ = vec1;
  fCenterK = vec1;

  vol = scalar;

}
//constructor -- initialize state vector with dummy variables
geomSlice::geomSlice( const int &li, const int &lj, const int &lk, const int &pblk, const int &pStartI, const int &pEndI, 
		      const int &pStartJ, const int &pEndJ, const int &pStartK, const int &pEndK ){
  // li -- size of direction i (cell)
  // lj -- size of direction j (cell)
  // lk -- size of direction k (cell)
  // pblk -- parent block that slice is coming from
  // pStartI -- starting i-index from parent block (face)
  // pEndI -- ending i-index from parent block (face)
  // pStartJ -- starting j-index from parent block (face)
  // pEndJ -- ending j-index from parent block (face)
  // pStartK -- starting k-index from parent block (face)
  // pEndK -- ending k-index from parent block (face)

  numI = li;
  numJ = lj;
  numK = lk;

  numCells = li * lj * lk;

  parBlock = pblk;
  parBlockStartI = pStartI;
  parBlockEndI = pEndI;
  parBlockStartJ = pStartJ;
  parBlockEndJ = pEndJ;
  parBlockStartK = pStartK;
  parBlockEndK = pEndK;

  int numIFaces = (numI+1)*(numJ)*(numK);  
  int numJFaces = (numI)*(numJ+1)*(numK);  
  int numKFaces = (numI)*(numJ)*(numK+1);  

  vector<vector3d<double> > vecIFaces(numIFaces);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vecJFaces(numJFaces);                 //dummy vector variable length of number of faces
  vector<vector3d<double> > vecKFaces(numKFaces);                 //dummy vector variable length of number of faces

  vector<vector3d<double> > vecCells(numCells);                   //dummy vector variable length of number of cells
  vector<double> scalar(numCells);                                //dummy scalar variable length of number of cells

  center = vecCells;
  fAreaI = vecIFaces;
  fAreaJ = vecJFaces;
  fAreaK = vecKFaces;
  fCenterI = vecIFaces;
  fCenterJ = vecJFaces;
  fCenterK = vecKFaces;

  vol = scalar;

}

//constructors for stateSlice class
stateSlice::stateSlice(){
  numCells = 1;
  numI = 1;
  numJ = 1;
  numK = 1;
  parBlock = 0;
  parBlockStartI = 0;
  parBlockEndI = 0;
  parBlockStartJ = 0;
  parBlockEndJ = 0;
  parBlockStartK = 0;
  parBlockEndK = 0;

  vector<primVars> prims(numCells);                          //dummy primVars variable length of number of cells

  state = prims;

}
//constructor -- initialize state vector with dummy variables
stateSlice::stateSlice( const int &li, const int &lj, const int &lk, const int &pblk, const int &pStartI, const int &pEndI, 
		      const int &pStartJ, const int &pEndJ, const int &pStartK, const int &pEndK ){
  // li -- size of direction i
  // lj -- size of direction j
  // lk -- size of direction k
  // pblk -- parent block that slice is coming from
  // pStartI -- starting i-index from parent block (face)
  // pEndI -- ending i-index from parent block (face)
  // pStartJ -- starting j-index from parent block (face)
  // pEndJ -- ending j-index from parent block (face)
  // pStartK -- starting k-index from parent block (face)
  // pEndK -- ending k-index from parent block (face)

  numI = li;
  numJ = lj;
  numK = lk;

  numCells = li * lj * lk;

  parBlock = pblk;
  parBlockStartI = pStartI;
  parBlockEndI = pEndI;
  parBlockStartJ = pStartJ;
  parBlockEndJ = pEndJ;
  parBlockStartK = pStartK;
  parBlockEndK = pEndK;

  vector<primVars> prims(numCells);                                //dummy primVars variable length of number of cells

  state = prims;

}


void stateSlice::PackSwapUnpackMPI( const interblock &inter, const MPI_Datatype &MPI_cellData, const int &rank ) {


  //swap with mpi_send_recv_replace
  //pack data into buffer, but first get size
  int bufSize = 0;
  int tempSize = 0;
  MPI_Pack_size((*this).NumCells(), MPI_cellData, MPI_COMM_WORLD, &tempSize); //add size for states
  bufSize += tempSize;
  MPI_Pack_size(11, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for ints in class stateSlice
  bufSize += tempSize;

  char *buffer = new char[bufSize]; //allocate buffer to pack data into

  //pack data into buffer
  int position = 0;
  MPI_Pack(&(*this).state[0], (*this).NumCells(), MPI_cellData, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numCells, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlock, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockStartI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockEndI, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockStartJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockEndJ, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockStartK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).parBlockEndK, 1, MPI_INT, buffer, bufSize, &position, MPI_COMM_WORLD);

  MPI_Status status;
  if ( rank == inter.RankFirst() ){ //send/recv with second entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankSecond(), 1, inter.RankSecond(), 1, MPI_COMM_WORLD, &status);
  }
  else{ //send/recv with first entry in interblock
    MPI_Sendrecv_replace(buffer, bufSize, MPI_PACKED, inter.RankFirst(), 1, inter.RankFirst(), 1, MPI_COMM_WORLD, &status);
  }

  //put slice back into stateSlice
  position = 0;
  MPI_Unpack(buffer, bufSize, &position, &(*this).state[0], (*this).NumCells(), MPI_cellData, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numCells, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).numK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlock, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockStartI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockEndI, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockStartJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockEndJ, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockStartK, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, bufSize, &position, &(*this).parBlockEndK, 1, MPI_INT, MPI_COMM_WORLD);

  delete [] buffer;

}

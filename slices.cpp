/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

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
geomSlice::geomSlice( const int &li, const int &lj, const int &lk, const int &pblk ){
  // li -- size of direction i (cell)
  // lj -- size of direction j (cell)
  // lk -- size of direction k (cell)
  // pblk -- parent block that slice is coming from

  numI = li;
  numJ = lj;
  numK = lk;

  numCells = li * lj * lk;

  parBlock = pblk;

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

  vector<primVars> prims(numCells);                          //dummy primVars variable length of number of cells

  state = prims;

}
//constructor -- initialize state vector with dummy variables
stateSlice::stateSlice( const int &li, const int &lj, const int &lk, const int &pblk ){
  // li -- size of direction i
  // lj -- size of direction j
  // lk -- size of direction k
  // pblk -- parent block that slice is coming from

  numI = li;
  numJ = lj;
  numK = lk;

  numCells = li * lj * lk;

  parBlock = pblk;

  vector<primVars> prims(numCells);                                //dummy primVars variable length of number of cells

  state = prims;

}

/*Member function to pack a stateslice into a buffer, swap it with its interblock partner, and then unpack it into a stateslice.*/
void stateSlice::PackSwapUnpackMPI( const interblock &inter, const MPI_Datatype &MPI_cellData, const int &rank ) {
  // inter -- interblock boundary for the swap
  // MPI_cellData -- MPI datatype to pass cell state data
  // rank -- processor rank

  //swap with mpi_send_recv_replace
  //pack data into buffer, but first get size
  int bufSize = 0;
  int tempSize = 0;
  MPI_Pack_size((*this).NumCells(), MPI_cellData, MPI_COMM_WORLD, &tempSize); //add size for states
  bufSize += tempSize;
  MPI_Pack_size(5, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for ints in class stateSlice
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

  delete [] buffer;

}

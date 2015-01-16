#include "geomSlice.h"
#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::max;
using std::min;

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


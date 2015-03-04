#include "boundaryConditions.h"
#include <math.h>       //sqrt
#include <iostream>     //cout
#include "vector3d.h"     //vector3d

using std::cout;
using std::endl;
using std::cerr;
using std::swap;

//constructor when passed no arguements
boundaryConditions::boundaryConditions(){
  numSurfI = 2;
  numSurfJ = 2;
  numSurfK = 2;
  int length = numSurfI + numSurfJ + numSurfK;

  boundarySurface bcSurf;
  vector<boundarySurface> dumVec(length,bcSurf);
  surfs = dumVec;
}

//constructor when passed number of i, j, k surfaces
boundaryConditions::boundaryConditions( const int &i, const int &j, const int &k){
  numSurfI = i;
  numSurfJ = j;
  numSurfK = k;
  int length = numSurfI + numSurfJ + numSurfK;

  boundarySurface bcSurf;
  vector<boundarySurface> dumVec(length,bcSurf);
  surfs = dumVec;
}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const boundaryConditions &bc){

  os << "Number of surfaces (I, J, K): " << bc.numSurfI << ", " << bc.numSurfJ << ", " << bc.numSurfK << endl;

  for ( int ii = 0; ii < bc.NumSurfaces(); ii++ ){
    os << bc.surfs[ii] << endl;
  }

  return os;
}

//operator to resize all of the vector components of the boundary conditions class
void boundaryConditions::ResizeVecs( const int &a){
  surfs.resize(a);
}

//operator to resize all of the vector components of the boundary conditions class
void boundaryConditions::ResizeVecs( const int &i, const int &j, const int &k){

  numSurfI = i;
  numSurfJ = j;
  numSurfK = k;

  surfs.resize(i + j + k);
}

//member function to return the boundary condition type given the i,j,k face coordinates and the surface type
string boundaryConditions::GetBCName(const int i, const int j, const int k, const string& surf)const{

  string bcName;
  int iStart = 0;
  int iEnd = 0;

  if (surf == "il" || surf == "iu"){ //i-surfaces search between 0 and number of i-surfaces
    iStart = 0;
    iEnd = (*this).NumSurfI();
  }
  else if (surf == "jl" || surf == "ju"){ //j-surfaces search between end of i-surfaces and end of j-surfaces
    iStart = (*this).NumSurfI();
    iEnd = iStart + (*this).NumSurfJ();
  }
  else if (surf == "kl" || surf == "ku"){ //k-surfaces search between end of j-surfaces and end of k-surfaces
    iStart = (*this).NumSurfI() + (*this).NumSurfJ();
    iEnd = iStart + (*this).NumSurfK();
  }
  else {
    cerr << "ERROR: Surface type " << surf << " is not recognized!" << endl;
  }

  //Determine which boundary condition should be applied
  for ( int nn = iStart; nn < iEnd; nn++ ){

    //Boundary mins and maxes start at 1 instead of 0, so 1 is subtracted
    //determine which boundary given i, j, k coordinates apply to
    if ( (i >= (*this).GetIMin(nn)-1 && i <= (*this).GetIMax(nn)-1 && j >= (*this).GetJMin(nn)-1 && j <= (*this).GetJMax(nn)-1
	  && k >= (*this).GetKMin(nn)-1 && k <= (*this).GetKMax(nn)-1) ){
      bcName = (*this).GetBCTypes(nn);
      break;
    }
  }

  return bcName;

}

//member function to fill one "row" of the vectors with data that has been read in from the input file. This function is called from
//input::ReadInput(). It is necessary so that the private data can be altered from another class's member function.
void boundaryConditions::AssignFromInput(const int &surfCounter, const vector<string> &tokens){

  boundarySurface bcSurf(tokens[0], atoi(tokens[1].c_str()), atoi(tokens[2].c_str()), atoi(tokens[3].c_str()), atoi(tokens[4].c_str()), 
			 atoi(tokens[5].c_str()), atoi(tokens[6].c_str()), atoi(tokens[7].c_str()) );
  surfs[surfCounter] = bcSurf;
}


void boundaryConditions::BordersInterblock(const int &ii, bool (&border)[4])const{
  // ii -- index of surface to test for border matches
  // border -- array of bools to show if boundarySurface is bordered by an interblock on any of its 4 sides

  boundarySurface surf = (*this).GetSurface(ii);

  //check that given boundarySurface is interblock
  if ( surf.BCType() != "interblock" ){
    cerr << "ERROR: Error in boundaryConditions::BordersInterblock(). Given index does not point to an interblock boundarySurface!" << endl;
    cerr << surf << endl;
    exit(0);
  }

  border[0] = false;
  border[1] = false;
  border[2] = false;
  border[3] = false;

  for ( int jj = 0; jj < (*this).NumSurfaces(); jj++ ){
    boundarySurface possibleBorder = (*this).GetSurface(jj);
    //if possible border is an interblock and of same surface type, test for border match
    if ( possibleBorder.BCType() == "interblock" && possibleBorder.SurfaceType() == surf.SurfaceType() ){
      if ( surf.Min1() == possibleBorder.Max1() ){
	border[0] = true;
      }
      if ( surf.Max1() == possibleBorder.Min1() ){
	border[1] = true;
      }
      if ( surf.Min2() == possibleBorder.Max2() ){
	border[2] = true;
      }
      if( surf.Max2() == possibleBorder.Min2() ){
	border[3] = true;
      }
    }
  }

}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const interblock &bc){

  os << "Ranks: " << bc.RankFirst() << ", " << bc.RankSecond() << endl;
  os << "Blocks: " << bc.BlockFirst() << ", " << bc.BlockSecond() << endl;
  os << "Local Blocks: " << bc.LocalBlockFirst() << ", " << bc.LocalBlockSecond() << endl;
  os << "Boundaries: " << bc.BoundaryFirst() << ", " << bc.BoundarySecond() << endl;
  os << "Direction 1 Starts: " << bc.Dir1StartFirst() << ", " << bc.Dir1StartSecond() << endl;
  os << "Direction 1 Ends: " << bc.Dir1EndFirst() << ", " << bc.Dir1EndSecond() << endl;
  os << "Direction 2 Starts: " << bc.Dir2StartFirst() << ", " << bc.Dir2StartSecond() << endl;
  os << "Direction 2 Ends: " << bc.Dir2EndFirst() << ", " << bc.Dir2EndSecond() << endl;
  os << "Direction 3 Constant Surface: " << bc.ConstSurfaceFirst() << ", " << bc.ConstSurfaceSecond() << endl;
  os << "Direction 1 Start Borders Interblock: " << bc.Dir1StartInterBorderFirst() << ", " << bc.Dir1StartInterBorderSecond() << endl;
  os << "Direction 1 End Borders Interblock: " << bc.Dir1EndInterBorderFirst() << ", " << bc.Dir1EndInterBorderSecond() << endl;
  os << "Direction 2 Start Borders Interblock: " << bc.Dir2StartInterBorderFirst() << ", " << bc.Dir2StartInterBorderSecond() << endl;
  os << "Direction 2 End Borders Interblock: " << bc.Dir2EndInterBorderFirst() << ", " << bc.Dir2EndInterBorderSecond() << endl;
  os << "Orientation: " << bc.Orientation() << endl;

  return os;
}

//constructor to take in two patches and fill an interblock. The orientation is left at the default value 0.
interblock::interblock(const patch &p1, const patch &p2){
  // p1 -- patch 1
  // p2 -- patch 2

  //fill interblock
  rank[0] = p1.Rank();
  rank[1] = p2.Rank();

  block[0] = p1.Block();
  block[1] =p2.Block();

  localBlock[0] = p1.LocalBlock();
  localBlock[1] = p2.LocalBlock();

  boundary[0] = p1.Boundary();
  boundary[1] = p2.Boundary();

  d1Start[0] = p1.Dir1Start();
  d1Start[1] = p2.Dir1Start();

  d1End[0] = p1.Dir1End();
  d1End[1] = p2.Dir1End();

  d2Start[0] = p1.Dir2Start();
  d2Start[1] = p2.Dir2Start();

  d2End[0] = p1.Dir2End();
  d2End[1] = p2.Dir2End();

  constSurf[0] = p1.ConstSurface();
  constSurf[1] = p2.ConstSurface();

  interblockBorder[0] = p1.Dir1StartInterBorder();
  interblockBorder[1] = p1.Dir1EndInterBorder();
  interblockBorder[2] = p1.Dir2StartInterBorder();
  interblockBorder[3] = p1.Dir2EndInterBorder();
  interblockBorder[4] = p2.Dir1StartInterBorder();
  interblockBorder[5] = p2.Dir1EndInterBorder();
  interblockBorder[6] = p2.Dir2StartInterBorder();
  interblockBorder[7] = p2.Dir2EndInterBorder();

  orientation = 0; //default value (real values 1-8)
}

//function to swap the order of an interblock so the 2nd entry in the pair will be the first, and vice versa
void interblock::SwapOrder(){

  swap(rank[0], rank[1]);
  swap(block[0], block[1]);
  swap(localBlock[0], localBlock[1]);
  swap(boundary[0], boundary[1]);
  swap(d1Start[0], d1Start[1]);
  swap(d1End[0], d1End[1]);
  swap(d2Start[0], d2Start[1]);
  swap(d2End[0], d2End[1]);
  swap(constSurf[0], constSurf[1]);

  swap(interblockBorder[0], interblockBorder[4]);
  swap(interblockBorder[1], interblockBorder[5]);
  swap(interblockBorder[2], interblockBorder[6]);
  swap(interblockBorder[3], interblockBorder[7]);

  //if orientation is 4 or 5, needs to be swapped because direction 1/2 are swapped and only one direction is reversed
  if (orientation == 4){
    orientation = 5;
  }
  else if (orientation == 5){
    orientation = 4;
  }
}

/* Function to go through the boundary conditions and pair the interblock BCs together and determine their orientation.
*/
vector<interblock> GetInterblockBCs( const vector<boundaryConditions> &bc, const vector<plot3dBlock> &grid, const vector<vector3d<int> > &rankParPos ){
  // bc -- vector of boundaryConditions for all blocks
  // grid -- vector of plot3Dblocks for entire computational mesh
  // rankParPos -- rank of blocks in grid, parent block in grid, local position on processor

  //isolate only the interblock BCs and their associated data from all of the BCs
  vector<boundarySurface> isolatedInterblocks; //outer vector for each interblock BC, inner vector for information about interblock
  vector<vector3d<int> > numRankPos;
  vector<int> surfaceNums;
  for ( unsigned int ii = 0; ii < bc.size(); ii++ ){ //loop over all blocks
    for ( int jj = 0; jj < bc[ii].NumSurfaces(); jj++ ){ //loop over number of surfaces in block

      if ( bc[ii].GetBCTypes(jj) == "interblock" ){ //if boundary condition is interblock, store data
	vector3d<int> temp(ii, rankParPos[ii][0], rankParPos[ii][2]);
	numRankPos.push_back(temp);                              //block number of bc, rank, local position
	isolatedInterblocks.push_back(bc[ii].GetSurface(jj));    //boundarySurface of bc
	surfaceNums.push_back(jj);
      }

    }
  }

  //----------------------------------------------------------------------------------------------------
  //intialize vector of interblocks to return
  //size is halved because each interblock pairs with another
  vector<interblock> connections(isolatedInterblocks.size()/2);

  //loop over isolated interblocks
  //ii counts by two because after a pair is found, that data is swapped to ii+1. This allows the next search to avoid the matched pair
  for ( unsigned int ii = 0; ii < isolatedInterblocks.size(); ii+=2 ){
    for ( unsigned int jj = ii+1; jj < isolatedInterblocks.size(); jj++ ){ //loop over possible matches

      //blocks and boundary surfaces between interblocks match
      if ( isolatedInterblocks[ii].PartnerBlock() == numRankPos[jj][0] && isolatedInterblocks[ii].PartnerSurface() == isolatedInterblocks[jj].SurfaceType() ){ //blocks between interblock BCs match

	bool border[4] = {false, false, false, false};
	bc[numRankPos[ii][0]].BordersInterblock(surfaceNums[ii], border);

	//get current patch
	patch cPatch( isolatedInterblocks[ii], grid[numRankPos[ii][0]], numRankPos[ii][0], border, numRankPos[ii][1], numRankPos[ii][2] );

	bc[numRankPos[jj][0]].BordersInterblock(surfaceNums[jj], border);
	//get new patch (possible match)
	patch nPatch( isolatedInterblocks[jj], grid[numRankPos[jj][0]], numRankPos[jj][0], border, numRankPos[jj][1], numRankPos[jj][2] );

	//test for match
	interblock match(cPatch, nPatch);
	if ( match.TestPatchMatch(cPatch, nPatch) ){ //match found
	  connections[ii/2] = match; //store interblock pair
	  swap(isolatedInterblocks[jj], isolatedInterblocks[ii+1]); //swap matched interblock BC to top portion of vector so it is not searched again
	  swap(numRankPos[jj], numRankPos[ii+1]); 
	  swap(surfaceNums[jj], surfaceNums[ii+1]);
	  break; //exit innermost loop and search for next interblock match
	}

      }
    }
  }

  return connections;
}

/* Function to take in two patches and return if they are matched. If there is a match it uses the patches to modify the given
interblock to contain the information on this match.

Each patch is on a constant i, j, or k surface and is a 4 sided rectangle. The match is tested for by determining if the vertexes
on one patch, match up with the vertexes on the other patch. Only 3 vertexes need to match, so only 3 are tested for. If there is a
match the patches can be oriented in 8 different ways with respect to each other. The orientation is stored in the interblock that
is modified.

                         Patch 1                             Patch 2                   Description
                     __________________                 __________________         
                    |C1             C12|	       |C1             C12|
                    |                  |	       |                  |
Orientation 1:      |                  |	       |                  |           Same orientation
                  D1|O               C2|	     D1|O               C2|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | -->D2              
                     __________________                 __________________         
                    |C1             C12|	       |C2             C12|
                    |                  |	       |                  |
Orientation 2:      |                  |	       |                  |           D1/D2 swapped
                  D1|O               C2|	     D2|O               C1|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | -->D1              
                     __________________                 __________________         
                    |C1             C12|	       |O               C2|
                    |                  |	       |                  |
Orientation 3:      |                  |	       |                  |           D1 reversed
                  D1|O               C2|	    -D1|C1             C12|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | -->D2              
                     __________________                 __________________         
                    |C1             C12|	       |C12             C2|
                    |                  |	       |                  |
Orientation 4:      |                  |	       |                  |           D1/D2 swapped, D1 reversed
                  D1|O               C2|	     D2|C1               O|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | D1<--                
                     __________________                 __________________         
                    |C1             C12|	       |O               C1|
                    |                  |	       |                  |
Orientation 5:      |                  |	       |                  |           D1/D2 swapped, D2 reversed
                  D1|O               C2|	    -D2|C2             C12|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | -->D1              
                     __________________                 __________________         
                    |C1             C12|	       |C12             C1|
                    |                  |	       |                  |
Orientation 6:      |                  |	       |                  |           D2 reversed
                  D1|O               C2|	     D1|C2               O|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | D2<--              
                     __________________                 __________________         
                    |C1             C12|	       |C1               O|
                    |                  |	       |                  |
Orientation 7:      |                  |	       |                  |           D1/D2 swapped, D1/D2 reversed
                  D1|O               C2|	    -D2|C12             C2|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | D1<--               
                     __________________                 __________________         
                    |C1             C12|	       |C2               O|
                    |                  |	       |                  |
Orientation 8:      |                  |	       |                  |           D1/D2 reversed
                  D1|O               C2|	    -D1|C12             C1|
                   ^|__________________|	      ^|__________________|
                   | -->D2              	      | D2<--              

The above diagrams show how patch 2 would have to be moved to match up with patch 1. D1 and D2 are the local patch directions. They are
cyclic, so on a constant i-patch, D1 is j, and D2 is k. On a constant j-patch, D1 is k, and D2 is i, etc. O is the origin which is always 
at the minimim of D1 and D2 on the patch. C1 is the corner where D1 is at a max, and D2 is zero. C2 is the corner where D2 is at a max, and
D1 is zero. C12 is the corner where both D1 and D2 are at a max.
*/
bool interblock::TestPatchMatch( const patch &p1, const patch &p2 ){
  // p1 -- first patch
  // p2 -- second patch

  bool match = false; //initialize match to false

  //determine if there is a potential match by comparing origins
  if ( p1.Origin() == p2.Origin() ){ //origins match -----------------------------------------------------------------------------------

    //if origin matches origin, corner 1 can only be at corner 1 or 2
    if ( p1.Corner1() == p2.Corner1() ){ //corner 1s match

      //if all 3 corners match, same orientation
      if ( p1.Corner2() == p2.Corner2() ){ //corner 2s match
	(*this).orientation = 1;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner2() ){ //corner 1 matches corner 2

      //if origins match and 1 matches 2, 2 must match 1
      if ( p1.Corner2() == p2.Corner1() ){ //corner 2 matches corner 1
	(*this).orientation = 2;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else{ //no match
      return match;
    }

  }
  else if ( p1.Origin() == p2.Corner1() ){ //origin matches corner 1 -----------------------------------------------------------------

    //if origin matches corner1, corner 1 can only be at corner 12 or origin
    if ( p1.Corner1() == p2.Origin() ){ 

      //corner 2 must match 12 for match
      if ( p1.Corner2() == p2.Corner12() ){ 
	(*this).orientation = 3;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner12() ){ 

      //corner 2 must match origin for match
      if ( p1.Corner2() == p2.Origin() ){ 
	(*this).orientation = 4;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else{ //no match
      return match;
    }

  }
  else if ( p1.Origin() == p2.Corner2() ){ //origin matches corner 2 -----------------------------------------------------------------

    //if origin matches corner2, corner 1 can only be at corner 12 or origin
    if ( p1.Corner1() == p2.Origin() ){ 

      //corner 2 must match 12 for match
      if ( p1.Corner2() == p2.Corner12() ){ 
	(*this).orientation = 5;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner12() ){ 

      //corner 2 must match origin for match
      if ( p1.Corner2() == p2.Origin() ){ 
	(*this).orientation = 6;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else{ //no match
      return match;
    }

  }
  else if ( p1.Origin() == p2.Corner12() ){ //origin matches opposite corner ---------------------------------------------------------

    //if origin matches corner12, corner 1 can only be at corner 1 or corner 2
    if ( p1.Corner1() == p2.Corner1() ){ 

      //corner 2 must match 2 for match
      if ( p1.Corner2() == p2.Corner2() ){ 
	(*this).orientation = 7;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner2() ){ 

      //corner 2 must match corner 1 for match
      if ( p1.Corner2() == p2.Corner2() ){ 
	(*this).orientation = 8;
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else{ //no match
      return match;
    }

  }
  else{ //no match
    return match;
  }

  return match;
}

/* Member function to adjust the interblock for use with a geomSlice or stateSlice
*/
void interblock::AdjustForSlice( const bool &blkFirst, const int &numG ){
  // blkFirst -- boolean that is true if block to insert into is first
  // numG -- number of ghost cells in block

  if (!blkFirst){ //block to insert into is second, swap order
    (*this).SwapOrder(); //have block be first entry, slice second
  }

  //if at an upper surface, start block at upper boundary (after including ghosts), if at lower surface, start block at 0
  int blkStart = ((*this).BoundaryFirst() % 2 == 0) ? (*this).ConstSurfaceFirst() + numG : 0;
  (*this).constSurf[1] = 0; //slice always starts at 0
  (*this).constSurf[0] = blkStart;
  //adjust direction 1 start and end for ghost cells
  (*this).d1End[1] = (*this).Dir1EndSecond() - (*this).Dir1StartSecond() + 2 * numG;
  (*this).d1End[0] = (*this).Dir1EndFirst() + 2 * numG;
  (*this).d1Start[1] = 0; //slice always starts at 0
  //adjust direction 2 start and end for ghost cells
  (*this).d2End[1] = (*this).Dir2EndSecond() - (*this).Dir2StartSecond() + 2 * numG;
  (*this).d2End[0] = (*this).Dir2EndFirst() + 2 * numG;
  (*this).d2Start[1] = 0; //slice always starts at 0

}

//Member function to get the addresses of an interblock to create an MPI_Datatype
void interblock::GetAddressesMPI(MPI_Aint (&disp)[11])const{

  //get addresses of each field
  MPI_Get_address(&(*this).rank[0],              &disp[0]);
  MPI_Get_address(&(*this).block[0],             &disp[1]);
  MPI_Get_address(&(*this).localBlock[0],        &disp[2]);
  MPI_Get_address(&(*this).boundary[0],          &disp[3]);
  MPI_Get_address(&(*this).d1Start[0],           &disp[4]);
  MPI_Get_address(&(*this).d1End[0],             &disp[5]);
  MPI_Get_address(&(*this).d2Start[0],           &disp[6]);
  MPI_Get_address(&(*this).d2End[0],             &disp[7]);
  MPI_Get_address(&(*this).constSurf[0],         &disp[8]);
  MPI_Get_address(&(*this).interblockBorder[0],  &disp[9]);
  MPI_Get_address(&(*this).orientation,          &disp[10]);

}

/* Member function to split boundary conditions along a given direction at a given index. The calling instance retains the lower portion of the split,
and the returned instance is the upper portion
*/
boundaryConditions boundaryConditions::Split(const string &dir, const int &ind, const int &numBlk, const int &newBlkNum, vector<boundarySurface> &aSurf){

  int indNG = ind + 1; //+1 because boundaries start at 1, not 0

  boundaryConditions bound1 = (*this);
  boundaryConditions bound2 = (*this);

  vector<boundarySurface> alteredSurf;

  vector<int> del1, del2;
  del1.reserve((*this).NumSurfaces()); //reserved for maximum number of deletions
  del2.reserve((*this).NumSurfaces()); //reserved for maximum number of deletions

  if ( dir == "i" ){ //split along i-plane

    int del1J = 0;
    int del1K = 0;
    int del2J = 0;
    int del2K = 0;

    for ( int ii = 0; ii < (*this).NumSurfaces(); ii++ ){
      if ( ii < (*this).NumSurfI() ){ //i-surface
	if ( (*this).GetIMax(ii) == 1){ //lower i surface
	  //no change to lower bc at lower i surface

	  //at lower i surface, upper bc is now interface
	  int tag = 2000 + numBlk; //lower surface matches with upper surface
	  bound2.surfs[ii].bcType = "interblock";          //bcType
	  bound2.surfs[ii].data[0] = (*this).GetIMin(ii);  //imin
	  bound2.surfs[ii].data[1] = (*this).GetIMax(ii);  //imax
	  bound2.surfs[ii].data[6] = tag;                  //tag
	}
	else{ //upper surface
	  //at upper i surface, lower bc is now interface
	  int tag = 1000 + newBlkNum; //upper surface matches with lower surface
	  bound1.surfs[ii].bcType = "interblock"; //bcType
	  bound1.surfs[ii].data[0] = indNG;       //imin
	  bound1.surfs[ii].data[1] = indNG;       //imax
	  bound1.surfs[ii].data[6] = tag;         //tag

	  //at upper i surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.surfs[ii].data[0] = (*this).GetIMax(ii) - indNG + 1;      //imin
	  bound2.surfs[ii].data[1] = (*this).GetIMax(ii) - indNG + 1;      //imax

	  //at upper i surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	  if ( (*this).GetBCTypes(ii) == "interblock" ){
	    alteredSurf.push_back((*this).GetSurface(ii));
	  }

	}
      }
      else { //j-surface or k-surface

	//at j/k surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	if ( (*this).GetIMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.surfs[ii].data[0] = (*this).GetIMax(ii) - indNG + 1;    //imin
	  bound2.surfs[ii].data[1] = (*this).GetIMax(ii) - indNG + 1;    //imax
	  if ( ii >= (*this).NumSurfI() && ii < (*this).NumSurfI() + (*this).NumSurfJ() ){ //j-surface
	    del1J++;
	  }
	  else{ //k-surface
	    del1K++;
	  }
	}
	else if ( (*this).GetIMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.surfs[ii].data[1] = indNG;                               //imax
	  bound2.surfs[ii].data[0] = 1;                                   //imin
	  bound2.surfs[ii].data[1] = (*this).GetIMax(ii) - indNG + 1;     //imax
	}
	else{ //this surface is only present in the lower split
	  del2.push_back(ii);
	  if ( ii >= (*this).NumSurfI() && ii < (*this).NumSurfI() + (*this).NumSurfJ() ){ //j-surface
	    del2J++;
	  }
	  else{ //k-surface
	    del2K++;
	  }
	}
      }
    }

    //delete unnecessary boundaries
    for ( unsigned int ii = 0; ii < del1.size(); ii++ ){
      bound1.surfs.erase(bound1.surfs.begin() + del1[ii]);
      bound1.numSurfJ -= del1J;
      bound1.numSurfK -= del1K;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.surfs.erase(bound2.surfs.begin() + del2[ii]);
      bound2.numSurfJ -= del2J;
      bound2.numSurfK -= del2K;
    }

  }
  else if ( dir == "j" ){ //split along j-plane

    int del1I = 0;
    int del1K = 0;
    int del2I = 0;
    int del2K = 0;

    for ( int ii = 0; ii < (*this).NumSurfaces(); ii++ ){
      if ( ii >= (*this).NumSurfI() && ii < (*this).NumSurfI() + (*this).NumSurfJ() ){ //j-surface
	if ( (*this).GetJMax(ii) == 1){ //lower j surface
	  //no change to lower bc at lower j surface

	  //at lower j surface, upper bc is now interface
	  int tag = 4000 + numBlk; //lower surface matches with upper surface
	  bound2.surfs[ii].bcType = "interblock";           //bctype
	  bound2.surfs[ii].data[2] = (*this).GetJMin(ii);   //jmin
	  bound2.surfs[ii].data[3] = (*this).GetJMax(ii);   //jmax
	  bound2.surfs[ii].data[6] = tag;                   //tag
	}
	else {
	  //at upper j surface, lower bc is now interface
	  int tag = 3000 + newBlkNum; //upper surface matches with lower surface
	  bound1.surfs[ii].bcType = "interblock";         //bctype
	  bound1.surfs[ii].data[2] = indNG;               //jmin
	  bound1.surfs[ii].data[3] = indNG;               //jmax
	  bound1.surfs[ii].data[6] = tag;                 //tag

	  //at upper j surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.surfs[ii].data[2] = (*this).GetJMax(ii) - indNG + 1;  //jmin
	  bound2.surfs[ii].data[3] = (*this).GetJMax(ii) - indNG + 1;  //jmax

	  //at upper j surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	  if ( (*this).GetBCTypes(ii) == "interblock" ){
	    alteredSurf.push_back((*this).GetSurface(ii));
	  }

	}
      }
      else { //i-surface or k-surface

	//at i/k surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	if ( (*this).GetJMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.surfs[ii].data[2] = (*this).GetJMax(ii) - indNG + 1;   //jmin
	  bound2.surfs[ii].data[3] = (*this).GetJMax(ii) - indNG + 1;   //jmax
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del1I++;
	  }
	  else{ //k-surface
	    del1K++;
	  }
	}
	else if ( (*this).GetJMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.surfs[ii].data[3] = indNG;                              //jmax
	  bound2.surfs[ii].data[2] = 1;                                  //jmin
	  bound2.surfs[ii].data[3] = (*this).GetJMax(ii) - indNG + 1;    //jmax
	}
	else{ //this surface is only present in the lower split
	  del2.push_back(ii);
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del2I++;
	  }
	  else{ //k-surface
	    del2K++;
	  }
	}
      }
    }

    //delete unnecessary boundaries
    for ( unsigned int ii = 0; ii < del1.size(); ii++ ){
      bound1.surfs.erase(bound1.surfs.begin() + del1[ii]);
      bound1.numSurfI -= del1I;
      bound1.numSurfK -= del1K;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.surfs.erase(bound2.surfs.begin() + del2[ii]);
      bound2.numSurfI -= del2I;
      bound2.numSurfK -= del2K;
    }

  }
  else if ( dir == "k" ){ //split along k-plane

    int del1I = 0;
    int del1J = 0;
    int del2I = 0;
    int del2J = 0;

    for ( int ii = 0; ii < (*this).NumSurfaces(); ii++ ){
      if ( ii >= (*this).NumSurfI() + (*this).NumSurfJ() ){ //k-surface
	if ( (*this).GetKMax(ii) == 1){ //lower k surface
	  //no change to lower bc at lower k surface

	  //at lower k surface, upper bc is now interface
	  int tag = 6000 + numBlk; //lower surface matches with upper surface
	  bound2.surfs[ii].bcType = "interblock";           //bctype
	  bound2.surfs[ii].data[4] = (*this).GetKMin(ii);   //kmin
	  bound2.surfs[ii].data[5] = (*this).GetKMax(ii);   //kmax
	  bound2.surfs[ii].data[6] = tag;                   //tag
	}
	else {
	  //at upper k surface, lower bc is now interface
	  int tag = 5000 + newBlkNum; //upper surface matches with lower surface
	  bound1.surfs[ii].bcType = "interblock";           //bctype
	  bound1.surfs[ii].data[4] = indNG;                 //kmin
	  bound1.surfs[ii].data[5] = indNG;                 //kmax
	  bound1.surfs[ii].data[6] = tag;                   //tag

	  //at upper k surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.surfs[ii].data[4] = (*this).GetKMax(ii) - indNG + 1;   //kmin
	  bound2.surfs[ii].data[5] = (*this).GetKMax(ii) - indNG + 1;   //kmax

	  //at upper k surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	  if ( (*this).GetBCTypes(ii) == "interblock" ){
	    alteredSurf.push_back((*this).GetSurface(ii));
	  }

	}
      }
      else { //i-surface or j-surface

	//at i/j surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	if ( (*this).GetKMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.surfs[ii].data[4] = (*this).GetKMax(ii) - indNG + 1;      //kmin
	  bound2.surfs[ii].data[5] = (*this).GetKMax(ii) - indNG + 1;      //kmax
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del1I++;
	  }
	  else{ //j-surface
	    del1J++;
	  }
	}
	else if ( (*this).GetKMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.surfs[ii].data[5] = indNG;                              //kmax
	  bound2.surfs[ii].data[4] = 1;                                  //kmin
	  bound2.surfs[ii].data[5] = (*this).GetKMax(ii) - indNG + 1;    //kmax
	}
	else{ //this surface is only present in the lower split
	  del2.push_back(ii);
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del2I++;
	  }
	  else{ //j-surface
	    del2J++;
	  }
	}
      }
    }

    //delete unnecessary boundaries
    for ( unsigned int ii = 0; ii < del1.size(); ii++ ){
      bound1.surfs.erase(bound1.surfs.begin() + del1[ii]);
      bound1.numSurfI -= del1I;
      bound1.numSurfJ -= del1J;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.surfs.erase(bound2.surfs.begin() + del2[ii]);
      bound2.numSurfI -= del2I;
      bound2.numSurfJ -= del2J;
    }

  }
  else{
    cerr << "ERROR: Error in boundaryCondition::Split(). Direction " << dir << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }

  (*this) = bound1;
  aSurf = alteredSurf;
  return bound2;

}


void boundaryConditions::DependentSplit(const boundarySurface &surf, const plot3dBlock &part, const plot3dBlock &self, const int &sblk, 
					const string &dir, const int &ind, const int &lblk, const int &ublk ){

  // surf -- boundarySurface of partner block
  // part -- plot3dBlock that surf is assigned to
  // self -- plot3dBlock that (*this) is assigned to
  // sblk -- block number of self
  // dir -- direction that partner split was in
  // ind -- index of split
  // lblk -- lower block number in partner split
  // ublk -- upper block number in partner split

  bool border[4] = {false, false, false, false};
  patch partner(surf, part, lblk, border);

  for ( int ii = 0; ii < (*this).NumSurfaces(); ii++ ){

    patch candidate((*this).GetSurface(ii), self, sblk, border);

    interblock match(candidate, partner);
    if ( match.TestPatchMatch(candidate, partner) ){ //match found

      boundarySurface lowSurf = (*this).GetSurface(ii);

      string candDir;
      int candInd;
      if (match.Orientation() == 1){ //same orientation

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}
	candInd = ind;

      }
      else if (match.Orientation() == 2){ //D1/D2 swapped

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}
	candInd = ind;

      }
      else if (match.Orientation() == 3){ //D1 reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = surf.Max1() -1 - ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }
      else if (match.Orientation() == 4){ //D1/D2 swapped, D1 reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = surf.Max1() - 1 - ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }
      else if (match.Orientation() == 5){ //D1/D2 swapped, D2 reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = surf.Max2() - 1 - ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }
      else if (match.Orientation() == 6){ //D2 reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = surf.Max2() - 1 - ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }
      else if (match.Orientation() == 7){ //D1/D2 swapped and reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = surf.Max1() - 1 - ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = surf.Max2() - 1 - ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }
      else{ //D1/D2 reversed

	if ( surf.Direction1() == dir ){ //split was in direction 1 of partner, needs to be direction 1 of candidate
	  candDir = lowSurf.Direction1();
	  candInd = surf.Max1() - 1 - ind;
	}
	else if ( surf.Direction2() == dir ){ //split was in direction 2 of partner, needs to be direction 2 of candidate
	  candDir = lowSurf.Direction2();
	  candInd = surf.Max2() - 1 - ind;
	}
	else if ( surf.Direction3() == dir ){ //split was in direction 3 of partner, needs to be direction 3 of candidate
	  candDir = lowSurf.Direction3();
	  candInd = ind;
	}
	else{
	  cerr<< "ERROR: Error in boundaryConditions::DependentSplit(). Direction " << dir << " is not recognized." << endl;
	  cerr << "Please choose i, j, or k." << endl;
	  exit(0);
	}

      }

      cout << "unsplit surface: " << lowSurf << endl;
      cout << "splitting with dir, ind, orientation: " << candDir << ", " << candInd << ", " << match.Orientation() << endl;

      //split matched surface
      bool split = false;
      boundarySurface upSurf = lowSurf.Split(candDir, candInd, lblk, ublk, split, match.Orientation()); 

      cout << "split surfaces" << endl;
      cout << lowSurf << endl;
      cout << upSurf << endl;


      //assign boundarySurface back into boundaryConditions, if surface wasn't split partner block was updated
      (*this).surfs[ii] = lowSurf;

      //if surface was split, insert it into the vector of boundarySurfaces and adjust the surface numbers
      if (split){ //boundary surface was split, insert new surface into boundary conditions
	(*this).surfs.insert( (*this).surfs.begin() + ii, upSurf );
	if ( upSurf.SurfaceType() <= 2 ) { //i-surface
	  (*this).numSurfI++;
	}
	else if ( upSurf.SurfaceType() <= 4 ){ //j-surface
	  (*this).numSurfJ++;
	}
	else{
	  (*this).numSurfK++;
	}
      }

      break;
    }
  }
}


/* Member function to join 2 boundaryConditions. It assumes that the calling instance is the "lower" boundary condition and the input instance
is the "upper" boundary condition.
*/
void boundaryConditions::Join( const boundaryConditions &bc, const string &dir, vector<boundarySurface> &aSurf ){

  vector<boundarySurface> alteredSurf;

  if ( dir == "i" ){ //split along i-plane

    //total number of i surfaces in joined block will be equal to the lower i surfaces from lower bc plus the upper i surfaces from the upper bc
    int numI = 0;
    for ( int ii = 0; ii < (*this).NumSurfI(); ii++ ){
      if ( (*this).GetIMax(ii) == 1 ){ //lower i surface
	numI++;
      }
    }
    for ( int ii = 0; ii < bc.NumSurfI(); ii++ ){
      if ( bc.GetIMax(ii) != 1 ){ //upper i surface
	numI++;
      }
    }

    //total number of j surfaces in joined block will be equal to all j surfaces from lower bc plus the j surfaces from the upper bc that only reside in the upper bc
    int numJ = (*this).NumSurfJ() + bc.NumSurfJ();

    //total number of k surfaces in joined block will be equal to all k surfaces from lower bc plus the k surfaces from the upper bc that only reside in the upper bc
    int numK = (*this).NumSurfK() + bc.NumSurfK();

    //initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    int cc = 0; //boundary condition counter

    //insert all i lower surfaces from lower bc
    int lowDimI = 0;
    for ( int ii = 0; ii < (*this).NumSurfI(); ii++ ){
      if ( (*this).GetIMax(ii) == 1 ){ //lower i surface
	newBC.surfs[cc] = (*this).surfs[ii];
	cc++;
      }
      else{
	lowDimI = (*this).GetIMax(ii); //maximum i-dimension for lower bc
      }
    }
    //insert all i upper surfaces from upper bc
    for ( int ii = 0; ii < bc.NumSurfI(); ii++ ){
      if ( bc.GetIMax(ii) != 1 ){ //upper i surface
	//at upper i surface, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	//adjust i coordinates for join
	boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii) + lowDimI - 1, bc.GetIMax(ii) + lowDimI - 1, bc.GetJMin(ii), bc.GetJMax(ii), 
			       bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

	newBC.surfs[cc] = bcSurf;
	cc++;
      }
    }

    //insert all j surfaces from lower and upper bcs
    for ( int ii = (*this).NumSurfI(); ii < (*this).NumSurfI() + (*this).NumSurfJ(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++ ){
      //at j surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust i coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii) + lowDimI - 1, bc.GetIMax(ii) + lowDimI - 1, bc.GetJMin(ii), bc.GetJMax(ii), 
			     bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    //insert all k surfaces from lower and upper bcs
    for ( int ii = (*this).NumSurfI() + (*this).NumSurfJ(); ii < (*this).NumSurfaces(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++ ){
      //at k surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust i coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii) + lowDimI - 1, bc.GetIMax(ii) + lowDimI - 1, bc.GetJMin(ii), bc.GetJMax(ii), 
			     bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;
  }
  else if ( dir == "j" ){ //split along j-plane

    //total number of j surfaces in joined block will be equal to the lower j surfaces from lower bc plus the upper j surfaces from the upper bc
    int numJ = 0;
    for ( int ii = (*this).NumSurfI(); ii < (*this).NumSurfI() + (*this).NumSurfJ(); ii++ ){
      if ( (*this).GetJMax(ii) == 1 ){ //lower j surface
	numJ++;
      }
    }
    for ( int ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++ ){
      if ( bc.GetJMax(ii) != 1 ){ //upper j surface
	numJ++;
      }
    }

    //total number of i surfaces in joined block will be equal to all i surfaces from lower bc plus the i surfaces from the upper bc that only reside in the upper bc
    int numI = (*this).NumSurfI() + bc.NumSurfI();

    //total number of k surfaces in joined block will be equal to all k surfaces from lower bc plus the k surfaces from the upper bc that only reside in the upper bc
    int numK = (*this).NumSurfK() + bc.NumSurfK();

    //initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    int cc = numI; //boundary condition counter

    //insert all j lower surfaces from lower bc
    int lowDimJ = 0;
    for ( int ii = (*this).NumSurfI(); ii < (*this).NumSurfI() + (*this).NumSurfJ(); ii++ ){
      if ( (*this).GetJMax(ii) == 1 ){ //lower j surface
	newBC.surfs[cc] = (*this).surfs[ii];
	cc++;
      }
      else{
	lowDimJ = (*this).GetJMax(ii); //maximum j-dimension for lower bc
      }
    }
    //insert all j upper surfaces from upper bc
    for ( int ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++ ){
      if ( bc.GetJMax(ii) != 1 ){ //upper j surface
	//at j upper surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	//adjust j coordinates for join
	boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii) + lowDimJ - 1, bc.GetJMax(ii) + lowDimJ - 1, 
			       bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

	newBC.surfs[cc] = bcSurf;
	cc++;
      }
    }

    cc = 0;
    //insert all i surfaces from lower and upper bcs
    for ( int ii = 0; ii < (*this).NumSurfI(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = 0; ii < bc.NumSurfI(); ii++ ){
      //at i surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust j coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii) + lowDimJ - 1, bc.GetJMax(ii) + lowDimJ - 1, 
			     bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    cc = numI + numJ;
    //insert all k surfaces from lower and upper bcs
    for ( int ii = (*this).NumSurfI() + (*this).NumSurfJ(); ii < (*this).NumSurfaces(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++ ){
      //at k surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust j coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii) + lowDimJ - 1, bc.GetJMax(ii) + lowDimJ - 1, 
			     bc.GetKMin(ii), bc.GetKMax(ii), bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;

  }
  else if ( dir == "k" ){ //split along k-plane

    //total number of k surfaces in joined block will be equal to the lower k surfaces from lower bc plus the upper k surfaces from the upper bc
    int numK = 0;
    for ( int ii = (*this).NumSurfI() + (*this).NumSurfJ(); ii < (*this).NumSurfaces(); ii++ ){
      if ( (*this).GetKMax(ii) == 1 ){ //lower k surface
	numK++;
      }
    }
    for ( int ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++ ){
      if ( bc.GetKMax(ii) != 1 ){ //upper k surface
	numK++;
      }
    }

    //total number of i surfaces in joined block will be equal to all i surfaces from lower bc plus the i surfaces from the upper bc that only reside in the upper bc
    int numI = (*this).NumSurfI() + bc.NumSurfI();

    //total number of j surfaces in joined block will be equal to all j surfaces from lower bc plus the j surfaces from the upper bc that only reside in the upper bc
    int numJ = (*this).NumSurfJ() + bc.NumSurfJ();

    //initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    int cc = numI + numJ; //boundary condition counter

    //insert all k lower surfaces from lower bc
    int lowDimK = 0;
    for ( int ii = (*this).NumSurfI() + (*this).NumSurfJ(); ii < (*this).NumSurfaces(); ii++ ){
      if ( (*this).GetKMax(ii) == 1 ){ //lower k surface
	newBC.surfs[cc] = (*this).surfs[ii];
	cc++;
      }
      else{
	lowDimK = (*this).GetKMax(ii); //maximum k-dimension for lower bc
      }
    }
    //insert all k upper surfaces from upper bc
    for ( int ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++ ){
      if ( bc.GetKMax(ii) != 1 ){ //upper k surface
	//at upper k surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
	if ( (*this).GetBCTypes(ii) == "interblock" ){
	  alteredSurf.push_back((*this).GetSurface(ii));
	}

	//adjust k coordinates for join
	boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii), bc.GetJMax(ii), 
			       bc.GetKMin(ii) + lowDimK - 1, bc.GetKMax(ii) + lowDimK - 1, bc.GetTag(ii) );

	newBC.surfs[cc] = bcSurf;
	cc++;
      }
    }

    cc = 0;
    //insert all i surfaces from lower and upper bcs
    for ( int ii = 0; ii < (*this).NumSurfI(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = 0; ii < bc.NumSurfI(); ii++ ){
      //at i surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust k coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii), bc.GetJMax(ii), 
			     bc.GetKMin(ii) + lowDimK - 1, bc.GetKMax(ii) + lowDimK - 1, bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    cc = numI;
    //insert all j surfaces from lower and upper bcs
    for ( int ii = (*this).NumSurfI(); ii < (*this).NumSurfI() + (*this).NumSurfJ(); ii++ ){
      newBC.surfs[cc] = (*this).surfs[ii];
      cc++;
    }
    for ( int ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++ ){
      //at j surface for upper block, if bc is interblock, store boundarySurface because partner block BC will need to be updated
      if ( (*this).GetBCTypes(ii) == "interblock" ){
	alteredSurf.push_back((*this).GetSurface(ii));
      }

      //adjust k coordinates for join
      boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii), bc.GetIMax(ii), bc.GetJMin(ii), bc.GetJMax(ii), 
			     bc.GetKMin(ii) + lowDimK - 1, bc.GetKMax(ii) + lowDimK - 1, bc.GetTag(ii) );

      newBC.surfs[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;

  }
  else{
    cerr << "ERROR: Error in procBlock::Join(). Direction " << dir << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }

  aSurf = alteredSurf;
}


//constructor when passed no arguements
patch::patch(){
  //initialize all variables to zero
  vector3d<double> zero(0.0, 0.0, 0.0);
  origin = zero;
  corner1 = zero;
  corner2 = zero;
  corner12 = zero;

  block = 0;
  boundary = 0;
  d1Start = 0;
  d1End = 0;
  d2Start = 0;
  d2End = 0;
  constSurf = 0;
  rank = 0;
  localBlock = 0;
  interblockBorder[0] = false;
  interblockBorder[1] = false;
  interblockBorder[2] = false;
  interblockBorder[3] = false;
}

//constructor with arguements passed
patch::patch( const int &bound, const int &b, const int &d1s, const int &d1e, const int &d2s, const int &d2e, const int &d3s, 
	      const int &d3e, const plot3dBlock &blk, const int &r, const int &l, const bool (&border)[4]){
  // bound -- boundary number which patch is on (1-6)
  // b -- parent block number
  // d1s -- direction 1 starting index
  // d1e -- direction 1 ending index
  // d2s -- direction 2 starting index
  // d2e -- direction 2 ending index
  // d3s -- direction 3 surface index (constant surface that patch is on)

  boundary = bound;
  block = b;
  rank = r;
  localBlock = l;
  interblockBorder[0] = border[0];
  interblockBorder[1] = border[1];
  interblockBorder[2] = border[2];
  interblockBorder[3] = border[3];

  if ( bound == 1 || bound == 2 ){ //patch on i-surface - dir1 = j, dir2 = k
    d1Start = d2s;
    d1End = d2e;
    d2Start = d3s;
    d2End = d3e;
    constSurf = d1s;

    //get corner points
    //origin at jmin, kmin
    int loc = GetLoc1D(constSurf, d1Start, d2Start, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at jmax, kmin
    loc = GetLoc1D(constSurf, d1End, d2Start, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at jmin, kmax
    loc = GetLoc1D(constSurf, d1Start, d2End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at jmax, kmax
    loc = GetLoc1D(constSurf, d1End, d2End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else if ( bound == 3 || bound == 4 ){ //patch on j-surface - dir1 = k, dir2 = i
    d1Start = d3s;
    d1End = d3e;
    d2Start = d1s;
    d2End = d1e;
    constSurf = d2s;

    //get corner points
    //origin at kmin, imin
    int loc = GetLoc1D(d2Start, constSurf, d1Start, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at kmax, imin
    loc = GetLoc1D(d2Start, constSurf, d1End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at kmin, imax
    loc = GetLoc1D(d2End, constSurf, d1Start, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at kmax, imax
    loc = GetLoc1D(d2End, constSurf, d1End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else if ( bound == 5 || bound == 6 ){ //patch on k-surface - dir1 = i, dir2 = j
    d1Start = d1s;
    d1End = d1e;
    d2Start = d2s;
    d2End = d2e;
    constSurf = d3s;

    //get corner points
    //origin at imin, jmin
    int loc = GetLoc1D(d1Start, d2Start, constSurf, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at imax, jmin
    loc = GetLoc1D(d1End, d2Start, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at imin, jmax
    loc = GetLoc1D(d1Start, d2End, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at imax, jmax
    loc = GetLoc1D(d1End, d2End, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else{
    cerr << "ERROR: Error in patch::patch(). Boundary surface " << bound << " is not recognized!" << endl;
    cerr << "Choose an integer between 1-6." << endl;
    exit(0);
  }

}

//constructor with arguements passed
patch::patch( const boundarySurface &surf, const plot3dBlock &blk, const int &bNum, const bool (&border)[4], int r, int l){

  boundary = surf.SurfaceType();
  block = bNum;
  rank = r;
  localBlock = l;
  interblockBorder[0] = border[0];
  interblockBorder[1] = border[1];
  interblockBorder[2] = border[2];
  interblockBorder[3] = border[3];

  if ( boundary == 1 || boundary == 2 ){ //patch on i-surface - dir1 = j, dir2 = k
    d1Start = surf.JMin() - 1;
    d1End = surf.JMax() - 1;
    d2Start = surf.KMin() - 1;
    d2End = surf.KMax() - 1;
    constSurf = surf.IMax() - 1;

    //get corner points
    //origin at jmin, kmin
    int loc = GetLoc1D(constSurf, d1Start, d2Start, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at jmax, kmin
    loc = GetLoc1D(constSurf, d1End, d2Start, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at jmin, kmax
    loc = GetLoc1D(constSurf, d1Start, d2End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at jmax, kmax
    loc = GetLoc1D(constSurf, d1End, d2End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else if ( boundary == 3 || boundary == 4 ){ //patch on j-surface - dir1 = k, dir2 = i
    d1Start = surf.KMin() - 1;
    d1End = surf.KMax() - 1;
    d2Start = surf.IMin() - 1;
    d2End = surf.IMax() - 1;
    constSurf = surf.JMax() - 1;

    //get corner points
    //origin at kmin, imin
    int loc = GetLoc1D(d2Start, constSurf, d1Start, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at kmax, imin
    loc = GetLoc1D(d2Start, constSurf, d1End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at kmin, imax
    loc = GetLoc1D(d2End, constSurf, d1Start, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at kmax, imax
    loc = GetLoc1D(d2End, constSurf, d1End, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else if ( boundary == 5 || boundary == 6 ){ //patch on k-surface - dir1 = i, dir2 = j
    d1Start = surf.IMin();
    d1End = surf.IMax();
    d2Start = surf.JMin();
    d2End = surf.JMax();
    constSurf = surf.KMax();

    //get corner points
    //origin at imin, jmin
    int loc = GetLoc1D(d1Start, d2Start, constSurf, blk.NumI(), blk.NumJ()); 
    vector3d<double> temp(blk.XLoc(loc),
			  blk.YLoc(loc),
			  blk.ZLoc(loc));
    origin = temp;

    //corner1 at imax, jmin
    loc = GetLoc1D(d1End, d2Start, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner1 = temp;

    //corner2 at imin, jmax
    loc = GetLoc1D(d1Start, d2End, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner2 = temp;

    //corner12 at imax, jmax
    loc = GetLoc1D(d1End, d2End, constSurf, blk.NumI(), blk.NumJ()); 
    temp.SetX(blk.XLoc(loc));
    temp.SetY(blk.YLoc(loc));
    temp.SetZ(blk.ZLoc(loc));
    corner12 = temp;

  }
  else{
    cerr << "ERROR: Error in patch::patch(). Boundary surface " << boundary << " is not recognized!" << endl;
    cerr << "Choose an integer between 1-6." << endl;
    exit(0);
  }

}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const patch &p){

  os << "Boundary: " << p.Boundary() << endl;
  os << "Block: " << p.Block() << endl;
  os << "Direction 1 Start: " << p.Dir1Start() << endl;
  os << "Direction 1 End: " << p.Dir1End() << endl;
  os << "Direction 2 Start: " << p.Dir2Start() << endl;
  os << "Direction 2 End: " << p.Dir2End() << endl;
  os << "Constant Surface: " << p.ConstSurface() << endl;
  os << "Rank: " << p.Rank() << endl;
  os << "Borders Interblock: " << p.Dir1StartInterBorder() << ", " << p.Dir1EndInterBorder() << ", " << p.Dir2StartInterBorder() << ", " << p.Dir2EndInterBorder() << endl;
  os << "Origin: " << p.origin << endl;
  os << "Corner 1: " << p.corner1 << endl;
  os << "Corner 2: " << p.corner2 << endl;
  os << "Corner 12: " << p.corner12 << endl;

  return os;
}

void boundaryConditions::PackBC( char *(&sendBuffer), const int &sendBufSize, int &position)const{


  //get string lengths for each boundary condition to be sent, so processors unpacking know how much data to unpack for each string
  vector<int> strLength((*this).NumSurfaces());
  for ( unsigned int jj = 0; jj < strLength.size(); jj++ ){
    strLength[jj] = (*this).GetBCTypes(jj).size() + 1; //+1 for c_str end character
  }

  MPI_Pack(&(*this).numSurfI, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numSurfJ, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numSurfK, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  for ( int jj = 0; jj < (*this).NumSurfaces(); jj++ ){
    MPI_Pack(&(*this).surfs[jj].data[0], 7, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  }
  MPI_Pack(&strLength[0], strLength.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  for ( int jj = 0; jj < (*this).NumSurfaces(); jj++ ){
    MPI_Pack((*this).surfs[jj].bcType.c_str(), (*this).surfs[jj].bcType.size()+1, MPI_CHAR, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD); //+1 for c_str end character
  }

}

void boundaryConditions::UnpackBC( char *(&recvBuffer), const int &recvBufSize, int &position){


  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfI, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of i-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfJ, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of j-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfK, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of k-surfaces

  (*this).ResizeVecs((*this).NumSurfaces());
  vector<int> strLength((*this).NumSurfaces());

  //unpack boundary condition data into appropriate vectors
  for ( int jj = 0; jj < (*this).NumSurfaces(); jj++ ){
    MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).surfs[jj].data[0], 7, MPI_INT, MPI_COMM_WORLD); //unpack bc surfaces
  }
  MPI_Unpack(recvBuffer, recvBufSize, &position, &strLength[0], strLength.size(), MPI_INT, MPI_COMM_WORLD); //unpack string sizes
  //unpack boundary condition names
  for ( unsigned int jj = 0; jj < strLength.size(); jj++ ){
    char *nameBuf = new char[strLength[jj]]; //allocate buffer to store BC name
    MPI_Unpack(recvBuffer, recvBufSize, &position, &nameBuf[0], strLength[jj], MPI_CHAR, MPI_COMM_WORLD); //unpack bc types
    string bcName(nameBuf, strLength[jj] - 1); //create string of bc name (-1 to exclude c_str end character)
    (*this).surfs[jj].bcType = bcName;
    delete [] nameBuf; //deallocate bc name buffer
  }

}

//constructor when passed no arguements
boundarySurface::boundarySurface(){
  bcType = "undefined";
  data[0] = 0;
  data[1] = 0;
  data[2] = 0;
  data[3] = 0;
  data[4] = 0;
  data[5] = 0;
  data[6] = 0;
}

boundarySurface::boundarySurface(const string &name, const int &imin, const int &imax, const int &jmin, const int &jmax, const int &kmin, const int &kmax, const int &tag){
  bcType = name;
  data[0] = imin;
  data[1] = imax;
  data[2] = jmin;
  data[3] = jmax;
  data[4] = kmin;
  data[5] = kmax;
  data[6] = tag;
}

int boundarySurface::SurfaceType() const {

  int surf = 0;

  if ( data[0] == data[1] ){ //i-surface
    if ( data[1] == 1 ){ //lower surface
      surf = 1;
    }
    else{ //upper surface
      surf = 2;
    }
  }
  else if ( data[2] == data[3] ){ //j-surface
    if ( data[3] == 1){ //lower surface
      surf = 3;
    }
    else{ //upper surface
      surf = 4;
    }
  }
  else if ( data[4] == data[5] ){ //k-surface
    if ( data[5] == 1){ //lower surface
      surf = 5;
    }
    else{ //upper surface
      surf = 6;
    }
  }
  else{
    cerr << "ERROR: Error in boundarySurface::SurfaceType(). Surface is defined incorrectly, it is neither an i, j, or k surface." << endl;
    cerr << (*this) << endl;
    exit(0);
  }

  return surf;
}

int boundarySurface::PartnerBlock() const {

  if ( (*this).bcType != "interblock" ){ 
    cerr << "ERROR: Partner blocks are only associated with interblock boundaries. Current boundary is " << (*this).bcType << endl;
    exit(0);
  }

  int subtract = (*this).PartnerSurface() * 1000;
  return (*this).Tag() - subtract;
}

int boundarySurface::PartnerSurface() const {

  if ( (*this).bcType != "interblock" ){ 
    cerr << "ERROR: Partner blocks are only associated with interblock boundaries. Current boundary is " << (*this).bcType << endl;
    exit(0);
  }

  int surf = 0;

  if ( (*this).Tag() < 2000 ){
    surf = 1; //i-lower surface
  }
  else if ( (*this).Tag() < 3000 ){
    surf = 2; //i-upper surface
  }
  else if ( (*this).Tag() < 4000 ){
    surf = 3; //j-lower surface
  }
  else if ( (*this).Tag() < 5000 ){
    surf = 4; //j-upper surface
  }
  else if ( (*this).Tag() < 6000 ){
    surf = 5; //k-lower surface
  }
  else if ( (*this).Tag() < 7000 ){
    surf = 6; //k-upper surface
  }
  else {
    cerr << "ERROR: Error in boundarySurface::PartnerSurface(). Tag does not fit in range. Tag must be between 1000 and 6999." << endl;
    cerr << (*this) << endl;
    exit(0);
  }


  return surf;
}

string boundarySurface::Direction1()const{

  string dir;
  if ( (*this).SurfaceType() <= 2 ){ //dir 3 is i, dir 1 is j, dir 2 is k
    dir = "j";
  }
  else if ( (*this).SurfaceType() <= 4){ //dir 3 is j, dir 1 is k, dir 2 is i
    dir = "k";
  }
  else{ //dir 3 is k, dir 1 is i, dir 2 is j
    dir = "i";
  }

  return dir;
}

string boundarySurface::Direction2()const{

  string dir;
  if ( (*this).SurfaceType() <= 2 ){ //dir 3 is i, dir 1 is j, dir 2 is k
    dir = "k";
  }
  else if ( (*this).SurfaceType() <= 4){ //dir 3 is j, dir 1 is k, dir 2 is i
    dir = "i";
  }
  else{ //dir 3 is k, dir 1 is i, dir 2 is j
    dir = "j";
  }

  return dir;
}

string boundarySurface::Direction3()const{

  string dir;
  if ( (*this).SurfaceType() <= 2 ){ //dir 3 is i, dir 1 is j, dir 2 is k
    dir = "i";
  }
  else if ( (*this).SurfaceType() <= 4){ //dir 3 is j, dir 1 is k, dir 2 is i
    dir = "j";
  }
  else{ //dir 3 is k, dir 1 is i, dir 2 is j
    dir = "k";
  }

  return dir;
}


int boundarySurface::Max1()const{
  int m = 0;
  if ( (*this).Direction1() == "i" ){
    m = (*this).IMax();
  }
  else if ( (*this).Direction1() == "j" ){
    m = (*this).JMax();
  }
  else{
    m = (*this).KMax();
  }
  return m;
}

int boundarySurface::Min1()const{
  int m = 0;
  if ( (*this).Direction1() == "i" ){
    m = (*this).IMin();
  }
  else if ( (*this).Direction1() == "j" ){
    m = (*this).JMin();
  }
  else{
    m = (*this).KMin();
  }
  return m;
}

int boundarySurface::Max2()const{
  int m = 0;
  if ( (*this).Direction2() == "i" ){
    m = (*this).IMax();
  }
  else if ( (*this).Direction2() == "j" ){
    m = (*this).JMax();
  }
  else{
    m = (*this).KMax();
  }
  return m;
}

int boundarySurface::Min2()const{
  int m = 0;
  if ( (*this).Direction2() == "i" ){
    m = (*this).IMin();
  }
  else if ( (*this).Direction2() == "j" ){
    m = (*this).JMin();
  }
  else{
    m = (*this).KMin();
  }
  return m;
}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const boundarySurface &bcSurf){

  os << bcSurf.bcType << "   " << bcSurf.IMin() << "   " << bcSurf.IMax() << "   " << bcSurf.JMin() << "   " << bcSurf.JMax() << "   "
     << bcSurf.KMin() << "   " << bcSurf.KMax() << "   " << bcSurf.Tag();

  return os;
}

void boundarySurface::UpdateTagForSplitJoin(const int &nBlk){
  (*this).data[6] = (*this).PartnerSurface() * 1000 + nBlk;
}

boundarySurface boundarySurface::Split(const string &dir, const int &ind, const int &lBlk, const int &uBlk, bool &split, int orientation){

  int indNG = ind + 1; //+1 because boundaries start at 1, not 0

  boundarySurface surf1 = (*this);
  boundarySurface surf2 = (*this);

  split = true;
  bool isReversed = surf1.SplitDirectionIsReversed(dir, orientation); //surf1 and surf2 have same orientation, so if reversed for 1, reversed for 2

  if ( dir == "i" ){ //split along i-plane

    if ( (*this).SurfaceType() == 1 || (*this).SurfaceType() == 2 ){ //cannot split an i-surface along i-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    }
    else{ //j or k surface

      if ( (*this).IMin() >= indNG ){ //this surface is only present in the upper split
	surf1.data[0] = (*this).IMax() - indNG + 1;    //imin
	surf1.data[1] = (*this).IMax() - indNG + 1;    //imax
	isReversed ? surf1.UpdateTagForSplitJoin(lBlk) : surf1.UpdateTagForSplitJoin(uBlk);
	split = false;
      }
      else if ( (*this).IMax() >= indNG ){ //this surface straddles the split
	surf2.data[0] = indNG;                            //imin
	isReversed ? surf2.UpdateTagForSplitJoin(lBlk) : surf2.UpdateTagForSplitJoin(uBlk);

	surf1.data[1] = indNG;                            //imax
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
      }
      else{ //this surface is only present in the lower split
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
	split = false;
      }
    }

  }
  else if ( dir == "j" ){ //split along j-plane

    if ( (*this).SurfaceType() == 3 || (*this).SurfaceType() == 4 ){ //cannot split a j-surface along j-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    }
    else{ //i or k surface

      if ( (*this).JMin() >= indNG ){ //this surface is only present in the upper split
	surf1.data[2] = (*this).JMax() - indNG + 1;    //jmin
	surf1.data[3] = (*this).JMax() - indNG + 1;    //jmax
	isReversed ? surf1.UpdateTagForSplitJoin(lBlk) : surf1.UpdateTagForSplitJoin(uBlk);
	split = false;
      }
      else if ( (*this).JMax() >= indNG ){ //this surface straddles the split
	surf2.data[2] = indNG;                            //jmin
	isReversed ? surf2.UpdateTagForSplitJoin(lBlk) : surf2.UpdateTagForSplitJoin(uBlk);

	surf1.data[3] = indNG;                            //jmax
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
      }
      else{ //this surface is only present in the lower split
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
	split = false;
      }
    }

  }
  else if ( dir == "k" ){ //split along k-plane

    if ( (*this).SurfaceType() == 5 || (*this).SurfaceType() == 6 ){ //cannot split a k-surface along k-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    }
    else{ //i or j surface

      if ( (*this).KMin() >= indNG ){ //this surface is only present in the upper split
	surf1.data[4] = (*this).KMax() - indNG + 1;    //kmin
	surf1.data[5] = (*this).KMax() - indNG + 1;    //kmax
	isReversed ? surf1.UpdateTagForSplitJoin(lBlk) : surf1.UpdateTagForSplitJoin(uBlk);
	split = false;
      }
      else if ( (*this).KMax() >= indNG ){ //this surface straddles the split
	surf2.data[4] = indNG;                            //kmin
	isReversed ? surf2.UpdateTagForSplitJoin(lBlk) : surf2.UpdateTagForSplitJoin(uBlk);

	surf1.data[5] = indNG;                            //kmax
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
      }
      else{ //this surface is only present in the lower split
	isReversed ? surf1.UpdateTagForSplitJoin(uBlk) : surf1.UpdateTagForSplitJoin(lBlk);
	split = false;
      }
    }

  }
  else{
    cerr << "ERROR: Error in boundarySurface::Split(). Direction " << dir << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }

  (*this) = surf1;
  return surf2;

}

bool boundarySurface::SplitDirectionIsReversed( const string &dir, const int &orientation) const {

  bool isReversed;

  //find out if split direction is 1, 2, or 3
  if ( (*this).Direction1() == dir ){ //split direction is direction 1 - reverse if dir 1 is reversed (relative to partner, taking into account D1/D2 swap)
    isReversed = ( orientation == 3 || orientation == 5 || orientation == 7 || orientation == 8 ) ? true : false;
  }
  else if ( (*this).Direction2() == dir ){ //split direction is direction 2 - reverse if dir 2 is reversed (relative to partner, taking into account D1/D2 swap)
    isReversed = ( orientation == 4 || orientation == 6 || orientation == 7 || orientation == 8 ) ? true : false;
  }
  else if ( (*this).Direction3() == dir ){ //split direction is direction 3 - no need to reverse
    isReversed = false;
  }
  else {
    cerr << "ERROR: Error in boundarySurface::SplitDirectionIsReversed(). Direction " << dir << " does not match i, j, or k!" << endl;
    exit(0);
  }

  return isReversed;
}

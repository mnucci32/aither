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

  vector<int> dumVecInt(length,0);
  vector<string> dumVecStr(length,"undefined");

  //assign dummy vectors to class variables
  bcTypes = dumVecStr;
  iMin = dumVecInt;
  iMax = dumVecInt;
  jMin = dumVecInt;
  jMax = dumVecInt;
  kMin = dumVecInt;
  kMax = dumVecInt;
  tag = dumVecInt;

}

//constructor when passed number of i, j, k surfaces
boundaryConditions::boundaryConditions( const int &i, const int &j, const int &k){
  numSurfI = i;
  numSurfJ = j;
  numSurfK = k;
  int length = numSurfI + numSurfJ + numSurfK;

  vector<int> dumVecInt(length,0);
  vector<string> dumVecStr(length,"undefined");

  //assign dummy vectors to class variables
  bcTypes = dumVecStr;
  iMin = dumVecInt;
  iMax = dumVecInt;
  jMin = dumVecInt;
  jMax = dumVecInt;
  kMin = dumVecInt;
  kMax = dumVecInt;
  tag = dumVecInt;

}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const boundaryConditions &bc){

  os << "Number of surfaces (I, J, K): " << bc.numSurfI << ", " << bc.numSurfJ << ", " << bc.numSurfK << endl;

  for ( unsigned int ii = 0; ii < bc.bcTypes.size(); ii++ ){
    os << bc.bcTypes[ii] << "   " << bc.iMin[ii] << "   " << bc.iMax[ii] << "   " << bc.jMin[ii] << "   " << bc.jMax[ii] << "   "
       << bc.kMin[ii] << "   " << bc.kMax[ii] << "   " << bc.tag[ii] << endl;

  }

  return os;
}

//operator to resize all of the vector components of the boundary conditions class
void boundaryConditions::ResizeVecs( const int &a){

  bcTypes.resize(a);
  iMin.resize(a);
  iMax.resize(a);
  jMin.resize(a);
  jMax.resize(a);
  kMin.resize(a);
  kMax.resize(a);
  tag.resize(a);

}

//operator to resize all of the vector components of the boundary conditions class
void boundaryConditions::ResizeVecs( const int &i, const int &j, const int &k){

  numSurfI = i;
  numSurfJ = j;
  numSurfK = k;

  bcTypes.resize(i + j + k);
  iMin.resize(i + j + k);
  iMax.resize(i + j + k);
  jMin.resize(i + j + k);
  jMax.resize(i + j + k);
  kMin.resize(i + j + k);
  kMax.resize(i + j + k);
  tag.resize(i + j + k);

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

  bcTypes[surfCounter] = tokens[0];
  iMin[surfCounter] = atoi(tokens[1].c_str());
  iMax[surfCounter] = atoi(tokens[2].c_str());
  jMin[surfCounter] = atoi(tokens[3].c_str());
  jMax[surfCounter] = atoi(tokens[4].c_str());
  kMin[surfCounter] = atoi(tokens[5].c_str());
  kMax[surfCounter] = atoi(tokens[6].c_str());
  tag[surfCounter] = atoi(tokens[7].c_str());

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
  os << "Orientation: " << bc.Orientation() << endl;

  return os;
}

//constructor to take in two patches and fill an interblock. The orientation is left at the default value 0.
interblock::interblock(const patch &p1, const patch &p2){
  // p1 -- patch 1
  // p2 -- patch 2

  //fill interblock
  rank[0] = 0; //default value is 0
  rank[1] = 0;

  block[0] = p1.Block();
  block[1] =p2.Block();

  localBlock[0] = 0; //default value is 0
  localBlock[1] = 0;

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

  orientation = 0; //default value (real values 1-6)
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
vector<interblock> GetInterblockBCs( const vector<boundaryConditions> &bc, const vector<plot3dBlock> &grid ){
  // bc -- vector of boundaryConditions for all blocks
  // grid -- vector of plot3Dblocks for entire computational mesh

  //isolate only the interblock BCs and their associated data from all of the BCs
  vector<vector<int> > isolatedInterblocks; //outer vector for each interblock BC, inner vector for information about interblock
  for ( unsigned int ii = 0; ii < bc.size(); ii++ ){ //loop over all blocks
    int numSurf = bc[ii].NumSurfI() + bc[ii].NumSurfJ() + bc[ii].NumSurfK(); //number of surfaces in block

    for ( int jj = 0; jj < numSurf; jj++ ){ //loop over number of surfaces in block

      if ( bc[ii].GetBCTypes(jj) == "interblock" ){ //if boundary condition is interblock, store data
	vector<int> temp (10,0); 
	temp[0] = ii;                                  //block number of bc
	
	//boundary number of bc (1-6)
	if ( jj < bc[ii].NumSurfI() ){ //i-surface
	  if ( bc[ii].GetIMin(jj) == 1) { //lower surface
	    temp[1] = 1;
	  }
	  else{ //upper surface
	    temp[1] = 2;
	  }
	}
	else if ( jj < bc[ii].NumSurfI() + bc[ii].NumSurfJ() ){ //j-surface
	  if ( bc[ii].GetJMin(jj) == 1) { //lower surface
	    temp[1] = 3;
	  }
	  else{ //upper surface
	    temp[1] = 4;
	  }
	}
	else{ //k-surface
	  if ( bc[ii].GetKMin(jj) == 1) { //lower surface
	    temp[1] = 5;
	  }
	  else{ //upper surface
	    temp[1] = 6;
	  }
	}

	//1 subtracted from indices because vectors start at 0.
	//these are grid point indices
	temp[2] = bc[ii].GetIMin(jj) - 1;              //i min of bc patch
	temp[3] = bc[ii].GetIMax(jj) - 1;              //i max of bc patch
	temp[4] = bc[ii].GetJMin(jj) - 1;              //j min of bc patch
	temp[5] = bc[ii].GetJMax(jj) - 1;              //j max of bc patch
	temp[6] = bc[ii].GetKMin(jj) - 1;              //k min of bc patch
	temp[7] = bc[ii].GetKMax(jj) - 1;              //k max of bc patch

	//determine block/boundary that interblock BC is supposed to match to
	int tag = bc[ii].GetTag(jj);
	int bound = -1;
	int blk = -1;
	if (tag >= 1000  && tag < 2000){ //at il boundary
	  bound = 1;
	  blk = tag - 1000;
	} 
	else if ( tag >= 2000 && tag < 3000){ //at iu boundary
	  bound = 2;
	  blk = tag - 2000;
	}
	else if ( tag >= 3000 && tag < 4000){ //at jl boundary
	  bound = 3;
	  blk = tag - 3000;
	}
	else if ( tag >= 4000 && tag < 5000){ //at ju boundary
	  bound = 4;
	  blk = tag - 4000;
	}
	else if ( tag >= 5000 && tag < 6000){ //at kl boundary
	  bound = 5;
	  blk = tag - 5000;
	}
	else if ( tag >= 6000 && tag < 7000){ //at ku boundary
	  bound = 6;
	  blk = tag - 6000;
	}
	else{
	  cerr << "ERROR: Error in boundaryConditions.cpp:GetInterblockBCs(). Not sure what to do with tag "
	       << tag << " from block " << ii << " at boundary " << jj << "." << endl;
	}

	temp[8] = bound;                                  //boundary of match
	temp[9] = blk;                                    //block of match

	isolatedInterblocks.push_back(temp); //add data to end of vector

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

      if ( isolatedInterblocks[ii][9] == isolatedInterblocks[jj][0] ) { //blocks between interblock BCs match
	if ( isolatedInterblocks[ii][8] == isolatedInterblocks[jj][1] ) { //boundary surfaces between interblock BCs match

	  //get current patch
	  patch cPatch(isolatedInterblocks[ii][1], isolatedInterblocks[ii][0], isolatedInterblocks[ii][2], isolatedInterblocks[ii][3], isolatedInterblocks[ii][4],
		       isolatedInterblocks[ii][5], isolatedInterblocks[ii][6], isolatedInterblocks[ii][7], grid[isolatedInterblocks[ii][0]]);

	  //get new patch (possible match)
	  patch nPatch(isolatedInterblocks[jj][1], isolatedInterblocks[jj][0], isolatedInterblocks[jj][2], isolatedInterblocks[jj][3], isolatedInterblocks[jj][4],
		       isolatedInterblocks[jj][5], isolatedInterblocks[jj][6], isolatedInterblocks[jj][7], grid[isolatedInterblocks[jj][0]]);

	  //test for match
	  interblock match(cPatch, nPatch);
	  if ( match.TestPatchMatch(cPatch, nPatch) ){ //match found
	    connections[ii/2] = match; //store interblock pair
	    swap(isolatedInterblocks[jj], isolatedInterblocks[ii+1]); //swap matched interblock BC to top portion of vector so it is not searched again
	  }

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
void interblock::GetAddressesMPI(MPI_Aint (&disp)[10])const{

  //get addresses of each field
  MPI_Get_address(&(*this).rank[0],       &disp[0]);
  MPI_Get_address(&(*this).block[0],      &disp[1]);
  MPI_Get_address(&(*this).localBlock[0], &disp[2]);
  MPI_Get_address(&(*this).boundary[0],   &disp[3]);
  MPI_Get_address(&(*this).d1Start[0],    &disp[4]);
  MPI_Get_address(&(*this).d1End[0],      &disp[5]);
  MPI_Get_address(&(*this).d2Start[0],    &disp[6]);
  MPI_Get_address(&(*this).d2End[0],      &disp[7]);
  MPI_Get_address(&(*this).constSurf[0],  &disp[8]);
  MPI_Get_address(&(*this).orientation,   &disp[9]);

}

/* Member function to split boundary conditions along a given direction at a given index. The calling instance retains the lower portion of the split,
and the returned instance is the upper portion
*/
boundaryConditions boundaryConditions::Split(const string &dir, const int &ind, const int &numBlk){

  int indNG = ind + 1; //+1 because boundaries start at 1, not 0

  boundaryConditions bound1 = (*this);
  boundaryConditions bound2 = (*this);

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
	  bound2.bcTypes[ii] = "interblock";
	  bound2.iMin[ii] = (*this).GetIMin(ii);
	  bound2.iMax[ii] = (*this).GetIMax(ii);
	  bound2.tag[ii] = tag;
	}
	else{ //upper surface
	  //at upper i surface, lower bc is now interface
	  int tag = 1000 + numBlk; //upper surface matches with lower surface
	  bound1.bcTypes[ii] = "interblock";
	  bound1.iMin[ii] = indNG;
	  bound1.iMax[ii] = indNG;
	  bound1.tag[ii] = tag;

	  //at upper i surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.iMin[ii] = (*this).GetIMax(ii) - indNG + 1;
	  bound2.iMax[ii] = (*this).GetIMax(ii) - indNG + 1;
	}
      }
      else { //j-surface or k-surface
	if ( (*this).GetIMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.iMin[ii] = (*this).GetIMax(ii) - indNG + 1;
	  bound2.iMax[ii] = (*this).GetIMax(ii) - indNG + 1;
	  if ( ii >= (*this).NumSurfI() && ii < (*this).NumSurfI() + (*this).NumSurfJ() ){ //j-surface
	    del1J++;
	  }
	  else{ //k-surface
	    del1K++;
	  }
	}
	else if ( (*this).GetIMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.iMax[ii] = indNG;
	  bound2.iMin[ii] = 1;
	  bound2.iMax[ii] = (*this).GetIMax(ii) - indNG + 1;
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
      bound1.iMin.erase(bound1.iMin.begin() + del1[ii]);
      bound1.iMax.erase(bound1.iMax.begin() + del1[ii]);
      bound1.jMin.erase(bound1.jMin.begin() + del1[ii]);
      bound1.jMax.erase(bound1.jMax.begin() + del1[ii]);
      bound1.kMin.erase(bound1.kMin.begin() + del1[ii]);
      bound1.kMax.erase(bound1.kMax.begin() + del1[ii]);
      bound1.tag.erase(bound1.tag.begin() + del1[ii]);
      bound1.numSurfJ -= del1J;
      bound1.numSurfK -= del1K;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.iMin.erase(bound2.iMin.begin() + del2[ii]);
      bound2.iMax.erase(bound2.iMax.begin() + del2[ii]);
      bound2.jMin.erase(bound2.jMin.begin() + del2[ii]);
      bound2.jMax.erase(bound2.jMax.begin() + del2[ii]);
      bound2.kMin.erase(bound2.kMin.begin() + del2[ii]);
      bound2.kMax.erase(bound2.kMax.begin() + del2[ii]);
      bound2.tag.erase(bound2.tag.begin() + del2[ii]);
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
	  bound2.bcTypes[ii] = "interblock";
	  bound2.jMin[ii] = (*this).GetJMin(ii);
	  bound2.jMax[ii] = (*this).GetJMax(ii);
	  bound2.tag[ii] = tag;
	}
	else {
	  //at upper j surface, lower bc is now interface
	  int tag = 3000 + numBlk; //upper surface matches with lower surface
	  bound1.bcTypes[ii] = "interblock";
	  bound1.jMin[ii] = indNG;
	  bound1.jMax[ii] = indNG;
	  bound1.tag[ii] = tag;

	  //at upper j surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.jMin[ii] = (*this).GetJMax(ii) - indNG + 1;
	  bound2.jMax[ii] = (*this).GetJMax(ii) - indNG + 1;
	}
      }
      else { //i-surface or k-surface
	if ( (*this).GetJMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.jMin[ii] = (*this).GetJMax(ii) - indNG + 1;
	  bound2.jMax[ii] = (*this).GetJMax(ii) - indNG + 1;
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del1I++;
	  }
	  else{ //k-surface
	    del1K++;
	  }
	}
	else if ( (*this).GetJMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.jMax[ii] = indNG;
	  bound2.jMin[ii] = 1;
	  bound2.jMax[ii] = (*this).GetJMax(ii) - indNG + 1;
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
      bound1.iMin.erase(bound1.iMin.begin() + del1[ii]);
      bound1.iMax.erase(bound1.iMax.begin() + del1[ii]);
      bound1.jMin.erase(bound1.jMin.begin() + del1[ii]);
      bound1.jMax.erase(bound1.jMax.begin() + del1[ii]);
      bound1.kMin.erase(bound1.kMin.begin() + del1[ii]);
      bound1.kMax.erase(bound1.kMax.begin() + del1[ii]);
      bound1.tag.erase(bound1.tag.begin() + del1[ii]);
      bound1.numSurfI -= del1I;
      bound1.numSurfK -= del1K;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.iMin.erase(bound2.iMin.begin() + del2[ii]);
      bound2.iMax.erase(bound2.iMax.begin() + del2[ii]);
      bound2.jMin.erase(bound2.jMin.begin() + del2[ii]);
      bound2.jMax.erase(bound2.jMax.begin() + del2[ii]);
      bound2.kMin.erase(bound2.kMin.begin() + del2[ii]);
      bound2.kMax.erase(bound2.kMax.begin() + del2[ii]);
      bound2.tag.erase(bound2.tag.begin() + del2[ii]);
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
	  bound2.bcTypes[ii] = "interblock";
	  bound2.kMin[ii] = (*this).GetKMin(ii);
	  bound2.kMax[ii] = (*this).GetKMax(ii);
	  bound2.tag[ii] = tag;
	}
	else {
	  //at upper k surface, lower bc is now interface
	  int tag = 5000 + numBlk; //upper surface matches with lower surface
	  bound1.bcTypes[ii] = "interblock";
	  bound1.kMin[ii] = indNG;
	  bound1.kMax[ii] = indNG;
	  bound1.tag[ii] = tag;

	  //at upper k surface, upper bc is same as original, but indices are adjusted for new block size
	  bound2.kMin[ii] = (*this).GetKMax(ii) - indNG + 1;
	  bound2.kMax[ii] = (*this).GetKMax(ii) - indNG + 1;
	}
      }
      else { //i-surface or j-surface
	if ( (*this).GetKMin(ii) >= indNG ){ //this surface is only present in the upper split
	  del1.push_back(ii);
	  bound2.kMin[ii] = (*this).GetKMax(ii) - indNG + 1;
	  bound2.kMax[ii] = (*this).GetKMax(ii) - indNG + 1;
	  if ( ii < (*this).NumSurfI() ){ //i-surface
	    del1I++;
	  }
	  else{ //j-surface
	    del1J++;
	  }
	}
	else if ( (*this).GetKMax(ii) >= indNG ){ //this surface straddles the split
	  bound1.kMax[ii] = indNG;
	  bound2.kMin[ii] = 1;
	  bound2.kMax[ii] = (*this).GetKMax(ii) - indNG + 1;
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
      bound1.iMin.erase(bound1.iMin.begin() + del1[ii]);
      bound1.iMax.erase(bound1.iMax.begin() + del1[ii]);
      bound1.jMin.erase(bound1.jMin.begin() + del1[ii]);
      bound1.jMax.erase(bound1.jMax.begin() + del1[ii]);
      bound1.kMin.erase(bound1.kMin.begin() + del1[ii]);
      bound1.kMax.erase(bound1.kMax.begin() + del1[ii]);
      bound1.tag.erase(bound1.tag.begin() + del1[ii]);
      bound1.numSurfI -= del1I;
      bound1.numSurfJ -= del1J;
    }
    for ( unsigned int ii = 0; ii < del2.size(); ii++ ){
      bound2.iMin.erase(bound2.iMin.begin() + del2[ii]);
      bound2.iMax.erase(bound2.iMax.begin() + del2[ii]);
      bound2.jMin.erase(bound2.jMin.begin() + del2[ii]);
      bound2.jMax.erase(bound2.jMax.begin() + del2[ii]);
      bound2.kMin.erase(bound2.kMin.begin() + del2[ii]);
      bound2.kMax.erase(bound2.kMax.begin() + del2[ii]);
      bound2.tag.erase(bound2.tag.begin() + del2[ii]);
      bound2.numSurfI -= del2I;
      bound2.numSurfJ -= del2J;
    }

  }
  else{
    cerr << "ERROR: Error in procBlock::Split(). Direction " << dir << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }

  (*this) = bound1;
  return bound2;

}

/* Member function to join 2 boundaryConditions. It assumes that the calling instance is the "lower" boundary condition and the input instance
is the "upper" boundary condition.
*/
void boundaryConditions::Join( const boundaryConditions &bc, const string &dir ){

  if ( dir == "i" ){ //split along i-plane

    for ( int ii = 0; ii < (*this).NumSurfI(); ii++ ){
      if ( (*this).GetIMax(ii) != 1){ //upper i surface
	//at upper i surface, get BCs from upper BC
	(*this).bcTypes[ii] = bc.GetBCTypes(ii);
	(*this).iMin[ii] += bc.GetIMin(ii) - 1;
	(*this).iMax[ii] += bc.GetIMax(ii) - 1;
	(*this).tag[ii] = bc.GetTag(ii);
      }
    }

  }
  else if ( dir == "j" ){ //split along j-plane

    for ( int ii = (*this).NumSurfI(); ii < (*this).NumSurfI() + (*this).NumSurfJ(); ii++ ){
      if ( (*this).GetJMax(ii) != 1){ //upper j surface
	//at upper i surface, get BCs from upper BC
	(*this).bcTypes[ii] = bc.GetBCTypes(ii);
	(*this).jMin[ii] += bc.GetJMin(ii) - 1;
	(*this).jMax[ii] += bc.GetJMax(ii) - 1;
	(*this).tag[ii] = bc.GetTag(ii);
      }
    }

  }
  else if ( dir == "k" ){ //split along k-plane

    for ( int ii = (*this).NumSurfI() + (*this).NumSurfJ(); ii < (*this).NumSurfaces(); ii++ ){
      if ( (*this).GetKMax(ii) != 1){ //upper k surface
	//at upper i surface, get BCs from upper BC
	(*this).bcTypes[ii] = bc.GetBCTypes(ii);
	(*this).kMin[ii] += bc.GetKMin(ii) - 1;
	(*this).kMax[ii] += bc.GetKMax(ii) - 1;
	(*this).tag[ii] = bc.GetTag(ii);
      }
    }

  }
  else{
    cerr << "ERROR: Error in procBlock::Join(). Direction " << dir << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }

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
}

//constructor with arguements passed
patch::patch( const int &bound, const int &b, const int &d1s, const int &d1e, const int &d2s, const int &d2e, const int &d3s, 
	      const int &d3e, const plot3dBlock &blk){
  // bound -- boundary number which patch is on (1-6)
  // b -- parent block number
  // d1s -- direction 1 starting index
  // d1e -- direction 1 ending index
  // d2s -- direction 2 starting index
  // d2e -- direction 2 ending index
  // d3s -- direction 3 surface index (constant surface that patch is on)

  boundary = bound;
  block = b;

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

void boundaryConditions::PackBC( char *(&sendBuffer), const int &sendBufSize, int &position)const{


  //get string lengths for each boundary condition to be sent, so processors unpacking know how much data to unpack for each string
  vector<int> strLength((*this).numSurfI + (*this).numSurfJ + (*this).numSurfK );
  for ( unsigned int jj = 0; jj < strLength.size(); jj++ ){
    strLength[jj] = (*this).GetBCTypes(jj).size();
  }

  MPI_Pack(&(*this).numSurfI, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numSurfJ, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).numSurfK, 1, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).iMin[0], (*this).iMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).iMax[0], (*this).iMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).jMin[0], (*this).jMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).jMax[0], (*this).jMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).kMin[0], (*this).kMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).kMax[0], (*this).kMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&(*this).tag[0], (*this).tag.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  MPI_Pack(&strLength[0], strLength.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  for ( int jj = 0; jj < ((*this).numSurfI + (*this).numSurfJ + (*this).numSurfK); jj++ ){
    MPI_Pack((*this).bcTypes[jj].c_str(), (*this).bcTypes[jj].size(), MPI_CHAR, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
  }

}


void boundaryConditions::UnpackBC( char *(&recvBuffer), const int &recvBufSize, int &position){


  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfI, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of i-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfJ, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of j-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).numSurfK, 1, MPI_INT, MPI_COMM_WORLD); //unpack number of k-surfaces

  (*this).ResizeVecs((*this).numSurfI + (*this).numSurfJ + (*this).numSurfK);
  vector<int> strLength((*this).numSurfI + (*this).numSurfJ + (*this).numSurfK);

  //unpack boundary condition data into appropriate vectors
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).iMin[0], (*this).iMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack i min coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).iMax[0], (*this).iMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack i max coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).jMin[0], (*this).jMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack j min coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).jMax[0], (*this).jMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack j max coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).kMin[0], (*this).kMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack k min coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).kMax[0], (*this).kMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack k max coordinates
  MPI_Unpack(recvBuffer, recvBufSize, &position, &(*this).tag[0], (*this).tag.size(), MPI_INT, MPI_COMM_WORLD); //unpack tags
  MPI_Unpack(recvBuffer, recvBufSize, &position, &strLength[0], strLength.size(), MPI_INT, MPI_COMM_WORLD); //unpack string sizes
  //unpack boundary condition names
  for ( unsigned int jj = 0; jj < strLength.size(); jj++ ){
    char *nameBuf = new char[strLength[jj]]; //allocate buffer to store BC name
    MPI_Unpack(recvBuffer, recvBufSize, &position, &nameBuf[0], strLength[jj], MPI_CHAR, MPI_COMM_WORLD); //unpack bc types
    string bcName(nameBuf, strLength[jj]); //create string of bc name
    (*this).bcTypes[jj] = bcName;
    delete [] nameBuf; //deallocate bc name buffer
  }


}

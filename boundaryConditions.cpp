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
boundaryConditions::boundaryConditions( int i, int j, int k){
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

//constructor when passed no arguements
interblock::interblock(){
  //initialize all variables to zero
  pair<int,int> zero(0,0);
  block = zero;
  boundary = zero;
  d1Start = zero;
  d1End = zero;
  d2Start = zero;
  d2End = zero;
  orientation = 0;
}

//operator overload for << - allows use of cout, cerr, etc.
ostream & operator<< (ostream &os, const interblock &bc){

  os << "Blocks: " << bc.BlockFirst() << ", " << bc.BlockSecond() << endl;
  os << "Boundaries: " << bc.BoundaryFirst() << ", " << bc.BoundarySecond() << endl;
  os << "Direction 1 Starts: " << bc.Dir1StartFirst() << ", " << bc.Dir1StartSecond() << endl;
  os << "Direction 1 Ends: " << bc.Dir1EndFirst() << ", " << bc.Dir1EndSecond() << endl;
  os << "Direction 2 Starts: " << bc.Dir2StartFirst() << ", " << bc.Dir2StartSecond() << endl;
  os << "Direction 2 Ends: " << bc.Dir2EndFirst() << ", " << bc.Dir2EndSecond() << endl;
  os << "Orientation: " << bc.Orientation() << endl;

  return os;
}

void interblock::SetInterblock(const patch &p1, const patch &p2){

  (*this).SetBlockFirst(p1.Block());
  (*this).SetBlockSecond(p2.Block());

  (*this).SetBoundaryFirst(p1.Boundary());
  (*this).SetBoundarySecond(p2.Boundary());

  (*this).SetDir1StartFirst(p1.Dir1Start());
  (*this).SetDir1StartSecond(p2.Dir1Start());

  (*this).SetDir1EndFirst(p1.Dir1End());
  (*this).SetDir1EndSecond(p2.Dir1End());

  (*this).SetDir2StartFirst(p1.Dir2Start());
  (*this).SetDir2StartSecond(p2.Dir2Start());

  orientation = 0;
}


/* Function to go through the boundary conditions and pair the interblock BCs together and determine their orientation.

*/
vector<interblock> GetInterblockBCs( const vector<boundaryConditions> &bc, const vector<plot3dBlock> &grid ){

  cout << "In GetInterblockBCs()" << endl;

  vector<vector<int> > isolatedInterblocks;

  for ( unsigned int ii = 0; ii < bc.size(); ii++ ){ //loop over all blocks
    int numSurf = bc[ii].NumSurfI() + bc[ii].NumSurfJ() + bc[ii].NumSurfK();

    for ( int jj = 0; jj < numSurf; jj++ ){ //loop over number of surfaces in block

      if ( bc[ii].GetBCTypes(jj) == "interblock" ){
	vector<int> temp (10,0);
	temp[0] = ii;                                  //block number of bc
	
	//boundary number of bc
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
	temp[2] = bc[ii].GetIMin(jj) - 1;              //i min of bc patch
	temp[3] = bc[ii].GetIMax(jj) - 1;              //i max of bc patch
	temp[4] = bc[ii].GetJMin(jj) - 1;              //j min of bc patch
	temp[5] = bc[ii].GetJMax(jj) - 1;              //j max of bc patch
	temp[6] = bc[ii].GetKMin(jj) - 1;              //k min of bc patch
	temp[7] = bc[ii].GetKMax(jj) - 1;              //k max of bc patch

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

	temp[8] = bound;                                  //boundary of pair
	temp[9] = blk;                                    //block of pair

	isolatedInterblocks.push_back(temp);

      }

    }
  }

  for (unsigned int ii = 0; ii < isolatedInterblocks.size(); ii++ ){
    for ( unsigned int jj = 0; jj < isolatedInterblocks[0].size(); jj++ ){
      cout << isolatedInterblocks[ii][jj] << ", ";
    }
    cout << endl;
  }

  //intialize vector of interblocks to return
  vector<interblock> connections(isolatedInterblocks.size()/2);

  //loop over isolated interblocks
  for ( unsigned int ii = 0; ii < isolatedInterblocks.size(); ii+=2 ){
    for ( unsigned int jj = ii+1; jj < isolatedInterblocks.size(); jj++ ){

      if ( isolatedInterblocks[ii][9] == isolatedInterblocks[jj][0] ) { //blocks between interblock BCs match
	if ( isolatedInterblocks[ii][8] == isolatedInterblocks[jj][1] ) { //boundary surfaces between interblock BCs match
	  //get current patch
	  patch cPatch(isolatedInterblocks[ii][1], isolatedInterblocks[ii][0], isolatedInterblocks[ii][2], isolatedInterblocks[ii][3], isolatedInterblocks[ii][4],
		       isolatedInterblocks[ii][5], isolatedInterblocks[ii][6], isolatedInterblocks[ii][7], grid[isolatedInterblocks[ii][0]]);

	  //get new patch
	  patch nPatch(isolatedInterblocks[jj][1], isolatedInterblocks[jj][0], isolatedInterblocks[jj][2], isolatedInterblocks[jj][3], isolatedInterblocks[jj][4],
		       isolatedInterblocks[jj][5], isolatedInterblocks[jj][6], isolatedInterblocks[jj][7], grid[isolatedInterblocks[jj][0]]);

	  //test for match
	  interblock match;
	  if ( TestPatchMatch(cPatch, nPatch, match) ){ //match found
	    connections[ii/2] = match;
	    swap(isolatedInterblocks[jj], isolatedInterblocks[ii+1]); //swap matched interblock BC to top portion of vector so it is not searched again
	  }

	}
      }
    }
  }

  for (unsigned int ii = 0; ii < connections.size(); ii++ ){
    cout << "Interblock " << ii << endl;
    cout << connections[ii] << endl;
  }


  exit(0);

  return connections;
}

/* Function to take in the origin of the patch to be matched, as well as the origin locations of the possible match and determine if there
is a match.
*/
bool TestPatchMatch( const patch &p1, const patch &p2, interblock &inter ){

  bool match = false; //initialize match to false

  //determine if there is a potential match by comparing origins
  if ( p1.Origin() == p2.Origin() ){ //origins match -----------------------------------------------------------------------------------

    //if origin matches origin, corner 1 can only be at corner 1 or 2
    if ( p1.Corner1() == p2.Corner1() ){ //corner 1s match

      //if all 3 corners match, same orientation
      if ( p1.Corner2() == p2.Corner2() ){ //corner 2s match
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(1);
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner2() ){ //corner 1 matches corner 2

      //if origins match and 1 matches 2, 2 must match 1
      if ( p1.Corner2() == p2.Corner1() ){ //corner 2 matches corner 1
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(2);
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
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(3);
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner12() ){ 

      //corner 2 must match origin for match
      if ( p1.Corner2() == p2.Origin() ){ 
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(4);
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
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(5);
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner12() ){ 

      //corner 2 must match origin for match
      if ( p1.Corner2() == p2.Origin() ){ 
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(6);
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
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(7);
	match = true;
      }
      else{ //no match
	return match;
      }
    }
    else if ( p1.Corner1() == p2.Corner2() ){ 

      //corner 2 must match corner 1 for match
      if ( p1.Corner2() == p2.Corner2() ){ 
	inter.SetInterblock(p1,p2);
	inter.SetOrientation(8);
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

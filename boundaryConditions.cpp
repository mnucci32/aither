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

      if ( isolatedInterblocks[ii][9] == isolatedInterblocks[jj][0] ) { //block i'm looking for matches block i'm at

	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	if ( (isolatedInterblocks[ii][8] == 1 || isolatedInterblocks[ii][8] == 2 ) && 
	     isolatedInterblocks[jj][2] == isolatedInterblocks[jj][3] ){ //boundary i'm looking for matches boundary i'm at i-i boundary

	    int iIndl = 0;
	    if ( isolatedInterblocks[ii][8] == 1 ){ //boundary i'm looking for is a i-lower boundary
	      iIndl = 2;
	    }
	    else{
	      iIndl = 3;
	    }

	    //get location of origin for patch to be matched
	    //at lower/upper i surface -- [ii][iIndl]
	    //origin is at jmin, kmin of patch -- [ii][4], [ii][6]
	    int loc = GetLoc1D(isolatedInterblocks[ii][iIndl], isolatedInterblocks[ii][4], isolatedInterblocks[ii][6], 
			       grid[isolatedInterblocks[ii][0]].NumI(), grid[isolatedInterblocks[ii][0]].NumJ()); 
	    vector3d<double> oLook(grid[isolatedInterblocks[ii][0]].XLoc(loc),
				   grid[isolatedInterblocks[ii][0]].YLoc(loc),
				   grid[isolatedInterblocks[ii][0]].ZLoc(loc));

	  if ( isolatedInterblocks[jj][1] == 1 || isolatedInterblocks[jj][1] == 2){ //test for i-surface to i-surface match

	    int iInd = 0;
	    if ( isolatedInterblocks[jj][1] == 1 ){ //boundary i'm testing for a match is a i-lower boundary
	      iInd = 2;
	    }
	    else{
	      iInd = 3;
	    }

	    //get location of origin for patch that is potential match
	    //at lower/upper i surface -- [jj][iInd]
	    //origin is at jmin, kmin of patch -- [jj][4], [jj][6]
	    int loc1 = GetLoc1D(isolatedInterblocks[jj][iInd], isolatedInterblocks[jj][4], isolatedInterblocks[jj][6], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr1(grid[isolatedInterblocks[jj][0]].XLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc1));

	    //get location of origin for patch that is potential match
	    //at lower/upper i surface -- [jj][iInd]
	    //origin is at jmin, kmax of patch -- [jj][4], [jj][7]
	    int loc2 = GetLoc1D(isolatedInterblocks[jj][iInd], isolatedInterblocks[jj][4], isolatedInterblocks[jj][7], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr2(grid[isolatedInterblocks[jj][0]].XLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc2));

	    //get location of origin for patch that is potential match
	    //at lower/upper i surface -- [jj][iInd]
	    //origin is at jmax, kmin of patch -- [jj][5], [jj][6]
	    int loc3 = GetLoc1D(isolatedInterblocks[jj][iInd], isolatedInterblocks[jj][5], isolatedInterblocks[jj][6], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr3(grid[isolatedInterblocks[jj][0]].XLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc3));

	    //get location of origin for patch that is potential match
	    //at lower/upper i surface -- [jj][iInd]
	    //origin is at jmax, kmax of patch -- [jj][5], [jj][7]
	    int loc4 = GetLoc1D(isolatedInterblocks[jj][iInd], isolatedInterblocks[jj][5], isolatedInterblocks[jj][7], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr4(grid[isolatedInterblocks[jj][0]].XLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc4));

	    if ( oLook == oCurr1 || oLook == oCurr2 || oLook == oCurr3 || oLook == oCurr4 ) { //match found
	      interblock temp;
	      temp.SetBlockFirst(isolatedInterblocks[ii][0]);
	      temp.SetBlockSecond(isolatedInterblocks[jj][0]);
	      temp.SetBoundaryFirst(isolatedInterblocks[ii][1]);
	      temp.SetBoundarySecond(isolatedInterblocks[jj][1]);
	      temp.SetDir1StartFirst(isolatedInterblocks[ii][4]); //for i-surfaces dir1 is j
	      temp.SetDir1StartSecond(isolatedInterblocks[jj][4]);
	      temp.SetDir1EndFirst(isolatedInterblocks[ii][5]); 
	      temp.SetDir1EndSecond(isolatedInterblocks[jj][5]);
	      temp.SetDir2StartFirst(isolatedInterblocks[ii][6]); //for i-surfaces dir2 is k
	      temp.SetDir2StartSecond(isolatedInterblocks[jj][6]);
	      temp.SetDir2EndFirst(isolatedInterblocks[ii][7]); 
	      temp.SetDir2EndSecond(isolatedInterblocks[jj][7]);
	      if ( oLook == oCurr1){
		temp.SetOrientation(1);
	      }
	      else if ( oLook == oCurr2){
		temp.SetOrientation(2);
	      }
	      else if ( oLook == oCurr3){
		temp.SetOrientation(3);
	      }
	      else{
		temp.SetOrientation(4);
	      }

	      connections[ii/2] = temp;
	      swap(isolatedInterblocks[ii+1], isolatedInterblocks[jj]); //swap matched bc so it is not searched for again
	    }

	  }
	  else{ //test for i-surface to j/k surface match

	  }
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( (isolatedInterblocks[ii][8] == 3 || isolatedInterblocks[ii][8] == 4) && 
		  isolatedInterblocks[jj][4] == isolatedInterblocks[jj][5] ){ //boundary i'm looking for matches boundary i'm at j-j boundary

	  int jIndl = 0;
	  if ( isolatedInterblocks[ii][8] == 3 ){ //boundary i'm looking for is a j-lower boundary
	    jIndl = 4;
	  }
	  else{
	    jIndl = 5;
	  }

	  //get location of origin for patch to be matched
	  //at lower/upper j surface -- [ii][jIndl]
	  //origin is at kmin, imin of patch -- [ii][6], [ii][2]
	  int loc = GetLoc1D(isolatedInterblocks[ii][2], isolatedInterblocks[ii][jIndl], isolatedInterblocks[ii][6], 
			     grid[isolatedInterblocks[ii][0]].NumI(), grid[isolatedInterblocks[ii][0]].NumJ()); 
	  vector3d<double> oLook(grid[isolatedInterblocks[ii][0]].XLoc(loc),
				 grid[isolatedInterblocks[ii][0]].YLoc(loc),
				 grid[isolatedInterblocks[ii][0]].ZLoc(loc));

	  if ( isolatedInterblocks[jj][1] == 3 || isolatedInterblocks[jj][1] == 4){ //test for j-surface to j-surface match

	    int jInd = 0;
	    if ( isolatedInterblocks[jj][1] == 3 ){ //boundary i'm testing for a match is a j-lower boundary
	      jInd = 4;
	    }
	    else{
	      jInd = 5;
	    }

	    //get location of origin for patch that is potential match
	    //at lower/upper j surface -- [jj][jInd]
	    //origin is at kmin, imin of patch -- [jj][6], [jj][2]
	    int loc1 = GetLoc1D(isolatedInterblocks[jj][2], isolatedInterblocks[jj][jInd], isolatedInterblocks[jj][6],  
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr1(grid[isolatedInterblocks[jj][0]].XLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc1));

	    //get location of origin for patch that is potential match
	    //at lower/upper j surface -- [jj][jInd]
	    //origin is at kmin, imax of patch -- [jj][6], [jj][3]
	    int loc2 = GetLoc1D( isolatedInterblocks[jj][3], isolatedInterblocks[jj][jInd], isolatedInterblocks[jj][6], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr2(grid[isolatedInterblocks[jj][0]].XLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc2));

	    //get location of origin for patch that is potential match
	    //at lower/upper j surface -- [jj][jInd]
	    //origin is at kmax, imin of patch -- [jj][7], [jj][2]
	    int loc3 = GetLoc1D( isolatedInterblocks[jj][2], isolatedInterblocks[jj][jInd], isolatedInterblocks[jj][7], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr3(grid[isolatedInterblocks[jj][0]].XLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc3));

	    //get location of origin for patch that is potential match
	    //at lower/upper j surface -- [jj][jInd]
	    //origin is at kmax, imax of patch -- [jj][7], [jj][3]
	    int loc4 = GetLoc1D( isolatedInterblocks[jj][3], isolatedInterblocks[jj][jInd], isolatedInterblocks[jj][7],
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr4(grid[isolatedInterblocks[jj][0]].XLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc4));

	    if ( oLook == oCurr1 || oLook == oCurr2 || oLook == oCurr3 || oLook == oCurr4 ) { //match found
	      interblock temp;
	      temp.SetBlockFirst(isolatedInterblocks[ii][0]);
	      temp.SetBlockSecond(isolatedInterblocks[jj][0]);
	      temp.SetBoundaryFirst(isolatedInterblocks[ii][1]);
	      temp.SetBoundarySecond(isolatedInterblocks[jj][1]);
	      temp.SetDir1StartFirst(isolatedInterblocks[ii][6]); //for j-surfaces dir1 is k
	      temp.SetDir1StartSecond(isolatedInterblocks[jj][6]);
	      temp.SetDir1EndFirst(isolatedInterblocks[ii][7]); 
	      temp.SetDir1EndSecond(isolatedInterblocks[jj][7]);
	      temp.SetDir2StartFirst(isolatedInterblocks[ii][2]); //for j-surfaces dir2 is i
	      temp.SetDir2StartSecond(isolatedInterblocks[jj][2]);
	      temp.SetDir2EndFirst(isolatedInterblocks[ii][3]); 
	      temp.SetDir2EndSecond(isolatedInterblocks[jj][3]);
	      if ( oLook == oCurr1){
		temp.SetOrientation(1);
	      }
	      else if ( oLook == oCurr2){
		temp.SetOrientation(2);
	      }
	      else if ( oLook == oCurr3){
		temp.SetOrientation(3);
	      }
	      else{
		temp.SetOrientation(4);
	      }

	      connections[ii/2] = temp;
	      swap(isolatedInterblocks[ii+1], isolatedInterblocks[jj]); //swap matched bc so it is not searched for again
	    }

	  }
	  else{ //test for j-surface to i/k surface match

	  }
	}
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------
	else if ( (isolatedInterblocks[ii][8] == 5 || isolatedInterblocks[ii][8] == 6) && 
		  isolatedInterblocks[jj][6] == isolatedInterblocks[jj][7] ){ //boundary i'm looking for matches boundary i'm at k-k boundary

	  int kIndl = 0;
	  if ( isolatedInterblocks[ii][8] == 5 ){ //boundary i'm looking for is a k-lower boundary
	    kIndl = 6;
	  }
	  else{
	    kIndl = 7;
	  }

	  //get location of origin for patch to be matched
	  //at lower/upper k surface -- [ii][kIndl]
	  //origin is at imin, jmin of patch -- [ii][2], [ii][4]
	  int loc = GetLoc1D(isolatedInterblocks[ii][2], isolatedInterblocks[ii][4], isolatedInterblocks[ii][kIndl], 
			     grid[isolatedInterblocks[ii][0]].NumI(), grid[isolatedInterblocks[ii][0]].NumJ()); 
	  vector3d<double> oLook(grid[isolatedInterblocks[ii][0]].XLoc(loc),
				 grid[isolatedInterblocks[ii][0]].YLoc(loc),
				 grid[isolatedInterblocks[ii][0]].ZLoc(loc));

	  if ( isolatedInterblocks[jj][1] == 5 || isolatedInterblocks[jj][1] == 6){ //test for k-surface to k-surface match

	    int kInd = 0;
	    if ( isolatedInterblocks[jj][1] == 5 ){ //boundary i'm testing for a match is a k-lower boundary
	      kInd = 6;
	    }
	    else{
	      kInd = 7;
	    }

	    //get location of origin for patch that is potential match
	    //at lower/upper k surface -- [jj][kInd]
	    //origin is at imin, jmin of patch -- [jj][2], [jj][4]
	    int loc1 = GetLoc1D(isolatedInterblocks[jj][2], isolatedInterblocks[jj][4], isolatedInterblocks[jj][kInd], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr1(grid[isolatedInterblocks[jj][0]].XLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc1),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc1));

	    //get location of origin for patch that is potential match
	    //at lower/upper k surface -- [jj][kInd]
	    //origin is at imax, jmin of patch -- [jj][3], [jj][4]
	    int loc2 = GetLoc1D( isolatedInterblocks[jj][3], isolatedInterblocks[jj][4], isolatedInterblocks[jj][kInd], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr2(grid[isolatedInterblocks[jj][0]].XLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc2),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc2));

	    //get location of origin for patch that is potential match
	    //at lower/upper k surface -- [jj][kInd]
	    //origin is at imin, jmax of patch -- [jj][2], [jj][5]
	    int loc3 = GetLoc1D( isolatedInterblocks[jj][2], isolatedInterblocks[jj][5], isolatedInterblocks[jj][kInd], 
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr3(grid[isolatedInterblocks[jj][0]].XLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc3),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc3));

	    //get location of origin for patch that is potential match
	    //at lower/upper k surface -- [jj][kInd]
	    //origin is at imax, jmax of patch -- [jj][3], [jj][5]
	    int loc4 = GetLoc1D( isolatedInterblocks[jj][3], isolatedInterblocks[jj][5], isolatedInterblocks[jj][kInd],
				grid[isolatedInterblocks[jj][0]].NumI(), grid[isolatedInterblocks[jj][0]].NumJ()); 
	    vector3d<double> oCurr4(grid[isolatedInterblocks[jj][0]].XLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].YLoc(loc4),
				    grid[isolatedInterblocks[jj][0]].ZLoc(loc4));

	    if ( oLook == oCurr1 || oLook == oCurr2 || oLook == oCurr3 || oLook == oCurr4 ) { //match found
	      interblock temp;
	      temp.SetBlockFirst(isolatedInterblocks[ii][0]);
	      temp.SetBlockSecond(isolatedInterblocks[jj][0]);
	      temp.SetBoundaryFirst(isolatedInterblocks[ii][1]);
	      temp.SetBoundarySecond(isolatedInterblocks[jj][1]);
	      temp.SetDir1StartFirst(isolatedInterblocks[ii][2]); //for k-surfaces dir1 is i
	      temp.SetDir1StartSecond(isolatedInterblocks[jj][2]);
	      temp.SetDir1EndFirst(isolatedInterblocks[ii][3]); 
	      temp.SetDir1EndSecond(isolatedInterblocks[jj][3]);
	      temp.SetDir2StartFirst(isolatedInterblocks[ii][4]); //for k-surfaces dir2 is j
	      temp.SetDir2StartSecond(isolatedInterblocks[jj][4]);
	      temp.SetDir2EndFirst(isolatedInterblocks[ii][5]); 
	      temp.SetDir2EndSecond(isolatedInterblocks[jj][5]);
	      if ( oLook == oCurr1){
		temp.SetOrientation(1);
	      }
	      else if ( oLook == oCurr2){
		temp.SetOrientation(2);
	      }
	      else if ( oLook == oCurr3){
		temp.SetOrientation(3);
	      }
	      else{
		temp.SetOrientation(4);
	      }

	      connections[ii/2] = temp;
	      swap(isolatedInterblocks[ii+1], isolatedInterblocks[jj]); //swap matched bc so it is not searched for again
	    }

	  }
	  else{ //test for k-surface to i/j surface match

	  }
	}

      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------------

    }
  }

  for (unsigned int ii = 0; ii < connections.size(); ii++ ){
    cout << "Interblock " << ii << endl;
    cout << connections[ii] << endl;
  }


  exit(0);

  return connections;
}

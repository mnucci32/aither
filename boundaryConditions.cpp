#include "boundaryConditions.h"
#include <math.h>       //sqrt
#include <iostream>     //cout

using std::cout;
using std::endl;
using std::cerr;

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

  os << "Number of I surfaces: " << bc.numSurfI << endl;
  os << "Number of J surfaces: " << bc.numSurfJ << endl;
  os << "Number of K surfaces: " << bc.numSurfK << endl;

  unsigned int ii = 0; 
  for ( ii = 0; ii < bc.bcTypes.size(); ii++ ){
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
  int nn = 0;

  if (surf == "il"){ //lower i-surface
    iStart = 0;
    iEnd = (*this).NumSurfI();
  }
  else if (surf == "iu"){ //upper i-surface
    iStart = 0;
    iEnd = (*this).NumSurfI();
  }
  else if (surf == "jl"){ //lower j-surface
    iStart = (*this).NumSurfI();
    iEnd = iStart + (*this).NumSurfJ();
  }
  else if (surf == "ju"){ //upper j-surface
    iStart = (*this).NumSurfI();
    iEnd = iStart + (*this).NumSurfJ();
    //j = j+1;
  }
  else if (surf == "kl"){ //lower k-surface
    iStart = (*this).NumSurfI() + (*this).NumSurfJ();
    iEnd = iStart + (*this).NumSurfK();
  }
  else if (surf == "ku"){ //upper k-surface
    iStart = (*this).NumSurfI() + (*this).NumSurfJ();
    iEnd = iStart + (*this).NumSurfK();
    //k = k+1;
  }
  else {
    cerr << "ERROR: Surface type " << surf << " is not recognized!" << endl;
  }

  //Determine which boundary condition should be applied
  for ( nn = iStart; nn < iEnd; nn++ ){

    //Bounary mins and maxes start at 1 instead of 0, so 1 is subtracted
    if ( (i >= (*this).GetIMin(nn)-1 && i <= (*this).GetIMax(nn)-1 && j >= (*this).GetJMin(nn)-1 && j <= (*this).GetJMax(nn)-1
	  && k >= (*this).GetKMin(nn)-1 && k <= (*this).GetKMax(nn)-1) ){
      bcName = (*this).GetBCTypes(nn);
      break;
    }
  }

  return bcName;

}

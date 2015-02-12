#ifndef GEOMSLICEHEADERDEF             //only if the macro GEOMSLICEHEADERDEF is not defined execute these lines of code
#define GEOMSLICEHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "primVars.h" //primVars
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

class geomSlice {

  vector<vector3d<double> > center;        //coordinates of cell center
  vector<vector3d<double> > fAreaI;        //face area vector for i-faces
  vector<vector3d<double> > fAreaJ;        //face area vector for j-faces
  vector<vector3d<double> > fAreaK;        //face area vector for k-faces
  vector<vector3d<double> > fCenterI;      //coordinates of i-face centers
  vector<vector3d<double> > fCenterJ;      //coordinates of j-face centers
  vector<vector3d<double> > fCenterK;      //coordinates of k-face centers

  vector<double> vol ;                     //cell volume

  int numCells;                            //number of cells in block
  int numI;                                //i-dimension of block (cells)
  int numJ;                                //j-dimension of block (cells)
  int numK;                                //k-dimension of block (cells)
  int parBlock;                            //parent block number
  int parBlockStartI;                      //parent block starting index for i
  int parBlockEndI;                        //parent block ending index for i
  int parBlockStartJ;                      //parent block starting index for j
  int parBlockEndJ;                        //parent block ending index for j
  int parBlockStartK;                      //parent block starting index for k
  int parBlockEndK;                        //parent block ending index for k

 public:
  //constructors
  geomSlice();
  geomSlice( const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&,
	     const int&, const int& );

  //member functions
  void SetNumCells( const int &a){numCells = a;}
  int NumCells() const {return numCells;}
  void SetI( const int &a){numI = a;}
  int NumI() const {return numI;}
  void SetJ( const int &a){numJ = a;}
  int NumJ() const {return numJ;}
  void SetK( const int &a){numK = a;}
  int NumK() const {return numK;}

  void SetParentBlock( const int &a){parBlock = a;}
  int ParentBlock() const {return parBlock;}
  void SetParentBlockStartI( const int &a){parBlockStartI = a;}
  int ParentBlockStartI() const {return parBlockStartI;}
  void SetParentBlockEndI( const int &a){parBlockEndI = a;}
  int ParentBlockEndI() const {return parBlockEndI;}
  void SetParentBlockStartJ( const int &a){parBlockStartJ = a;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  void SetParentBlockEndJ( const int &a){parBlockEndJ = a;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  void SetParentBlockStartK( const int &a){parBlockStartK = a;}
  int ParentBlockStartK() const {return parBlockStartK;}
  void SetParentBlockEndK( const int &a){parBlockEndK = a;}
  int ParentBlockEndK() const {return parBlockEndK;}

  void SetVol( const double &a, const int &ind){vol[ind] = a;}
  double Vol(const int &ind) const {return vol[ind];}
  void SetCenter( const vector3d<double> &a, const int &ind){center[ind] = a;}
  vector3d<double> Center(const int &ind) const {return center[ind];}

  void SetFAreaI( const vector3d<double> &a, const int &ind){fAreaI[ind] = a;}
  vector3d<double> FAreaI(const int &ind) const {return fAreaI[ind];}
  void SetFAreaJ( const vector3d<double> &a, const int &ind){fAreaJ[ind] = a;}
  vector3d<double> FAreaJ(const int &ind) const {return fAreaJ[ind];}
  void SetFAreaK( const vector3d<double> &a, const int &ind){fAreaK[ind] = a;}
  vector3d<double> FAreaK(const int &ind) const {return fAreaK[ind];}
  void SetFCenterI( const vector3d<double> &a, const int &ind){fCenterI[ind] = a;}
  vector3d<double> FCenterI(const int &ind) const {return fCenterI[ind];}
  void SetFCenterJ( const vector3d<double> &a, const int &ind){fCenterJ[ind] = a;}
  vector3d<double> FCenterJ(const int &ind) const {return fCenterJ[ind];}
  void SetFCenterK( const vector3d<double> &a, const int &ind){fCenterK[ind] = a;}
  vector3d<double> FCenterK(const int &ind) const {return fCenterK[ind];}

  //destructor
  ~geomSlice() {}

};

class stateSlice {

 public:
  vector<primVars> state ;                     //cell states

  int numCells;                            //number of cells in block
  int numI;                                //i-dimension of block (cells)
  int numJ;                                //j-dimension of block (cells)
  int numK;                                //k-dimension of block (cells)
  int parBlock;                            //parent block number
  int parBlockStartI;                      //parent block starting index for i
  int parBlockEndI;                        //parent block ending index for i
  int parBlockStartJ;                      //parent block starting index for j
  int parBlockEndJ;                        //parent block ending index for j
  int parBlockStartK;                      //parent block starting index for k
  int parBlockEndK;                        //parent block ending index for k

  //constructors
  stateSlice();
  stateSlice( const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&,
	     const int&, const int& );

  //member functions
  void SetNumCells( const int &a){numCells = a;}
  int NumCells() const {return numCells;}
  void SetI( const int &a){numI = a;}
  int NumI() const {return numI;}
  void SetJ( const int &a){numJ = a;}
  int NumJ() const {return numJ;}
  void SetK( const int &a){numK = a;}
  int NumK() const {return numK;}

  void SetParentBlock( const int &a){parBlock = a;}
  int ParentBlock() const {return parBlock;}
  void SetParentBlockStartI( const int &a){parBlockStartI = a;}
  int ParentBlockStartI() const {return parBlockStartI;}
  void SetParentBlockEndI( const int &a){parBlockEndI = a;}
  int ParentBlockEndI() const {return parBlockEndI;}
  void SetParentBlockStartJ( const int &a){parBlockStartJ = a;}
  int ParentBlockStartJ() const {return parBlockStartJ;}
  void SetParentBlockEndJ( const int &a){parBlockEndJ = a;}
  int ParentBlockEndJ() const {return parBlockEndJ;}
  void SetParentBlockStartK( const int &a){parBlockStartK = a;}
  int ParentBlockStartK() const {return parBlockStartK;}
  void SetParentBlockEndK( const int &a){parBlockEndK = a;}
  int ParentBlockEndK() const {return parBlockEndK;}

  void SetState( const primVars &a, const int &ind){state[ind] = a;}
  primVars State(const int &ind) const {return state[ind];}

  //destructor
  ~stateSlice() {}

};


//function definitions

#endif

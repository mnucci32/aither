#ifndef PLOT3DHEADERDEF             //only if the macro PLOT3DHEADERDEF is not defined execute these lines of code

#define PLOT3DHEADERDEF             //define the macro

//This header file contains the functions and data structures associated with dealing with PLOT3D grid files

#include <vector>
#include <string>
#include "vector3d.h"

using std::vector;
using std::string;

//---------------------------------------------------------------------------------------------------------------//
//Class for an individual plot3d block
class plot3dBlock {
  //by default everything above the public: declaration is private
  vector<double> x;      //vector of x-coordinates
  vector<double> y;      //vector of y-coordinates
  vector<double> z;      //vector of z-coordinates
  int numi;              //number of points in i-direction
  int numj;              //number of points in j-direction
  int numk;              //number of points in k-direction

 public:
   //constructor -- create a plot3d block by passing the above quantities
   plot3dBlock( const int&, const int&, const int&, const vector<double>&, const vector<double>&, const vector<double>&);
   plot3dBlock( const int&, const int&, const int& );
   plot3dBlock();

   //member functions
   const vector<double> Volume() const;
   const vector<vector3d<double> > FaceAreaI() const;
   const vector<vector3d<double> > FaceAreaJ() const;
   const vector<vector3d<double> > FaceAreaK() const;
   const vector<vector3d<double> > Centroid() const;
   const vector<vector3d<double> > FaceCenterI() const;
   const vector<vector3d<double> > FaceCenterJ() const;
   const vector<vector3d<double> > FaceCenterK() const;
   int NumI() const {return numi;}
   int NumJ() const {return numj;}
   int NumK() const {return numk;}
   vector<double> const X(){return x;}
   vector<double> const Y(){return y;}
   vector<double> const Z(){return z;}
   double XLoc(const int &a)const{return x[a];}
   double YLoc(const int &a)const{return y[a];}
   double ZLoc(const int &a)const{return z[a];}

   plot3dBlock Split(const string&, const int&);
   void Join(const plot3dBlock&, const string&);

   //destructor
   ~plot3dBlock() {}

};

//------------------------------------------------------------------------------------------------------------------
//function declarations
vector<plot3dBlock> ReadP3dGrid(const string &gridName, double &numCells);
vector3d<int> GetIJK(const int &, const int &, const int &, const int &);

//input cell coordinates, get face coordinates
int GetUpperFaceI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetLowerFaceI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetUpperFaceJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetLowerFaceJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetUpperFaceK(const int &, const int &, const int &, const int &, const int &, int=1);
int GetLowerFaceK(const int &, const int &, const int &, const int &, const int &, int=1);

//input cell coordinates get neighbor cell coordinates
int GetNeighborUpI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetNeighborLowI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetNeighborUpJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetNeighborLowJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetNeighborUpK(const int &, const int &, const int &, const int &, const int &, int=1);
int GetNeighborLowK(const int &, const int &, const int &, const int &, const int &, int=1);

//input face coordinates, get cell coordinates
int GetCellFromFaceUpperI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetCellFromFaceLowerI(const int &, const int &, const int &, const int &, const int &, int=1);
int GetCellFromFaceUpperJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetCellFromFaceLowerJ(const int &, const int &, const int &, const int &, const int &, int=1);
int GetCellFromFaceUpperK(const int &, const int &, const int &, const int &, const int &, int=1);
int GetCellFromFaceLowerK(const int &, const int &, const int &, const int &, const int &, int=1);

//get location inside of 1D array
int GetLoc1D(const int &, const int &, const int &, const int &, const int &);

//find out if matrix should have data at the indicated cell
bool IsMatrixData(const int&, const int&, const int&, const int&, const int&, const int&, const string&);

//function to reorder block by hyperplanes
vector<vector3d<int> > HyperplaneReorder(const int &, const int &, const int &);

#endif

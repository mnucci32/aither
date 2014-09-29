#ifndef PLOT3DHEADERDEF             //only if the macro PLOT3DHEADERDEF is not defined execute these lines of code

#define PLOT3DHEADERDEF             //define the macro

//This header file contains the functions and data structures associated with dealing with PLOT3D grid and solution files

#include <vector>
#include <string>
#include "vector3d.h"


using std::vector;
using std::string;

//---------------------------------------------------------------------------------------------------------------//
//Class for an individual plot3d block
class plot3dBlock {
  //by default everything above the public: declaration is private
  int numi;              //number of points in i-direction
  int numj;              //number of points in j-direction
  int numk;              //number of points in k-direction
  vector<double> x;      //vector of x-coordinates
  vector<double> y;      //vector of y-coordinates
  vector<double> z;      //vector of z-coordinates

 public:
   //constructor -- create a plot3d block by passing the above quantities
   plot3dBlock( int, int, int, vector<double>&, vector<double>&, vector<double>&);
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
   void SetI(const int & dim){numi = dim;}
   void SetJ(const int & dim){numj = dim;}
   void SetK(const int & dim){numk = dim;}
   void SetX(const vector<double> & data){x = data;}
   void SetY(const vector<double> & data){y = data;}
   void SetZ(const vector<double> & data){z = data;}

   //destructor
   ~plot3dBlock() {}

};

//---------------------------------------------------------------------------------------------------------------//
//Class for an individual plot3d solution block
class plot3dQBlock {
  //by default everything above the public: declaration is private
  int numi;               //number of points in i-direction
  int numj;               //number of points in j-direction
  int numk;               //number of points in k-direction
  double mach;            //mach number
  double alpha;           //angle of attack
  double Re;              //Reynolds number
  double ti;              //time
  vector<double> q1;      //vector of first conserved variable
  vector<double> q2;      //vector of second conserved variable
  vector<double> q3;      //vector of third conserved variable
  vector<double> q4;      //vector of fourth conserved variable
  vector<double> q5;      //vector of fifth conserved variable

 public:
   //constructor -- create a plot3d block by passing the above quantities
   plot3dQBlock(const plot3dBlock&);
   plot3dQBlock();

   //member functions
   int NumI() const {return numi;}
   int NumJ() const {return numj;}
   int NumK() const {return numk;}
   double Mach() const {return mach;}
   double Alpha() const {return alpha;}
   double Reynolds() const {return Re;}
   double Ti() const {return ti;}
   vector<double> Q1() const {return q1;}
   vector<double> Q2() const {return q2;}
   vector<double> Q3() const {return q3;}
   vector<double> Q4() const {return q4;}
   vector<double> Q5() const {return q5;}
   void SetI(const int & dim){numi = dim;}
   void SetJ(const int & dim){numj = dim;}
   void SetK(const int & dim){numk = dim;}
   void SetMach(const double & dim){mach = dim;}
   void SetAlpha(const double & dim){alpha = dim;}
   void SetReynolds(const double & dim){Re = dim;}
   void SetTi(const double & dim){ti = dim;}
   void SetQ1(const vector<double> & data){q1 = data;}
   void SetQ2(const vector<double> & data){q2 = data;}
   void SetQ3(const vector<double> & data){q3 = data;}
   void SetQ4(const vector<double> & data){q4 = data;}
   void SetQ5(const vector<double> & data){q5 = data;}

   //destructor
   ~plot3dQBlock() {}

};

//---------------------------------------------------------------------------------------------------------------//
//Class for a plot3d mesh which consists of individual plot3d blocks
class plot3dMesh {
  //by default everything above public: is private
  vector<plot3dBlock> blocks;            //vector of plot3dBlocks

 public:
  //constructor -- create a plot3d mesh by passing a single plot3dBlock
  plot3dMesh( const plot3dBlock&);
  plot3dMesh();
  plot3dMesh(int);

  //member functions
  void AddP3dBlock(const plot3dBlock&); 
  void ReplaceBlock(int, const plot3dBlock&);
  plot3dBlock Blocks(const int &ii) const {return blocks[ii];}
  int NumBlocks() const {return blocks.size();}

  //destructor
  ~plot3dMesh() {}

};

//---------------------------------------------------------------------------------------------------------------//
//Class for a plot3d solution mesh which consists of individual plot3d solution blocks
class plot3dQMesh {
  //by default everything above public: is private
  vector<plot3dQBlock> qblocks;            //vector of plot3dQBlocks

 public:
  //constructor -- create a plot3d solution mesh by passing a plot3d mesh
  plot3dQMesh(const plot3dMesh&);

  //member functions
  vector<plot3dQBlock> QBlocks() const {return qblocks;}
  void WriteData() const;

  //destructor
  ~plot3dQMesh() {}

};


//function declarations
plot3dMesh ReadP3dGrid(const string &gridName, double &numCells);
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

int GetMatrixDiagUpperFromMainI(const int &);
int GetMatrixDiagLowerFromMainI(const int &);
int GetMatrixDiagUpperFromMainJ(const int &, const int &);
int GetMatrixDiagLowerFromMainJ(const int &, const int &);
int GetMatrixDiagUpperFromMainK(const int &, const int &, const int &);
int GetMatrixDiagLowerFromMainK(const int &, const int &, const int &);

int GetDiagPosUpperI(const int &);
int GetDiagPosLowerI(const int &);
int GetDiagPosUpperJ(const int &, const int &);
int GetDiagPosLowerJ(const int &, const int &);
int GetDiagPosUpperK(const int &, const int &, const int &);
int GetDiagPosLowerK(const int &, const int &, const int &);

//find out if matrix should have data at the indicated cell
bool IsMatrixData(const int&, const int&, const int&, const int&, const int&, const int&, const string&);

//function to reorder block by hyperplanes
vector<int> HyperplaneReorder(const int &, const int &, const int &);

#endif

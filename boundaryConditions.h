#ifndef BOUNDARYCONDITIONSHEADERDEF             //only if the macro BOUNDARYCONDITIONSHEADERDEF is not defined execute these lines of code
#define BOUNDARYCONDITIONSHEADERDEF             //define the macro

/* This header contains the class boundaryConditions.

This class stores the information needed to specify the boundary conditions for one block. */

#include <vector>  //vector
#include <string>  //string
#include <iostream> //ostream
#include "plot3d.h" //plot3dBlock
#include "mpi.h" //parallelism

using std::ostream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

// a class to store the necessary information for the boundary conditions of a block
class boundaryConditions {
  vector<string> bcTypes;                 //vector of boundary condition names for each surface
  vector<int> iMin;                       //vector of min i coordinate defining surface
  vector<int> iMax;                       //vector of max i coordinate defining surface
  vector<int> jMin;                       //vector of min j coordinate defining surface
  vector<int> jMax;                       //vector of max j coordinate defining surface
  vector<int> kMin;                       //vector of min k coordinate defining surface
  vector<int> kMax;                       //vector of max k coordinate defining surface
  vector<int> tag;                        //vector of boundary condition tags
  int numSurfI;                           //number of i-surfaces to define boundary on in block
  int numSurfJ;                           //number of j-surfaces to define boundary on in block
  int numSurfK;                           //number of k-surfaces to define boundary on in block

 public:
  //constructor
  boundaryConditions();
  boundaryConditions( const int&, const int&, const int& );

  //member functions
  int NumSurfI()const{return numSurfI;}
  int NumSurfJ()const{return numSurfJ;}
  int NumSurfK()const{return numSurfK;}
  int NumSurfaces()const{return numSurfI + numSurfJ + numSurfK;}

  string GetBCTypes( const int &a)const{return bcTypes[a];}
  int GetIMin( const int &a)const{return iMin[a];}
  int GetJMin( const int &a)const{return jMin[a];}
  int GetKMin( const int &a)const{return kMin[a];}
  int GetIMax( const int &a)const{return iMax[a];}
  int GetJMax( const int &a)const{return jMax[a];}
  int GetKMax( const int &a)const{return kMax[a];}
  int GetTag( const int &a)const{return tag[a];}

  void ResizeVecs( const int&);
  void ResizeVecs( const int&, const int&, const int& );

  friend ostream & operator<< (ostream &os, const boundaryConditions&);

  string GetBCName(const int, const int, const int, const string&)const;

  void AssignFromInput(const int&, const vector<string>& );

  boundaryConditions Split(const string&, const int&, const int&);
  void Join(const boundaryConditions&, const string&);

  void PackBC( char*(&), const int&, int&)const;
  void UnpackBC( char*(&), const int&, int&);

  //destructor
  ~boundaryConditions() {}

};

/* A class to store the necessary information for the boundary condition patches. A patch is a 2D surface on a block boundary that
is assigned the same boundary condition.
*/
class patch {
  vector3d<double> origin;              //coordinates of patch origin
  vector3d<double> corner1;             //coordinates of direction 1 max, direction 2 zero
  vector3d<double> corner2;             //coordinates of direction 1 zero, direction 2 max
  vector3d<double> corner12;            //coordinates of direction 1/2 max
  int boundary;                         //boundary number (1-6)
  int block;                            //parent block number
  int d1Start;                          //direction 1 start index
  int d1End;                            //direction 1 end index
  int d2Start;                          //direction 2 start index
  int d2End;                            //direction 2 end index
  int constSurf;                        //index of direction 3

 public:
  //constructor
  patch();
  patch(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const plot3dBlock&);

  //member functions
  vector3d<double> Origin()const{return origin;}
  vector3d<double> Corner1()const{return corner1;}
  vector3d<double> Corner2()const{return corner2;}
  vector3d<double> Corner12()const{return corner12;}

  int Boundary()const{return boundary;}
  int Block()const{return block;}
  int Dir1Start()const{return d1Start;}
  int Dir1End()const{return d1End;}
  int Dir2Start()const{return d2Start;}
  int Dir2End()const{return d2End;}
  int ConstSurface()const{return constSurf;}

  friend ostream & operator<< (ostream &os, const patch&);

  //destructor
  ~patch() {}

};

/* A class to store the necessary information for the interblock boundary conditions. The data is stored in pairs, where each
pair is patch on a boundary that is point matched.
*/
class interblock {

  int rank [2];             //processor location of boundaries
  int block [2];            //block numbers (global)
  int localBlock [2];       //local (on processor) block numbers
  int boundary [2];         //boundary numbers
  int d1Start [2];          //first direction start numbers for surface
  int d1End [2];            //first direction end numbers for surface
  int d2Start [2];          //second direction start numbers for surface
  int d2End [2];            //second direction end numbers for surface
  int constSurf [2];        //index of direction 3
  int orientation;          //defines how patches are oriented relative to one another (1-8)

 public:
  //constructor
 interblock() : rank{0,0}, block{0,0}, localBlock{0,0}, boundary{0,0}, d1Start{0,0}, d1End{0,0}, d2Start{0,0}, d2End{0,0}, constSurf{0,0}, orientation(0) {};
 interblock( const patch&, const patch& );

  //member functions
  int RankFirst()const{return rank[0];}
  void SetRankFirst(const int &a){rank[0] = a;} //setter needed to assign to processor during decomposition
  int RankSecond()const{return rank[1];}
  void SetRankSecond(const int &a){rank[1] = a;} //setter needed to assign to processor during decomposition

  int BlockFirst()const{return block[0];}
  int BlockSecond()const{return block[1];}

  int LocalBlockFirst()const{return localBlock[0];}
  int LocalBlockSecond()const{return localBlock[1];}

  int BoundaryFirst()const{return boundary[0];}
  int BoundarySecond()const{return boundary[1];}

  int Dir1StartFirst()const{return d1Start[0];}
  int Dir1StartSecond()const{return d1Start[1];}

  int Dir1EndFirst()const{return d1End[0];}
  int Dir1EndSecond()const{return d1End[1];}

  int Dir2StartFirst()const{return d2Start[0];}
  int Dir2StartSecond()const{return d2Start[1];}

  int Dir2EndFirst()const{return d2End[0];}
  int Dir2EndSecond()const{return d2End[1];}

  int ConstSurfaceFirst()const{return constSurf[0];}
  int ConstSurfaceSecond()const{return constSurf[1];}

  int Orientation()const{return orientation;}

  void SwapOrder();
  void AdjustForSlice( const bool&, const int& );
  bool TestPatchMatch(const patch&, const patch&);
  void GetAddressesMPI(MPI_Aint (&)[10])const;

  friend ostream & operator<< (ostream &os, const interblock&);

  //destructor
  ~interblock() {}

};


//function declarations
vector<interblock> GetInterblockBCs( const vector<boundaryConditions>&, const vector<plot3dBlock>& );



#endif

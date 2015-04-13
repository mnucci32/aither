/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef BOUNDARYCONDITIONSHEADERDEF             //only if the macro BOUNDARYCONDITIONSHEADERDEF is not defined execute these lines of code
#define BOUNDARYCONDITIONSHEADERDEF             //define the macro

/* This header contains the class boundaryConditions.

This class stores the information needed to specify the boundary conditions for one block. */

#include <vector>  //vector
#include <string>  //string
#include <iostream> //ostream
#include "plot3d.hpp" //plot3dBlock
#include "mpi.h" //parallelism

using std::ostream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

class boundarySurface {
  string bcType;                 //boundary condition name for surface
  int data [7];                  //data for boundary surface: imin, imax, jmin, jmax, kmin, kmax, tag

 public:
  //constructor
  boundarySurface();
  boundarySurface(const string&, const int&, const int&, const int&, const int&, const int&, const int&, const int&);

  friend class boundaryConditions;

  //member functions
  string BCType()const{return bcType;}
  int IMin()const{return data[0];}
  int IMax()const{return data[1];}
  int JMin()const{return data[2];}
  int JMax()const{return data[3];}
  int KMin()const{return data[4];}
  int KMax()const{return data[5];}
  int Tag()const{return data[6];}

  int SurfaceType()const;
  string Direction1()const;
  string Direction2()const;
  string Direction3()const;
  int Max1()const;
  int Max2()const;
  int Min1()const;
  int Min2()const;

  int PartnerBlock()const;
  int PartnerSurface()const;
  void UpdateTagForSplitJoin(const int&);
  boundarySurface Split(const string&, const int&, const int&, const int&, bool&, int=0);
  bool SplitDirectionIsReversed(const string&, const int&)const;

  friend ostream & operator<< (ostream &os, const boundarySurface&);

  //destructor
  ~boundarySurface() {}

};

/* A class to store the necessary information for the boundary condition patches. A patch is a 2D surface on a block boundary that
is assigned the same boundary condition.
*/
class patch {
  vector3d<double> origin;              //coordinates of patch origin
  vector3d<double> corner1;             //coordinates of direction 1 max, direction 2 zero
  vector3d<double> corner2;             //coordinates of direction 1 zero, direction 2 max
  vector3d<double> corner12;            //coordinates of direction 1/2 max
  bool interblockBorder[4];             //array of booleans for 4 sides of patch (1 if borders an interblock)
  int boundary;                         //boundary number (1-6)
  int block;                            //parent block number
  int d1Start;                          //direction 1 start index
  int d1End;                            //direction 1 end index
  int d2Start;                          //direction 2 start index
  int d2End;                            //direction 2 end index
  int constSurf;                        //index of direction 3
  int rank;                             //rank of block that patch belongs to
  int localBlock;                       //position of block on processor

 public:
  //constructor
  patch();
  patch(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const plot3dBlock&, const int&, const int&, const bool(&)[4]);
  patch(const boundarySurface&, const plot3dBlock&, const int&, const bool(&)[4], int=0, int=0);

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
  int Rank()const{return rank;}
  int LocalBlock()const{return localBlock;}
  bool Dir1StartInterBorder()const{return interblockBorder[0];}
  bool Dir1EndInterBorder()const{return interblockBorder[1];}
  bool Dir2StartInterBorder()const{return interblockBorder[2];}
  bool Dir2EndInterBorder()const{return interblockBorder[3];}

  friend ostream & operator<< (ostream &os, const patch&);

  //destructor
  ~patch() {}

};

// a class to store the necessary information for the boundary conditions of a block
class boundaryConditions {
  vector<boundarySurface> surfs;          //vector of boundary condition surfaces defining block
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

  string GetBCTypes( const int &a)const{return surfs[a].BCType();}
  int GetIMin( const int &a)const{return surfs[a].IMin();}
  int GetJMin( const int &a)const{return surfs[a].JMin();}
  int GetKMin( const int &a)const{return surfs[a].KMin();}
  int GetIMax( const int &a)const{return surfs[a].IMax();}
  int GetJMax( const int &a)const{return surfs[a].JMax();}
  int GetKMax( const int &a)const{return surfs[a].KMax();}
  int GetTag( const int &a)const{return surfs[a].Tag();}
  boundarySurface GetSurface( const int &a)const{return surfs[a];}

  int BlockDimI()const;
  int BlockDimJ()const;
  int BlockDimK()const;

  void ResizeVecs( const int&);
  void ResizeVecs( const int&, const int&, const int& );

  friend ostream & operator<< (ostream &os, const boundaryConditions&);

  string GetBCName(const int, const int, const int, const string&)const;

  void AssignFromInput(const int&, const vector<string>& );

  boundaryConditions Split(const string&, const int&, const int&, const int&, vector<boundarySurface>&);
  void DependentSplit(const boundarySurface&, const plot3dBlock&, const plot3dBlock&, const int&, const string&, const int&, const int&, const int&);
  void Join(const boundaryConditions&, const string&, vector<boundarySurface>&);

  void BordersInterblock(const int&, bool (&)[4])const;

  void PackBC( char*(&), const int&, int&)const;
  void UnpackBC( char*(&), const int&, int&);

  //destructor
  ~boundaryConditions() {}

};

/* A class to store the necessary information for the interblock boundary conditions. The data is stored in pairs, where each
pair is patch on a boundary that is point matched.
*/
class interblock {

  int rank [2];              //processor location of boundaries
  int block [2];             //block numbers (global)
  int localBlock [2];        //local (on processor) block numbers
  int boundary [2];          //boundary numbers
  int d1Start [2];           //first direction start numbers for surface
  int d1End [2];             //first direction end numbers for surface
  int d2Start [2];           //second direction start numbers for surface
  int d2End [2];             //second direction end numbers for surface
  int constSurf [2];         //index of direction 3
  bool interblockBorder [8]; //borders interblock on sides of patch
  int orientation;           //defines how patches are oriented relative to one another (1-8)

 public:
  //constructor
 interblock() : rank{0,0}, block{0,0}, localBlock{0,0}, boundary{0,0}, d1Start{0,0}, d1End{0,0}, d2Start{0,0}, d2End{0,0}, constSurf{0,0}, 
	        interblockBorder{false, false, false, false, false, false, false, false}, orientation(0) {};
 interblock( const patch&, const patch& );

  //member functions
  int RankFirst()const{return rank[0];}
  int RankSecond()const{return rank[1];}

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

  bool Dir1StartInterBorderFirst()const{return interblockBorder[0];}
  bool Dir1EndInterBorderFirst()const{return interblockBorder[1];}
  bool Dir2StartInterBorderFirst()const{return interblockBorder[2];}
  bool Dir2EndInterBorderFirst()const{return interblockBorder[3];}
  bool Dir1StartInterBorderSecond()const{return interblockBorder[4];}
  bool Dir1EndInterBorderSecond()const{return interblockBorder[5];}
  bool Dir2StartInterBorderSecond()const{return interblockBorder[6];}
  bool Dir2EndInterBorderSecond()const{return interblockBorder[7];}

  int Orientation()const{return orientation;}

  string Direction1First()const;
  string Direction2First()const;
  string Direction3First()const;
  string Direction1Second()const;
  string Direction2Second()const;
  string Direction3Second()const;

  void UpdateBorderFirst(const int&);
  void UpdateBorderSecond(const int&);
  void SwapOrder();
  void AdjustForSlice( const bool&, const int& );
  bool TestPatchMatch(const patch&, const patch&);
  void GetAddressesMPI(MPI_Aint (&)[11])const;

  friend ostream & operator<< (ostream &os, const interblock&);

  //destructor
  ~interblock() {}

};

class decomposition {
  vector<int> rank;                 //rank of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> parBlock;             //parent block of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> localPos;             //local position of each procBlock (vector size equals number of procBlocks after decomp)
  vector<int> splitHistBlkLow;      //lower block of split (vector size equals number of splits)
  vector<int> splitHistBlkUp;       //upper block of split (vector size equals number of splits)
  vector<int> splitHistIndex;       //index of split (vector size equals number of splits)
  vector<string> splitHistDir;      //direction of split (vector size equals number of splits)
  int numProcs;                     //number of processors

 public:
  //constructor
  decomposition();
  decomposition( const int&, const int& );

  //member functions
  int Rank( const int &a)const{return rank[a];}
  int ParentBlock( const int &a)const{return parBlock[a];}
  int LocalPosition( const int &a)const{return localPos[a];}
  double IdealLoad( const vector<plot3dBlock>& )const;
  double MaxLoad( const vector<plot3dBlock>&)const;
  double MinLoad( const vector<plot3dBlock>&)const;
  double ProcLoad( const vector<plot3dBlock>&, const int& )const;
  double LoadRatio( const vector<plot3dBlock>&, const int& )const;
  int MostOverloadedProc( const vector<plot3dBlock>&, double& )const;
  int MostUnderloadedProc( const vector<plot3dBlock>&, double& )const;
  int NumBlocksOnProc( const int&)const;
  vector<int> NumBlocksOnAllProc()const;
  void SendToProc( const int&, const int&, const int&);
  void Split( const int&, const int&, const string& );
  int SendWholeOrSplit( const vector<plot3dBlock>&, const int&, const int&, int&, string& )const;
  int Size()const{return (int)rank.size();}

  int NumSplits()const{return (int)splitHistDir.size();}
  int SplitHistBlkLower(const int &a)const{return splitHistBlkLow[a];}
  int SplitHistBlkUpper(const int &a)const{return splitHistBlkUp[a];}
  int SplitHistIndex(const int &a)const{return splitHistIndex[a];}
  string SplitHistDir(const int &a)const{return splitHistDir[a];}

  void PrintDiagnostics( const vector<plot3dBlock>& )const;

  friend ostream & operator<< (ostream &os, const decomposition&);

  //destructor
  ~decomposition() {}

};


//function declarations
vector<interblock> GetInterblockBCs( const vector<boundaryConditions>&, const vector<plot3dBlock>&, const decomposition& );



#endif

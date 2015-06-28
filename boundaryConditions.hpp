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

// Only if the macro BOUNDARYCONDITIONSHEADERDEF is not defined
// execute these lines of code
#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF   // define the macro

/* This header contains the class boundaryConditions.

This class stores the information needed to specify the boundary_ conditions for one block_. */

#include <vector>  // vector
#include <string>  // string
#include <iostream>  // ostream
#include "mpi.h"  // parallelism
#include "vector3d.hpp"

using std::ostream;
using std::vector;
using std::string;
using std::cout;
using std::endl;

// forward class declaration
class plot3dBlock;

class boundarySurface {
  string bcType_;    // boundary_ condition name for surface
  // data for boundary_ surface: imin, imax, jmin, jmax, kmin, kmax, tag
  int data_[7];

 public:
  // Constructor
  boundarySurface();
  boundarySurface(const string&, const int&, const int&, const int&,
                  const int&, const int&, const int&, const int&);

  friend class boundaryConditions;

  // Member functions
  string BCType() const {return bcType_;}
  int IMin() const {return data_[0];}
  int IMax() const {return data_[1];}
  int JMin() const {return data_[2];}
  int JMax() const {return data_[3];}
  int KMin() const {return data_[4];}
  int KMax() const {return data_[5];}
  int Tag() const {return data_[6];}

  int SurfaceType() const;
  string Direction1() const;
  string Direction2() const;
  string Direction3() const;
  int Max1() const;
  int Max2() const;
  int Min1() const;
  int Min2() const;

  int PartnerBlock() const;
  int PartnerSurface() const;
  void UpdateTagForSplitJoin(const int&);
  boundarySurface Split(const string&, const int&, const int&,
                        const int&, bool&, int = 0);
  bool SplitDirectionIsReversed(const string&, const int&) const;

  friend ostream & operator<< (ostream &os, const boundarySurface&);

  // Destructor
  ~boundarySurface() {}
};

/* A class to store the necessary information for the boundary_ condition patches.
   A patch is a 2D surface on a block_ boundary_ that is assigned the same
   boundary_ condition. */
class patch {
  vector3d<double> origin_;      // coordinates of patch origin
  // Coordinates of direction 1 max, direction 2 zero
  vector3d<double> corner1_;
  // Coordinates of direction 1 zero, direction 2 max
  vector3d<double> corner2_;
  vector3d<double> corner12_;    // coordinates of direction 1/2 max
  // Array of booleans for 4 sides of patch (1 if borders an interblock)
  bool interblockBorder_[4];
  int boundary_;                 // boundary number (1-6)
  int block_;                    // parent block number
  int d1Start_;                  // direction 1 start index
  int d1End_;                    // direction 1 end index
  int d2Start_;                  // direction 2 start index
  int d2End_;                    // direction 2 end index
  int constSurf_;                // index of direction 3
  int rank_;                     // rank of block that patch belongs to
  int localBlock_;               // position of block on processor

 public:
  // Constructor
  patch();
  patch(const int&, const int&, const int&, const int&, const int&, const int&,
        const int&, const int&, const plot3dBlock&, const int&, const int&,
        const bool(&)[4]);
  patch(const boundarySurface&, const plot3dBlock&, const int&,
        const bool(&)[4], int = 0, int = 0);

  // Member functions
  vector3d<double> Origin() const {return origin_;}
  vector3d<double> Corner1() const {return corner1_;}
  vector3d<double> Corner2() const {return corner2_;}
  vector3d<double> Corner12() const {return corner12_;}
  int Boundary() const {return boundary_;}
  int Block() const {return block_;}
  int Dir1Start() const {return d1Start_;}
  int Dir1End() const {return d1End_;}
  int Dir2Start() const {return d2Start_;}
  int Dir2End() const {return d2End_;}
  int ConstSurface() const {return constSurf_;}
  int Rank() const {return rank_;}
  int LocalBlock() const {return localBlock_;}
  bool Dir1StartInterBorder() const {return interblockBorder_[0];}
  bool Dir1EndInterBorder() const {return interblockBorder_[1];}
  bool Dir2StartInterBorder() const {return interblockBorder_[2];}
  bool Dir2EndInterBorder() const {return interblockBorder_[3];}

  friend ostream & operator<< (ostream &os, const patch&);

  // Destructor
  ~patch() {}
};

// A class to store the necessary information for the boundary
// conditions of a block
class boundaryConditions {
  // Vector of boundary condition surfaces defining block
  vector<boundarySurface> surfs_;
  int numSurfI_;        // number of i-surfaces to define boundary on block
  int numSurfJ_;        // number of j-surfaces to define boundary on block
  int numSurfK_;        // number of k-surfaces to define boundary on block

 public:
  // Constructor
  boundaryConditions();
  boundaryConditions(const int&, const int&, const int&);

  // Member functions
  int NumSurfI() const {return numSurfI_;}
  int NumSurfJ() const {return numSurfJ_;}
  int NumSurfK() const {return numSurfK_;}
  int NumSurfaces() const {return numSurfI_ + numSurfJ_ + numSurfK_;}

  string GetBCTypes(const int &a) const {return surfs_[a].BCType();}
  int GetIMin(const int &a) const {return surfs_[a].IMin();}
  int GetJMin(const int &a) const {return surfs_[a].JMin();}
  int GetKMin(const int &a) const {return surfs_[a].KMin();}
  int GetIMax(const int &a) const {return surfs_[a].IMax();}
  int GetJMax(const int &a) const {return surfs_[a].JMax();}
  int GetKMax(const int &a) const {return surfs_[a].KMax();}
  int GetTag(const int &a) const {return surfs_[a].Tag();}
  boundarySurface GetSurface(const int &a) const {return surfs_[a];}

  int BlockDimI() const;
  int BlockDimJ() const;
  int BlockDimK() const;

  void ResizeVecs(const int&);
  void ResizeVecs(const int&, const int&, const int&);

  friend ostream & operator<< (ostream &os, const boundaryConditions&);

  string GetBCName(const int, const int, const int, const string&) const;

  void AssignFromInput(const int&, const vector<string>&);

  boundaryConditions Split(const string&, const int&, const int&,
                           const int&, vector<boundarySurface>&);
  void DependentSplit(const boundarySurface&, const plot3dBlock&,
                      const plot3dBlock&, const int&, const string&,
                      const int&, const int&, const int&);
  void Join(const boundaryConditions&, const string&, vector<boundarySurface>&);

  void BordersInterblock(const int&, bool (&)[4]) const;

  void PackBC(char*(&), const int&, int&) const;
  void UnpackBC(char*(&), const int&, int&);

  // Destructor
  ~boundaryConditions() {}
};

/* A class to store the necessary information for the interblock boundary_ conditions.
   The data_ is stored in pairs, where each pair is patch on a boundary_ that is point matched. */
class interblock {
  int rank_[2];               // processor location of boundaries
  int block_[2];              // block_ numbers (global)
  int localBlock_[2];         // local (on processor) block_ numbers
  int boundary_[2];           // boundary_ numbers
  int d1Start_[2];            // first direction start numbers for surface
  int d1End_[2];              // first direction end numbers for surface
  int d2Start_[2];            // second direction start numbers for surface
  int d2End_[2];              // second direction end numbers for surface
  int constSurf_[2];          // index of direction 3
  bool interblockBorder_[8];  // borders interblock on sides of patch
  // Defines how patches are oriented relative to one another (1-8)
  int orientation_;

 public:
  // Constructor
  interblock() : rank_{0, 0}, block_{0, 0}, localBlock_{0, 0}, boundary_{0, 0},
    d1Start_{0, 0}, d1End_{0, 0}, d2Start_{0, 0}, d2End_{0, 0},
    constSurf_{0, 0}, interblockBorder_{false, false, false, false, false,
                            false, false, false}, orientation_(0) {}

  interblock(const patch&, const patch&);

  // Member functions
  int RankFirst() const {return rank_[0];}
  int RankSecond() const {return rank_[1];}

  int BlockFirst() const {return block_[0];}
  int BlockSecond() const {return block_[1];}

  int LocalBlockFirst() const {return localBlock_[0];}
  int LocalBlockSecond() const {return localBlock_[1];}

  int BoundaryFirst() const {return boundary_[0];}
  int BoundarySecond() const {return boundary_[1];}

  int Dir1StartFirst() const {return d1Start_[0];}
  int Dir1StartSecond() const {return d1Start_[1];}

  int Dir1EndFirst() const {return d1End_[0];}
  int Dir1EndSecond() const {return d1End_[1];}

  int Dir2StartFirst() const {return d2Start_[0];}
  int Dir2StartSecond() const {return d2Start_[1];}

  int Dir2EndFirst() const {return d2End_[0];}
  int Dir2EndSecond() const {return d2End_[1];}

  int ConstSurfaceFirst() const {return constSurf_[0];}
  int ConstSurfaceSecond() const {return constSurf_[1];}

  bool Dir1StartInterBorderFirst() const {return interblockBorder_[0];}
  bool Dir1EndInterBorderFirst() const {return interblockBorder_[1];}
  bool Dir2StartInterBorderFirst() const {return interblockBorder_[2];}
  bool Dir2EndInterBorderFirst() const {return interblockBorder_[3];}
  bool Dir1StartInterBorderSecond() const {return interblockBorder_[4];}
  bool Dir1EndInterBorderSecond() const {return interblockBorder_[5];}
  bool Dir2StartInterBorderSecond() const {return interblockBorder_[6];}
  bool Dir2EndInterBorderSecond() const {return interblockBorder_[7];}

  int Orientation() const {return orientation_;}

  string Direction1First() const;
  string Direction2First() const;
  string Direction3First() const;
  string Direction1Second() const;
  string Direction2Second() const;
  string Direction3Second() const;

  void UpdateBorderFirst(const int&);
  void UpdateBorderSecond(const int&);
  void SwapOrder();
  void AdjustForSlice(const bool&, const int&);
  bool TestPatchMatch(const patch&, const patch&);
  void GetAddressesMPI(MPI_Aint (&)[11])const;

  friend ostream & operator<< (ostream &os, const interblock&);

  // Destructor
  ~interblock() {}
};

class decomposition {
  // rank_ of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> rank_;
  // parent block_ of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> parBlock_;
  // local position of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> localPos_;
  // lower block_ of split (vector size equals number of splits)
  vector<int> splitHistBlkLow_;
  // upper block_ of split (vector size equals number of splits)
  vector<int> splitHistBlkUp_;
  // index of split (vector size equals number of splits)
  vector<int> splitHistIndex_;
  // direction of split (vector size equals number of splits)
  vector<string> splitHistDir_;
  int numProcs;                     // number of processors

 public:
  // Constructor
  decomposition();
  decomposition(const int&, const int&);

  // Member functions
  int Rank(const int &a) const {return rank_[a];}
  int ParentBlock(const int &a) const {return parBlock_[a];}
  int LocalPosition(const int &a) const {return localPos_[a];}
  double IdealLoad(const vector<plot3dBlock>&) const;
  double MaxLoad(const vector<plot3dBlock>&) const;
  double MinLoad(const vector<plot3dBlock>&) const;
  double ProcLoad(const vector<plot3dBlock>&, const int&) const;
  double LoadRatio(const vector<plot3dBlock>&, const int&) const;
  int MostOverloadedProc(const vector<plot3dBlock>&, double&) const;
  int MostUnderloadedProc(const vector<plot3dBlock>&, double&) const;
  int NumBlocksOnProc(const int&) const;
  vector<int> NumBlocksOnAllProc() const;
  void SendToProc(const int&, const int&, const int&);
  void Split(const int&, const int&, const string&);
  int SendWholeOrSplit(const vector<plot3dBlock>&, const int&,
                       const int&, int&, string&) const;
  int Size() const {return static_cast<int> (rank_.size());}

  int NumSplits() const {return static_cast<int> (splitHistDir_.size());}
  int SplitHistBlkLower(const int &a) const {return splitHistBlkLow_[a];}
  int SplitHistBlkUpper(const int &a) const {return splitHistBlkUp_[a];}
  int SplitHistIndex(const int &a) const {return splitHistIndex_[a];}
  string SplitHistDir(const int &a) const {return splitHistDir_[a];}

  void PrintDiagnostics(const vector<plot3dBlock>&) const;

  friend ostream & operator<< (ostream &os, const decomposition&);

  // Destructor
  ~decomposition() {}
};


// Function declarations
vector<interblock> GetInterblockBCs(const vector<boundaryConditions>&,
                                    const vector<plot3dBlock>&,
                                    const decomposition&);


#endif

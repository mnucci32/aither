/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
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

This class stores the information needed to specify the boundary conditions for
one block. */

#include <vector>  // vector
#include <array>   // array
#include <string>  // string
#include <iostream>  // ostream
#include <memory>  // unique_ptr
#include "mpi.h"  // parallelism
#include "vector3d.hpp"
#include "range.hpp"

using std::ostream;
using std::vector;
using std::array;
using std::string;
using std::cout;
using std::endl;
using std::shared_ptr;

// forward class declaration
class plot3dBlock;
class decomposition;
class inputState;
class input;

class boundarySurface {
  string bcType_;    // boundary condition name for surface
  // data for boundary surface: imin, imax, jmin, jmax, kmin, kmax, tag
  int data_[7];

 public:
  // Constructor
  boundarySurface(const string &name, const int &imin, const int &imax,
                  const int &jmin, const int &jmax, const int &kmin,
                  const int &kmax, const int &tag) :
      bcType_(name), data_{imin, imax, jmin, jmax, kmin, kmax, tag} {}
  boundarySurface() : boundarySurface("undefined", 0, 0, 0, 0, 0, 0, 0) {}

  // move constructor and assignment operator
  boundarySurface(boundarySurface&&) noexcept = default;
  boundarySurface& operator=(boundarySurface&&) noexcept = default;

  // copy constructor and assignment operator
  boundarySurface(const boundarySurface&) = default;
  boundarySurface& operator=(const boundarySurface&) = default;

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
  int NumFaces() const {return this->NumI() * this->NumJ() * this->NumK();}

  int NumI() const {return this->RangeI().Size();}
  int NumJ() const {return this->RangeJ().Size();}
  int NumK() const {return this->RangeK().Size();}

  range RangeI() const;
  range RangeJ() const;
  range RangeK() const;
  range RangeDir1() const;
  range RangeDir2() const;
  range RangeDir3() const;

  int PartnerBlock() const;
  int PartnerSurface() const;
  void UpdateTagForSplitJoin(const int&);
  boundarySurface Split(const string&, const int&, const int&,
                        const int&, bool&, int = 0);
  bool SplitDirectionIsReversed(const string&, const int&) const;
  bool IsConnection() const {
    return bcType_ == "interblock" || bcType_ == "periodic";
  }

  void PackBoundarySurface(char*(&), const int&, int&) const;
  void UnpackBoundarySurface(char*(&), const int&, int&);

  bool operator==(const boundarySurface &) const;

  // Destructor
  ~boundarySurface() noexcept {}
};

/* A class to store the necessary information for the boundary condition patches.
   A patch is a 2D surface on a block boundary that is assigned the same
   boundary condition. */
class patch {
  vector3d<double> origin_;      // coordinates of patch origin
  // Coordinates of direction 1 max, direction 2 zero
  vector3d<double> corner1_;
  // Coordinates of direction 1 zero, direction 2 max
  vector3d<double> corner2_;
  vector3d<double> corner12_;    // coordinates of direction 1/2 max
  // Array of booleans for 4 sides of patch (1 if borders another patch)
  array<bool, 4> patchBorder_;
  int boundary_;                 // boundary number (1-6)
  int block_;                    // parent block number
  int d1Start_;                  // direction 1 start index
  int d1End_;                    // direction 1 end index
  int d2Start_;                  // direction 2 start index
  int d2End_;                    // direction 2 end index
  int constSurf_;                // index of direction 3
  int rank_;                     // rank of block that patch belongs to
  int localBlock_;               // position of block on processor
  string bcName_;                // name of bc on patch

 public:
  // Constructor
  patch();
  patch(const int&, const int&, const int&, const int&, const int&, const int&,
        const int&, const int&, const plot3dBlock&, const int&, const int&,
        const array<bool, 4>&, const string&);
  patch(const boundarySurface &surf, const plot3dBlock &blk, const int &bNum,
        const array<bool, 4> &border, int r = 0, int l = 0) :
      patch(surf.SurfaceType(), bNum, surf.IMin(), surf.IMax(), surf.JMin(),
            surf.JMax(), surf.KMin(), surf.KMax(), blk, r, l, border,
            surf.BCType()) {}

  // move constructor and assignment operator
  patch(patch&&) noexcept = default;
  patch& operator=(patch&&) noexcept = default;

  // copy constructor and assignment operator
  patch(const patch&) = default;
  patch& operator=(const patch&) = default;

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
  bool Dir1StartInterBorder() const {return patchBorder_[0];}
  bool Dir1EndInterBorder() const {return patchBorder_[1];}
  bool Dir2StartInterBorder() const {return patchBorder_[2];}
  bool Dir2EndInterBorder() const {return patchBorder_[3];}
  string BCType() const {return bcName_;}
  void Transform(const shared_ptr<inputState> &);
  void Translate(const vector3d<double> &);
  void Rotate(const vector3d<double> &, const vector3d<double> &,
              const double &);

  // Destructor
  ~patch() noexcept {}
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
  boundaryConditions(const int&, const int&, const int&);
  boundaryConditions() : boundaryConditions(2, 2, 2) {}

  // move constructor and assignment operator
  boundaryConditions(boundaryConditions&&) noexcept = default;
  boundaryConditions& operator=(boundaryConditions&&) noexcept = default;

  // copy constructor and assignment operator
  boundaryConditions(const boundaryConditions&) = default;
  boundaryConditions& operator=(const boundaryConditions&) = default;

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
  int GetSurfaceType(const int &a) const {return surfs_[a].SurfaceType();}
  boundarySurface GetSurface(const int &a) const {return surfs_[a];}
  int NumViscousFaces() const;
  int NumViscousSurfaces() const;
  bool IsConnection(const int &a) const {return surfs_[a].IsConnection();}

  int BlockDimI() const;
  int BlockDimJ() const;
  int BlockDimK() const;

  range RangeI(const int &a) const {return surfs_[a].RangeI();}
  range RangeJ(const int &a) const {return surfs_[a].RangeJ();}
  range RangeK(const int &a) const {return surfs_[a].RangeK();}
  range RangeDir1(const int &a) const {return surfs_[a].RangeDir1();}
  range RangeDir2(const int &a) const {return surfs_[a].RangeDir2();}
  range RangeDir3(const int &a) const {return surfs_[a].RangeDir3();}

  string Direction1(const int &a) const {return surfs_[a].Direction1();}
  string Direction2(const int &a) const {return surfs_[a].Direction2();}
  string Direction3(const int &a) const {return surfs_[a].Direction3();}

  void ResizeVecs(const int&);
  void ResizeVecs(const int&, const int&, const int&);

  boundarySurface GetBCSurface(const int &, const int &, const int &,
                               const int &) const;
  string GetBCName(const int &ii, const int &jj, const int &kk,
                   const int &surf) const {
    return this->GetBCSurface(ii, jj, kk, surf).BCType();
  }
  bool BCIsConnection(const int &ii, const int &jj, const int &kk,
                      const int &surf) const {
    const auto name = this->GetBCName(ii, jj, kk, surf);
    return name == "interblock" || name == "periodic";
  }

  void AssignFromInput(const int&, const vector<string>&);

  boundaryConditions Split(const string&, const int&, const int&,
                           const int&, vector<boundarySurface>&);
  void DependentSplit(const boundarySurface&, const plot3dBlock&,
                      const plot3dBlock&, const int&, const string&,
                      const int&, const int&, const int&);
  void Join(const boundaryConditions&, const string&, vector<boundarySurface>&);

  void BordersSurface(const int&, array<bool, 4>&) const;

  void PackBC(char*(&), const int&, int&) const;
  void UnpackBC(char*(&), const int&, int&);

  // Destructor
  ~boundaryConditions() noexcept {}
};

/* A class to store the necessary information for the connection boundary conditions.
   The data is stored in pairs, where each pair is patch on a boundary that is point matched. */
class connection {
  int rank_[2];               // processor location of boundaries
  int block_[2];              // block numbers (global)
  int localBlock_[2];         // local (on processor) block numbers
  int boundary_[2];           // boundary numbers
  int d1Start_[2];            // first direction start numbers for surface
  int d1End_[2];              // first direction end numbers for surface
  int d2Start_[2];            // second direction start numbers for surface
  int d2End_[2];              // second direction end numbers for surface
  int constSurf_[2];          // index of direction 3
  array<bool, 8> patchBorder_;  // borders another patch on sides of patch
  // Defines how patches are oriented relative to one another (1-8)
  int orientation_;
  bool isInterblock_;

 public:
  // Constructor
  connection() : rank_{0, 0}, block_{0, 0}, localBlock_{0, 0}, boundary_{0, 0},
                 d1Start_{0, 0}, d1End_{0, 0}, d2Start_{0, 0}, d2End_{0, 0},
                 constSurf_{0, 0}, patchBorder_{false, false, false, false, false,
                                         false, false, false}, orientation_(0),
                 isInterblock_(true) {}

  connection(const patch&, const patch&);

  // move constructor and assignment operator
  connection(connection&&) noexcept = default;
  connection& operator=(connection&&) noexcept = default;

  // copy constructor and assignment operator
  connection(const connection&) = default;
  connection& operator=(const connection&) = default;

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

  int Dir1LenFirst() const {return d1End_[0] - d1Start_[0];}
  int Dir1LenSecond() const {return d1End_[1] - d1Start_[1];}

  int Dir2StartFirst() const {return d2Start_[0];}
  int Dir2StartSecond() const {return d2Start_[1];}

  int Dir2EndFirst() const {return d2End_[0];}
  int Dir2EndSecond() const {return d2End_[1];}

  int Dir2LenFirst() const {return d2End_[0] - d2Start_[0];}
  int Dir2LenSecond() const {return d2End_[1] - d2Start_[1];}

  range Dir1RangeFirst() const;
  range Dir1RangeSecond() const;
  range Dir2RangeFirst() const;
  range Dir2RangeSecond() const;

  int ConstSurfaceFirst() const {return constSurf_[0];}
  int ConstSurfaceSecond() const {return constSurf_[1];}

  bool IsLowerFirst() const {return constSurf_[0] == 0;}
  bool IsLowerSecond() const {return constSurf_[1] == 0;}

  bool Dir1StartInterBorderFirst() const {return patchBorder_[0];}
  bool Dir1EndInterBorderFirst() const {return patchBorder_[1];}
  bool Dir2StartInterBorderFirst() const {return patchBorder_[2];}
  bool Dir2EndInterBorderFirst() const {return patchBorder_[3];}
  bool Dir1StartInterBorderSecond() const {return patchBorder_[4];}
  bool Dir1EndInterBorderSecond() const {return patchBorder_[5];}
  bool Dir2StartInterBorderSecond() const {return patchBorder_[6];}
  bool Dir2EndInterBorderSecond() const {return patchBorder_[7];}

  int Orientation() const {return orientation_;}
  bool IsInterblock() const {return isInterblock_;}

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
  void GetAddressesMPI(MPI_Aint (&)[12])const;

  void FirstSliceIndices(int&, int&, int&, int&, int&, int&, const int&) const;
  void SecondSliceIndices(int&, int&, int&, int&, int&, int&, const int&) const;

  bool IsLowerLowerOrUpperUpper() const {
    return (boundary_[0] + boundary_[1]) % 2 == 0;
  }

  // Destructor
  ~connection() noexcept {}
};


// Function declarations
vector<connection> GetConnectionBCs(const vector<boundaryConditions>&,
                                    const vector<plot3dBlock>&,
                                    const decomposition&, const input&);

ostream & operator<< (ostream &os, const boundaryConditions&);
ostream & operator<< (ostream &os, const boundarySurface&);
ostream & operator<< (ostream &os, const patch&);
ostream & operator<< (ostream &os, const connection&);


array<int, 3> GetSwapLoc(const int &, const int &, const int &, const int &,
                         const connection &, const int &, const bool &);

#endif

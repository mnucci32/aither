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

#include <math.h>     // sqrt
#include <iostream>   // cout
#include <iomanip>    // setw
#include <algorithm>  // max
#include <string>     // string
#include <vector>     // vector
#include "boundaryConditions.hpp"
#include "vector3d.hpp"  // vector3d
#include "plot3d.hpp"  // plot3dBlock
#include "parallel.hpp"  // decomposition
#include "inputStates.hpp"  // inputState
#include "input.hpp"  // input

using std::cout;
using std::endl;
using std::cerr;
using std::swap;

// Constructor when passed number of i, j, k surfaces
boundaryConditions::boundaryConditions(const int &i, const int &j,
                                       const int &k) {
  // i -- number of i surfaces
  // j -- number of j surfaces
  // k -- number of k surfaces

  numSurfI_ = i;
  numSurfJ_ = j;
  numSurfK_ = k;
  const auto length = numSurfI_ + numSurfJ_ + numSurfK_;

  boundarySurface bcSurf_;
  vector<boundarySurface> dumVec(length, bcSurf_);
  surfs_ = dumVec;
}

// Operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const boundaryConditions &bc) {
  // bc -- boundary conditions to print
  // os -- ostream to print to

  os << "Number of surfaces (I, J, K): " << bc.NumSurfI() << ", " <<
      bc.NumSurfJ() << ", " << bc.NumSurfK() << endl;

  for (auto ii = 0; ii < bc.NumSurfaces(); ii++) {
    os << bc.GetSurface(ii) << endl;
  }
  return os;
}

// Operator to resize all of the vector components of the
// boundary conditions class
void boundaryConditions::ResizeVecs(const int &a) {
  // a -- total number of surfaces (i+j+k)
  surfs_.resize(a);
}

// Operator to resize all of the vector components of the
// boundary conditions class
void boundaryConditions::ResizeVecs(const int &i, const int &j, const int &k) {
  // i -- number of i surfaces
  // j -- number of j surfaces
  // k -- number of k surfaces

  numSurfI_ = i;
  numSurfJ_ = j;
  numSurfK_ = k;

  surfs_.resize(i + j + k);
}

// Member function to return the boundary condition type given the
// i,j,k face coordinates and the surface type
string boundaryConditions::GetBCName(const int &i, const int &j, const int &k,
                                     const int &surf) const {
  // ii -- i coordinate
  // jj -- j coordinate
  // kk -- k coordinate
  // surf -- boundary condition surface type [1-6]

  auto iStart = 0;
  auto iEnd = 0;

  // i-surfaces search between 0 and number of i-surfaces
  if (surf == 1 || surf == 2) {
    iStart = 0;
    iEnd = this->NumSurfI();
    // j-surfaces search between end of i-surfaces and end of j-surfaces
  } else if (surf == 3 || surf == 4) {
    iStart = this->NumSurfI();
    iEnd = iStart + this->NumSurfJ();
    // k-surfaces search between end of j-surfaces and end of k-surfaces
  } else if (surf == 5 || surf == 6) {
    iStart = this->NumSurfI() + this->NumSurfJ();
    iEnd = iStart + this->NumSurfK();
  } else {
    cerr << "ERROR: Surface type " << surf << " is not recognized!" << endl;
  }

  string bcName = "undefined";

  // Determine which boundary condition should be applied
  for (auto nn = iStart; nn < iEnd; nn++) {
    // Determine which boundary given i, j, k coordinates apply to
    if ((i >= this->GetIMin(nn) && i <= this->GetIMax(nn) &&
         j >= this->GetJMin(nn) && j <= this->GetJMax(nn) &&
         k >= this->GetKMin(nn) && k <= this->GetKMax(nn))) {
      bcName = this->GetBCTypes(nn);
      break;
    }
  }

  return bcName;
}


// Member function to return the boundary condition tag given the
// i,j,k face coordinates and the surface type
int boundaryConditions::GetBCTag(const int &i, const int &j, const int &k,
                                 const int &surf) const {
  // ii -- i coordinate
  // jj -- j coordinate
  // kk -- k coordinate
  // surf -- boundary condition surface type [1-6]

  auto iStart = 0;
  auto iEnd = 0;

  // i-surfaces search between 0 and number of i-surfaces
  if (surf == 1 || surf == 2) {
    iStart = 0;
    iEnd = this->NumSurfI();
    // j-surfaces search between end of i-surfaces and end of j-surfaces
  } else if (surf == 3 || surf == 4) {
    iStart = this->NumSurfI();
    iEnd = iStart + this->NumSurfJ();
    // k-surfaces search between end of j-surfaces and end of k-surfaces
  } else if (surf == 5 || surf == 6) {
    iStart = this->NumSurfI() + this->NumSurfJ();
    iEnd = iStart + this->NumSurfK();
  } else {
    cerr << "ERROR: Surface type " << surf << " is not recognized!" << endl;
  }

  auto bcTag = -1;

  // Determine which boundary condition should be applied
  for (auto nn = iStart; nn < iEnd; nn++) {
    // Determine which boundary given i, j, k coordinates apply to
    if ((i >= this->GetIMin(nn) && i <= this->GetIMax(nn) &&
         j >= this->GetJMin(nn) && j <= this->GetJMax(nn) &&
         k >= this->GetKMin(nn) && k <= this->GetKMax(nn))) {
      bcTag = this->GetTag(nn);
      break;
    }
  }

  return bcTag;
}

// Member function to fill one "row" of the vectors with data that has been
// read in from the input file. This function is called from
// input::ReadInput(). It is necessary so that the private data can be
// altered from another class's member function.
void boundaryConditions::AssignFromInput(const int &surfCounter,
                                         const vector<string> &tokens) {
  // surfCounter -- index at which to place data
  // tokens -- vector of strings read from input file

  boundarySurface bcSurf_(tokens[0], stoi(tokens[1]),
                          stoi(tokens[2]), stoi(tokens[3]),
                          stoi(tokens[4]), stoi(tokens[5]),
                          stoi(tokens[6]), stoi(tokens[7]));
  surfs_[surfCounter] = bcSurf_;
}

/* Member function to determine of what sides of a boundary condition surface
   border another boundarySurface. This is necessary to complete the ghost cell
   swap properly. The function alters an array of 4 bools, which return true if
   the boundary surface borders an connection on the side that they represent.
   The order of the 4 bools is as follows [direction 1 start, direction 1 end,
   direction 2 start, direction 2 end].*/
void boundaryConditions::BordersSurface(const int &ii,
                                        bool (&border)[4]) const {
  // ii -- index of surface to test for border matches
  // border -- array of bools to show if boundarySurface is bordered
  // by another boundarySurface on any of its 4 sides

  // Get surface to test for borders
  const auto surf = this->GetSurface(ii);

  // Check that given boundarySurface is connection
  if (!surf.IsConnection()) {
    cerr << "ERROR: Error in boundaryConditions::BordersSurface(). "
         << "Given index does not point to an connection boundarySurface!"
         << endl;
    cerr << surf << endl;
    exit(EXIT_FAILURE);
  }

  // Initialize array of bools to false (does not border surface)
  border[0] = false;
  border[1] = false;
  border[2] = false;
  border[3] = false;

  // Loop over all surfaces in boundary conditions
  for (auto jj = 0; jj < this->NumSurfaces(); jj++) {
    const auto possibleBorder = this->GetSurface(jj);
    // If possible border is of same surface type, test for border match
    if (possibleBorder.SurfaceType() == surf.SurfaceType()) {
      // borders on direction 1 start side
      if (surf.Min1() == possibleBorder.Max1()) {
        border[0] = true;
      }
      // borders on direction 1 end side
      if (surf.Max1() == possibleBorder.Min1()) {
        border[1] = true;
      }
      // borders on direction 2 start side
      if (surf.Min2() == possibleBorder.Max2()) {
        border[2] = true;
      }
      // borders on direction 2 end side
      if (surf.Max2() == possibleBorder.Min2()) {
        border[3] = true;
      }
    }
  }
}

// Operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const connection &bc) {
  // os -- ostream to print to
  // bc -- connection to print

  os << "Is Interblock: " << bc.IsInterblock() << endl;
  os << "Ranks: " << bc.RankFirst() << ", " << bc.RankSecond() << endl;
  os << "Blocks: " << bc.BlockFirst() << ", " << bc.BlockSecond() << endl;

  os << "Local Blocks: " << bc.LocalBlockFirst() << ", "
     << bc.LocalBlockSecond() << endl;

  os << "Boundaries: " << bc.BoundaryFirst() << ", " << bc.BoundarySecond()
     << endl;

  os << "Direction 1 Starts: " << bc.Dir1StartFirst() << ", "
     << bc.Dir1StartSecond() << endl;

  os << "Direction 1 Ends: " << bc.Dir1EndFirst() << ", " << bc.Dir1EndSecond()
     << endl;

  os << "Direction 2 Starts: " << bc.Dir2StartFirst() << ", "
     << bc.Dir2StartSecond() << endl;

  os << "Direction 2 Ends: " << bc.Dir2EndFirst() << ", " << bc.Dir2EndSecond()
     << endl;

  os << "Direction 3 Constant Surface: " << bc.ConstSurfaceFirst() << ", "
     << bc.ConstSurfaceSecond() << endl;

  os << "Direction 1 Start Borders Surface: "
     << bc.Dir1StartInterBorderFirst() << ", "
     << bc.Dir1StartInterBorderSecond() << endl;

  os << "Direction 1 End Borders Surface: " << bc.Dir1EndInterBorderFirst()
     << ", " << bc.Dir1EndInterBorderSecond() << endl;

  os << "Direction 2 Start Borders Surface: "
     << bc.Dir2StartInterBorderFirst() << ", "
     << bc.Dir2StartInterBorderSecond() << endl;

  os << "Direction 2 End Borders Surface: " << bc.Dir2EndInterBorderFirst()
     << ", " << bc.Dir2EndInterBorderSecond() << endl;

  os << "Orientation: " << bc.Orientation() << endl;

  return os;
}

// Constructor to take in two patches and fill a connection.
// The orientation is left at the default value 0.
connection::connection(const patch &p1, const patch &p2) {
  // p1 -- patch 1
  // p2 -- patch 2

  // Fill connection
  rank_[0] = p1.Rank();
  rank_[1] = p2.Rank();

  block_[0] = p1.Block();
  block_[1] = p2.Block();

  localBlock_[0] = p1.LocalBlock();
  localBlock_[1] = p2.LocalBlock();

  boundary_[0] = p1.Boundary();
  boundary_[1] = p2.Boundary();

  d1Start_[0] = p1.Dir1Start();
  d1Start_[1] = p2.Dir1Start();

  d1End_[0] = p1.Dir1End();
  d1End_[1] = p2.Dir1End();

  d2Start_[0] = p1.Dir2Start();
  d2Start_[1] = p2.Dir2Start();

  d2End_[0] = p1.Dir2End();
  d2End_[1] = p2.Dir2End();

  constSurf_[0] = p1.ConstSurface();
  constSurf_[1] = p2.ConstSurface();

  patchBorder_[0] = p1.Dir1StartInterBorder();
  patchBorder_[1] = p1.Dir1EndInterBorder();
  patchBorder_[2] = p1.Dir2StartInterBorder();
  patchBorder_[3] = p1.Dir2EndInterBorder();
  patchBorder_[4] = p2.Dir1StartInterBorder();
  patchBorder_[5] = p2.Dir1EndInterBorder();
  patchBorder_[6] = p2.Dir2StartInterBorder();
  patchBorder_[7] = p2.Dir2EndInterBorder();

  orientation_ = 0;  // default value (real values 1-8)
  isInterblock_ = (p1.BCType() == "interblock" && p2.BCType() == "interblock")
      ? true : false;
}

// Function to swap the order of an connection so the 2nd entry
// in the pair will be the first, and vice versa
void connection::SwapOrder() {
  swap(rank_[0], rank_[1]);
  swap(block_[0], block_[1]);
  swap(localBlock_[0], localBlock_[1]);
  swap(boundary_[0], boundary_[1]);
  swap(d1Start_[0], d1Start_[1]);
  swap(d1End_[0], d1End_[1]);
  swap(d2Start_[0], d2Start_[1]);
  swap(d2End_[0], d2End_[1]);
  swap(constSurf_[0], constSurf_[1]);

  swap(patchBorder_[0], patchBorder_[4]);
  swap(patchBorder_[1], patchBorder_[5]);
  swap(patchBorder_[2], patchBorder_[6]);
  swap(patchBorder_[3], patchBorder_[7]);

  // If orientation is 4 or 5, needs to be swapped because direction
  // 1/2 are swapped and only one direction is reversed
  if (orientation_ == 4) {
    orientation_ = 5;
  } else if (orientation_ == 5) {
    orientation_ = 4;
  }
}

range connection::Dir1RangeFirst() const {
  return {d1Start_[0], d1End_[0]};
}

range connection::Dir1RangeSecond() const {
  return {d1Start_[1], d1End_[1]};
}

range connection::Dir2RangeFirst() const {
  return {d2Start_[0], d2End_[0]};
}

range connection::Dir2RangeSecond() const {
  return {d2Start_[1], d2End_[1]};
}

/* Function to go through the boundary conditions and pair the connection
   BCs together and determine their orientation.*/
vector<connection> GetConnectionBCs(const vector<boundaryConditions> &bc,
                                    const vector<plot3dBlock> &grid,
                                    const decomposition &decomp,
                                    const input &inp) {
  // bc -- vector of boundaryConditions for all blocks
  // grid -- vector of plot3Dblocks for entire computational mesh
  // decomp -- decomposition of grid onto processors
  // inp -- input variables

  // Isolate only the connection BCs and their associated data
  // from all of the BCs
  // Outer vector for each connection BC, inner vector for
  // information about connection
  vector<boundarySurface> isolatedConnections;

  // Block number of bc, rank of block, local position on processor
  vector<array<int, 3>> numRankPos;
  vector<int> surfaceNums;  // surface number of connection

  // loop over all blocks
  for (auto ii = 0U; ii < bc.size(); ii++) {
    // Loop over number of surfaces in block
    for (auto jj = 0; jj < bc[ii].NumSurfaces(); jj++) {
      // If boundary condition is connection, store data
      if (bc[ii].IsConnection(jj)) {
        // block number of bc, rank, local position boundarySurface of bc
        const array<int, 3> temp = {static_cast<int>(ii), decomp.Rank(ii),
                                    decomp.LocalPosition(ii)};
        numRankPos.push_back(temp);
        isolatedConnections.push_back(bc[ii].GetSurface(jj));
        surfaceNums.push_back(jj);
      }
    }
  }

  // ----------------------------------------------------------------------
  // Intialize vector of connections to return
  // Size is halved because each connection pairs with another
  vector<connection> connections(isolatedConnections.size() / 2);

  // Loop over isolated connections
  // ii counts by two because after a pair is found, that data is swapped
  // to ii+1. This allows the next search to avoid the matched pair
  for (auto ii = 0U; ii < isolatedConnections.size(); ii += 2) {
    // Loop over possible matches
    for (auto jj = ii + 1U; jj < isolatedConnections.size(); jj++) {
      // Blocks and boundary surfaces between connections match
      // blocks between connection BCs match
      // or both are periodic
      if ((isolatedConnections[ii].PartnerBlock() == numRankPos[jj][0] &&
           isolatedConnections[ii].PartnerSurface() ==
           isolatedConnections[jj].SurfaceType()) ||
          (isolatedConnections[ii].BCType() == "periodic" &&
           isolatedConnections[jj].BCType() == "periodic")) {
        // Determine if surface borders any other surfaces
        bool border[4] = {false, false, false, false};
        bc[numRankPos[ii][0]].BordersSurface(surfaceNums[ii], border);

        // Get current patch
        patch cPatch(isolatedConnections[ii], grid[numRankPos[ii][0]],
                     numRankPos[ii][0], border, numRankPos[ii][1],
                     numRankPos[ii][2]);
        if (cPatch.BCType() == "periodic") {
          const auto &bcData = inp.BCData(isolatedConnections[ii].Tag());
          if (bcData->StartTag() == isolatedConnections[ii].Tag()) {
            // need to transform data if surface is startTag
            cPatch.Transform(bcData);
          }
        }

        // Determine if surface borders any other surfaces
        bc[numRankPos[jj][0]].BordersSurface(surfaceNums[jj], border);

        // Get new patch (possible match)
        patch nPatch(isolatedConnections[jj], grid[numRankPos[jj][0]],
                     numRankPos[jj][0], border, numRankPos[jj][1],
                     numRankPos[jj][2]);
        if (nPatch.BCType() == "periodic") {
          const auto &bcData = inp.BCData(isolatedConnections[jj].Tag());
          if (bcData->StartTag() == isolatedConnections[jj].Tag()) {
            // need to transform data if surface is startTag
            nPatch.Transform(bcData);
          }
        }

        // Test for match
        connection match(cPatch, nPatch);
        if (match.TestPatchMatch(cPatch, nPatch)) {  // match found
          connections[ii / 2] = match;               // store connection pair

          // Swap matched connection BC to top portion of vector so
          // it is not searched again
          swap(isolatedConnections[jj], isolatedConnections[ii + 1]);
          swap(numRankPos[jj], numRankPos[ii + 1]);
          swap(surfaceNums[jj], surfaceNums[ii + 1]);
          break;  // exit innermost loop and search for next connection match
        }
      }
    }
  }

  return connections;
}

/* Function to take in two patches and return if they are matched. If there is a
   match it uses the patches to modify the given connection to contain the
   information on this match.

   Each patch is on a constant i, j, or k surface and is a 4 sided rectangle.
   The match is tested for by determining if the vertexes on one patch, match
   up with the vertexes on the other patch. Only 3 vertexes need to match, so
   only 3 are tested for. If there is a match the patches can be oriented in
   8 different ways with respect to each other. The orientation is stored in the
   connection that is modified.

                         Patch 1                       Patch 2 Description
                     __________________           __________________
                    |C1             C12|         |C1             C12|
                    |                  |         |                  |
Orientation 1:      |                  |         |                  |  Same orientation
                  D1|O               C2|       D1|O               C2|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | -->D2
                     __________________           __________________
                    |C1             C12|         |C2             C12|
                    |                  |         |                  |
Orientation 2:      |                  |         |                  |  D1/D2 swapped
                  D1|O               C2|       D2|O               C1|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | -->D1
                     __________________           __________________
                    |C1             C12|         |O               C2|
                    |                  |         |                  |
Orientation 3:      |                  |         |                  |  D1 reversed
                  D1|O               C2|      -D1|C1             C12|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | -->D2
                     __________________           __________________
                    |C1             C12|         |C12             C2|
                    |                  |         |                  |
Orientation 4:      |                  |         |                  |  D1/D2 swapped, D1 reversed
                  D1|O               C2|       D2|C1               O|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | D1<--
                     __________________           __________________
                    |C1             C12|         |O               C1|
                    |                  |         |                  |
Orientation 5:      |                  |         |                  |  D1/D2 swapped, D2 reversed
                  D1|O               C2|      -D2|C2             C12|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | -->D1
                     __________________           __________________
                    |C1             C12|         |C12             C1|
                    |                  |         |                  |
Orientation 6:      |                  |         |                  |  D2 reversed
                  D1|O               C2|       D1|C2               O|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | D2<--
                     __________________           __________________
                    |C1             C12|         |C1               O|
                    |                  |         |                  |
Orientation 7:      |                  |         |                  |  D1/D2 swapped, D1/D2 reversed
                  D1|O               C2|      -D2|C12             C2|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | D1<--
                     __________________           __________________
                    |C1             C12|         |C2               O|
                    |                  |         |                  |
Orientation 8:      |                  |         |                  |  D1/D2 reversed
                  D1|O               C2|      -D1|C12             C1|
                   ^|__________________|        ^|__________________|
                   | -->D2                      | D2<--

The above diagrams show how patch 2 would have to be moved to match up with
patch 1.
D1 and D2 are the local patch directions. They are cyclic, so on a constant
i-patch,
D1 is j, and D2 is k. On a constant j-patch, D1 is k, and D2 is i, etc. O is the
origin which is always at the minimim of D1 and D2 on the patch. C1 is the corner
where D1 is at a max, and D2 is zero. C2 is the corner where D2 is at a max, and
D1 is zero. C12 is the corner where both D1 and D2 are at a max.*/
bool connection::TestPatchMatch(const patch &p1, const patch &p2) {
  // p1 -- first patch
  // p2 -- second patch

  auto match = false;  // initialize match to false

  // test if bc types are the same
  if (p1.BCType() != p2.BCType()) {
    return match;
  }

  // Determine if there is a potential match by comparing origins
  if (p1.Origin().CompareWithTol(p2.Origin())) {  // origins match ------------
    // If origin matches origin, corner 1 can only be at corner 1 or 2
    if (p1.Corner1().CompareWithTol(p2.Corner1())) {  // corner 1s match
      // If all 3 corners match, same orientation
      if (p1.Corner2().CompareWithTol(p2.Corner2())) {  // corner 2s match
        orientation_ = 1;
        match = true;
      } else {  // no match
        return match;
      }
    } else if (p1.Corner1().CompareWithTol(p2.Corner2())) {  // match 1/2
      // If origins match and 1 matches 2, 2 must match 1
      if (p1.Corner2().CompareWithTol(p2.Corner1())) {  // corner 2 matches corner 1
        orientation_ = 2;
        match = true;
      } else {  // no match
        return match;
      }
    } else {  // no match
      return match;
    }
  } else if (p1.Origin().CompareWithTol(p2.Corner1())) {  // origin match corner 1
    // If origin matches corner1, corner 1 can only be at corner 12 or origin
    if (p1.Corner1().CompareWithTol(p2.Origin())) {
      // Corner 2 must match 12 for match
      if (p1.Corner2().CompareWithTol(p2.Corner12())) {
        orientation_ = 3;
        match = true;
      } else {  // no match
        return match;
      }
    } else if (p1.Corner1().CompareWithTol(p2.Corner12())) {
      // Corner 2 must match origin for match
      if (p1.Corner2().CompareWithTol(p2.Origin())) {
        orientation_ = 4;
        match = true;
      } else {  // no match
        return match;
      }
    } else {  // no match
      return match;
    }
  } else if (p1.Origin().CompareWithTol(p2.Corner2())) {  // origin match corner 2
    // If origin matches corner2, corner 1 can only be at corner 12 or origin
    if (p1.Corner1().CompareWithTol(p2.Origin())) {
      // Corner 2 must match 12 for match
      if (p1.Corner2().CompareWithTol(p2.Corner12())) {
        orientation_ = 5;
        match = true;
      } else {  // no match
        return match;
      }
    } else if (p1.Corner1().CompareWithTol(p2.Corner12())) {
      // Corner 2 must match origin_ for match
      if (p1.Corner2().CompareWithTol(p2.Origin())) {
        orientation_ = 6;
        match = true;
      } else {  // no match
        return match;
      }
    } else {  // no match
      return match;
    }
  } else if (p1.Origin().CompareWithTol(p2.Corner12())) {  // origin match 12
    // If origin matches corner 12, corner 1 can only be at corner 1 or corner 2
    if (p1.Corner1().CompareWithTol(p2.Corner1())) {
      // Corner 2 must match 2 for match
      if (p1.Corner2().CompareWithTol(p2.Corner2())) {
        orientation_ = 7;
        match = true;
      } else {  // no match
        return match;
      }
    } else if (p1.Corner1().CompareWithTol(p2.Corner2())) {
      // Corner 2 must match corner 1 for match
      if (p1.Corner2().CompareWithTol(p2.Corner2())) {
        orientation_ = 8;
        match = true;
      } else {  // no match
        return match;
      }
    } else {  // no match
      return match;
    }
  } else {  // no match
    return match;
  }

  return match;
}

/* Member function to adjust the connection for use with a geomSlice */
void connection::AdjustForSlice(const bool &blkFirst, const int &numG) {
  // blkFirst -- boolean that is true if block to insert into is first
  // numG -- number of ghost cells in block

  if (!blkFirst) {      // block to insert into is second, swap order
    this->SwapOrder();  // have block be first entry, slice second
  }

  // If at an upper surface, start block at upper boundary
  // if at lower surface, start block at last ghost layer
  const auto blkStart = (this->BoundaryFirst() % 2 == 0)
      ? this->ConstSurfaceFirst() : -numG;

  constSurf_[1] = 0;  // slice always starts at 0
  constSurf_[0] = blkStart;
  // Adjust direction 1 start and end for ghost cells
  d1End_[1] = this->Dir1EndSecond() - this->Dir1StartSecond() + 2 * numG;
  d1End_[0] = this->Dir1EndFirst() + numG;
  d1Start_[1] = 0;  // slice always starts at 0
  d1Start_[0] = this->Dir1StartFirst() - numG;
  // Adjust direction 2 start and end for ghost cells
  d2End_[1] = this->Dir2EndSecond() - this->Dir2StartSecond() + 2 * numG;
  d2End_[0] = this->Dir2EndFirst() + numG;
  d2Start_[1] = 0;  // slice always starts at 0
  d2Start_[0] = this->Dir2StartFirst() - numG;
}

// Member function to get the addresses of an connection to create
// an MPI_Datatype
void connection::GetAddressesMPI(MPI_Aint (&disp)[12]) const {
  // Get addresses of each field
  MPI_Get_address(&rank_[0], &disp[0]);
  MPI_Get_address(&block_[0], &disp[1]);
  MPI_Get_address(&localBlock_[0], &disp[2]);
  MPI_Get_address(&boundary_[0], &disp[3]);
  MPI_Get_address(&d1Start_[0], &disp[4]);
  MPI_Get_address(&d1End_[0], &disp[5]);
  MPI_Get_address(&d2Start_[0], &disp[6]);
  MPI_Get_address(&d2End_[0], &disp[7]);
  MPI_Get_address(&constSurf_[0], &disp[8]);
  MPI_Get_address(&patchBorder_[0], &disp[9]);
  MPI_Get_address(&orientation_, &disp[10]);
  MPI_Get_address(&isInterblock_, &disp[11]);
}

// Function to return which direction (i,j,k) is direction 1 in the
// first partner in the connection
string connection::Direction1First() const {
  string dir = "";
  if (this->BoundaryFirst() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "j";
  } else if (this->BoundaryFirst() <= 4) {  // dir 3 - j, dir 1 - k, dir 2 - i
    dir = "k";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "i";
  }

  return dir;
}

// Function to return which direction (i,j,k) is direction 1 in the second
// partner in the connection
string connection::Direction1Second() const {
  string dir = "";
  if (this->BoundarySecond() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "j";
  } else if (this->BoundarySecond() <= 4) {  // dir 3-j, dir 1-k, dir 2-i
    dir = "k";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "i";
  }

  return dir;
}

// Function to return which direction (i,j,k) is direction 2 in the first
// partner in the connection
string connection::Direction2First() const {
  string dir = "";
  if (this->BoundaryFirst() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "k";
  } else if (this->BoundaryFirst() <= 4) {  // dir 3 - j, dir 1 - k, dir 2 - i
    dir = "i";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "j";
  }

  return dir;
}

// Function to return which direction (i,j,k) is direction 2 in the second
// partner in the connection
string connection::Direction2Second() const {
  string dir = "";
  if (this->BoundarySecond() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "k";
  } else if (this->BoundarySecond() <= 4) {  // dir 3-j, dir 1-k, dir 2-i
    dir = "i";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "j";
  }

  return dir;
}

// function to return which direction (i,j,k) is direction 3 in the first
// partner in the connection
string connection::Direction3First() const {
  string dir = "";
  if (this->BoundaryFirst() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "i";
  } else if (this->BoundaryFirst() <= 4) {  // dir 3 - j, dir 1 - k, dir 2 - i
    dir = "j";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "k";
  }

  return dir;
}

// Function to return which direction (i,j,k) is direction 3 in the second
// partner in the connection
string connection::Direction3Second() const {
  string dir = "";
  if (this->BoundarySecond() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "i";
  } else if (this->BoundarySecond() <= 4) {  // dir 3-j, dir 1-k, dir 2-i
    dir = "j";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "k";
  }

  return dir;
}

/* Member function to update the connection border matching portion of the first
   partner in an connection. An index between 1-4 representing the 4 sides of
   the
   connection surface are passed to the function, and it ensures that the
   connectionBorder variable at that position reads true. This is used during
   the
   first PutGeomSlice ghost cell exchange to locate a "t" intersection of blocks
   and ensure that they are treated properly. */
void connection::UpdateBorderFirst(const int &a) {
  // a -- position in patchBorder_ to update (0-3) -
  // (dir 1 start, dir1 end, dir2 start, dir2 end)

  if (a >= 0 && a <= 3) {
    patchBorder_[a] = true;
  } else {
    cerr << "ERROR: Error in connection::UpdateBorderFirst(). "
         << "Position to update is out of range. Choose between 0-3. "
         << "Position input was " << a << endl;
    exit(EXIT_FAILURE);
  }
}

/* Member function to update the connection border matching portion of the
   second
   partner in an connection. An index between 1-4 representing the 4 sides of
   the
   connection surface are passed to the function, and it ensures that the
   connectionBorder variable at that position reads true. This is used during
   the first
   PutGeomSlice ghost cell exchange to locate a "t" intersection of blocks and
   ensure
   that they are treated properly. */
void connection::UpdateBorderSecond(const int &a) {
  // a -- position in patchBorder_ to update (0-3) -
  // (dir 1 start, dir1 end, dir2 start, dir2 end)

  if (a >= 0 && a <= 3) {
    patchBorder_[a + 4] = true;
  } else {
    cerr << "ERROR: Error in connection::UpdateBorderSecond(). "
         << "Position to update is out of range. Choose between 0-3. "
         << "Position input was " << a << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to get the indices of a slice for purposes of swapping
// for the first block in an connection
void connection::FirstSliceIndices(int &is1, int &ie1, int &js1, int &je1,
                                   int &ks1, int &ke1,
                                   const int &numGhosts1) const {
  // is1 -- starting i index for first slice
  // ie1 -- ending i index for first slice
  // js1 -- starting j index for first slice
  // je1 -- ending j index for first slice
  // ks1 -- starting k index for first slice
  // ke1 -- ending k index for first slice
  // numGhosts1 -- number of ghost cells in first slice

  // if at lower boundary no need to adjust for ghost cells as constant surface
  // is already at the interior cells when acounting for ghost cells
  // if at the upper boundary adjust the constant surface by the number of ghost
  // cells to get to the nth interior cell
  const auto upLowFac = (this->BoundaryFirst() % 2 == 0) ? -numGhosts1 : 0;

  if (this->BoundaryFirst() == 1 ||
      this->BoundaryFirst() == 2) {  // direction 3 is i
    // extend min/maxes to cover ghost cells
    is1 = this->ConstSurfaceFirst() + upLowFac;
    ie1 = is1 + numGhosts1;

    // direction 1 is j
    js1 = this->Dir1StartFirst() - numGhosts1;
    je1 = this->Dir1EndFirst() + numGhosts1;

    // direction 2 is k
    ks1 = this->Dir2StartFirst() - numGhosts1;
    ke1 = this->Dir2EndFirst() + numGhosts1;

  } else if (this->BoundaryFirst() == 3 ||
             this->BoundaryFirst() == 4) {  // direction 3 is j
    // extend min/maxes to cover ghost cells
    js1 = this->ConstSurfaceFirst() + upLowFac;
    je1 = js1 + numGhosts1;

    // direction 1 is k
    ks1 = this->Dir1StartFirst() - numGhosts1;
    ke1 = this->Dir1EndFirst() + numGhosts1;

    // direction 2 is i
    is1 = this->Dir2StartFirst() - numGhosts1;
    ie1 = this->Dir2EndFirst() + numGhosts1;

  } else if (this->BoundaryFirst() == 5 ||
             this->BoundaryFirst() == 6) {  // direction 3 is k
    // extend min/maxes to cover ghost cells
    ks1 = this->ConstSurfaceFirst() + upLowFac;
    ke1 = ks1 + numGhosts1;

    // direction 1 is i
    is1 = this->Dir1StartFirst() - numGhosts1;
    ie1 = this->Dir1EndFirst() + numGhosts1;

    // direction 2 is j
    js1 = this->Dir2StartFirst() - numGhosts1;
    je1 = this->Dir2EndFirst() + numGhosts1;

  } else {
    cerr << "ERROR: Error in connection::FirstSliceIndices(). Surface boundary "
         << this->BoundaryFirst() << " is not recognized!" << endl;
    cerr << "connection is: " << (*this) << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to get the indices of a slice for purposes of swapping
// for the first block in an connection
void connection::SecondSliceIndices(int &is2, int &ie2, int &js2, int &je2,
                                    int &ks2, int &ke2,
                                    const int &numGhosts2) const {
  // is2 -- starting i index for second slice
  // ie2 -- ending i index for second slice
  // js2 -- starting j index for second slice
  // je2 -- ending j index for second slice
  // ks2 -- starting k index for second slice
  // ke2 -- ending k index for second slice
  // numGhosts2 -- number of ghost cells in second slice

  // if at lower boundary no need to adjust for ghost cells as constant surface
  // is already at the interior cells when acounting for ghost cells
  // if at the upper boundary adjust the constant surface by the number of ghost
  // cells to get to the nth interior cell
  const auto upLowFac = (this->BoundarySecond() % 2 == 0) ? -numGhosts2 : 0;

  if (this->BoundarySecond() == 1 ||
      this->BoundarySecond() == 2) {  // direction 3 is i
    // extend min/maxes to cover ghost cells
    is2 = this->ConstSurfaceSecond() + upLowFac;
    ie2 = is2 + numGhosts2;

    // direction 1 is j
    js2 = this->Dir1StartSecond() - numGhosts2;
    je2 = this->Dir1EndSecond() + numGhosts2;

    // direction 2 is k
    ks2 = this->Dir2StartSecond() - numGhosts2;
    ke2 = this->Dir2EndSecond() + numGhosts2;

  } else if (this->BoundarySecond() == 3 ||
             this->BoundarySecond() == 4) {  // direction 3 is j
    // extend min/maxes to cover ghost cells
    js2 = this->ConstSurfaceSecond() + upLowFac;
    je2 = js2 + numGhosts2;

    // direction 1 is k
    ks2 = this->Dir1StartSecond() - numGhosts2;
    ke2 = this->Dir1EndSecond() + numGhosts2;

    // direction 2 is i
    is2 = this->Dir2StartSecond() - numGhosts2;
    ie2 = this->Dir2EndSecond() + numGhosts2;

  } else if (this->BoundarySecond() == 5 ||
             this->BoundarySecond() == 6) {  // direction 3 is k
    // extend min/maxes to cover ghost cells
    ks2 = this->ConstSurfaceSecond() + upLowFac;
    ke2 = ks2 + numGhosts2;

    // direction 1 is i
    is2 = this->Dir1StartSecond() - numGhosts2;
    ie2 = this->Dir1EndSecond() + numGhosts2;

    // direction 2 is j
    js2 = this->Dir2StartSecond() - numGhosts2;
    je2 = this->Dir2EndSecond() + numGhosts2;

  } else {
    cerr << "ERROR: Error in connection::SecondSliceIndices(). " <<
        "Surface boundary " << this->BoundarySecond() <<
        " is not recognized!" << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to determine the number of faces with a viscous wall BC
int boundaryConditions::NumViscousFaces() const {
  auto nFaces = 0;
  for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
    if (this->GetBCTypes(ii) == "viscousWall") {
      nFaces += surfs_[ii].NumFaces();
    }
  }
  return nFaces;
}

int boundaryConditions::BlockDimI() const {
  auto dim = 0;
  for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
    dim = std::max(dim, this->GetIMax(ii));
  }
  return dim;
}

int boundaryConditions::BlockDimJ() const {
  auto dim = 0;
  for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
    dim = std::max(dim, this->GetJMax(ii));
  }
  return dim;
}

int boundaryConditions::BlockDimK() const {
  auto dim = 0;
  for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
    dim = std::max(dim, this->GetKMax(ii));
  }
  return dim;
}

/* Member function to split boundary conditions along a given direction at a
   given index. The calling instance retains the lower portion of the split,
   and the returned instance is the upper portion. */
boundaryConditions boundaryConditions::Split(const string &dir, const int &ind,
                                             const int &numBlk,
                                             const int &newBlkNum,
                                             vector<boundarySurface> &aSurf) {
  // dir -- direction to split it (i, j, k)
  // ind -- index of cell to split at
  //        (this index is the last cell that remains in the lower split)
  // numBlk -- block number that (*this) is assocatied with
  // newBlkNum -- block number for upper split
  // aSurf -- vector of any interblocks that are split,
  //          because their partners will need to be altered for the split as
  //          well

  const auto indNG = ind;

  auto bound1 = (*this);  // lower boundaryConditions
  auto bound2 = (*this);  // upper boundaryConditions

  // Initialize vector of surfaces that will be altered by this split
  vector<boundarySurface> alteredSurf;

  // If a surface in entirely on one side of the split, it should be removed
  // from the other side - these variables keep track of the surfaces to remove
  vector<int> del1, del2;
  // Reserved for maximum number of deletions
  del1.reserve(this->NumSurfaces());
  // Reserved for maximum number of deletions
  del2.reserve(this->NumSurfaces());

  if (dir == "i") {  // split along i-plane
    // Initialize deletion numbers to 0
    auto del1I = 0;
    auto del1J = 0;
    auto del1K = 0;
    auto del2I = 0;
    auto del2J = 0;
    auto del2K = 0;

    auto numInterL = 0;
    auto numInterU = 0;

    // Loop over all surfaces
    for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
      if (ii < this->NumSurfI()) {                      // i-surface
        if (this->GetSurface(ii).SurfaceType() == 1) {  // lower i surface
          // No change to lower bc at lower i surface

          // At lower i surface, upper bc is now interface
          // lower surface matches with upper surface
          const auto tag = 2000 + numBlk;
          bound2.surfs_[ii].bcType_ = "interblock";          // bcType
          bound2.surfs_[ii].data_[0] = this->GetIMin(ii);  // imin
          bound2.surfs_[ii].data_[1] = this->GetIMax(ii);  // imax
          bound2.surfs_[ii].data_[6] = tag;                  // tag

          // There should only be one surface between the split blocks
          bound2.surfs_[ii].data_[2] = 0;                    // jmin
          bound2.surfs_[ii].data_[3] = this->BlockDimJ();  // jmax
          bound2.surfs_[ii].data_[4] = 0;                    // kmin
          bound2.surfs_[ii].data_[5] = this->BlockDimK();  // kmax

          if (numInterL > 0) {
            del2.push_back(ii);
            del2I++;
          }
          numInterL++;
        } else {  // upper surface
          // At upper i surface, lower bc is now interface
          // upper surface matches with lower surface
          const auto tag = 1000 + newBlkNum;
          bound1.surfs_[ii].bcType_ = "interblock";  // bcType
          bound1.surfs_[ii].data_[0] = indNG;        // imin
          bound1.surfs_[ii].data_[1] = indNG;        // imax
          bound1.surfs_[ii].data_[6] = tag;          // tag

          // There should only be one surface between the split blocks
          bound1.surfs_[ii].data_[2] = 0;                    // jmin
          bound1.surfs_[ii].data_[3] = this->BlockDimJ();  // jmax
          bound1.surfs_[ii].data_[4] = 0;                    // kmin
          bound1.surfs_[ii].data_[5] = this->BlockDimK();  // kmax

          if (numInterU > 0) {
            del1.push_back(ii);
            del1I++;
          }

          numInterU++;

          // At upper i surface, upper bc is same as original,
          // but indices are adjusted for new block_ size
          bound2.surfs_[ii].data_[0] = this->GetIMax(ii) - indNG;  // imin
          bound2.surfs_[ii].data_[1] = this->GetIMax(ii) - indNG;  // imax

          // At upper i surface, if bc is interblock, store boundarySurface
          // because partner block BC will need to be updated
          if (this->GetBCTypes(ii) == "interblock") {
            alteredSurf.push_back(this->GetSurface(ii));
          }
        }
      } else {  // j-surface or k-surface
        // At j/k surface, if bc is interblock, store boundarySurface
        // because partner block BC will need to be updated
        if (this->GetBCTypes(ii) == "interblock") {
          alteredSurf.push_back(this->GetSurface(ii));
        }

        // This surface is only present in the upper split
        if (this->GetIMin(ii) >= indNG) {
          bound2.surfs_[ii].data_[0] = this->GetIMin(ii) - indNG;  // imin
          bound2.surfs_[ii].data_[1] = this->GetIMax(ii) - indNG;  // imax
          // Can delete it from lower split
          del1.push_back(ii);
          if (ii >= this->NumSurfI() &&
              ii < this->NumSurfI() + this->NumSurfJ()) {  // j-surface
            del1J++;
          } else {  // k-surface
            del1K++;
          }
        } else if (this->GetIMax(ii) > indNG) {  // surf straddles the split
          bound1.surfs_[ii].data_[1] = indNG;       // imax
          bound2.surfs_[ii].data_[0] = 0;           // imin
          bound2.surfs_[ii].data_[1] = this->GetIMax(ii) - indNG;  // imax
        } else {  // surf only present in the lower split - delete from upper
          del2.push_back(ii);
          if (ii >= this->NumSurfI() &&
              ii < this->NumSurfI() + this->NumSurfJ()) {  // j-surface
            del2J++;
          } else {  // k-surface
            del2K++;
          }
        }
      }
    }

    // Delete unnecessary boundaries and change number of surfaces in i,j,k
    // to appropriate number

    // need to delete from bottom of vector so indices are preserved
    // need to cast to int because value must be negative for termination
    for (auto ii = static_cast<int>(del1.size()) - 1; ii >= 0; --ii) {
      bound1.surfs_.erase(bound1.surfs_.begin() + del1[ii]);
    }
    bound1.numSurfI_ -= del1I;
    bound1.numSurfJ_ -= del1J;
    bound1.numSurfK_ -= del1K;

    for (auto ii = static_cast<int>(del2.size()) - 1; ii >= 0; --ii) {
      bound2.surfs_.erase(bound2.surfs_.begin() + del2[ii]);
    }
    bound2.numSurfI_ -= del2I;
    bound2.numSurfJ_ -= del2J;
    bound2.numSurfK_ -= del2K;

  } else if (dir == "j") {  // split along j-plane
    // Initialize deletion numbers to 0
    auto del1I = 0;
    auto del1J = 0;
    auto del1K = 0;
    auto del2I = 0;
    auto del2J = 0;
    auto del2K = 0;

    auto numInterL = 0;
    auto numInterU = 0;

    // Loop over all surfaces
    for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
      if (ii >= this->NumSurfI() &&
          ii < this->NumSurfI() + this->NumSurfJ()) {  // j-surface
        if (this->GetSurface(ii).SurfaceType() == 3) {   // lower j surface
          // No change to lower bc at lower j surface

          // At lower j surface, upper bc is now interface
          // lower surface matches with upper surface
          const auto tag = 4000 + numBlk;
          bound2.surfs_[ii].bcType_ = "interblock";          // bctype
          bound2.surfs_[ii].data_[2] = this->GetJMin(ii);  // jmin
          bound2.surfs_[ii].data_[3] = this->GetJMax(ii);  // jmax
          bound2.surfs_[ii].data_[6] = tag;                  // tag

          // There should only be one surface between the split blocks
          bound2.surfs_[ii].data_[0] = 0;                    // imin
          bound2.surfs_[ii].data_[1] = this->BlockDimI();  // imax
          bound2.surfs_[ii].data_[4] = 0;                    // kmin
          bound2.surfs_[ii].data_[5] = this->BlockDimK();  // kmax

          if (numInterL > 0) {
            del2.push_back(ii);
            del2J++;
          }
          numInterL++;
        } else {
          // At upper j surface, lower bc is now interface
          // upper surface matches with lower
          const auto tag = 3000 + newBlkNum;
          bound1.surfs_[ii].bcType_ = "interblock";  // bctype
          bound1.surfs_[ii].data_[2] = indNG;        // jmin
          bound1.surfs_[ii].data_[3] = indNG;        // jmax
          bound1.surfs_[ii].data_[6] = tag;          // tag

          // There should only be one surface between the split blocks
          bound1.surfs_[ii].data_[0] = 0;                    // imin
          bound1.surfs_[ii].data_[1] = this->BlockDimI();  // imax
          bound1.surfs_[ii].data_[4] = 0;                    // kmin
          bound1.surfs_[ii].data_[5] = this->BlockDimK();  // kmax

          if (numInterU > 0) {
            del1.push_back(ii);
            del1J++;
          }
          numInterU++;

          // At upper j surface, upper bc is same as original, but indices
          // are adjusted for new block size
          bound2.surfs_[ii].data_[2] = this->GetJMax(ii) - indNG;  // jmin
          bound2.surfs_[ii].data_[3] = this->GetJMax(ii) - indNG;  // jmax

          // At upper j surface, if bc is interblock, store boundarySurface
          // because partner block BC will need to be updated
          if (this->GetBCTypes(ii) == "interblock") {
            alteredSurf.push_back(this->GetSurface(ii));
          }
        }
      } else {  // i-surface or k-surface
        // At i/k surface, if bc is interblock, store boundarySurface
        // because partner block BC will need to be updated
        if (this->GetBCTypes(ii) == "interblock") {
          alteredSurf.push_back(this->GetSurface(ii));
        }

        if (this->GetJMin(ii) >= indNG) {  // surf is only in the upper split
          bound2.surfs_[ii].data_[2] = this->GetJMin(ii) - indNG;  // jmin
          bound2.surfs_[ii].data_[3] = this->GetJMax(ii) - indNG;  // jmax
          // can delete from lower split
          del1.push_back(ii);
          if (ii < this->NumSurfI()) {  // i-surface
            del1I++;
          } else {  // k-surface
            del1K++;
          }
        } else if (this->GetJMax(ii) > indNG) {  // surf straddles the split
          bound1.surfs_[ii].data_[3] = indNG;       // jmax
          bound2.surfs_[ii].data_[2] = 0;           // jmin
          bound2.surfs_[ii].data_[3] = this->GetJMax(ii) - indNG;  // jmax
        } else {  // surface is only in the lower split - can delete from upper
          del2.push_back(ii);
          if (ii < this->NumSurfI()) {  // i-surface
            del2I++;
          } else {  // k-surface
            del2K++;
          }
        }
      }
    }

    // Delete unnecessary boundaries - and set number of surfaces (i, j, k)
    // to appropriate number

    // need to delete from bottom of vector so indices are preserved
    // need to cast to int because value must be negative for termination
    for (auto ii = static_cast<int>(del1.size()) - 1; ii >= 0; --ii) {
      bound1.surfs_.erase(bound1.surfs_.begin() + del1[ii]);
    }
    bound1.numSurfI_ -= del1I;
    bound1.numSurfJ_ -= del1J;
    bound1.numSurfK_ -= del1K;

    for (auto ii = static_cast<int>(del2.size()) - 1; ii >= 0; --ii) {
      bound2.surfs_.erase(bound2.surfs_.begin() + del2[ii]);
    }
    bound2.numSurfI_ -= del2I;
    bound2.numSurfJ_ -= del2J;
    bound2.numSurfK_ -= del2K;

  } else if (dir == "k") {  // split along k-plane
    // Initialize deletion numbers to 0
    auto del1I = 0;
    auto del1J = 0;
    auto del1K = 0;
    auto del2I = 0;
    auto del2J = 0;
    auto del2K = 0;

    auto numInterL = 0;
    auto numInterU = 0;

    // Loop over all surfaces
    for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
      if (ii >= this->NumSurfI() + this->NumSurfJ()) {  // k-surface
        if (this->GetSurface(ii).SurfaceType() == 5) {    // lower k surface
          // No change to lower bc at lower k surface

          // At lower k surface, upper bc is now interface
          // lower surface matches with upper
          const auto tag = 6000 + numBlk;
          bound2.surfs_[ii].bcType_ = "interblock";          // bctype
          bound2.surfs_[ii].data_[4] = this->GetKMin(ii);  // kmin
          bound2.surfs_[ii].data_[5] = this->GetKMax(ii);  // kmax
          bound2.surfs_[ii].data_[6] = tag;                  // tag

          // There should only be one surface between the split blocks
          bound2.surfs_[ii].data_[0] = 0;                    // imin
          bound2.surfs_[ii].data_[1] = this->BlockDimI();  // imax
          bound2.surfs_[ii].data_[2] = 0;                    // jmin
          bound2.surfs_[ii].data_[3] = this->BlockDimJ();  // jmax

          if (numInterL > 0) {
            del2.push_back(ii);
            del2K++;
          }
          numInterL++;
        } else {
          // At upper k surface, lower bc is now interface
          // upper surface matches with lower
          const auto tag = 5000 + newBlkNum;
          bound1.surfs_[ii].bcType_ = "interblock";  // bctype
          bound1.surfs_[ii].data_[4] = indNG;        // kmin
          bound1.surfs_[ii].data_[5] = indNG;        // kmax
          bound1.surfs_[ii].data_[6] = tag;          // tag

          // There should only be one surface between the split blocks
          bound1.surfs_[ii].data_[0] = 0;                    // imin
          bound1.surfs_[ii].data_[1] = this->BlockDimI();  // imax
          bound1.surfs_[ii].data_[2] = 0;                    // jmin
          bound1.surfs_[ii].data_[3] = this->BlockDimJ();  // jmax

          if (numInterU > 0) {
            del1.push_back(ii);
            del1K++;
          }
          numInterU++;

          // At upper k surface, upper bc is same as original,
          // but indices are adjusted for new block size
          bound2.surfs_[ii].data_[4] = this->GetKMax(ii) - indNG;  // kmin
          bound2.surfs_[ii].data_[5] = this->GetKMax(ii) - indNG;  // kmax

          // At upper k surface, if bc is interblock, store boundarySurface
          // because partner block BC will need to be updated
          if (this->GetBCTypes(ii) == "interblock") {
            alteredSurf.push_back(this->GetSurface(ii));
          }
        }
      } else {  // i-surface or j-surface
        // At i/j surface, if bc is interblock, store boundarySurface because
        // partner block BC will need to be updated
        if (this->GetBCTypes(ii) == "interblock") {
          alteredSurf.push_back(this->GetSurface(ii));
        }

        if (this->GetKMin(ii) >= indNG) {  // surface is only in upper split
          bound2.surfs_[ii].data_[4] = this->GetKMin(ii) - indNG;  // kmin
          bound2.surfs_[ii].data_[5] = this->GetKMax(ii) - indNG;  // kmax
          // Can delete from lower surface
          del1.push_back(ii);
          if (ii < this->NumSurfI()) {  // i-surface
            del1I++;
          } else {  // j-surface
            del1J++;
          }
        } else if (this->GetKMax(ii) > indNG) {  // surf straddles the split
          bound1.surfs_[ii].data_[5] = indNG;       // kmax
          bound2.surfs_[ii].data_[4] = 0;           // kmin
          bound2.surfs_[ii].data_[5] = this->GetKMax(ii) - indNG;  // kmax
        } else {  // surface is only in the lower split - can delete from upper
          del2.push_back(ii);
          if (ii < this->NumSurfI()) {  // i-surface
            del2I++;
          } else {  // j-surface
            del2J++;
          }
        }
      }
    }

    // Delete unnecessary boundaries and set surface numbers (i,j,k) to
    // appropriate value

    // need to delete from bottom of vector so indices are preserved
    // need to cast to int because value must be negative for termination
    for (auto ii = static_cast<int>(del1.size()) - 1; ii >= 0; --ii) {
      bound1.surfs_.erase(bound1.surfs_.begin() + del1[ii]);
    }
    bound1.numSurfI_ -= del1I;
    bound1.numSurfJ_ -= del1J;
    bound1.numSurfK_ -= del1K;

    for (auto ii = static_cast<int>(del2.size()) - 1; ii >= 0; --ii) {
      bound2.surfs_.erase(bound2.surfs_.begin() + del2[ii]);
    }
    bound2.numSurfI_ -= del2I;
    bound2.numSurfJ_ -= del2J;
    bound2.numSurfK_ -= del2K;

  } else {
    cerr << "ERROR: Error in boundaryCondition::Split(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(EXIT_FAILURE);
  }

  (*this) = bound1;     // assign lower split to (*this)
  aSurf = alteredSurf;  // assign vector of altered surfaces
  return bound2;
}

/* Member function to split the surfaces of a boundaryCondtions accordingly when
   one of its connection partners has been altered. The connection partner may
   have been split, its block number updated, or both. In order to correctly
   match up the dependents of the connection must be updated for the split.*/
void boundaryConditions::DependentSplit(const boundarySurface &surf,
                                        const plot3dBlock &part,
                                        const plot3dBlock &self,
                                        const int &sblk, const string &dir,
                                        const int &ind, const int &lblk,
                                        const int &ublk) {
  // surf -- boundarySurface of partner block
  // part -- plot3dBlock that surf is assigned to
  // self -- plot3dBlock that (*this) is assigned to
  // sblk -- block number of self
  // dir -- direction that partner split was in
  // ind -- index of split
  // lblk -- lower block number in partner split
  // ublk -- upper block number in partner split

  // dummy value used because connection is only used to test for match
  bool border[4] = {false, false, false, false};

  const patch partner(surf, part, lblk, border);  // create patch for partner

  // loop over all surfaces
  for (auto ii = 0; ii < this->NumSurfaces(); ii++) {
    // create patch for candidate match
    const patch candidate(this->GetSurface(ii), self, sblk, border);

    connection match(candidate, partner);
    if (match.TestPatchMatch(candidate, partner)) {  // match found
      auto lowSurf = this->GetSurface(ii);  // get match surface

      // determine direction and index to split surface
      string candDir = "";
      auto candInd = 0;
      if (match.Orientation() == 1) {  // same orientation
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;  // candInd doesn't matter for dir 3 b/c block cannot
                          // be split, only partner block updated
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 2) {  // D1/D2 swapped
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 3) {  // D1 reversed
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = surf.Max1() - ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 4) {  // D1/D2 swapped, D1 reversed
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = surf.Max1() - ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 5) {  // D1/D2 swapped, D2 reversed
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = surf.Max2() - ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 6) {  // D2 reversed
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = surf.Max2() - ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else if (match.Orientation() == 7) {  // D1/D2 swapped and reversed
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = surf.Max1() - ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = surf.Max2() - ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }

      } else {  // D1/D2 reversed (orientation 8)
        if (surf.Direction1() == dir) {  // split was in direction 1 of partner,
                                         // needs to be direction 1 of candidate
          candDir = lowSurf.Direction1();
          candInd = surf.Max1() - ind - surf.Min1() + lowSurf.Min1();
        } else if (surf.Direction2() == dir) {  // split was in direction 2 of
                                                // partner, needs to be
                                                // direction 2 of candidate
          candDir = lowSurf.Direction2();
          candInd = surf.Max2() - ind - surf.Min2() + lowSurf.Min2();
        } else if (surf.Direction3() == dir) {  // split was in direction 3 of
                                                // partner, needs to be
                                                // direction 3 of candidate
          candDir = lowSurf.Direction3();
          candInd = ind;
        } else {
          cerr << "ERROR: Error in boundaryConditions::DependentSplit(). "
                  "Direction " << dir << " is not recognized." << endl;
          cerr << "Please choose i, j, or k." << endl;
          exit(EXIT_FAILURE);
        }
      }

      // split matched surface
      auto split = false;  // flag to tell if surface was split (or just if
                           // block number updated)
      const auto upSurf = lowSurf.Split(candDir, candInd, lblk, ublk, split,
                                        match.Orientation());

      // assign boundarySurface back into boundaryConditions, if surface wasn't
      // split partner block was updated
      surfs_[ii] = lowSurf;

      // if surface was split, insert it into the vector of boundarySurfaces and
      // adjust the surface numbers
      if (split) {  // boundary_ surface was split, insert new surface into
                    // boundary_ conditions
        surfs_.insert(surfs_.begin() + ii, upSurf);
        if (upSurf.SurfaceType() <= 2) {  // i-surface
          numSurfI_++;
        } else if (upSurf.SurfaceType() <= 4) {  // j-surface
          numSurfJ_++;
        } else {
          numSurfK_++;
        }
      }

      break;
    }
  }
}

/* Member function to join 2 boundaryConditions. It assumes that the calling
instance is the "lower" boundary condition and the input instance
is the "upper" boundary condition.
*/
void boundaryConditions::Join(const boundaryConditions &bc, const string &dir,
                              vector<boundarySurface> &aSurf) {
  // bc -- boundary_ condition (upper) to join
  // dir -- direction of join plane
  // aSurf -- vector of connections whose partners will need to be altered by
  // the join

  vector<boundarySurface> alteredSurf;  // initialize vector of boundary
                                        // surfaces whose partners will need to
                                        // be altered by the join

  if (dir == "i") {  // split along i-plane
    // total number of i surfaces in joined block will be equal to the lower i
    // surfaces from lower bc plus the upper i surfaces from the upper bc
    auto numI = 0;
    for (auto ii = 0; ii < this->NumSurfI(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 1) {  // lower i surface
        numI++;
      }
    }
    for (auto ii = 0; ii < bc.NumSurfI(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 2) {  // upper i surface
        numI++;
      }
    }

    // total number of j surfaces in joined block will be equal to all j
    // surfaces from lower bc plus the j surfaces from the upper bc that only
    // reside in the upper bc
    const auto numJ = this->NumSurfJ() + bc.NumSurfJ();

    // total number of k surfaces in joined block_ will be equal to all k
    // surfaces from lower bc plus the k surfaces from the upper bc that only
    // reside in the upper bc
    const auto numK = this->NumSurfK() + bc.NumSurfK();

    // initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    auto cc = 0;  // boundary condition counter

    // insert all i lower surfaces from lower bc
    auto lowDimI = 0;
    for (auto ii = 0; ii < this->NumSurfI(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 1) {  // lower i surface
        newBC.surfs_[cc] = surfs_[ii];
        cc++;
      } else {
        lowDimI = this->GetIMax(ii);  // maximum i-dimension for lower bc
      }
    }
    // insert all i upper surfaces from upper bc
    for (auto ii = 0; ii < bc.NumSurfI(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 2) {  // upper i surface
        // at upper i surface, if bc is connection, store boundarySurface
        // because partner block BC will need to be updated
        if (bc.IsConnection(ii)) {
          alteredSurf.push_back(bc.GetSurface(ii));
        }

        // adjust i coordinates for join
        const boundarySurface bcSurf(bc.GetBCTypes(ii),
                                     bc.GetIMin(ii) + lowDimI,
                                     bc.GetIMax(ii) + lowDimI,
                                     bc.GetJMin(ii),
                                     bc.GetJMax(ii), bc.GetKMin(ii),
                                     bc.GetKMax(ii),
                                     bc.GetTag(ii));

        newBC.surfs_[cc] = bcSurf;
        cc++;
      }
    }

    // insert all j surfaces from lower and upper bcs
    for (auto ii = this->NumSurfI();
         ii < this->NumSurfI() + this->NumSurfJ(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++) {
      // at j surface for upper block, if bc is connection, store
      // boundarySurface because partner block BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust i coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii),
                                   bc.GetIMin(ii) + lowDimI,
                                   bc.GetIMax(ii) + lowDimI,
                                   bc.GetJMin(ii),
                                   bc.GetJMax(ii), bc.GetKMin(ii),
                                   bc.GetKMax(ii),
                                   bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    // insert all k surfaces from lower and upper bcs
    for (auto ii = this->NumSurfI() + this->NumSurfJ();
         ii < this->NumSurfaces(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++) {
      // at k surface for upper block, if bc is connection, store
      // boundarySurface because partner block BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust i coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii),
                                   bc.GetIMin(ii) + lowDimI,
                                   bc.GetIMax(ii) + lowDimI,
                                   bc.GetJMin(ii),
                                   bc.GetJMax(ii), bc.GetKMin(ii),
                                   bc.GetKMax(ii),
                                   bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;
  } else if (dir == "j") {  // split along j-plane
    // total number of j surfaces in joined block_ will be equal to the lower j
    // surfaces from lower bc plus the upper j surfaces from the upper bc
    auto numJ = 0;
    for (auto ii = this->NumSurfI();
         ii < this->NumSurfI() + this->NumSurfJ(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 3) {  // lower j surface
        numJ++;
      }
    }
    for (auto ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 4) {  // upper j surface
        numJ++;
      }
    }

    // total number of i surfaces in joined block_ will be equal to all i
    // surfaces from lower bc plus the i surfaces from the upper bc that only
    // reside in the upper bc
    const auto numI = this->NumSurfI() + bc.NumSurfI();

    // total number of k surfaces in joined block_ will be equal to all k
    // surfaces from lower bc plus the k surfaces from the upper bc that only
    // reside in the upper bc
    const auto numK = this->NumSurfK() + bc.NumSurfK();

    // initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    auto cc = numI;  // boundary_ condition counter

    // insert all j lower surfaces from lower bc
    auto lowDimJ = 0;
    for (auto ii = this->NumSurfI();
         ii < this->NumSurfI() + this->NumSurfJ(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 3) {  // lower j surface
        newBC.surfs_[cc] = surfs_[ii];
        cc++;
      } else {
        lowDimJ = this->GetJMax(ii);  // maximum j-dimension for lower bc
      }
    }
    // insert all j upper surfaces from upper bc
    for (auto ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 4) {  // upper j surface
        // at j upper surface for upper block, if bc is connection, store
        // boundarySurface because partner block BC will need to be updated
        if (bc.IsConnection(ii)) {
          alteredSurf.push_back(bc.GetSurface(ii));
        }

        // adjust j coordinates for join
        const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                     bc.GetIMax(ii),
                                     bc.GetJMin(ii) + lowDimJ,
                                     bc.GetJMax(ii) + lowDimJ,
                                     bc.GetKMin(ii),
                                     bc.GetKMax(ii), bc.GetTag(ii));

        newBC.surfs_[cc] = bcSurf;
        cc++;
      }
    }

    cc = 0;
    // insert all i surfaces from lower and upper bcs
    for (auto ii = 0; ii < this->NumSurfI(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = 0; ii < bc.NumSurfI(); ii++) {
      // at i surface for upper block, if bc is connection, store
      // boundarySurface because partner block BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust j coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                   bc.GetIMax(ii),
                                   bc.GetJMin(ii) + lowDimJ,
                                   bc.GetJMax(ii) + lowDimJ, bc.GetKMin(ii),
                                   bc.GetKMax(ii), bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    cc = numI + numJ;
    // insert all k surfaces from lower and upper bcs
    for (auto ii = this->NumSurfI() + this->NumSurfJ();
         ii < this->NumSurfaces(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++) {
      // at k surface for upper block, if bc is connection, store
      // boundarySurface because partner block BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust j coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                    bc.GetIMax(ii),
                                    bc.GetJMin(ii) + lowDimJ,
                                    bc.GetJMax(ii) + lowDimJ,
                                    bc.GetKMin(ii),
                                    bc.GetKMax(ii), bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;

  } else if (dir == "k") {  // split along k-plane
    // total number of k surfaces in joined block_ will be equal to the lower k
    // surfaces from lower bc plus the upper k surfaces from the upper bc
    auto numK = 0;
    for (auto ii = this->NumSurfI() + this->NumSurfJ();
         ii < this->NumSurfaces(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 5) {  // lower k surface
        numK++;
      }
    }
    for (auto ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 6) {  // upper k surface
        numK++;
      }
    }

    // total number of i surfaces in joined block_ will be equal to all i
    // surfaces from lower bc plus the i surfaces from the upper bc that only
    // reside in the upper bc
    const auto numI = this->NumSurfI() + bc.NumSurfI();

    // total number of j surfaces in joined block_ will be equal to all j
    // surfaces from lower bc plus the j surfaces from the upper bc that only
    // reside in the upper bc
    const auto numJ = this->NumSurfJ() + bc.NumSurfJ();

    // initialze bc with new surface numbers
    boundaryConditions newBC(numI, numJ, numK);
    auto cc = numI + numJ;  // boundary_ condition counter

    // insert all k lower surfaces from lower bc
    auto lowDimK = 0;
    for (auto ii = this->NumSurfI() + this->NumSurfJ();
         ii < this->NumSurfaces(); ii++) {
      if (this->GetSurface(ii).SurfaceType() == 5) {  // lower k surface
        newBC.surfs_[cc] = surfs_[ii];
        cc++;
      } else {
        lowDimK = this->GetKMax(ii);  // maximum k-dimension for lower bc
      }
    }
    // insert all k upper surfaces from upper bc
    for (auto ii = bc.NumSurfI() + bc.NumSurfJ(); ii < bc.NumSurfaces(); ii++) {
      if (bc.GetSurface(ii).SurfaceType() == 6) {  // upper k surface
        // at upper k surface for upper block, if bc is connection, store
        // boundarySurface because partner block BC will need to be updated
        if (bc.IsConnection(ii)) {
          alteredSurf.push_back(bc.GetSurface(ii));
        }

        // adjust k coordinates for join
        const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                     bc.GetIMax(ii), bc.GetJMin(ii),
                                     bc.GetJMax(ii),
                                     bc.GetKMin(ii) + lowDimK,
                                     bc.GetKMax(ii) + lowDimK,
                                     bc.GetTag(ii));

        newBC.surfs_[cc] = bcSurf;
        cc++;
      }
    }

    cc = 0;
    // insert all i surfaces from lower and upper bcs
    for (auto ii = 0; ii < this->NumSurfI(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = 0; ii < bc.NumSurfI(); ii++) {
      // at i surface for upper block_, if bc is connection, store
      // boundarySurface because partner block_ BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust k coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                   bc.GetIMax(ii),
                                   bc.GetJMin(ii), bc.GetJMax(ii),
                                   bc.GetKMin(ii) + lowDimK,
                                   bc.GetKMax(ii) + lowDimK, bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    cc = numI;
    // insert all j surfaces from lower and upper bcs
    for (auto ii = this->NumSurfI();
         ii < this->NumSurfI() + this->NumSurfJ(); ii++) {
      newBC.surfs_[cc] = surfs_[ii];
      cc++;
    }
    for (auto ii = bc.NumSurfI(); ii < bc.NumSurfI() + bc.NumSurfJ(); ii++) {
      // at j surface for upper block, if bc is connection, store
      // boundarySurface because partner block BC will need to be updated
      if (bc.IsConnection(ii)) {
        alteredSurf.push_back(bc.GetSurface(ii));
      }

      // adjust k coordinates for join
      const boundarySurface bcSurf(bc.GetBCTypes(ii), bc.GetIMin(ii),
                                   bc.GetIMax(ii),
                                   bc.GetJMin(ii), bc.GetJMax(ii),
                                   bc.GetKMin(ii) + lowDimK,
                                   bc.GetKMax(ii) + lowDimK, bc.GetTag(ii));

      newBC.surfs_[cc] = bcSurf;
      cc++;
    }

    (*this) = newBC;

  } else {
    cerr << "ERROR: Error in procBlock::Join(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(EXIT_FAILURE);
  }

  aSurf = alteredSurf;
}

// constructor when passed no arguements
patch::patch() {
  // initialize all variables to zero
  vector3d<double> zero(0.0, 0.0, 0.0);
  origin_ = zero;
  corner1_ = zero;
  corner2_ = zero;
  corner12_ = zero;

  block_ = 0;
  boundary_ = 0;
  d1Start_ = 0;
  d1End_ = 0;
  d2Start_ = 0;
  d2End_ = 0;
  constSurf_ = 0;
  rank_ = 0;
  localBlock_ = 0;
  patchBorder_[0] = false;
  patchBorder_[1] = false;
  patchBorder_[2] = false;
  patchBorder_[3] = false;
  bcName_ = "undefined";
}

// constructor with arguements passed
patch::patch(const int &bound, const int &b, const int &d1s, const int &d1e,
             const int &d2s, const int &d2e, const int &d3s, const int &d3e,
             const plot3dBlock &blk, const int &r, const int &l,
             const bool (&border)[4], const string &name) {
  // bound -- boundary number which patch is on (1-6)
  // b -- parent block number
  // d1s -- direction 1 starting index
  // d1e -- direction 1 ending index
  // d2s -- direction 2 starting index
  // d2e -- direction 2 ending index
  // d3s -- direction 3 surface index (constant surface that patch is on)
  // blk -- plot3dBlock that patch is on
  // r -- rank of block
  // l -- local position of block
  // border -- flags indicating if patch borders an connection bc on sides 1/2
  // name -- bc name

  boundary_ = bound;
  block_ = b;
  rank_ = r;
  localBlock_ = l;
  patchBorder_[0] = border[0];
  patchBorder_[1] = border[1];
  patchBorder_[2] = border[2];
  patchBorder_[3] = border[3];
  bcName_ = name;

  if (bound == 1 || bound == 2) {  // patch on i-surface - dir1 = j, dir2 = k
    d1Start_ = d2s;
    d1End_ = d2e;
    d2Start_ = d3s;
    d2End_ = d3e;
    constSurf_ = d1s;

    // get corner points
    // origin_ at jmin, kmin
    origin_ = vector3d<double>(blk.Coords(constSurf_, d1Start_, d2Start_));

    // corner1_ at jmax, kmin
    corner1_ = vector3d<double>(blk.Coords(constSurf_, d1End_, d2Start_));

    // corner2_ at jmin, kmax
    corner2_ = vector3d<double>(blk.Coords(constSurf_, d1Start_, d2End_));

    // corner12_ at jmax, kmax
    corner12_ = vector3d<double>(blk.Coords(constSurf_, d1End_, d2End_));

  } else if (bound == 3 ||
             bound == 4) {  // patch on j-surface - dir1 = k, dir2 = i
    d1Start_ = d3s;
    d1End_ = d3e;
    d2Start_ = d1s;
    d2End_ = d1e;
    constSurf_ = d2s;

    // get corner points
    // origin_ at kmin, imin
    origin_ = vector3d<double>(blk.Coords(d2Start_, constSurf_, d1Start_));

    // corner1_ at kmax, imin
    corner1_ = vector3d<double>(blk.Coords(d2Start_, constSurf_, d1End_));

    // corner2_ at kmin, imax
    corner2_ = vector3d<double>(blk.Coords(d2End_, constSurf_, d1Start_));

    // corner12_ at kmax, imax
    corner12_ = vector3d<double>(blk.Coords(d2End_, constSurf_, d1End_));

  } else if (bound == 5 ||
             bound == 6) {  // patch on k-surface - dir1 = i, dir2 = j
    d1Start_ = d1s;
    d1End_ = d1e;
    d2Start_ = d2s;
    d2End_ = d2e;
    constSurf_ = d3s;

    // get corner points
    // origin_ at imin, jmin
    origin_ = vector3d<double>(blk.Coords(d1Start_, d2Start_, constSurf_));

    // corner1_ at imax, jmin
    corner1_ = vector3d<double>(blk.Coords(d1End_, d2Start_, constSurf_));

    // corner2_ at imin, jmax
    corner2_ = vector3d<double>(blk.Coords(d1Start_, d2End_, constSurf_));

    // corner12_ at imax, jmax
    corner12_ = vector3d<double>(blk.Coords(d1End_, d2End_, constSurf_));

  } else {
    cerr << "ERROR: Error in patch::patch(). Boundary surface " << bound
         << " is not recognized!" << endl;
    cerr << "Choose an integer between 1-6." << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to transform patch based on given bc data
// this is used for periodic boundaries
void patch::Transform(const unique_ptr<inputState> &bcData) {
  if (bcData->IsTranslation()) {
    const auto translation = bcData->Translation();
    this->Translate(translation);
  } else if (bcData->IsRotation()) {
    const auto axis = bcData->Axis();
    const auto rotation = bcData->Rotation();
    this->Rotate(axis, rotation);
  } else {
    cerr << "ERROR. BC data for transformation is not translation or rotation!"
         << endl;
    exit(EXIT_FAILURE);
  }
}

// member function to transform patch based on given translation
// this is used for periodic boundaries
void patch::Translate(const vector3d<double> &translate) {
  origin_ += translate;
  corner1_ += translate;
  corner2_ += translate;
  corner12_ += translate;
}

// member function to transform patch based on given rotation
// this is used for periodic boundaries
void patch::Rotate(const vector3d<double> &axis, const double &rotation) {
  cerr << "Rotate transformation is not currently supported!" << endl;
  exit(EXIT_FAILURE);
}


// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const patch &p) {
  // os -- ostream to print to
  // p -- patch to print

  os << "BC Type: " << p.BCType() << endl;
  os << "Boundary: " << p.Boundary() << endl;
  os << "Block: " << p.Block() << endl;
  os << "Direction 1 Start: " << p.Dir1Start() << endl;
  os << "Direction 1 End: " << p.Dir1End() << endl;
  os << "Direction 2 Start: " << p.Dir2Start() << endl;
  os << "Direction 2 End: " << p.Dir2End() << endl;
  os << "Constant Surface: " << p.ConstSurface() << endl;
  os << "Rank: " << p.Rank() << endl;
  os << "Borders Connection: " << p.Dir1StartInterBorder() << ", "
     << p.Dir1EndInterBorder() << ", " << p.Dir2StartInterBorder() << ", "
     << p.Dir2EndInterBorder() << endl;
  os << "Origin: " << p.Origin() << endl;
  os << "Corner 1: " << p.Corner1() << endl;
  os << "Corner 2: " << p.Corner2() << endl;
  os << "Corner 1-2: " << p.Corner12() << endl;

  return os;
}

/*Member function to pack a boundaryConditions into a buffer so that in can be
 * sent with MPI.*/
void boundaryConditions::PackBC(char *(&sendBuffer), const int &sendBufSize,
                                int &position) const {
  // sendBuffer -- buffer to pack data into
  // sendBufSize -- size of buffer
  // position -- location within buffer

  // get string lengths for each boundary_ condition to be sent, so processors
  // unpacking know how much data to unpack for each string
  vector<int> strLength(this->NumSurfaces());
  for (auto jj = 0U; jj < strLength.size(); jj++) {
    // +1 for c_str end character
    strLength[jj] = this->GetBCTypes(jj).size() + 1;
  }

  // pack surface numbers
  MPI_Pack(&numSurfI_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numSurfJ_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  MPI_Pack(&numSurfK_, 1, MPI_INT, sendBuffer, sendBufSize, &position,
           MPI_COMM_WORLD);
  // pack non-string data from boundarySurface
  for (auto jj = 0; jj < this->NumSurfaces(); jj++) {
    MPI_Pack(&surfs_[jj].data_[0], 7, MPI_INT, sendBuffer, sendBufSize,
             &position, MPI_COMM_WORLD);
  }
  // pack string size
  MPI_Pack(&strLength[0], strLength.size(), MPI_INT, sendBuffer, sendBufSize,
           &position, MPI_COMM_WORLD);
  // pack boundary condtion names
  for (auto jj = 0; jj < this->NumSurfaces(); jj++) {
    MPI_Pack(surfs_[jj].bcType_.c_str(),
             surfs_[jj].bcType_.size() + 1, MPI_CHAR, sendBuffer,
             sendBufSize, &position,
             MPI_COMM_WORLD);  // +1 for c_str end character
  }
}

/*Member function to unpack data from a buffer into a boundaryConditions. Used
 * with MPI receive*/
void boundaryConditions::UnpackBC(char *(&recvBuffer), const int &recvBufSize,
                                  int &position) {
  // recvBuffer -- buffer to unpack data from
  // recvBufSize -- size of buffer
  // position -- location within buffer

  MPI_Unpack(recvBuffer, recvBufSize, &position, &numSurfI_, 1, MPI_INT,
             MPI_COMM_WORLD);  // unpack number of i-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &numSurfJ_, 1, MPI_INT,
             MPI_COMM_WORLD);  // unpack number of j-surfaces
  MPI_Unpack(recvBuffer, recvBufSize, &position, &numSurfK_, 1, MPI_INT,
             MPI_COMM_WORLD);  // unpack number of k-surfaces

  // resize boundaryCondition
  this->ResizeVecs(this->NumSurfaces());
  // allocate vector for string lengths
  vector<int> strLength(this->NumSurfaces());

  // unpack boundary condition surface data (non-string) into appropriate
  // vectors
  for (auto jj = 0; jj < this->NumSurfaces(); jj++) {
    MPI_Unpack(recvBuffer, recvBufSize, &position, &surfs_[jj].data_[0],
               7, MPI_INT, MPI_COMM_WORLD);  // unpack bc surfaces
  }
  MPI_Unpack(recvBuffer, recvBufSize, &position, &strLength[0],
             strLength.size(), MPI_INT, MPI_COMM_WORLD);  // unpack string sizes
  // unpack boundary condition names
  for (auto jj = 0U; jj < strLength.size(); jj++) {
    auto *nameBuf = new char[strLength[jj]];  // allocate buffer to store BC
                                              // name
    MPI_Unpack(recvBuffer, recvBufSize, &position, &nameBuf[0], strLength[jj],
               MPI_CHAR, MPI_COMM_WORLD);  // unpack bc types
    string bcName(nameBuf, strLength[jj] - 1);  // create string of bc name (-1
                                                // to exclude c_str end
                                                // character)
    surfs_[jj].bcType_ = bcName;
    delete[] nameBuf;  // deallocate bc name buffer
  }
}

// member function to retrn the suface type of a boundarySurface. The surface
// type is an integer from 1 to 6 representing the imin, imax, jmin, jmax, kmin,
// or kmax surfaces respectively.
int boundarySurface::SurfaceType() const {
  auto surf = 0;

  if (data_[0] == data_[1]) {  // i-surface
    if (data_[1] == 0) {  // lower surface
      surf = 1;
    } else {  // upper surface
      surf = 2;
    }
  } else if (data_[2] == data_[3]) {  // j-surface
    if (data_[3] == 0) {  // lower surface
      surf = 3;
    } else {  // upper surface
      surf = 4;
    }
  } else if (data_[4] == data_[5]) {  // k-surface
    if (data_[5] == 0) {  // lower surface
      surf = 5;
    } else {  // upper surface
      surf = 6;
    }
  } else {
    cerr << "ERROR: Error in boundarySurface::SurfaceType(). Surface is "
            "defined incorrectly, it is neither an i, j, or k surface." << endl;
    cerr << (*this) << endl;
    exit(EXIT_FAILURE);
  }

  return surf;
}

// member function to determine the partner block of an interblock
// boundarySurface. The partner block is where the surface gets its ghost cells
// from.
int boundarySurface::PartnerBlock() const {
  if (bcType_ != "interblock") {
    // partner block only defined for interblock
    return -1;
  } else {
    const auto subtract = this->PartnerSurface() * 1000;
    return this->Tag() - subtract;
  }
}

// member function to determine the partner surface of an interblock
// boundarySurface. The partner surface is where the surface gets its ghost
// cells from.
int boundarySurface::PartnerSurface() const {
  if (bcType_ != "interblock") {
    // partner surface only defined for interblock
    return -1;
  } else {
    auto surf = 0;

    if (this->Tag() < 2000) {
      surf = 1;  // i-lower surface
    } else if (this->Tag() < 3000) {
      surf = 2;  // i-upper surface
    } else if (this->Tag() < 4000) {
      surf = 3;  // j-lower surface
    } else if (this->Tag() < 5000) {
      surf = 4;  // j-upper surface
    } else if (this->Tag() < 6000) {
      surf = 5;  // k-lower surface
    } else if (this->Tag() < 7000) {
      surf = 6;  // k-upper surface
    } else {
      cerr << "ERROR: Error in boundarySurface::PartnerSurface(). Tag does not "
          "fit in range. Tag must be between 1000 and 6999." << endl;
      cerr << (*this) << endl;
      exit(EXIT_FAILURE);
    }

    return surf;
  }
}

// member function to return a string corresponding to which direction (i,j,k)
// is direction 1 of the boundary_ surface
string boundarySurface::Direction1() const {
  string dir = "";
  if (this->SurfaceType() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "j";
  } else if (this->SurfaceType() <= 4) {  // dir 3 is j, dir 1 is k, dir 2 is
                                            // i
    dir = "k";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "i";
  }

  return dir;
}

// member function to return a string corresponding to which direction (i,j,k)
// is direction 2 of the boundary_ surface
string boundarySurface::Direction2() const {
  string dir = "";
  if (this->SurfaceType() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "k";
  } else if (this->SurfaceType() <= 4) {  // dir 3 is j, dir 1 is k, dir 2 is
                                            // i
    dir = "i";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "j";
  }

  return dir;
}

// member function to return a string corresponding to which direction (i,j,k)
// is direction 3 of the boundary_ surface.
string boundarySurface::Direction3() const {
  string dir = "";
  if (this->SurfaceType() <= 2) {  // dir 3 is i, dir 1 is j, dir 2 is k
    dir = "i";
  } else if (this->SurfaceType() <= 4) {  // dir 3 is j, dir 1 is k, dir 2 is
                                            // i
    dir = "j";
  } else {  // dir 3 is k, dir 1 is i, dir 2 is j
    dir = "k";
  }

  return dir;
}

// member function to return the maximum dimension of the surface in direction
// 1.
int boundarySurface::Max1() const {
  auto m = 0;
  if (this->Direction1() == "i") {
    m = this->IMax();
  } else if (this->Direction1() == "j") {
    m = this->JMax();
  } else {
    m = this->KMax();
  }
  return m;
}

// member function to return the minimum dimension of the surface in direction
// 1.
int boundarySurface::Min1() const {
  auto m = 0;
  if (this->Direction1() == "i") {
    m = this->IMin();
  } else if (this->Direction1() == "j") {
    m = this->JMin();
  } else {
    m = this->KMin();
  }
  return m;
}

// member function to return the maximum dimension of the surface in direction
// 2.
int boundarySurface::Max2() const {
  auto m = 0;
  if (this->Direction2() == "i") {
    m = this->IMax();
  } else if (this->Direction2() == "j") {
    m = this->JMax();
  } else {
    m = this->KMax();
  }
  return m;
}

// member function to return the minimum dimension of the surface in direction
// 2.
int boundarySurface::Min2() const {
  auto m = 0;
  if (this->Direction2() == "i") {
    m = this->IMin();
  } else if (this->Direction2() == "j") {
    m = this->JMin();
  } else {
    m = this->KMin();
  }
  return m;
}

// member function to return the number of faces this boundary surface is
// applied to
int boundarySurface::NumFaces() const {
  auto nFaces = 0;
  if (this->SurfaceType() <= 2) {  // i-surface
    nFaces = (this->JMax() - this->JMin()) * (this->KMax() - this->KMin());
  } else if (this->SurfaceType() <= 4) {  // j-surface
    nFaces = (this->IMax() - this->IMin()) * (this->KMax() - this->KMin());
  } else {  // k-surface
    nFaces = (this->IMax() - this->IMin()) * (this->JMax() - this->JMin());
  }
  return nFaces;
}

range boundarySurface::RangeI() const {
  return (data_[0] == data_[1]) ? data_[0] : range{data_[0], data_[1]};
}

range boundarySurface::RangeJ() const {
  return (data_[2] == data_[3]) ? data_[2] : range{data_[2], data_[3]};
}

range boundarySurface::RangeK() const {
  return (data_[4] == data_[5]) ? data_[4] : range{data_[4], data_[5]};
}

range boundarySurface::RangeDir1() const {
  if (this->Direction1() == "i") {
    return this->RangeI();
  } else if (this->Direction1() == "j") {
    return this->RangeJ();
  } else {
    return this->RangeK();
  }
}

range boundarySurface::RangeDir2() const {
  if (this->Direction2() == "i") {
    return this->RangeI();
  } else if (this->Direction2() == "j") {
    return this->RangeJ();
  } else {
    return this->RangeK();
  }
}

range boundarySurface::RangeDir3() const {
  if (this->Direction3() == "i") {
    return this->RangeI();
  } else if (this->Direction3() == "j") {
    return this->RangeJ();
  } else {
    return this->RangeK();
  }
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const boundarySurface &bcSurf) {
  // os -- ostream to print to
  // bcSurf -- boundarySurface to print

  os << std::left << std::setw(24) << bcSurf.BCType() << std::right
     << std::setw(5) << bcSurf.IMin() << std::setw(5) << bcSurf.IMax()
     << std::setw(5) << bcSurf.JMin() << std::setw(5) << bcSurf.JMax()
     << std::setw(5) << bcSurf.KMin() << std::setw(5) << bcSurf.KMax()
     << std::setw(5) << bcSurf.Tag();

  return os;
}

// member funtion to updat the tag of an interblock to partner with the block
// input to the function
void boundarySurface::UpdateTagForSplitJoin(const int &nBlk) {
  if (bcType_ != "interblock") {
    cerr << "ERROR: Only tags associated with interblock boundaries need to be "
            "updated. Current boundary is " << bcType_ << endl;
    exit(EXIT_FAILURE);
  }

  data_[6] = this->PartnerSurface() * 1000 + nBlk;
}

// member function to split a boundary_ surface. The calling instance retains
// the lower portion of the split, and the returned instance retains the upper
// portion
boundarySurface boundarySurface::Split(const string &dir, const int &ind,
                                       const int &lBlk, const int &uBlk,
                                       bool &split, int orientation) {
  // dir -- direction to split the surface in
  // ind -- index at which to split the surface
  // lBlk -- lower block number of split
  // uBlk -- upper block number of split
  // split -- flag to determine whether block was split or just tag updated, if
  // no split, upper surface returned is meaningless
  // orientation -- if called from DependentSplit, orientation of partner split

  auto surf1 = (*this);  // lower surface
  auto surf2 = (*this);  // upper surface

  split = true;  // initialize split flag
  // flag to determine if the split direction is reversed - used to determine
  // which block should match lower/upper surfaces
  // surf1 and surf2 have same orientation, so if reversed for 1, reversed for 2
  const auto isReversed = surf1.SplitDirectionIsReversed(dir, orientation);

  if (dir == "i") {  // split along i-plane
    if (this->SurfaceType() == 1 || this->SurfaceType() == 2 || ind <= 0) {
      // cannot split an i-surface along i-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    } else {  // j or k surface
      if (this->IMin() >= ind) {
        // this surface is only present in the upper split
        surf1.data_[0] = this->IMin() - ind + 1 + this->IMin();  // imin
        surf1.data_[1] = this->IMax() - ind + 1 + this->IMin();  // imax
        isReversed ? surf1.UpdateTagForSplitJoin(lBlk)
                   : surf1.UpdateTagForSplitJoin(uBlk);
        split = false;
      } else if (this->IMax() > ind) {  // this surface straddles the split
        surf2.data_[0] = ind;  // imin
        isReversed ? surf2.UpdateTagForSplitJoin(lBlk)
                   : surf2.UpdateTagForSplitJoin(uBlk);

        surf1.data_[1] = ind;  // imax
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
      } else {  // this surface is only present in the lower split
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
        split = false;
      }
    }

  } else if (dir == "j") {  // split along j-plane
    if (this->SurfaceType() == 3 || this->SurfaceType() == 4 || ind <= 0) {
      // cannot split a j-surface along j-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    } else {  // i or k surface
      if (this->JMin() >= ind) {
        // this surface is only present in the upper split
        surf1.data_[2] = this->JMin() - ind + 1 + this->JMin();  // jmin
        surf1.data_[3] = this->JMax() - ind + 1 + this->JMin();  // jmax
        isReversed ? surf1.UpdateTagForSplitJoin(lBlk)
                   : surf1.UpdateTagForSplitJoin(uBlk);
        split = false;
      } else if (this->JMax() > ind) {  // this surface straddles the split
        surf2.data_[2] = ind;  // jmin
        isReversed ? surf2.UpdateTagForSplitJoin(lBlk)
                   : surf2.UpdateTagForSplitJoin(uBlk);

        surf1.data_[3] = ind;  // jmax
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
      } else {  // this surface is only present in the lower split
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
        split = false;
      }
    }

  } else if (dir == "k") {  // split along k-plane
    if (this->SurfaceType() == 5 || this->SurfaceType() == 6 || ind <= 0) {
      // cannot split a k-surface along k-plane, just update block
      surf1.UpdateTagForSplitJoin(uBlk);
      split = false;
    } else {  // i or j surface
      if (this->KMin() >= ind) {  // surface only present in the upper split
        surf1.data_[4] = this->KMin() - ind + 1 + this->KMin();  // kmin
        surf1.data_[5] = this->KMax() - ind + 1 + this->KMin();  // kmax
        isReversed ? surf1.UpdateTagForSplitJoin(lBlk)
                   : surf1.UpdateTagForSplitJoin(uBlk);
        split = false;
      } else if (this->KMax() > ind) {  // this surface straddles the split
        surf2.data_[4] = ind;  // kmin
        isReversed ? surf2.UpdateTagForSplitJoin(lBlk)
                   : surf2.UpdateTagForSplitJoin(uBlk);

        surf1.data_[5] = ind;  // kmax
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
      } else {  // this surface is only present in the lower split
        isReversed ? surf1.UpdateTagForSplitJoin(uBlk)
                   : surf1.UpdateTagForSplitJoin(lBlk);
        split = false;
      }
    }

  } else {
    cerr << "ERROR: Error in boundarySurface::Split(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(EXIT_FAILURE);
  }

  (*this) = surf1;  // return lower surface as (*this)
  return surf2;
}

/*member function to determine if the split direction is reversed. A "reversed"
split direction occurs when the two patches that make an connection have their
split direction
running in opposite directions.
          ______       ______
         |      |   | |      |
         |   1  |  \/ |   2  |
split->  |______|     |______|
         |      |     |      |
         |      |^    |      |
         |______||    |______|

As you can see, the index for the split would be different in block_ 1 vs block_
2 because the direction of the increasing index is reversed. This situation is a
"reversed" split
direction. This function finds cases such as these.
 */
bool boundarySurface::SplitDirectionIsReversed(const string &dir,
                                               const int &orientation) const {
  // dir -- direction of split
  // orientation -- orientation of one patch to its partner

  auto isReversed = false;

  // find out if split direction is 1, 2, or 3
  if (this->Direction1() ==
      dir) {  // split direction is direction 1 - reverse if dir 1 is reversed
              // (relative to partner, taking into account D1/D2 swap)
    isReversed = (orientation == 3 || orientation == 5 || orientation == 7 ||
                  orientation == 8)
                     ? true
                     : false;
  } else if (this->Direction2() == dir) {  // split direction is direction 2 -
                                             // reverse if dir 2 is reversed
                                             // (relative to partner, taking
                                             // into account D1/D2 swap)
    isReversed = (orientation == 4 || orientation == 6 || orientation == 7 ||
                  orientation == 8)
                     ? true
                     : false;
  } else if (this->Direction3() == dir) {
    // split direction is direction 3 - no need to reverse
    isReversed = false;
  } else {
    cerr << "ERROR: Error in boundarySurface::SplitDirectionIsReversed(). "
            "Direction " << dir << " does not match i, j, or k!" << endl;
    exit(EXIT_FAILURE);
  }

  return isReversed;
}

/* Function to return an array of location indicies for ghost cells at an
connection boundary. The array is formatted as shown below:

  array = [i j k]

The array will contain 3 entries corresponding to the i, j, and k locations of
either the first or second pair in the connection, depending on what is
specified in the 'first' variable. The indices returned will correspond to cell
locations and will take into account the orientation of the patches that
comprise the connection with relation to each other.
*/
array<int, 3> GetSwapLoc(const int &l1, const int &l2, const int &l3,
                         const int &numGhosts, const connection &inter,
                         const int &d3, const bool &first) {
  // l1 -- index of direction 1 within slice to insert
  // l2 -- index of direction 2 within slice to insert
  // l3 -- index of direction 3 within slice to insert
  // numGhosts -- number of layers of ghost cells
  // inter -- connection boundary condition
  // d3 -- length of normal direction of connection
  // first -- flag for first or second block in connection match

  // preallocate array to return
  array<int, 3> loc = {0, 0, 0};

  if (first) {  // working on first in pair ------------------------------
    // first patch in pair is calculated using orientation 1
    if (inter.Direction3First() == "i") {  // i-patch
      // get direction 1 length
      loc[1] = inter.Dir1StartFirst() + l1;  // direction 1 is j
      loc[2] = inter.Dir2StartFirst() + l2;  // direction 2 is k
      loc[0] = inter.IsLowerFirst() ? l3 - numGhosts : inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else if (inter.Direction3First() == "j") {  // j-patch
      // get direction 1 length
      loc[2] = inter.Dir1StartFirst() + l1;  // direction 1 is k
      loc[0] = inter.Dir2StartFirst() + l2;  // direction 2 is i
      loc[1] = inter.IsLowerFirst() ? l3 - numGhosts : inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else if (inter.Direction3First() == "k") {  // k-patch
      // get direction 1 length
      loc[0] = inter.Dir1StartFirst() + l1;  // direction 1 is i
      loc[1] = inter.Dir2StartFirst() + l2;  // direction 2 is j
      loc[2] = inter.IsLowerFirst() ? l3 - numGhosts : inter.ConstSurfaceFirst() +
               l3;  // add l3 to get to ghost cells (cell index instead of face)
    } else {
      cerr << "ERROR: Error in procBlock:GetSwapLoc(). Boundary direction "
           << inter.Direction3First() << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  //--------------------------------------------------------------------------
  // need to use orientation for second in pair
  } else {  // working on second in pair ---------------------------------------
    if (inter.Direction3Second() == "i") {  // i-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is j (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[2] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is k (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[1] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is j -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[1] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 2 is k -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[2] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      // add l3 to get to ghost cells
      // if lower/lower or upper/upper, need to reverse direction 3
      if (inter.IsLowerLowerOrUpperUpper()) {
        loc[0] = inter.IsLowerSecond() ? d3 - l3 - 1:
            inter.ConstSurfaceSecond() + d3 - l3 - 1;
      } else {
        loc[0] = inter.IsLowerSecond() ? l3 - numGhosts :
            inter.ConstSurfaceSecond() + l3;
      }

    //-------------------------------------------------------------------------
    } else if (inter.Direction3Second() == "j") {  // j-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is k (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[0] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is i (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[2] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is k -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[2] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 2 is i -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[0] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      // add l3 to get to ghost cells
      // if lower/lower or upper/upper, need to reverse direction 3
      if (inter.IsLowerLowerOrUpperUpper()) {
        loc[1] = inter.IsLowerSecond() ? d3 - l3 - 1:
            inter.ConstSurfaceSecond() + d3 - l3 - 1;
      } else {
        loc[1] = inter.IsLowerSecond() ? l3 - numGhosts :
            inter.ConstSurfaceSecond() + l3;
      }


    //------------------------------------------------------------------------
    } else if (inter.Direction3Second() == "k") {  // k-patch
      if (inter.Orientation() == 2 || inter.Orientation() == 4 ||
          inter.Orientation() == 5 ||
          inter.Orientation() == 7) {  // swap dir 1 and 2
        // direction 1 is i (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[1] = (inter.Orientation() == 5 || inter.Orientation() == 7)
                     ? inter.Dir2EndSecond() - 1 - l1
                     : inter.Dir2StartSecond() + l1;

        // direction 2 is j (swapped) -- if true direction reversed -- subtract
        // 1 from End to get to cell index
        loc[0] = (inter.Orientation() == 4 || inter.Orientation() == 7)
                     ? inter.Dir1EndSecond() - 1 - l2
                     : inter.Dir1StartSecond() + l2;
      } else {  // no direction swap
        // direction 1 is i -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[0] = (inter.Orientation() == 3 || inter.Orientation() == 8)
                     ? inter.Dir1EndSecond() - 1 - l1
                     : inter.Dir1StartSecond() + l1;

        // direction 2 is j -- if true direction reversed -- subtract 1 from End
        // to get to cell index
        loc[1] = (inter.Orientation() == 6 || inter.Orientation() == 8)
                     ? inter.Dir2EndSecond() - 1 - l2
                     : inter.Dir2StartSecond() + l2;
      }

      // calculate index for all ghost layers
      // add l3 to get to ghost cells
      // if lower/lower or upper/upper, need to reverse direction 3
      if (inter.IsLowerLowerOrUpperUpper()) {
        loc[2] = inter.IsLowerSecond() ? d3 - l3 - 1:
            inter.ConstSurfaceSecond() + d3 - l3 - 1;
      } else {
        loc[2] = inter.IsLowerSecond() ? l3 - numGhosts :
            inter.ConstSurfaceSecond() + l3;
      }


    //--------------------------------------------------------------------------
    } else {
      cerr << "ERROR: Error in procBlock.cpp:GetSwapLoc(). Boundary surface of "
           << inter.Direction3Second() << " is not recognized!" << endl;
      exit(EXIT_FAILURE);
    }
  }

  return loc;
}

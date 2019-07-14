/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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

#include <algorithm>  // max_element
#include <iterator>   // distance
#include <cstring>
#include <vector>  // vector
#include <string>  // string
#include <memory>  // make_unique
#include "parallel.hpp"
#include "output.hpp"
#include "vector3d.hpp"            // vector3d
#include "plot3d.hpp"              // plot3d
#include "primitive.hpp"            // primitive
#include "procBlock.hpp"           // procBlock
#include "boundaryConditions.hpp"  // connection
#include "resid.hpp"               // resid
#include "gridLevel.hpp"
#include "macros.hpp"

using std::max_element;
using std::min_element;
using std::distance;
using std::unique_ptr;
using std::make_unique;

/* Function to return processor list for manual decomposition. Manual
decomposition assumes that each block will reside on it's own processor.
The processor list tells how many procBlocks a processor will have.
*/
decomposition ManualDecomposition(vector<plot3dBlock> &grid,
                                  vector<boundaryConditions> &bcs,
                                  const int &numProc) {
  // grid -- vector of procBlocks (no need to split procBlocks or combine them
  // with manual decomposition)
  // bcs -- vector of boundary conditions for each block in grid
  // numProc -- number of processors in run

  if (static_cast<int>(grid.size()) != numProc) {
    cerr << "ERROR: Error in parallel.cpp:ManualDecomposition(). Manual "
            "decomposition assumes that the number of processors used is equal "
            "to the "
         << "number of blocks in the grid. This grid has " << grid.size()
         << " blocks and the simulation is using " << numProc << " processors."
         << endl;
    exit(EXIT_FAILURE);
  }

  cout << "--------------------------------------------------------------------"
          "------------" << endl;
  cout << "Using manual grid decomposition." << endl;

  decomposition decomp(grid.size(), numProc);
  for (auto ii = 1U; ii < grid.size(); ii++) {
    decomp.SendToProc(ii, ROOTP, ii);  // send block from ROOT to processor
  }

  decomp.PrintDiagnostics(grid);
  cout << endl;
  cout << "Ideal Load: " << decomp.IdealLoad(grid) << endl;
  cout << "Max Load: " << decomp.MaxLoad(grid) << endl;

  auto loaded = 0.0;
  auto ol = decomp.MostOverloadedProc(grid, loaded);
  cout << "Most overloaded processor is " << ol << "; overloaded by " << loaded
       << endl;
  auto ul = decomp.MostUnderloadedProc(grid, loaded);
  cout << "Most underloaded processor is " << ul << "; underloaded by "
       << loaded << endl;

  cout << "Ratio of most loaded processor to average processor is : "
       << decomp.MaxLoad(grid) / decomp.IdealLoad(grid) << endl;
  cout << "--------------------------------------------------------------------"
          "------------" << endl << endl;

  return decomp;
}

/* Function to return processor list for cubic decomposition.
The processor list tells how many procBlocks a processor will have.
*/
decomposition CubicDecomposition(vector<plot3dBlock> &grid,
                                 vector<boundaryConditions> &bcs,
                                 const int &numProc) {
  // grid -- vector of procBlocks (no need to split procBlocks or combine them
  // with manual decomposition)
  // bcs -- vector of boundary conditions for all blocks
  // numProc -- number of processors in run

  cout << "--------------------------------------------------------------------"
          "------------" << endl;
  cout << "Using cubic grid decomposition." << endl;

  decomposition decomp(grid.size(), numProc);
  // average number of cells per processor
  auto idealLoad = decomp.IdealLoad(grid);
  auto count = 0;

  const auto maxSplits = numProc * 10;
  while (decomp.MaxLoad(grid) / idealLoad > 1.1 && count < maxSplits) {
    auto loaded = 0.0;
    auto ol = decomp.MostOverloadedProc(grid, loaded);
    auto ul = decomp.MostUnderloadedProc(grid, loaded);

    string dir = "";
    auto blk = 0;
    auto ind = decomp.SendWholeOrSplit(grid, ol, ul, blk, dir);

    if (ind < 0) {  // send whole
      decomp.SendToProc(blk, ol, ul);
    } else {  // split/send
      auto newBlk = static_cast<int>(grid.size());
      // find all interblocks that could be altered by this split along with
      // their partners and orientation
      auto affectedConnections = GetBlockInterConnBCs(bcs, grid, blk);

      // split grid
      const auto uBlk = grid[blk].Split(dir, ind);
      grid.push_back(uBlk);

      // split bcs
      vector<boundarySurface> altSurf;
      auto newBcs = bcs[blk].Split(dir, ind, blk, newBlk, altSurf);
      bcs.push_back(newBcs);

      // update interblock partners affected by split
      for (auto &alt : altSurf) {
        bcs[alt.PartnerBlock()].DependentSplit(
            alt, affectedConnections.at(alt).first,
            affectedConnections.at(alt).second, alt.PartnerBlock(), dir, ind,
            blk, newBlk);
      }

      decomp.Split(blk, ind, dir);
      decomp.SendToProc(blk, ol, ul);
    }

    count++;
  }

  if (count == maxSplits - 1) {
    cout << "WARNING: Maximum number of splits in decomposition has been "
            "reached." << endl;
  }

  decomp.PrintDiagnostics(grid);
  cout << endl;
  cout << "Ideal Load: " << idealLoad << endl;
  cout << "Max Load: " << decomp.MaxLoad(grid) << endl;

  auto loaded = 0.0;
  auto ol = decomp.MostOverloadedProc(grid, loaded);
  cout << "Most overloaded processor is " << ol << "; overloaded by " << loaded
       << endl;
  auto ul = decomp.MostUnderloadedProc(grid, loaded);
  cout << "Most underloaded processor is " << ul << "; underloaded by "
       << loaded << endl;

  cout << "Ratio of most loaded processor to average processor is : "
       << decomp.MaxLoad(grid) / idealLoad << endl;
  cout << "--------------------------------------------------------------------"
          "------------" << endl << endl;

  return decomp;
}

// function to send each processor the number of procBlocks that it should
// contain
void SendNumProcBlocks(const vector<int> &loadBal, int &numProcBlock) {
  MPI_Scatter(&loadBal[0], 1, MPI_INT, &numProcBlock, 1, MPI_INT, ROOTP,
              MPI_COMM_WORLD);
}

/* Function to set custom MPI datatypes to allow for easier data transmission */
void SetDataTypesMPI(MPI_Datatype &MPI_vec3d,
                     MPI_Datatype &MPI_procBlockInts,
                     MPI_Datatype &MPI_connection,
                     MPI_Datatype &MPI_DOUBLE_5INT,
                     MPI_Datatype &MPI_vec3dMag,
                     MPI_Datatype &MPI_uncoupledScalar,
                     MPI_Datatype &MPI_tensorDouble) {
  // MPI_vec3d -- output MPI_Datatype for a vector3d<double>
  // MPI_procBlockInts -- output MPI_Datatype for 14 INTs (14 INTs in procBlock
  //                      class)
  // MPI_connection -- output MPI_Datatype to send connection class
  // MPI_DOUBLE_5INT -- output MPI_Datatype for a double followed by 5 ints
  // MPI_vec3dMag -- output MPI_Datatype for a unitVect3dMag<double>
  // MPI_uncoupledScalar -- output MPI_Datatype for a uncoupledScalar
  // MPI_tensorDouble -- output MPI_Datatype for a tensor<double>

  // create vector3d<double> MPI datatype
  MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_vec3d);
  MPI_Type_commit(&MPI_vec3d);

  // create vector3d<double> MPI datatype
  MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_vec3dMag);
  MPI_Type_commit(&MPI_vec3dMag);

  // create uncoupledScalar MPI datatype
  MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_uncoupledScalar);
  MPI_Type_commit(&MPI_uncoupledScalar);

  // create tensor<double> MPI datatype
  MPI_Type_contiguous(9, MPI_DOUBLE, &MPI_tensorDouble);
  MPI_Type_commit(&MPI_tensorDouble);

  // create MPI datatype for all the integers in the procBlock class
  MPI_Type_contiguous(15, MPI_INT, &MPI_procBlockInts);
  MPI_Type_commit(&MPI_procBlockInts);

  // create MPI datatype for a double followed by 5 ints
  int fieldCounts[2] = {1, 5};  // number of entries per field
  MPI_Datatype fieldTypes[2] = {MPI_DOUBLE, MPI_INT};  // field types
  MPI_Aint displacement[2], lBound, ext;
  resid res;  // dummy resid to get layout of class
  res.GetAddressesMPI(displacement);

  // make addresses relative to first field
  for (auto ii = 1; ii >= 0; ii--) {
    displacement[ii] -= displacement[0];
  }
  MPI_Type_create_struct(2, fieldCounts, displacement, fieldTypes,
                         &MPI_DOUBLE_5INT);

  // check that datatype has the correct extent, if it doesn't change the extent
  // this is necessary to portably send an array of this type
  MPI_Type_get_extent(MPI_DOUBLE_5INT, &lBound, &ext);
  if (ext != sizeof(res)) {
    MPI_Datatype temp = MPI_DOUBLE_5INT;
    MPI_Type_create_resized(temp, 0, sizeof(res), &MPI_DOUBLE_5INT);
    MPI_Type_free(&temp);
  }

  MPI_Type_commit(&MPI_DOUBLE_5INT);

  // create MPI datatype for connection class
  constexpr auto nFields = 12;
  int counts[nFields] = {2, 2, 2, 2, 2, 2,
                         2, 2, 2, 8, 1, 1};  // number of entries per field
  // field types
  MPI_Datatype types[nFields] = {MPI_INT, MPI_INT,    MPI_INT, MPI_INT,
                                 MPI_INT, MPI_INT,    MPI_INT, MPI_INT,
                                 MPI_INT, MPI_C_BOOL, MPI_INT, MPI_C_BOOL};
  MPI_Aint disp[nFields], lowerBound, extent;
  connection inter;  // dummy connection to get layout of class
  inter.GetAddressesMPI(disp);

  // make addresses relative to first field
  for (auto ii = nFields - 1; ii >= 0; ii--) {
    disp[ii] -= disp[0];
  }
  MPI_Type_create_struct(11, counts, disp, types, &MPI_connection);

  // check that datatype has the correct extent, if it doesn't change the extent
  // this is necessary to portably send an array of this type
  MPI_Type_get_extent(MPI_connection, &lowerBound, &extent);
  if (extent != sizeof(inter)) {
    auto temp = MPI_connection;
    MPI_Type_create_resized(temp, 0, sizeof(inter), &MPI_connection);
    MPI_Type_free(&temp);
  }

  MPI_Type_commit(&MPI_connection);
}

/* Function to free custom MPI datatypesn */
void FreeDataTypesMPI(MPI_Datatype &MPI_vec3d,
                      MPI_Datatype &MPI_procBlockInts,
                      MPI_Datatype &MPI_connection,
                      MPI_Datatype &MPI_DOUBLE_5INT,
                      MPI_Datatype &MPI_vec3dMag,
                      MPI_Datatype &MPI_uncoupledScalar,
                      MPI_Datatype &MPI_tensorDouble) {
  // MPI_vec3d -- MPI_Datatype for a vector3d<double>
  // MPI_procBlockInts -- MPI_Datatype for 14 INTs (14 INTs in procBlock class)
  // MPI_connection -- MPI_Datatype to send connection class
  // MPI_DOUBLE_5INT -- MPI_Datatype for a double followed by 5 ints
  // MPI_vec3dMag -- MPI_Datatype for a unitVect3dMag<double>
  // MPI_uncoupledScalar -- MPI_Datatype for a uncoupledScalar
  // MPI_tensorDouble -- MPI_Datatype for a tensor<double>

  // free vector3d<double> MPI datatype
  MPI_Type_free(&MPI_vec3d);

  // free unitVect3dMag<double> MPI datatype
  MPI_Type_free(&MPI_vec3dMag);

  // free MPI datatype for uncoupledScalar
  MPI_Type_free(&MPI_uncoupledScalar);

  // free MPI datatype for all the integers in the procBlock class
  MPI_Type_free(&MPI_procBlockInts);

  // free MPI datatype for a double followed by 5 ints
  MPI_Type_free(&MPI_DOUBLE_5INT);

  // free MPI datatype for connection class
  MPI_Type_free(&MPI_connection);

  // free MPI datatype for tensor<double> class
  MPI_Type_free(&MPI_tensorDouble);
}

/* Function to send procBlocks to the root processor. In this function, the
non-ROOT processors pack the procBlocks and send them to the ROOT processor.
The ROOT processor receives and unpacks the data from the non-ROOT processors.
This is used to get all the data on the ROOT processor to write out results.
*/
void GetProcBlocks(vector<procBlock> &blocks,
                   const vector<procBlock> &localBlocks, const int &rank,
                   const MPI_Datatype &MPI_uncoupledScalar,
                   const MPI_Datatype &MPI_vec3d,
                   const MPI_Datatype &MPI_tensorDouble, const input &inp) {
  // blocks -- full vector of all procBlocks. This is only used on ROOT
  //           processor, all other processors just need a dummy variable to
  //           call the function
  // localBlocks -- procBlocks local to each processor. These are sent to ROOT
  // rank -- processor rank. Used to determine if process should send or
  //         receive
  // MPI_uncoupledScalar -- MPI_Datatype used for uncoupledScalar transmission
  // MPI_vec3d -- MPI_Datatype used for vector3d<double> transmission
  // MPI_tensorDouble -- MPI_Datatype used for tensor<double> transmission
  // input -- input variables

  //--------------------------------------------------------------------------
  //                                      ROOT
  //--------------------------------------------------------------------------
  if (rank == ROOTP) {  // may have to recv and unpack data
    // loop over ALL blocks
    auto count = 0;
    for (auto &blk : blocks) {
      if (blk.Rank() == ROOTP) {  // data already on ROOT processor
        // assign local block to global block in order of local vector
        blk = localBlocks[blk.LocalPosition()];
      } else {  // recv data from sending processors
        blk.RecvUnpackSolMPI(MPI_uncoupledScalar, MPI_vec3d, MPI_tensorDouble,
                             inp);
      }
      count++;
    }
    //-------------------------------------------------------------------------
    //                                   NON - ROOT
    //-------------------------------------------------------------------------
  } else {  // pack and send data (non-root)
    // get vector of local positions
    vector<std::pair<int, int>> localPos;
    localPos.reserve(localBlocks.size());
    for (auto &lb : localBlocks) {
      localPos.push_back(std::make_pair(lb.LocalPosition(), lb.GlobalPos()));
    }
    // sort by global position
    // need to send data in order of global position, not local position to
    // prevent deadlock
    std::sort(
        std::begin(localPos), std::end(localPos),
        [](const auto &d1, const auto &d2) { return d1.second < d2.second; });

    for (auto &lp : localPos) {
      localBlocks[lp.first].PackSendSolMPI(MPI_uncoupledScalar, MPI_vec3d,
                                           MPI_tensorDouble);
    }
  }
}

/*function to broadcast a string from ROOT to all processors. This is needed
because it is not garunteed in the MPI standard that the commmand
line arguments will be on any processor but ROOT.  */
void BroadcastString(string &str) {
  // str -- string to broadcast to all processors

  // get size of string (+1 for c_str end character)
  auto strSize = static_cast<int>(str.size() + 1);
  MPI_Bcast(&strSize, 1, MPI_INT, ROOTP,
            MPI_COMM_WORLD);  // broadcast string size
  
  // allocate a char buffer of string size
  auto buf = std::make_unique<char[]>(strSize);
  snprintf(buf.get(), strSize, "%s", str.c_str());  // copy string into buffer
  MPI_Bcast(buf.get(), strSize, MPI_CHAR, ROOTP,
            MPI_COMM_WORLD);  // broadcast string as char

  // create new string and assign to old string
  string newStr(buf.get(), strSize - 1);  // -1 to not include c_str end character
  str = newStr;
}

// constructor with arguements
decomposition::decomposition(const int &num, const int &nProcs) {
  // num -- number of grid blocks

  // default configuration is all blocks on rank 0
  // default value for rank is 0
  vector<int> temp(num, 0);
  rank_ = temp;
  parBlock_ = temp;
  localPos_ = temp;
  for (auto ii = 0; ii < num; ii++) {
    parBlock_[ii] = ii;
    localPos_[ii] = ii;
  }

  // no splits for default
  vector<int> temp2;
  splitHistBlkLow_ = temp2;
  splitHistBlkUp_ = temp2;
  splitHistIndex_ = temp2;

  // no splits for default
  vector<string> temp3;
  splitHistDir_ = temp3;

  numProcs_ = nProcs;
}

/*Member function to determine the ideal load given the mesh. The ideal load the
 * the total number of cells divided by the number of processors.*/
double decomposition::IdealLoad(const vector<plot3dBlock> &grid) const {
  // grid -- vector of plot3dBlocks containing entire grid

  auto totalCells = 0;
  for (auto ii = 0U; ii < grid.size(); ii++) {
    totalCells += grid[ii].NumCells();
  }

  return static_cast<double>(totalCells) /
         static_cast<double>(numProcs_);
}

/*Member function to determine the maximum load (number of cells) on a
 * processor.*/
double decomposition::MaxLoad(const vector<plot3dBlock> &grid) const {
  // grid -- vector of plot3dBlocks containing entire grid (split for
  // decomposition)

  vector<int> load(numProcs_, 0);
  for (auto ii = 0U; ii < grid.size(); ii++) {
    load[rank_[ii]] += grid[ii].NumCells();
  }

  return static_cast<double>(*max_element(load.begin(), load.end()));
}

/*Member function to determine the minimum load (number of cells) on a
 * processor.*/
double decomposition::MinLoad(const vector<plot3dBlock> &grid) const {
  // grid -- vector of plot3dBlocks containing entire grid (split for
  // decomposition)

  vector<int> load(numProcs_, 0);
  for (auto ii = 0U; ii < grid.size(); ii++) {
    load[rank_[ii]] += grid[ii].NumCells();
  }

  return static_cast<double>(*min_element(load.begin(), load.end()));
}

/*Member function to determine the index of the maximum loaded (number of cells)
 * processor.*/
int decomposition::MostOverloadedProc(const vector<plot3dBlock> &grid,
                                      double &overload) const {
  // grid -- vector of plot3dBlocks containing entire grid (split for
  // decomposition)
  // overload -- how much the processor is overloaded by

  vector<int> load(numProcs_, 0);
  for (auto ii = 0U; ii < grid.size(); ii++) {
    load[rank_[ii]] += grid[ii].NumCells();
  }

  overload = this->MaxLoad(grid) - this->IdealLoad(grid);
  return distance(load.begin(), max_element(load.begin(), load.end()));
}

/*Member function to determine the index of the minimum loaded (number of cells)
 * processor.*/
int decomposition::MostUnderloadedProc(const vector<plot3dBlock> &grid,
                                       double &underload) const {
  // grid -- vector of plot3dBlocks containing entire grid (split for
  // decomposition)
  // underload -- how much processor is underloaded by

  vector<int> load(numProcs_, 0);
  for (auto ii = 0U; ii < grid.size(); ii++) {
    load[rank_[ii]] += grid[ii].NumCells();
  }

  underload = this->IdealLoad(grid) - this->MinLoad(grid);
  return distance(load.begin(), min_element(load.begin(), load.end()));
}

/*Member function to return the number of blocks on a given processor.*/
int decomposition::NumBlocksOnProc(const int &a) const {
  // a -- processor rank to find number of blocks on

  auto num = 0;
  for (auto &rank : rank_) {
    if (rank == a) {
      num++;
    }
  }
  return num;
}

/*Member function to return the number of blocks on all processors.*/
vector<int> decomposition::NumBlocksOnAllProc() const {
  vector<int> num(numProcs_, 0);
  for (auto &rank : rank_) {
    num[rank]++;
  }
  return num;
}


/*Member function to send a block to a given processor*/
void decomposition::SendToProc(const int &blk, const int &fromProc,
                               const int &toProc) {
  // blk -- sending block index
  // fromProc -- processor to send from
  // toProc -- processor to send to

  // only fields to change are local position on processor and processor rank

  auto oldPos = localPos_[blk];
  // local position is now equal to the number of blocks on given processor
  // (equal b/c indexing starts at 0)
  localPos_[blk] = this->NumBlocksOnProc(toProc);

  // change rank of procBlock
  rank_[blk] = toProc;

  // all procBlocks on same processor with a local position higher than oldPos
  // should be moved down one
  for (auto ii = 0U; ii < localPos_.size(); ii++) {
    if (rank_[ii] == fromProc &&
        localPos_[ii] > oldPos) {  // procBlock is on given processor
                                   // and in "higher" position than the
                                   // block to be moved
      localPos_[ii]--;
    }
  }
}

/*Member function to add data for a split*/
void decomposition::Split(const int &low, const int &ind, const string &dir) {
  // low -- index of lower block in split
  // dir -- direction of split
  // ind -- index of split

  splitHistBlkLow_.push_back(low);  // assign lower block index in split (given)
  // assign upper block index in split (one more than current max index)
  splitHistBlkUp_.push_back(rank_.size());
  splitHistIndex_.push_back(ind);  // assign index of split
  splitHistDir_.push_back(dir);    // assign split direction (given)

  rank_.push_back(rank_[low]);  // rank of upper portion of
                                // split block is same as lower
  // parent block of upper portion of split block is same as lower
  parBlock_.push_back(parBlock_[low]);

  // local position of upper portion of split block is equal to number of blocks
  // on processor (b/c indexing starts at 0)
  // -1 b/c indexing starts at 0
  localPos_.push_back(this->NumBlocksOnProc(rank_[low]) - 1);
}

int decomposition::GlobalPos(const int &rank, const int &localPos) const {
  auto globalPos = -1;
  for (auto ii = 0U; ii < rank_.size(); ++ii) {
    if (rank_[ii] == rank && localPos_[ii] == localPos) {
        globalPos = ii;
        break;
      }
  }
  MSG_ASSERT(globalPos >= 0, "could not find global pos");
  return globalPos;
}

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const decomposition &d) {
  // os -- stream to print to
  // d -- decomposition to print

  os << "Decomposition for " << d.NumProcs() << " processors" << endl;
  for (auto ii = 0; ii < d.Size(); ii++) {
    os << "Block: " << ii << "; Rank: " << d.Rank(ii)
       << ", Parent Block: " << d.ParentBlock(ii)
       << ", Local Position: " << d.LocalPosition(ii) << endl;
  }
  os << "Split History" << endl;
  for (auto ii = 0; ii < d.NumSplits(); ii++) {
    os << "Split Number: " << ii << "; Lower Index: " << d.SplitHistBlkLower(ii)
       << ", Upper Index: " << d.SplitHistBlkUpper(ii)
       << ", Direction: " << d.SplitHistDir(ii)
       << ", Split Index: " << d.SplitHistIndex(ii) << endl;
  }
  return os;
}

double decomposition::ProcLoad(const vector<plot3dBlock> &grid,
                               const int &proc) const {
  // grid -- vector of plot3dBlocks making up entire grid
  // proc -- rank of processor to calculate load for

  auto load = 0;
  for (auto ii = 0U; ii < grid.size(); ii++) {
    if (rank_[ii] == proc) {
      load += grid[ii].NumCells();
    }
  }
  return static_cast<double>(load);
}

double decomposition::LoadRatio(const vector<plot3dBlock> &grid,
                                const int &proc) const {
  // grid -- vector of plot3dBlocks making up entire grid
  // proc -- rank of processor to calculate ratio for

  auto ideal = this->IdealLoad(grid);
  auto load = this->ProcLoad(grid, proc);
  return fabs(1.0 - load / ideal);
}

/* Member function to determine whether to send a whole block or a split block
from a given processor. If it is determined to send a split block, the index of
the split is returned, and the direction string is changed to the appropriate
value. If a whole block is to be sent, the index returned is -1.
*/
int decomposition::SendWholeOrSplit(const vector<plot3dBlock> &grid,
                                    const int &send, const int &recv, int &blk,
                                    string &dir) const {
  // grid -- vector of plot3dBlocks making up entire grid
  // send -- rank of processor to sending block
  // recv -- rank of process to receive block
  // blk -- block to split or send
  // dir -- direction of split

  auto ind = -1;
  dir = "none";

  auto sendRatio = this->LoadRatio(grid, send);
  auto recvRatio = this->LoadRatio(grid, recv);
  auto ideal = this->IdealLoad(grid);
  auto sendLoad = this->ProcLoad(grid, send);
  auto recvLoad = this->ProcLoad(grid, recv);

  // find out if there is a block that can be sent that would improve both
  // ratios
  for (auto ii = 0U; ii < grid.size(); ii++) {
    if (rank_[ii] == send) {
      auto newSendRatio = fabs(
          1.0 - (sendLoad - static_cast<double>(grid[ii].NumCells())) / ideal);
      auto newRecvRatio = fabs(
          1.0 - (recvLoad + static_cast<double>(grid[ii].NumCells())) / ideal);
      if (newSendRatio < sendRatio &&
          newRecvRatio < recvRatio) {  // can send whole block
        blk = ii;
        return ind;
      }
    }
  }

  // find out which block to split - largest
  auto bSize = 0;
  for (auto ii = 0U; ii < grid.size(); ii++) {
    if (rank_[ii] == send) {
      if (grid[ii].NumCells() > bSize) {
        blk = ii;
        bSize = grid[ii].NumCells();
      }
    }
  }

  // get block split direction, plane size, and possible split length
  auto planeSize = 0;
  auto splitLen = 0;
  if (grid[blk].NumK() >= grid[blk].NumJ() &&
      grid[blk].NumK() >= grid[blk].NumI()) {
    dir = "k";
    // -1 to get cell sizes
    planeSize = grid[blk].NumCellsJ() * grid[blk].NumCellsI();
    splitLen = grid[blk].NumK();
  } else if (grid[blk].NumJ() >= grid[blk].NumI()) {
    dir = "j";
    planeSize = grid[blk].NumCellsK() * grid[blk].NumCellsI();
    splitLen = grid[blk].NumJ();
  } else {
    dir = "i";
    planeSize = grid[blk].NumCellsJ() * grid[blk].NumCellsK();
    splitLen = grid[blk].NumI();
  }

  // get split index
  // splitLen - 2 bc index is last cell that is kept in lower portion of split,
  // if entire length is kept, no split, need to keep at least 2 cells thick for
  // ghost cell passing
  // starting index at 2 so split is at least 2 cells thick for ghost cell
  // passing - ind is the face index to split at
  for (auto ii = 2; ii < splitLen - 2; ii++) {
    auto newSendRatio = fabs(
        1.0 - (sendLoad - static_cast<double>(planeSize * ii)) / ideal);
    auto newRecvRatio = fabs(
        1.0 - (recvLoad + static_cast<double>(planeSize * ii)) / ideal);
    if (newSendRatio < sendRatio &&
        newRecvRatio < recvRatio) {  // can send block at index
      sendRatio = newSendRatio;
      recvRatio = newRecvRatio;
      ind = ii;
    }
  }

  return ind;
}


void decomposition::PrintDiagnostics(const vector<plot3dBlock> &grid) const {
  cout << "Decomposition for " << numProcs_ << " processors" << endl;
  for (auto ii = 0U; ii < rank_.size(); ii++) {
    cout << "Block: " << ii << "; Rank: " << rank_[ii]
         << ", Parent Block: " << parBlock_[ii]
         << ", Local Position: " << localPos_[ii]
         << ", NumI: " << grid[ii].NumI() - 1
         << ", NumJ: " << grid[ii].NumJ() - 1
         << ", NumK: " << grid[ii].NumK() - 1
         << ", Num Cells: " << grid[ii].NumCells() << endl;
  }
  cout << "Split History" << endl;
  for (auto ii = 0U; ii < splitHistBlkLow_.size(); ii++) {
    cout << "Split Number: " << ii
         << "; Lower Index: " << splitHistBlkLow_[ii]
         << ", Upper Index: " << splitHistBlkUp_[ii]
         << ", Direction: " << splitHistDir_[ii]
         << ", Split Index: " << splitHistIndex_[ii] << endl;
  }
}

void decomposition::Broadcast() { 
  // broadcast rank
  auto vecSize = rank_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  rank_.resize(vecSize);
  MPI_Bcast(&(*std::begin(rank_)), vecSize, MPI_INT, ROOTP, MPI_COMM_WORLD);

  // broadcast parent block
  vecSize = parBlock_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  parBlock_.resize(vecSize);
  MPI_Bcast(&(*std::begin(parBlock_)), vecSize, MPI_INT, ROOTP, MPI_COMM_WORLD);

  // broadcast local position
  vecSize = localPos_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  localPos_.resize(vecSize);
  MPI_Bcast(&(*std::begin(localPos_)), vecSize, MPI_INT, ROOTP, MPI_COMM_WORLD);

  // broadcast lower block split history
  vecSize = splitHistBlkLow_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  splitHistBlkLow_.resize(vecSize);
  MPI_Bcast(&(*std::begin(splitHistBlkLow_)), vecSize, MPI_INT, ROOTP,
            MPI_COMM_WORLD);

  // broadcast upper block split history
  vecSize = splitHistBlkUp_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  splitHistBlkUp_.resize(vecSize);
  MPI_Bcast(&(*std::begin(splitHistBlkUp_)), vecSize, MPI_INT, ROOTP,
            MPI_COMM_WORLD);

  // broadcast index split history
  vecSize = splitHistIndex_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  splitHistIndex_.resize(vecSize);
  MPI_Bcast(&(*std::begin(splitHistIndex_)), vecSize, MPI_INT, ROOTP,
            MPI_COMM_WORLD);

  // broadcast direction split history
  vecSize = splitHistDir_.size(); 
  MPI_Bcast(&vecSize, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
  splitHistDir_.resize(vecSize);
  for (auto &dir : splitHistDir_) {
    BroadcastString(dir);
  }

  // broadcast number of processors
  MPI_Bcast(&numProcs_, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);
}

void BroadcastViscFaces(const MPI_Datatype &MPI_vec3d,
                        vector<vector3d<double>> &viscFaces) {
  // first determine the number of viscous faces and send that to all processors
  auto nFaces = static_cast<int>(viscFaces.size());
  MPI_Bcast(&nFaces, 1, MPI_INT, ROOTP, MPI_COMM_WORLD);

  viscFaces.resize(nFaces);  // allocate space to receive the viscous faces

  // broadcast all viscous faces to all processors
  MPI_Bcast(&viscFaces[0], viscFaces.size(), MPI_vec3d, ROOTP,
            MPI_COMM_WORLD);
}

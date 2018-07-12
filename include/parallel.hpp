/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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

#ifndef PARALLEL_HEADER_DEF
#define PARALLEL_HEADER_DEF

/* This header contains the function declarations for many of the parallel
 * functions in the code */

#include <iostream>
#include <vector>                  // vector
#include <string>                  // string
#include "mpi.h"                   // parallelism
#include "vector3d.hpp"
#include "blkMultiArray3d.hpp"

using std::vector;
using std::string;
using std::ios;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

// forward class declarations
class boundaryConditions;
class procBlock;
class plot3dBlock;
class connection;
class resid;
class genArray;

class decomposition {
  // rank of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> rank_;
  // parent block of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> parBlock_;
  // local position of each procBlock
  // (vector size equals number of procBlocks after decomp)
  vector<int> localPos_;
  // lower block of split (vector size equals number of splits)
  vector<int> splitHistBlkLow_;
  // upper block of split (vector size equals number of splits)
  vector<int> splitHistBlkUp_;
  // index of split (vector size equals number of splits)
  vector<int> splitHistIndex_;
  // direction of split (vector size equals number of splits)
  vector<string> splitHistDir_;
  int numProcs_;                     // number of processors

 public:
  // Constructor
  decomposition(const int&, const int&);
  decomposition() : decomposition(1, 1) {}

  // move constructor and assignment operator
  decomposition(decomposition&&) noexcept = default;
  decomposition& operator=(decomposition&&) noexcept = default;

  // copy constructor and assignment operator
  decomposition(const decomposition&) = default;
  decomposition& operator=(const decomposition&) = default;

  // Member functions
  int Rank(const int &a) const {return rank_[a];}
  int ParentBlock(const int &a) const {return parBlock_[a];}
  int LocalPosition(const int &a) const {return localPos_[a];}
  int NumProcs() const {return numProcs_;}
  double IdealLoad(const vector<plot3dBlock>&) const;
  double MaxLoad(const vector<plot3dBlock>&) const;
  double MinLoad(const vector<plot3dBlock>&) const;
  double ProcLoad(const vector<plot3dBlock>&, const int&) const;
  double LoadRatio(const vector<plot3dBlock>&, const int&) const;
  int MostOverloadedProc(const vector<plot3dBlock>&, double&) const;
  int MostUnderloadedProc(const vector<plot3dBlock>&, double&) const;
  int NumBlocksOnProc(const int&) const;
  vector<int> NumBlocksOnAllProc() const;
  int NumBlocks() const {return rank_.size();}
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
  template <typename T>
  void DecompArray(vector<blkMultiArray3d<T>> &) const;
  void PrintDiagnostics(const vector<plot3dBlock>&) const;

  // Destructor
  ~decomposition() noexcept {}
};

// function definitions
ostream & operator<< (ostream &os, const decomposition&);

decomposition ManualDecomposition(vector<plot3dBlock>&,
                                  vector<boundaryConditions>&, const int&);
decomposition CubicDecomposition(vector<plot3dBlock>&,
                                 vector<boundaryConditions>&, const int&);

void SendNumProcBlocks(const vector<int>&, int&);

void SetDataTypesMPI(MPI_Datatype &, MPI_Datatype &, MPI_Datatype &,
                     MPI_Datatype &, MPI_Datatype &, MPI_Datatype &,
                     MPI_Datatype &);
void FreeDataTypesMPI(MPI_Datatype &, MPI_Datatype &, MPI_Datatype &,
                      MPI_Datatype &, MPI_Datatype &, MPI_Datatype &,
                      MPI_Datatype &);

void MaxLinf(resid*, resid*, int*, MPI_Datatype*);

void BroadcastString(string& str);

void BroadcastViscFaces(const MPI_Datatype&, vector<vector3d<double>> &);


template <typename T>
void decomposition::DecompArray(vector<blkMultiArray3d<T>> &arr) const {
  // resize vector for split blocks
  arr.resize(this->NumBlocks());

  // split array by decomposition
  for (auto ii = 0; ii < this->NumSplits(); ++ii) {
    auto ind = this->SplitHistIndex(ii);
    auto lower = this->SplitHistBlkLower(ii);
    auto upper = this->SplitHistBlkUpper(ii);
    auto dir = this->SplitHistDir(ii);
    auto lowerArray = arr[lower].Slice(dir, {arr[lower].Start(dir), ind + 1});
    auto upperArray = arr[lower].Slice(dir, {ind, arr[lower].End(dir)});

    arr[lower] = lowerArray;
    arr[upper] = upperArray;
  }
}

#endif

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

#ifndef OUTPUTHEADERDEF  // only if the macro OUTPUTHEADERDEF is not defined
                         // execute these lines of code
#define OUTPUTHEADERDEF  // define the macro

/* This header contains the function declarations for file output.

It contains function headers to write out the grid at the cell centers in Plot3D
format, as well as the Plot3D function files.
It also writes out a master file in Ensight format to name the Plot3D functions.
*/

#include <fstream>
#include <iostream>
#include <vector>        // vector
#include <string>        // string
#include "multiArray3d.hpp"
#include "blkMultiArray3d.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

// forward class declarations
class procBlock;
class conserved;
class decomposition;
class physics;
class resid;
class input;
class primitive;
class residual;
class gridLevel;

// function definitions
template<typename T>
void WriteBlockDims(ofstream &, const vector<T> &, int = 0);

void WriteNodes(const string &, const vector<plot3dBlock>&);
void WriteCellCenter(const string &, const vector<procBlock> &,
                     const decomposition &, const input &);
void WriteWallFaceCenter(const string &, const vector<procBlock> &,
                         const double &);
void WriteOutput(const vector<procBlock> &, const physics &, const int &,
                 const decomposition &, const input &);
void WriteFunFile(const vector<procBlock> &, const vector<procBlock> &,
                  const physics &, const decomposition &, const string &,
                  const input &);
void WriteCenterFun(const vector<procBlock> &, const vector<procBlock> &,
                    const physics &, const int &, const decomposition &,
                    const input &);
void WriteNodeFun(const vector<procBlock> &, vector<procBlock> &,
                  const physics &phys, const int &, const decomposition &,
                  const input &);
void WriteWallFun(const vector<procBlock> &, const physics &phys, const int &,
                  const input &);
void WriteMeta(const input &, const int &, const bool &);
void WriteWallMeta(const input &, const int &);

void WriteRestart(const vector<procBlock> &, const physics &, const int &,
                  const decomposition &, const input &, const residual &);
void ReadRestart(gridLevel &, const string &, const decomposition &,
                 input &, const physics &, residual &,
                 const vector<vector3d<int>> &);

blkMultiArray3d<primitive> ReadSolFromRestart(ifstream &, const input &,
                                              const physics &,
                                              const vector<string> &,
                                              const int &, const int &,
                                              const int &, const int &);
blkMultiArray3d<conserved> ReadSolNm1FromRestart(ifstream &, const input &,
                                                 const physics &,
                                                 const vector<string> &,
                                                 const int &, const int &,
                                                 const int &, const int &);

void PrintResiduals(const input &, residual &, const residual &, const resid &,
                    const double &, const int &, const int &, ostream &);
void PrintHeaders(const input &, ostream &);

vector<procBlock> Recombine(const vector<procBlock> &, const decomposition &);
int SplitBlockNumber(const vector<procBlock> &, const decomposition &,
                     const int &, const int &, const int &, const int &);

// ---------------------------------------------------------------------------
// function definitions
template<typename T>
void WriteBlockDims(ofstream &outFile, const vector<T> &vars,
                    int numVars) {
  // write number of blocks to file
  auto numBlks = static_cast<int>(vars.size());
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  // loop over all blocks and write out imax, jmax, kmax, numVars
  for (auto &blk : vars) {
    auto dumInt = blk.NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = blk.NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = blk.NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    
    if (numVars > 0) {
      outFile.write(reinterpret_cast<char *>(&numVars), sizeof(numVars));
    }
  }
}



#endif

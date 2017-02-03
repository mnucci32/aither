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
#include <memory>        // unique_ptr

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;
using std::unique_ptr;
using std::ostream;

// forward class declarations
class procBlock;
class genArray;
class decomposition;
class idealGas;
class sutherland;
class resid;
class input;
class turbModel;

// function definitions
void WriteBlockDims(ofstream &, const vector<procBlock> &, int = 0);

void WriteCellCenter(const string &, const vector<procBlock> &,
                     const decomposition &, const double &);
void WriteFun(const vector<procBlock> &, const idealGas &,
              const sutherland &, const int &, const decomposition &,
              const input &, const unique_ptr<turbModel> &);
void WriteRes(const input &, const int &);
void WriteMeta(const input &, const int &);

void WriteRestart(const vector<procBlock> &, const idealGas &,
                  const sutherland &, const int &, const decomposition &,
                  const input &, const genArray &);
void ReadRestart(vector<procBlock> &, const string &, const input &,
                 const idealGas &, const sutherland &,
                 const unique_ptr<turbModel> &, genArray &);

void WriteResiduals(const input &, genArray &, const genArray &, const resid &,
                    const double &, const int &, const int &, ostream &);
void PrintResiduals(const input &, genArray &, const genArray &, const resid &,
                    const double &, const int &, const int &, ostream &);
void PrintHeaders(const input &, ostream &);

vector<procBlock> Recombine(const vector<procBlock> &, const decomposition &);
int SplitBlockNumber(const vector<procBlock> &, const decomposition &,
                     const int &, const int &, const int &, const int &);

#endif

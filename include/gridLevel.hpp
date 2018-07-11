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

#ifndef GRIDLEVEL_HEADER_DEF
#define GRIDLEVEL_HEADER_DEF

#include <vector>
#include <string>
#include "procBlock.hpp"
#include "boundaryConditions.hpp"
#include "vector3d.hpp"

using std::string;
using std::vector;

// forward class declaration
class plot3dBlock;
class input;
class decomposition;
class input;
class physics;
class residual;

class gridLevel {
  vector<procBlock> blocks_;
  vector<connection> connections_;

 public:
  // Constructor
  gridLevel(const vector<plot3dBlock>& mesh,
            const vector<boundaryConditions>& bcs, const decomposition& decomp,
            const physics& phys, const vector<vector3d<int>>& origGridSizes,
            const string& restartFile, input& inp, residual& first);

  // move constructor and assignment operator
  gridLevel(gridLevel&&) noexcept = default;
  gridLevel& operator=(gridLevel&&) noexcept = default;

  // copy constructor and assignment operator
  gridLevel(const gridLevel&) = default;
  gridLevel& operator=(const gridLevel&) = default;

  // Member functions
  int NumBlocks() const { return blocks_.size(); }
  const vector<procBlock>& Blocks() const { return blocks_; }
  const procBlock& Block(const int &ii) const { return blocks_[ii]; }

  int NumConnections() const { return connections_.size(); }
  const vector<connection>& Connections() const { return connections_; }
  const connection& Connection(const int& ii) const { return connections_[ii]; }

  // Destructor
  ~gridLevel() noexcept {}
};

#endif

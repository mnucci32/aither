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

#include <iostream>     // cout
#include <cstdlib>      // exit()
#include <vector>
#include <string>
#include "gridLevel.hpp"
#include "utility.hpp"
#include "parallel.hpp"
#include "input.hpp"
#include "boundaryConditions.hpp"
#include "plot3d.hpp"
#include "physicsModels.hpp"
#include "output.hpp"
#include "resid.hpp"
#include "vector3d.hpp"
#include "macros.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// constructor
gridLevel::gridLevel(const vector<plot3dBlock>& mesh,
                     const vector<boundaryConditions>& bcs,
                     const decomposition& decomp, const physics& phys,
                     const vector<vector3d<int>>& origGridSizes,
                     const string& restartFile, input& inp, residual& first) {
  auto connections_ = GetConnectionBCs(bcs, mesh, decomp, inp);
  blocks_.reserve(mesh.size());
  for (auto ll = 0U; ll < mesh.size(); ++ll) {
    blocks_.emplace_back(mesh[ll], decomp.ParentBlock(ll), bcs[ll], ll,
                         decomp.Rank(ll), decomp.LocalPosition(ll), inp);
    blocks_.back().InitializeStates(inp, phys);
    blocks_.back().AssignGhostCellsGeom();
  }

  // if restart, get data from restart file
  if (inp.IsRestart()) {
    ReadRestart(blocks_, restartFile, decomp, inp, phys, first, origGridSizes);
  }

  // Swap geometry for connection BCs
  for (auto& conn : connections_) {
    SwapGeomSlice(conn, blocks_[conn.BlockFirst()],
                  blocks_[conn.BlockSecond()]);
  }
  // Get ghost cell edge data
  for (auto& block : blocks_) {
    block.AssignGhostCellsGeomEdge();
  }
}

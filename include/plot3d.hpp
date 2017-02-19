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

#ifndef PLOT3DHEADERDEF  // only if the macro PLOT3DHEADERDEF is not defined
                         // execute these lines of code

#define PLOT3DHEADERDEF  // define the macro

// This header file contains the functions and data structures associated with
// dealing with PLOT3D grid files

#include <vector>
#include <string>
#include "vector3d.hpp"
#include "multiArray3d.hpp"

using std::vector;
using std::string;

//-------------------------------------------------------------------------
// Class for an individual plot3d block
class plot3dBlock {
  // by default everything above the public: declaration is private
  multiArray3d<vector3d<double>> coords_;  // coordinates of nodes in block

 public:
  // constructor -- create a plot3d block by passing the above quantities
  explicit plot3dBlock(const multiArray3d<vector3d<double>> &coordinates) :
      coords_(coordinates) {}
  plot3dBlock(const int &ii, const int &jj, const int &kk) :
      coords_(ii, jj, kk, 0) {}
  plot3dBlock() : coords_(1, 1, 1, 0) {}

  // move constructor and assignment operator
  plot3dBlock(plot3dBlock&&) noexcept = default;
  plot3dBlock& operator=(plot3dBlock&&) noexcept = default;

  // copy constructor and assignment operator
  plot3dBlock(const plot3dBlock&) = default;
  plot3dBlock& operator=(const plot3dBlock&) = default;

  // member functions
  multiArray3d<double> Volume() const;
  multiArray3d<unitVec3dMag<double>> FaceAreaI() const;
  multiArray3d<unitVec3dMag<double>> FaceAreaJ() const;
  multiArray3d<unitVec3dMag<double>> FaceAreaK() const;
  multiArray3d<vector3d<double>> Centroid() const;
  multiArray3d<vector3d<double>> FaceCenterI() const;
  multiArray3d<vector3d<double>> FaceCenterJ() const;
  multiArray3d<vector3d<double>> FaceCenterK() const;

  int NumI() const { return coords_.NumI(); }
  int NumJ() const { return coords_.NumJ(); }
  int NumK() const { return coords_.NumK(); }
  int NumCellsI() const { return coords_.NumI() - 1; }
  int NumCellsJ() const { return coords_.NumJ() - 1; }
  int NumCellsK() const { return coords_.NumK() - 1; }
  int NumCells() const {
    return this->NumCellsI() * this->NumCellsJ() * this->NumCellsK();
  }
  double X(const int &ii, const int &jj, const int &kk) const {
    return coords_(ii, jj, kk)[0];
  }
  double Y(const int &ii, const int &jj, const int &kk) const {
    return coords_(ii, jj, kk)[1];
  }
  double Z(const int &ii, const int &jj, const int &kk) const {
    return coords_(ii, jj, kk)[2];
  }
  vector3d<double> Coords(const int &ii, const int &jj, const int &kk) const {
    return coords_(ii, jj, kk);
  }

  void Split(const string &, const int &, plot3dBlock &, plot3dBlock &) const;
  void Join(const plot3dBlock &, const string &);

  // destructor
  ~plot3dBlock() noexcept {}
};

//-------------------------------------------------------------------------
// function declarations
vector<plot3dBlock> ReadP3dGrid(const string &, const double &, double &);


#endif

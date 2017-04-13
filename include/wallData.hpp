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

#ifndef WALLDATAHEADERDEF
#define WALLDATAHEADERDEF

/* This header file contains the wallData class
*/

#include <vector>
#include "multiArray3d.hpp"
#include "inputStates.hpp"
#include "boundaryConditions.hpp"
#include "range.hpp"

using std::vector;

class wallData {
  double inviscidForce_;
  double viscousForce_;
  viscousWall bcData_;
  boundarySurface surf_;
  struct wallVars;  // forward declaration
  multiArray3d<wallVars> data_;

 public:
  // constructor
  wallData(const boundarySurface &surf, const viscousWall &bc)
      : inviscidForce_(0.0),
        viscousForce_(0.0),
        bcData_(bc),
        surf_(surf),
        data_(surf.NumI(), surf.NumJ(), surf.NumK(), 0) {}

  // move constructor and assignment operator
  wallData(wallData &&) = default;
  wallData &operator=(wallData &&) = default;

  // copy constructor and assignment operator
  wallData(const wallData &) = default;
  wallData &operator=(const wallData &) = default;

  // member functions
  double InviscidForce() const { return inviscidForce_; }
  double ViscousForce() const { return viscousForce_; }
  vector3d<double> WallShearStress(const int &ii, const int &jj,
                                   const int &kk) const;
  double WallHeatFlux(const int &ii, const int &jj, const int &kk) const;
  double Yplus(const int &ii, const int &jj, const int &kk) const;
  double WallTemperature(const int &ii, const int &jj, const int &kk) const;
  double WallEddyViscosity(const int &ii, const int &jj, const int &kk) const;
  double WallViscosity(const int &ii, const int &jj, const int &kk) const;
  double WallDensity(const int &ii, const int &jj, const int &kk) const;
  double WallFrictionVelocity(const int &ii, const int &jj,
                              const int &kk) const;

  // operator overloads
  wallVars &operator()(const int &ii, const int &jj, const int &kk) {
    return data_(ii - surf_.IMin(), jj - surf_.JMin(), kk - surf_.KMin());
   }
   const wallVars &operator()(const int &ii, const int &jj,
                              const int &kk) const {
     return data_(ii - surf_.IMin(), jj - surf_.JMin(), kk - surf_.KMin());
   }

   // destructor
   ~wallData() noexcept {}
};


#endif

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

#include "wallData.hpp"
#include "vector3d.hpp"

// member functions
vector3d<double> wallData::WallShearStress(const int &ii, const int &jj,
                                           const int &kk) const {
  return (*this)(ii, jj, kk).shearStress_;
}

double wallData::WallHeatFlux(const int &ii, const int &jj,
                              const int &kk) const {
  return (*this)(ii, jj, kk).heatFlux_;
}

double wallData::Yplus(const int &ii, const int &jj, const int &kk) const {
  return (*this)(ii, jj, kk).yplus_;
}

double wallData::WallTemperature(const int &ii, const int &jj,
                                 const int &kk) const {
  return (*this)(ii, jj, kk).temperature_;
}

double wallData::WallEddyViscosity(const int &ii, const int &jj,
                                   const int &kk) const {
  return (*this)(ii, jj, kk).turbEddyVisc_;
}

double wallData::WallViscosity(const int &ii, const int &jj,
                               const int &kk) const {
  return (*this)(ii, jj, kk).viscosity_;
}

double wallData::WallDensity(const int &ii, const int &jj,
                             const int &kk) const {
  return (*this)(ii, jj, kk).density_;
}

double wallData::WallFrictionVelocity(const int &ii, const int &jj,
                                      const int &kk) const {
  return (*this)(ii, jj, kk).frictionVelocity_;
}

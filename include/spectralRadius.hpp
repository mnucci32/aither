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

#ifndef SPECTRALRADIUSHEADERDEF
#define SPECTRALRADIUSHEADERDEF

/* This header contains the functions to calculate the inviscid and viscous
 * spectral radii
 */

#include <string>                  // string
#include <memory>                  // unique_ptr
#include "vector3d.hpp"            // vector3d

using std::unique_ptr;

// forward class declarations
class eos;
class thermodynamic;
class transport;
class primitive;

double InvCellSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                             const unitVec3dMag<double> &,
                             const unique_ptr<thermodynamic> &,
                             const unique_ptr<eos> &);
double InvFaceSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                             const unique_ptr<thermodynamic> &,
                             const unique_ptr<eos> &);

double ViscCellSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                              const unitVec3dMag<double> &,
                              const unique_ptr<thermodynamic> &,
                              const unique_ptr<eos> &,
                              const unique_ptr<transport> &, const double &,
                              const double &, const double &,
                              const unique_ptr<turbModel> &);
double ViscFaceSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                              const unique_ptr<thermodynamic> &,
                              const unique_ptr<eos> &,
                              const unique_ptr<transport> &, const double &,
                              const double &, const double &,
                              const unique_ptr<turbModel> &);

double CellSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                          const unitVec3dMag<double> &,
                          const unique_ptr<thermodynamic> &,
                          const unique_ptr<eos> &,
                          const unique_ptr<transport> &, const double &,
                          const double &, const double &,
                          const unique_ptr<turbModel> &, const bool &);
double FaceSpectralRadius(const primitive &, const unitVec3dMag<double> &,
                          const unique_ptr<thermodynamic> &,
                          const unique_ptr<eos> &,
                          const unique_ptr<transport> &, const double &,
                          const double &, const double &,
                          const unique_ptr<turbModel> &, const bool &);

#endif

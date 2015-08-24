/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef GRADIENTSHEADERDEF  // only if the macro GRADIENTSHEADERDEF is not
                            // defined execute these lines of code
#define GRADIENTSHEADERDEF  // define the macro

/* This header contains the gradients class.

   The gradients class stores the gradient terms for the inviscid fluxes
   viscous fluxes, and source terms. */

#include <vector>  // vector
#include <iostream>
#include "vector3d.hpp"
#include "tensor.hpp"
#include "multiArray3d.hpp"

using std::cout;
using std::endl;
using std::cerr;

// forward class declarations
class procBlock;
class idealGas;

class gradients {
  // velocity gradients at cell faces
  multiArray3d<tensor<double> > velocityI_;
  multiArray3d<tensor<double> > velocityJ_;
  multiArray3d<tensor<double> > velocityK_;

  // temperature gradients at cell faces
  multiArray3d<vector3d<double> > temperatureI_;
  multiArray3d<vector3d<double> > temperatureJ_;
  multiArray3d<vector3d<double> > temperatureK_;

  // tke gradients at cell faces
  multiArray3d<vector3d<double> > tkeI_;
  multiArray3d<vector3d<double> > tkeJ_;
  multiArray3d<vector3d<double> > tkeK_;

  // omega gradients at cell faces
  multiArray3d<vector3d<double> > omegaI_;
  multiArray3d<vector3d<double> > omegaJ_;
  multiArray3d<vector3d<double> > omegaK_;

 public:
  // constructors
  gradients();
  gradients(const bool &, const procBlock &, const idealGas &);

  // member functions
  tensor<double> VelGradI(const int &a) const { return velocityI_[a]; }
  tensor<double> VelGradJ(const int &a) const { return velocityJ_[a]; }
  tensor<double> VelGradK(const int &a) const { return velocityK_[a]; }
  vector3d<double> TempGradI(const int &a) const { return temperatureI_[a]; }
  vector3d<double> TempGradJ(const int &a) const { return temperatureJ_[a]; }
  vector3d<double> TempGradK(const int &a) const { return temperatureK_[a]; }
  vector3d<double> TkeGradI(const int &a) const { return tkeI_[a]; }
  vector3d<double> TkeGradJ(const int &a) const { return tkeJ_[a]; }
  vector3d<double> TkeGradK(const int &a) const { return tkeK_[a]; }
  vector3d<double> OmegaGradI(const int &a) const { return omegaI_[a]; }
  vector3d<double> OmegaGradJ(const int &a) const { return omegaJ_[a]; }
  vector3d<double> OmegaGradK(const int &a) const { return omegaK_[a]; }

  // functions to get number of cells in block
  int NumI() const { return velocityI_.NumI() - 1; }
  int NumJ() const { return velocityI_.NumJ(); }
  int NumK() const { return velocityI_.NumK(); }

  tensor<double> VelGradCell(const int &, const int &, const int &) const;
  vector3d<double> TempGradCell(const int &, const int &, const int &) const;
  vector3d<double> TkeGradCell(const int &, const int &, const int &) const;
  vector3d<double> OmegaGradCell(const int &, const int &, const int &) const;

  // destructor
  ~gradients() {}
};

// function definitions

#endif

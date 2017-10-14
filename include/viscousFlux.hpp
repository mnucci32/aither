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

#ifndef VISCFLUXHEADERDEF  // only if the macro VISCFLUXHEADERDEF is not defined
                           // execute these lines of code
#define VISCFLUXHEADERDEF  // define the macro

#include <iostream>        // cout
#include <memory>          // unique_ptr
#include "vector3d.hpp"    // vector3d
#include "tensor.hpp"      // tensor
#include "varArray.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declarations
class primitive;
class eos;
class transport;
class thermodynamic;
class turbModel;
class squareMatrix;
struct wallVars;

class viscousFlux : public varArray {
 public:
  // constructors
  viscousFlux() : varArray() {}
  viscousFlux(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}

  // move constructor and assignment operator
  viscousFlux(viscousFlux&&) noexcept = default;
  viscousFlux& operator=(viscousFlux&&) noexcept = default;

  // copy constructor and assignment operator
  viscousFlux(const viscousFlux&) = default;
  viscousFlux& operator=(const viscousFlux&) = default;

  // member functions
  void CalcFlux(const tensor<double> &, const unique_ptr<transport> &,
                const unique_ptr<thermodynamic> &,
                const unique_ptr<eos> &eqnState, const vector3d<double> &,
                const vector3d<double> &, const vector3d<double> &,
                const vector3d<double> &, const unique_ptr<turbModel> &,
                const primitive &, const double &, const double &,
                const double &);
  wallVars CalcWallFlux(const tensor<double> &, const unique_ptr<transport> &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<eos> &, const vector3d<double> &,
                        const vector3d<double> &, const vector3d<double> &,
                        const vector3d<double> &, const unique_ptr<turbModel> &,
                        const primitive &, const double &, const double &,
                        const double &);
  void CalcWallLawFlux(const vector3d<double> &, const double &, const double &,
                       const double &, const vector3d<double> &,
                       const vector3d<double> &, const vector3d<double> &,
                       const vector3d<double> &, const unique_ptr<turbModel> &);

  // destructor
  ~viscousFlux() noexcept {}
};

// function definitions
ostream &operator<<(ostream &os, const viscousFlux &);

#endif

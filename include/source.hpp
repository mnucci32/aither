/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef SOURCEHEADERDEF  // only if the macro SOURCEHEADERDEF is not defined
                         // execute these lines of code
#define SOURCEHEADERDEF  // define the macro

/* This header contains the source class.

   The source class stores the source terms for the Euler and Navier-Stokes
   equations. */

#include <iostream>
#include <memory>
#include "varArray.hpp"
#include "arrayView.hpp"
#include "vector3d.hpp"
#include "tensor.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declaration
class turbModel;
class transport;
class squareMatrix;
class physics;

class source : public residual {
 public:
  // constructor
  source(const int &numEqns, const int &numSpecies)
      : residual(numEqns, numSpecies) {}

  // move constructor and assignment operator
  source(source&&) noexcept = default;
  source& operator=(source&&) noexcept = default;

  // copy constructor and assignment operator
  source(const source&) = default;
  source& operator=(const source&) = default;

  // member functions
  squareMatrix CalcTurbSrc(const unique_ptr<turbModel> &, const primitiveView &,
                           const tensor<double> &, const vector3d<double> &,
                           const vector3d<double> &, const vector3d<double> &,
                           const unique_ptr<transport> &, const double &,
                           const double &, const double &, const double &,
                           const double &);
  squareMatrix CalcChemSrc(const physics &, const primitiveView &,
                           const double &);

  // destructor
  ~source() noexcept {}
};

ostream &operator<<(ostream &os, const source &);


// function definitions -------------------------------------


#endif

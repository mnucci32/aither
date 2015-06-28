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

#ifndef INVFLUXHEADERDEF  // only if the macro INVFLUXHEADERDEF is not defined
                          // execute these lines of code
#define INVFLUXHEADERDEF  // define the macro

#include <vector>        // vector
#include <string>        // string
#include <iostream>      // cout
#include "vector3d.hpp"  // vector3d
#include "macros.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

// forward class declaration
class idealGas;
class primVars;
class genArray;
class squareMatrix;

class inviscidFlux {
  double data_[NUMVARS];  // rho dot velocity vector
  // rho dot velocity vector * u-velocity + pressure * i-dir-vector
  // rho dot velocity vector * v-velocity + pressure * j-dir-vector
  // rho dot velocity vector * w-velocity + pressure * k-dir-vector
  // rho dot velocity vector * enthalpy

 public:
  // constructors
  inviscidFlux() : data_{0.0} {}
  inviscidFlux(const primVars&, const idealGas&, const vector3d<double>&);
  inviscidFlux(const genArray&, const idealGas&, const vector3d<double>&);

  // member functions
  double RhoVel() const { return data_[0]; }
  double RhoVelU() const { return data_[1]; }
  double RhoVelV() const { return data_[2]; }
  double RhoVelW() const { return data_[3]; }
  double RhoVelH() const { return data_[4]; }
  double RhoVelK() const { return data_[5]; }
  double RhoVelO() const { return data_[6]; }

  void RoeFlux(const inviscidFlux&, const genArray&);

  inviscidFlux operator*(const double&);
  inviscidFlux operator/(const double&);

  inviscidFlux operator+(const inviscidFlux&) const;
  inviscidFlux operator-(const inviscidFlux&) const;

  friend ostream& operator<<(ostream& os, inviscidFlux&);
  friend inviscidFlux operator*(const double&, const inviscidFlux&);
  friend inviscidFlux operator/(const double&, const inviscidFlux&);

  genArray ConvertToGenArray() const;

  // destructor
  ~inviscidFlux() {}
};

// function definitions
// function to calculate Roe flux with entropy fix
inviscidFlux RoeFlux(const primVars&, const primVars&, const idealGas&,
                     const vector3d<double>&, double&);

// function to calculate Roe flux with entropy fix for implicit methods
void ApproxRoeFluxJacobian(const primVars&, const primVars&, const idealGas&,
                           const vector3d<double>&, double&, squareMatrix&,
                           squareMatrix&);

genArray ConvectiveFluxUpdate(const primVars&, const idealGas&,
                              const vector3d<double>&, const genArray&);

#endif

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

#ifndef SOURCEHEADERDEF  // only if the macro SOURCEHEADERDEF is not defined
                         // execute these lines of code
#define SOURCEHEADERDEF  // define the macro

/* This header contains the source class.

   The source class stores the source terms for the Euler and Navier-Stokes
   equations. */

#include <vector>  // vector
#include <string>  // string
#include <iostream>
#include "macros.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

// forward class declaration
class primVars;
class turbModel;
class sutherland;
class idealGas;
class gradients;

class source {
  double data_[NUMVARS];  // source variables at cell center

 public:
  // constructors
  source() : data_{0.0} {}

  // member functions
  double SrcMass() const { return data_[0]; }
  double SrcMomX() const { return data_[1]; }
  double SrcMomY() const { return data_[2]; }
  double SrcMomZ() const { return data_[3]; }
  double SrcEngy() const { return data_[4]; }
  double SrcTke() const { return data_[5]; }
  double SrcOmg() const { return data_[6]; }

  void CalcTurbSrc(const turbModel *, const primVars &, const gradients &,
                   const sutherland &, const idealGas &, const int &,
                   const int &, const int &);

  // operator overloads for addition and subtraction of states
  source operator+(const source &) const;
  source operator-(const source &) const;
  source operator*(const source &) const;
  source operator/(const source &) const;

  source operator+(const double &) const;
  source operator-(const double &) const;
  source operator*(const double &) const;
  source operator/(const double &) const;

  friend source operator+(const double &, const source &);
  friend source operator-(const double &, const source &);
  friend source operator*(const double &, const source &);
  friend source operator/(const double &, const source &);
  friend ostream &operator<<(ostream &os, const source &);

  // destructor
  ~source() {}
};

// function definitions

#endif

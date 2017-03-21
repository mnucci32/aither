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

#ifndef WALLLAWHEADERDEF
#define WALLLAWHEADERDEF

/* This defines the classes for the turbulent wall law implementation */

#include "vector3d.hpp"
#include "tensor.hpp"

// forward class declaration
class primVars;
class idealGas;
class sutherland;

class wallLaw {
  const double vonKarmen_;
  const double wallConst_;

  // private member functions
  double Gamma(const idealGas &, const double &, const double &) const;
  double Beta() const {return 0.0;}

 public:
  // constructor
  wallLaw(const double &k, const double &c) : vonKarmen_(k), wallConst_(c) {}
  wallLaw() : wallLaw(0.41, 5.1) {}
  
  // move constructor and assignment operator
  wallLaw(wallLaw&&) = default;
  wallLaw& operator=(wallLaw&&) = default;

  // copy constructor and assignment operator
  wallLaw(const wallLaw&) = default;
  wallLaw& operator=(const wallLaw&) = default;

  // member functions
  double VonKarmen() const {return vonKarmen_;}
  double WallConstant() const {return wallConst_;}
  double WallShearStress(const primVars &, const vector3d<double> &,
                         const tensor<double> &, const idealGas &,
                         const double &, const double &, const double &,
                         const double &, const double &,
                         const sutherland &) const;

  // destructor
  ~wallLaw() noexcept {}
};


// function declarations



#endif

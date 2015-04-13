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

#ifndef TURBHEADERDEF             //only if the macro TURBHEADERDEF is not defined execute these lines of code
#define TURBHEADERDEF             //define the macro

/* This defines the classes for the turbulence  models implemented in the code. */

#include <vector>  //vector
#include <string>  //string
#include <math.h>       //sqrt

using std::vector;
using std::string;

class turbModel {
  string eddyViscMethod;

 public:
  //constructor
  turbModel() {
    string temp = "boussinesq";
    eddyViscMethod = temp;
  }
 turbModel(const string &meth) : eddyViscMethod(meth) {}

  //member functions
  virtual double BoussinesqEddyVisc()const=0;

  //destructor
  virtual ~turbModel() {}

};


class turbNone : public turbModel {

 public:
  //constructor
 turbNone() : turbModel() {}

  //member functions
  virtual double BoussinesqEddyVisc()const{return 0.0;};

  //destructor
  virtual ~turbNone() {}

};


class turbKWWilcox : public turbModel {

 public:
  //constructor
 turbKWWilcox() : turbModel() {}

  //member functions
  virtual double BoussinesqEddyVisc()const{return 0.0;};

  //destructor
  virtual ~turbKWWilcox() {}

};

//function declarations


#endif

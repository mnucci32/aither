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

#ifndef GRADIENTSHEADERDEF             //only if the macro GRADIENTSHEADERDEF is not defined execute these lines of code
#define GRADIENTSHEADERDEF             //define the macro

/* This header contains the gradients class.

   The gradients class stores the gradients needed for the calculation of the viscous flux and source terms. */

#include <vector>  //vector
#include <string>  //string
#include "vector3d.hpp" //vector3d
#include "tensor.hpp" //tensor
#include "plot3d.hpp" //vector3d
#include "primVars.hpp"
#include <iostream>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class gradients {
  vector<tensor<double> > velocity;              //velocity gradients at cell face
  vector<vector3d<double> > temperature;         //temperature gradients at cell face
  vector<vector3d<double> > tke;                 //tke gradients at cell face
  vector<vector3d<double> > omega;               //omega gradients at cell face

 public:
  //constructors
  gradients(const bool&, const int&, const int&, const int&, const vector<primVars>&);

  //member functions
  tensor<double> VelGrad(const int &a) const {return velocity[a];}
  vector3d<double> TempGrad(const int &a) const {return temperature[a];}
  vector3d<double> TkeGrad(const int &a) const {return tke[a];}
  vector3d<double> OmegaGrad(const int &a) const {return omega[a];}

  //destructor
  ~gradients() {}

};

//function definitions

#endif

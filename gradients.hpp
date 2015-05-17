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
#include "procBlock.hpp"
#include <iostream>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class gradients {
  vector<tensor<double> > velocityI;              //velocity gradients at cell i-face
  vector<tensor<double> > velocityJ;              //velocity gradients at cell j-face
  vector<tensor<double> > velocityK;              //velocity gradients at cell k-face
  vector<vector3d<double> > temperatureI;         //temperature gradients at cell i-face
  vector<vector3d<double> > temperatureJ;         //temperature gradients at cell j-face
  vector<vector3d<double> > temperatureK;         //temperature gradients at cell k-face
  vector<vector3d<double> > tkeI;                 //tke gradients at cell i-face
  vector<vector3d<double> > tkeJ;                 //tke gradients at cell j-face
  vector<vector3d<double> > tkeK;                 //tke gradients at cell k-face
  vector<vector3d<double> > omegaI;               //omega gradients at cell i-face
  vector<vector3d<double> > omegaJ;               //omega gradients at cell j-face
  vector<vector3d<double> > omegaK;               //omega gradients at cell k-face

 public:
  //constructors
  gradients();
  gradients(const bool&, const procBlock&);

  //member functions
  tensor<double> VelGradI(const int &a) const {return velocityI[a];}
  tensor<double> VelGradJ(const int &a) const {return velocityJ[a];}
  tensor<double> VelGradK(const int &a) const {return velocityK[a];}
  vector3d<double> TempGradI(const int &a) const {return temperatureI[a];}
  vector3d<double> TempGradJ(const int &a) const {return temperatureJ[a];}
  vector3d<double> TempGradK(const int &a) const {return temperatureK[a];}
  vector3d<double> TkeGradI(const int &a) const {return tkeI[a];}
  vector3d<double> TkeGradJ(const int &a) const {return tkeJ[a];}
  vector3d<double> TkeGradK(const int &a) const {return tkeK[a];}
  vector3d<double> OmegaGradI(const int &a) const {return omegaI[a];}
  vector3d<double> OmegaGradJ(const int &a) const {return omegaJ[a];}
  vector3d<double> OmegaGradK(const int &a) const {return omegaK[a];}

  //destructor
  ~gradients() {}

};

//function definitions

#endif

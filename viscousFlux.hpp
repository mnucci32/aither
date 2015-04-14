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

#ifndef VISCFLUXHEADERDEF             //only if the macro VISCFLUXHEADERDEF is not defined execute these lines of code
#define VISCFLUXHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.hpp" //vector3d
#include "tensor.hpp" //tensor
#include "eos.hpp"  //idealGas
#include "primVars.hpp" //primVars
#include "matrix.hpp" //squareMatrix
#include "input.hpp" //input
#include "turbulence.hpp" //turbModel
#include <iostream> //cout

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;

class viscousFlux {
  double data[NUMVARS-1];   //viscous flux for x-momentum equation
                            //viscous flux for y-momentum equation
                            //viscous flux for z-momentum equation
                            //viscous flux for energy equation

 public:
  //constructors
  viscousFlux() : data{0.0, 0.0, 0.0, 0.0, 0.0, 0.0} {}
  viscousFlux( const tensor<double>&, const vector3d<double>&, const double&, const double&, const sutherland&, const idealGas&, const vector3d<double>&, const vector3d<double>&, const double& );

  //member functions
  double MomX() const {return data[0];}
  double MomY() const {return data[1];}
  double MomZ() const {return data[2];}
  double Engy() const {return data[3];}

  viscousFlux operator * (const double&);
  viscousFlux operator / (const double&);

  friend ostream & operator<< (ostream &os, viscousFlux&);
  friend viscousFlux operator * (const double&, const viscousFlux&);
  friend viscousFlux operator / (const double&, const viscousFlux&);

  //destructor
  ~viscousFlux() {}

};

//function definitions
void CalcTSLFluxJac( const double&, const double&, const idealGas&, const vector3d<double>&, const primVars&, const primVars&, 
		     const double&, squareMatrix&, squareMatrix&, const sutherland&, const double&);

tensor<double> CalcVelGradTSL(const primVars&, const primVars&, const vector3d<double>&, const double&);

#endif

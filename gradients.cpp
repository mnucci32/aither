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

#include "gradients.hpp"
#include <iostream>
#include <vector>
#include <string>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

//constructor
gradients::gradients(){

  //default to 2 faces
  velocityI.resize(2);
  velocityJ.resize(2);
  velocityK.resize(2);

  temperatureI.resize(2);
  temperatureJ.resize(2);
  temperatureK.resize(2);

  tkeI.resize(2);
  tkeJ.resize(2);
  tkeK.resize(2);

  omegaI.resize(2);
  omegaJ.resize(2);
  omegaK.resize(2);

}

//alternate constructor to calculate gradients from state variables
gradients::gradients( const bool &turbFlag, const procBlock &blk ){
  // turbFlag -- flag to determine if simulation is turbulent and therefore requires tke and omega gradients
  // blk -- procBlock to calculate gradients for

  velocityI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
  velocityJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
  velocityK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

  temperatureI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
  temperatureJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
  temperatureK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

  if (turbFlag){
    tkeI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
    tkeJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
    tkeK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

    omegaI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
    omegaJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
    omegaK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );
  }
  else{
    tkeI.resize(2);
    tkeJ.resize(2);
    tkeK.resize(2);

    omegaI.resize(2);
    omegaJ.resize(2);
    omegaK.resize(2);
  }

}

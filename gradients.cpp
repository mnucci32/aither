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

#include "procBlock.hpp"
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
gradients::gradients( const bool &turbFlag, const procBlock &blk, const idealGas &eos ){
  // turbFlag -- flag to determine if simulation is turbulent and therefore requires tke and omega gradients
  // blk -- procBlock to calculate gradients for
  // eos -- equation of state

  //size vectors for input procBlock
  velocityI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
  velocityJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
  velocityK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

  temperatureI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
  temperatureJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
  temperatureK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

  //if flow is not turbulent, don't need tke and omega gradients so they can be set to a minimum size
  if (turbFlag){
    tkeI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
    tkeJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
    tkeK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );

    omegaI.resize( (blk.NumI() + 1) * blk.NumJ() * blk.NumK() );
    omegaJ.resize( blk.NumI() * (blk.NumJ() + 1) * blk.NumK() );
    omegaK.resize( blk.NumI() * blk.NumJ() * (blk.NumK() + 1) );
  }
  else{
    tkeI.resize(1);
    tkeJ.resize(1);
    tkeK.resize(1);

    omegaI.resize(1);
    omegaJ.resize(1);
    omegaK.resize(1);
  }

  //loop over i-faces and calculate gradients
  //loop over all physical faces
  for ( int kk = 0; kk < blk.NumK(); kk++){   
    for ( int jj = 0; jj < blk.NumJ(); jj++){    
      for ( int ii = 0; ii < blk.NumI() + 1; ii++){      

	//get face location
	int loc = GetLoc1D(ii, jj, kk, blk.NumI() + 1, blk.NumJ());

	if (turbFlag){
	  blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityI[loc], temperatureI[loc], tkeI[loc], omegaI[loc]);
	}
	else{
	  vector3d<double> dum1, dum2;
	  blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityI[loc], temperatureI[loc], dum1, dum2);
	}

      }
    }
  }


  //loop over j-faces and calculate gradients
  //loop over all physical faces
  for ( int kk = 0; kk < blk.NumK(); kk++){   
    for ( int jj = 0; jj < blk.NumJ() + 1; jj++){    
      for ( int ii = 0; ii < blk.NumI(); ii++){      

	//get face location
	int loc = GetLoc1D(ii, jj, kk, blk.NumI(), blk.NumJ() + 1);

	if (turbFlag){
	  blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityJ[loc], temperatureJ[loc], tkeJ[loc], omegaJ[loc]);
	}
	else{
	  vector3d<double> dum1, dum2;
	  blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityJ[loc], temperatureJ[loc], dum1, dum2);
	}

      }
    }
  }


  //loop over k-faces and calculate gradients
  //loop over all physical faces
  for ( int kk = 0; kk < blk.NumK() + 1; kk++){   
    for ( int jj = 0; jj < blk.NumJ(); jj++){    
      for ( int ii = 0; ii < blk.NumI(); ii++){      

	//get face location
	int loc = GetLoc1D(ii, jj, kk, blk.NumI(), blk.NumJ());

	if (turbFlag){
	  blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityK[loc], temperatureK[loc], tkeK[loc], omegaK[loc]);
	}
	else{
	  vector3d<double> dum1, dum2;
	  blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(), kk + blk.NumGhosts(), eos, turbFlag, velocityK[loc], temperatureK[loc], dum1, dum2);
	}

      }
    }
  }

}

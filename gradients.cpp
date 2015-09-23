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

#include <iostream>
#include <vector>
#include <string>
#include "gradients.hpp"
#include "procBlock.hpp"
#include "eos.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::string;

// constructor
gradients::gradients() {
  // default to 1 faces
  velocityI_ = multiArray3d<tensor<double> >(1, 1, 1);
  velocityJ_ = multiArray3d<tensor<double> >(1, 1, 1);
  velocityK_ = multiArray3d<tensor<double> >(1, 1, 1);

  temperatureI_ = multiArray3d<vector3d<double> >(1, 1, 1);
  temperatureJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
  temperatureK_ = multiArray3d<vector3d<double> >(1, 1, 1);

  tkeI_ = multiArray3d<vector3d<double> >(1, 1, 1);
  tkeJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
  tkeK_ = multiArray3d<vector3d<double> >(1, 1, 1);

  omegaI_ = multiArray3d<vector3d<double> >(1, 1, 1);
  omegaJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
  omegaK_ = multiArray3d<vector3d<double> >(1, 1, 1);
}

// alternate constructor to calculate gradients from state variables
gradients::gradients(const bool &turbFlag, const procBlock &blk,
                     const idealGas &eos) {
  // turbFlag -- flag to determine if simulation is turbulent and therefore
  // requires tke and omega gradients
  // blk -- procBlock to calculate gradients for
  // eos -- equation of state

  // size vectors for input procBlock
  velocityI_ = multiArray3d<tensor<double> >(blk.NumI() + 1, blk.NumJ(),
                                             blk.NumK());
  velocityJ_ = multiArray3d<tensor<double> >(blk.NumI(), blk.NumJ() + 1,
                                             blk.NumK());
  velocityK_ = multiArray3d<tensor<double> >(blk.NumI(), blk.NumJ(),
                                             blk.NumK() + 1);

  temperatureI_ = multiArray3d<vector3d<double> >(blk.NumI() + 1,
                                                  blk.NumJ(), blk.NumK());
  temperatureJ_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ() + 1,
                                                  blk.NumK());
  temperatureK_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ(),
                                                  blk.NumK() + 1);

  // if flow is not turbulent, don't need tke and omega gradients so they can be
  // set to a minimum size
  if (turbFlag) {
    tkeI_ = multiArray3d<vector3d<double> >(blk.NumI() + 1, blk.NumJ(),
                                            blk.NumK());
    tkeJ_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ() + 1,
                                            blk.NumK());
    tkeK_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ(),
                                            blk.NumK() + 1);

    omegaI_ = multiArray3d<vector3d<double> >(blk.NumI() + 1, blk.NumJ(),
                                              blk.NumK());
    omegaJ_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ() + 1,
                                              blk.NumK());
    omegaK_ = multiArray3d<vector3d<double> >(blk.NumI(), blk.NumJ(),
                                              blk.NumK() + 1);
  } else {
    tkeI_ = multiArray3d<vector3d<double> >(1, 1, 1);
    tkeJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
    tkeK_ = multiArray3d<vector3d<double> >(1, 1, 1);

    omegaI_ = multiArray3d<vector3d<double> >(1, 1, 1);
    omegaJ_ = multiArray3d<vector3d<double> >(1, 1, 1);
    omegaK_ = multiArray3d<vector3d<double> >(1, 1, 1);
  }

  // loop over i-faces and calculate gradients
  // loop over all physical faces
  for (int kk = 0; kk < blk.NumK(); kk++) {
    for (int jj = 0; jj < blk.NumJ(); jj++) {
      for (int ii = 0; ii < blk.NumI() + 1; ii++) {
        if (turbFlag) {
          blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityI_(ii, jj, kk), temperatureI_(ii, jj, kk),
                         tkeI_(ii, jj, kk), omegaI_(ii, jj, kk));
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityI_(ii, jj, kk), temperatureI_(ii, jj, kk),
                         dum1, dum2);
        }
      }
    }
  }

  // loop over j-faces and calculate gradients
  // loop over all physical faces
  for (int kk = 0; kk < blk.NumK(); kk++) {
    for (int jj = 0; jj < blk.NumJ() + 1; jj++) {
      for (int ii = 0; ii < blk.NumI(); ii++) {
        // get face location
        if (turbFlag) {
          blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityJ_(ii, jj, kk), temperatureJ_(ii, jj, kk),
                         tkeJ_(ii, jj, kk), omegaJ_(ii, jj, kk));
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityJ_(ii, jj, kk), temperatureJ_(ii, jj, kk),
                         dum1, dum2);
        }
      }
    }
  }

  // loop over k-faces and calculate gradients
  // loop over all physical faces
  for (int kk = 0; kk < blk.NumK() + 1; kk++) {
    for (int jj = 0; jj < blk.NumJ(); jj++) {
      for (int ii = 0; ii < blk.NumI(); ii++) {
        if (turbFlag) {
          blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityK_(ii, jj, kk), temperatureK_(ii, jj, kk),
                         tkeK_(ii, jj, kk), omegaK_(ii, jj, kk));
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag,
                         velocityK_(ii, jj, kk), temperatureK_(ii, jj, kk),
                         dum1, dum2);
        }
      }
    }
  }
}

// member function to calculate velocity gradient at cell center
tensor<double> gradients::VelGradCell(const int &ii, const int &jj,
                                      const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (velocityI_(ii, jj, kk) + velocityI_(ii + 1, jj, kk) +
                        velocityJ_(ii, jj, kk) + velocityJ_(ii, jj + 1, kk) +
                        velocityK_(ii, jj, kk) + velocityK_(ii, jj, kk + 1));
}

// member function to calculate temperature gradient at cell center
vector3d<double> gradients::TempGradCell(const int &ii, const int &jj,
                                         const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (temperatureI_(ii, jj, kk) +
                        temperatureI_(ii + 1, jj, kk) +
                        temperatureJ_(ii, jj, kk) +
                        temperatureJ_(ii, jj + 1, kk) +
                        temperatureK_(ii, jj, kk) +
                        temperatureK_(ii, jj, kk + 1));
}

// member function to calculate tke gradient at cell center
vector3d<double> gradients::TkeGradCell(const int &ii, const int &jj,
                                        const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (tkeI_(ii, jj, kk) + tkeI_(ii + 1, jj, kk) +
                        tkeJ_(ii, jj, kk) + tkeJ_(ii, jj + 1, kk) +
                        tkeK_(ii, jj, kk) + tkeK_(ii, jj, kk + 1));
}

// member function to calculate omega gradient at cell center
vector3d<double> gradients::OmegaGradCell(const int &ii, const int &jj,
                                          const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (omegaI_(ii, jj, kk) + omegaI_(ii + 1, jj, kk) +
                        omegaJ_(ii, jj, kk) + omegaJ_(ii, jj + 1, kk) +
                        omegaK_(ii, jj, kk) + omegaK_(ii, jj, kk + 1));
}

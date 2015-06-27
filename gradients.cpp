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

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

// constructor
gradients::gradients() {
  // default to 2 faces
  velocityI_.resize(2);
  velocityJ_.resize(2);
  velocityK_.resize(2);

  temperatureI_.resize(2);
  temperatureJ_.resize(2);
  temperatureK_.resize(2);

  tkeI_.resize(2);
  tkeJ_.resize(2);
  tkeK_.resize(2);

  omegaI_.resize(2);
  omegaJ_.resize(2);
  omegaK_.resize(2);

  imax_ = 1;
  jmax_ = 1;
  kmax_ = 1;
}

// alternate constructor to calculate gradients from state variables
gradients::gradients(const bool &turbFlag, const procBlock &blk,
                     const idealGas &eos) {
  // turbFlag -- flag to determine if simulation is turbulent and therefore
  // requires tke and omega gradients
  // blk -- procBlock to calculate gradients for
  // eos -- equation of state

  // set dimensions (cell-size)
  imax_ = blk.NumI();
  jmax_ = blk.NumJ();
  kmax_ = blk.NumK();

  // size vectors for input procBlock
  velocityI_.resize((blk.NumI() + 1) * blk.NumJ() * blk.NumK());
  velocityJ_.resize(blk.NumI() * (blk.NumJ() + 1) * blk.NumK());
  velocityK_.resize(blk.NumI() * blk.NumJ() * (blk.NumK() + 1));

  temperatureI_.resize((blk.NumI() + 1) * blk.NumJ() * blk.NumK());
  temperatureJ_.resize(blk.NumI() * (blk.NumJ() + 1) * blk.NumK());
  temperatureK_.resize(blk.NumI() * blk.NumJ() * (blk.NumK() + 1));

  // if flow is not turbulent, don't need tke and omega gradients so they can be
  // set to a minimum size
  if (turbFlag) {
    tkeI_.resize((blk.NumI() + 1) * blk.NumJ() * blk.NumK());
    tkeJ_.resize(blk.NumI() * (blk.NumJ() + 1) * blk.NumK());
    tkeK_.resize(blk.NumI() * blk.NumJ() * (blk.NumK() + 1));

    omegaI_.resize((blk.NumI() + 1) * blk.NumJ() * blk.NumK());
    omegaJ_.resize(blk.NumI() * (blk.NumJ() + 1) * blk.NumK());
    omegaK_.resize(blk.NumI() * blk.NumJ() * (blk.NumK() + 1));
  } else {
    tkeI_.resize(1);
    tkeJ_.resize(1);
    tkeK_.resize(1);

    omegaI_.resize(1);
    omegaJ_.resize(1);
    omegaK_.resize(1);
  }

  // loop over i-faces and calculate gradients
  // loop over all physical faces
  for (int kk = 0; kk < blk.NumK(); kk++) {
    for (int jj = 0; jj < blk.NumJ(); jj++) {
      for (int ii = 0; ii < blk.NumI() + 1; ii++) {
        // get face location
        int loc = GetLoc1D(ii, jj, kk, blk.NumI() + 1, blk.NumJ());

        if (turbFlag) {
          blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityI_[loc],
                         temperatureI_[loc], tkeI_[loc], omegaI_[loc]);
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsI(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityI_[loc],
                         temperatureI_[loc], dum1, dum2);
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
        int loc = GetLoc1D(ii, jj, kk, blk.NumI(), blk.NumJ() + 1);

        if (turbFlag) {
          blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityJ_[loc],
                         temperatureJ_[loc], tkeJ_[loc], omegaJ_[loc]);
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsJ(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityJ_[loc],
                         temperatureJ_[loc], dum1, dum2);
        }
      }
    }
  }

  // loop over k-faces and calculate gradients
  // loop over all physical faces
  for (int kk = 0; kk < blk.NumK() + 1; kk++) {
    for (int jj = 0; jj < blk.NumJ(); jj++) {
      for (int ii = 0; ii < blk.NumI(); ii++) {
        // get face location
        int loc = GetLoc1D(ii, jj, kk, blk.NumI(), blk.NumJ());

        if (turbFlag) {
          blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityK_[loc],
                         temperatureK_[loc], tkeK_[loc], omegaK_[loc]);
        } else {
          vector3d<double> dum1, dum2;
          blk.CalcGradsK(ii + blk.NumGhosts(), jj + blk.NumGhosts(),
                         kk + blk.NumGhosts(), eos, turbFlag, velocityK_[loc],
                         temperatureK_[loc], dum1, dum2);
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

  return (1.0 / 6.0) * (velocityI_[GetUpperFaceI(ii, jj, kk, imax_, jmax_)] +
                        velocityI_[GetLowerFaceI(ii, jj, kk, imax_, jmax_)] +
                        velocityJ_[GetUpperFaceJ(ii, jj, kk, imax_, jmax_)] +
                        velocityJ_[GetLowerFaceJ(ii, jj, kk, imax_, jmax_)] +
                        velocityK_[GetUpperFaceK(ii, jj, kk, imax_, jmax_)] +
                        velocityK_[GetLowerFaceK(ii, jj, kk, imax_, jmax_)]);
}

// member function to calculate temperature gradient at cell center
vector3d<double> gradients::TempGradCell(const int &ii, const int &jj,
                                         const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (temperatureI_[GetUpperFaceI(ii, jj, kk, imax_, jmax_)] +
                        temperatureI_[GetLowerFaceI(ii, jj, kk, imax_, jmax_)] +
                        temperatureJ_[GetUpperFaceJ(ii, jj, kk, imax_, jmax_)] +
                        temperatureJ_[GetLowerFaceJ(ii, jj, kk, imax_, jmax_)] +
                        temperatureK_[GetUpperFaceK(ii, jj, kk, imax_, jmax_)] +
                        temperatureK_[GetLowerFaceK(ii, jj, kk, imax_, jmax_)]);
}

// member function to calculate tke gradient at cell center
vector3d<double> gradients::TkeGradCell(const int &ii, const int &jj,
                                        const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (tkeI_[GetUpperFaceI(ii, jj, kk, imax_, jmax_)] +
                        tkeI_[GetLowerFaceI(ii, jj, kk, imax_, jmax_)] +
                        tkeJ_[GetUpperFaceJ(ii, jj, kk, imax_, jmax_)] +
                        tkeJ_[GetLowerFaceJ(ii, jj, kk, imax_, jmax_)] +
                        tkeK_[GetUpperFaceK(ii, jj, kk, imax_, jmax_)] +
                        tkeK_[GetLowerFaceK(ii, jj, kk, imax_, jmax_)]);
}

// member function to calculate omega gradient at cell center
vector3d<double> gradients::OmegaGradCell(const int &ii, const int &jj,
                                          const int &kk) const {
  // ii - i dimension of cell to get gradient at
  // jj - j dimension of cell to get gradient at
  // kk - k dimension of cell to get gradient at

  return (1.0 / 6.0) * (omegaI_[GetUpperFaceI(ii, jj, kk, imax_, jmax_)] +
                        omegaI_[GetLowerFaceI(ii, jj, kk, imax_, jmax_)] +
                        omegaJ_[GetUpperFaceJ(ii, jj, kk, imax_, jmax_)] +
                        omegaJ_[GetLowerFaceJ(ii, jj, kk, imax_, jmax_)] +
                        omegaK_[GetUpperFaceK(ii, jj, kk, imax_, jmax_)] +
                        omegaK_[GetLowerFaceK(ii, jj, kk, imax_, jmax_)]);
}

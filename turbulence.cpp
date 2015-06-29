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

#include <cstdlib>    // exit()
#include <iostream>   // cout
#include <algorithm>  // max
#include "turbulence.hpp"
#include "primVars.hpp"  // primVars

using std::cout;
using std::endl;
using std::cerr;

// member functions for k-omega Wilcox 2006 model
// member function to return the eddy viscosity calculated with the stress
// limiter
double turbKWWilcox::EddyVisc(const primVars &state,
                              const tensor<double> &vGrad) const {
  return state.Rho() * state.Tke() / (*this).OmegaTilda(state, vGrad);
}

// member functionto return the eddy viscosity calculated without the stress
// limiter
double turbKWWilcox::EddyViscNoLim(const primVars &state) const {
  return state.Rho() * state.Tke() / state.Omega();
}

double turbKWWilcox::SigmaD(const vector3d<double> &kGrad,
                            const vector3d<double> &wGrad) const {
  // kGrad -- tke gradient
  // wGrad -- omega gradient

  if (kGrad.DotProd(wGrad) <= 0.0) {
    return 0.0;
  } else {
    return (*this).SigmaD0();
  }
}

double turbKWWilcox::Beta(const primVars &state,
                          const tensor<double> &velGrad) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  return (*this).Beta0() * (*this).FBeta(state, velGrad);
}

double turbKWWilcox::FBeta(const primVars &state,
                           const tensor<double> &velGrad) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  double xw = (*this).Xw(state, velGrad);
  return (1.0 + 85.0 * xw) / (1.0 + 100.0 * xw);
}

double turbKWWilcox::Xw(const primVars &state,
                        const tensor<double> &velGrad) const {
  // state -- primative variables
  // velGrad -- velocity gradient

  tensor<double> vorticity = 0.5 * (velGrad - velGrad.Transpose());

  // using DoubleDotTrans for speed
  // both tensors are symmetric so result is the same
  // vorticity is asymmetric but vorticity * vorticity is symmetric
  return fabs( (vorticity * vorticity).DoubleDotTrans((*this).StrainKI(velGrad))
               / pow((*this).BetaStar() * state.Omega(), 3.0) );
}

tensor<double> turbKWWilcox::StrainKI(const tensor<double> &velGrad) const {
  // velGrad -- velocity gradient
  tensor<double> I;
  I.Identity();
  return 0.5 * (velGrad + velGrad.Transpose()) - 0.5 * velGrad.Trace() * I;
}

double turbKWWilcox::OmegaTilda(const primVars &state,
                                const tensor<double> &velGrad) const {
  // state -- primative variables
  // velGrad -- velocity gradient
  tensor<double> I;
  I.Identity();
  tensor<double> sHat = 0.5 * (velGrad + velGrad.Transpose()) - 1.0 / 3.0 *
      velGrad.Trace() * I;

  // using DoubleDotTrans instead of DoubleDot for speed
  // since tensors are symmetric, result is the same
  return std::max(state.Omega(), (*this).CLim() *
                  sqrt(2.0 * sHat.DoubleDotTrans(sHat) / (*this).BetaStar()));
}

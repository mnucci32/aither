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
#include <cmath>
#include "source.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;

// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const source &src) {
  for (int ii = 0; ii < NUMVARS; ii++) {
    os << src.data_[ii];
    if (ii != NUMVARS - 1) {
      os << ", ";
    }
  }
  return os;
}

// operator overload for addition
source source::operator+(const source &src2) const {
  source src1 = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] += src2.data_[ii];
  }
  return src1;
}

// operator overload for addition with a scalar
source operator+(const double &scalar, const source &src2) {
  source src1;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] = src2.data_[ii] + scalar;
  }
  return src1;
}

// operator overload for subtraction
source source::operator-(const source &src2) const {
  source src1 = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] -= src2.data_[ii];
  }
  return src1;
}

// operator overload for subtraction with a scalar
source operator-(const double &scalar, const source &src2) {
  source src1;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] = scalar - src2.data_[ii];
  }
  return src1;
}

// operator overload for elementwise multiplication
source source::operator*(const source &src2) const {
  source src1 = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] *= src2.data_[ii];
  }
  return src1;
}

// member function for scalar multiplication
source source::operator*(const double &scalar) {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] *= scalar;
  }
  return temp;
}

// member function for scalar addition
source source::operator+(const double &scalar) {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] += scalar;
  }
  return temp;
}

// member function for scalar subtraction
source source::operator-(const double &scalar) {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] -= scalar;
  }
  return temp;
}

// member function for scalar division
source source::operator/(const double &scalar) {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] /= scalar;
  }
  return temp;
}

// operator overload for multiplication with a scalar
source operator*(const double &scalar, const source &src2) {
  source src1;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] = src2.data_[ii] * scalar;
  }
  return src1;
}

// operator overload for elementwise division
source source::operator/(const source &src2) const {
  source src1 = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] /= src2.data_[ii];
  }
  return src1;
}

// operator overload for division with a scalar
source operator/(const double &scalar, const source &src2) {
  source src1;
  for (int ii = 0; ii < NUMVARS; ii++) {
    src1.data_[ii] = scalar / src2.data_[ii];
  }
  return src1;
}

// Member function to calculate the source terms for the turbulence equations
void source::CalcTurbSrc(const turbModel *turb, const primVars &state,
                         const gradients &grads, const sutherland &suth,
                         const idealGas &eqnState, const int &ii, const int &jj,
                         const int &kk, const input &inp) {
  // turb -- turbulence model
  // state -- primative variables
  // grads -- gradients

  // calculate wall shear stress (double dot) velocity gradient
  // wall shear stress
  double mut = turb->EddyVisc(state) * suth.NondimScaling();
  double lambda = suth.Lambda(mut);

  tensor<double> tau =
      lambda * grads.VelGradCell(ii, jj, kk).Trace() + mut *
      (grads.VelGradCell(ii, jj, kk) +
       grads.VelGradCell(ii, jj, kk).Transpose());

  // Using DoubleDotTrans instead of DoubleDot for speed
  // Since tensors are symmetric result is the same
  data_[5] = state.Rho() * tau.DoubleDotTrans(grads.VelGradCell(ii, jj, kk))
      - turb->Eqn1DissipationCoeff() * state.Rho() * state.Omega() *
      state.Tke();

  data_[6] = turb->Eqn2ProductionCoeff() * state.Omega() / state.Tke() *
      state.Rho() * tau.DoubleDotTrans(grads.VelGradCell(ii, jj, kk))
      - turb->Eqn2DissipationCoeff(state, grads.VelGradCell(ii, jj, kk))
      * state.Rho() * state.Omega() * state.Omega();
}

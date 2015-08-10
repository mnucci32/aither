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
#include "turbulence.hpp"
#include "primVars.hpp"
#include "gradients.hpp"

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
source source::operator*(const double &scalar) const {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] *= scalar;
  }
  return temp;
}

// member function for scalar addition
source source::operator+(const double &scalar) const {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] += scalar;
  }
  return temp;
}

// member function for scalar subtraction
source source::operator-(const double &scalar) const {
  source temp = *this;
  for (int ii = 0; ii < NUMVARS; ii++) {
    temp.data_[ii] -= scalar;
  }
  return temp;
}

// member function for scalar division
source source::operator/(const double &scalar) const {
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
                         const int &ii, const int &jj, const int &kk) {
  // turb -- turbulence model
  // state -- primative variables
  // grads -- gradients
  // suth -- sutherland's law for viscosity
  // eqnState -- equation of state
  // ii -- cell i-location to calculate source terms at
  // jj -- cell j-location to calculate source terms at
  // kk -- cell k-location to calculate source terms at

  // get cell gradients
  tensor<double> vGrad = grads.VelGradCell(ii, jj, kk);
  vector3d<double> kGrad = grads.TkeGradCell(ii, jj, kk);
  vector3d<double> wGrad = grads.OmegaGradCell(ii, jj, kk);

  // calculate turbulent source terms
  data_[5] = suth.NondimScaling() * turb->Production1(state, vGrad, suth)
      - suth.InvNondimScaling() * turb->Destruction1(state);

  data_[6] = suth.NondimScaling() * turb->Production2(state, vGrad, suth)
      - suth.InvNondimScaling() * turb->Destruction2(state, vGrad, suth)
      + suth.NondimScaling() * turb->CrossDiff2(state, kGrad, wGrad);
}

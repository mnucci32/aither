/*  This file is part of aither.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    Aither is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Aither is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

#include <cmath>  // sqrt
#include <string>
#include <memory>
#include "fluxJacobian.hpp"
#include "turbulence.hpp"  // turbModel
#include "input.hpp"       // input
#include "primVars.hpp"    // primVars
#include "genArray.hpp"    // genArray

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std:: unique_ptr;

// constructor
fluxJacobian::fluxJacobian(const primVars &state,
                           const unitVec3dMag<double> &fAreaL,
                           const unitVec3dMag<double> &fAreaR,
                           const idealGas &eos, const sutherland &suth,
                           const double &vol,
                           const unique_ptr<turbModel> &turb,
                           const input &inp, const bool &mainDiagonal) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model
  // inp -- input variables
  // mainDiagonal -- flag to determine if jacobian is for main diagonal of
  //                 implicit matrix

  // initialize jacobians
  flowJacobian_ = 0.0;
  turbJacobian_ = 0.0;

  const auto isTurbulent = inp.IsTurbulent();

  this->AddInviscidJacobian(state, fAreaL, fAreaR, eos, turb, isTurbulent);

  // if viscous, add viscous jacobians
  if (inp.IsViscous()) {
    this->AddViscousJacobian(state, fAreaL, fAreaR, eos, suth, vol,
                             turb, isTurbulent);
  }

  // source term is only added to main diagonal,
  // and only for turbulence equations
  if (mainDiagonal && isTurbulent) {
    this->AddTurbSourceJacobian(state, suth, vol, turb);
  }
}

// member functions
// member function to add inviscid jacobians
void fluxJacobian::AddInviscidJacobian(const primVars &state,
                                       const unitVec3dMag<double> &fAreaL,
                                       const unitVec3dMag<double> &fAreaR,
                                       const idealGas &eos,
                                       const unique_ptr<turbModel> &turb,
                                       const bool &isTurbulent) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // turb -- turbulence model
  // isTurbulent -- flag to determine if simulation is turbulent

  flowJacobian_ += state.CellSpectralRadius(fAreaL, fAreaR, eos);

  if (isTurbulent) {
    turbJacobian_ += turb->InviscidSpecRad(state, fAreaL, fAreaR);
  }
}

// member function to add viscous jacobians
void fluxJacobian::AddViscousJacobian(const primVars &state,
                                      const unitVec3dMag<double> &fAreaL,
                                      const unitVec3dMag<double> &fAreaR,
                                      const idealGas &eos,
                                      const sutherland &suth,
                                      const double &vol,
                                      const unique_ptr<turbModel> &turb,
                                      const bool &isTurbulent) {
  // state -- primative variables
  // fAreaL -- face area for left face
  // fAreaR -- face area for right face
  // eos -- equation of state
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model
  // isTurbulent -- flag to determine if simulation is turbulent

  // factor of 2 because viscous spectral radius is not halved (Blazek 6.53)
  flowJacobian_ += 2.0 * state.ViscCellSpectralRadius(fAreaL, fAreaR, eos, suth,
                                                      vol, turb);

  if (isTurbulent) {
    // factor of 2 because viscous spectral radius is not halved (Blazek 6.53)
    turbJacobian_ += 2.0 * turb->ViscSpecRad(state, fAreaL, fAreaR, eos, suth,
                                             vol);
  }
}

// member function to add source jacobians
// this should only be used for the main diagonal
void fluxJacobian::AddTurbSourceJacobian(const primVars &state,
                                         const sutherland &suth,
                                         const double &vol,
                                         const unique_ptr<turbModel> &turb) {
  // state -- primative variables
  // suth -- sutherland's law for viscosity
  // vol -- cell volume
  // turb -- turbulence model

  turbJacobian_ -= turb->SrcSpecRad(state, suth) * vol;
}

// member function to multiply the flux jacobians with a genArray
genArray fluxJacobian::ArrayMult(genArray arr) const {
  arr[0] *= flowJacobian_;
  arr[1] *= flowJacobian_;
  arr[2] *= flowJacobian_;
  arr[3] *= flowJacobian_;
  arr[4] *= flowJacobian_;

  arr[5] *= turbJacobian_;
  arr[6] *= turbJacobian_;

  return arr;
}

// member function to take the inverse of a flux jacobian
fluxJacobian fluxJacobian::Inverse(const bool &isTurbulent) const {
  auto inv = *this;
  inv.flowJacobian_ = 1.0 / inv.flowJacobian_;

  if (isTurbulent) {
    inv.turbJacobian_ = 1.0 / inv.turbJacobian_;
  }
  return inv;
}

// non-member functions
// ----------------------------------------------------------------------------
// operator overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, fluxJacobian &jacobian) {
  os << jacobian.FlowJacobian() << endl;
  os << jacobian.TurbulenceJacobian() << endl;
  return os;
}


/*  This file is part of aither.
    Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef PRIMITIVEHEADERDEF
#define PRIMITIVEHEADERDEF

/* This header contains the primitive class.

The primitive class stores the primitive variables for the Euler and
Navier-Stokes equations [rho, u, v, w, P]. It contains several member functions
to manipulate these primitive varibles. It also contains member functions to
extrapolate the primitive variables from the cell center to the cell
face using constant and MUSCL reconstruction. It is has a member function to
supply a ghost state given a boundary condition and boundary cell.  */

#include <fstream>
#include <iostream>
#include <vector>                  // vector
#include <string>                  // string
#include <memory>                  // unique_ptr
#include "vector3d.hpp"            // vector3d
#include "tensor.hpp"              // tensor
#include "eos.hpp"                 // equation of state
#include "transport.hpp"           // transport model
#include "thermodynamic.hpp"       // thermodynamic model
#include "multiArray3d.hpp"        // multiArray3d
#include "varArray.hpp"
#include "conserved.hpp"

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;
using std::unique_ptr;

// forward class declarations
class input;
class turbModel;
struct wallVars;

class primitive : public varArray {
 public:
  // constructors
  primitive(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}
  primitive(const int &numEqns, const int &numSpecies, const double &val)
      : varArray(numEqns, numSpecies, val) {}
  primitive(const conserved &, const unique_ptr<eos> &,
            const unique_ptr<thermodynamic> &, const unique_ptr<turbModel> &);
  primitive(const vector<const double>::iterator &b,
            const vector<const double>::iterator &e, const int &numSpecies)
      : varArray(b, e, numSpecies) {}

  // move constructor and assignment operator
  primitive(primitive&&) noexcept = default;
  primitive& operator=(primitive&&) noexcept = default;

  // copy constructor and assignment operator
  primitive(const primitive&) = default;
  primitive& operator=(const primitive&) = default;

  // member functions
  const double & RhoN(const int &ii) const { return this->SpeciesN(ii); }
  double Rho() const { return this->SpeciesSum(); }
  double MassFractionN(const int &ii) const {
    return this->RhoN(ii) / this->Rho();
  }
  const double & U() const { return this->MomentumX(); }
  const double & V() const { return this->MomentumY(); }
  const double & W() const { return this->MomentumZ(); }
  const double & P() const { return this->varArray::Energy(); }
  const double & Tke() const { return this->TurbulenceN(0); }
  const double & Omega() const { return this->TurbulenceN(1); }
  const double & TurbN(const int &ii) const { return this->TurbulenceN(ii); }

  void NondimensionalInitialize(const unique_ptr<eos> &, const input &,
                                const unique_ptr<transport> &, const int &,
                                const unique_ptr<turbModel> &);

  primitive Abs() const;

  vector3d<double> Velocity() const {
    return {this->U(), this->V(), this->W()};
  }

  double Energy(const unique_ptr<eos> &eqnState,
                const unique_ptr<thermodynamic> &thermo) const {
    const auto t = this->Temperature(eqnState);
    return eqnState->Energy(eqnState->SpecEnergy(thermo, t),
                            this->Velocity().Mag());
  }
  double Enthalpy(const unique_ptr<eos> &eqnState,
                  const unique_ptr<thermodynamic> &thermo) const {
    const auto t = this->Temperature(eqnState);
    return eqnState->Enthalpy(thermo, t, this->Velocity().Mag());
  }
  double Temperature(const unique_ptr<eos> &eqnState) const {
    return eqnState->Temperature(this->P(), this->Rho());
  }
  double SoS(const unique_ptr<thermodynamic> &thermo,
                    const unique_ptr<eos> &eqnState) const {
    return sqrt(thermo->Gamma(this->Temperature(eqnState)) * this->P() /
                this->Rho());
  }

  inline conserved ConsVars(const unique_ptr<eos> &,
                            const unique_ptr<thermodynamic> &) const;
  primitive UpdateWithConsVars(const unique_ptr<eos> &,
                               const unique_ptr<thermodynamic> &,
                               const conserved &,
                               const unique_ptr<turbModel> &) const;

  void ApplyFarfieldTurbBC(const vector3d<double> &, const double &,
                           const double &, const unique_ptr<transport> &,
                           const unique_ptr<eos> &,
                           const unique_ptr<turbModel> &);
  void LimitTurb(const unique_ptr<turbModel> &);

  // destructor
  ~primitive() noexcept {}
};

// function definitions
// member function to calculate conserved variables from primitive variables
conserved primitive::ConsVars(const unique_ptr<eos> &eqnState,
                            const unique_ptr<thermodynamic> &thermo) const {
  conserved cv(this->Size(), this->NumSpecies());
  for (auto ii = 0; ii < cv.NumSpecies(); ++ii) {
    cv[ii] = (*this)[ii];
  }
  const auto rho = this->Rho();
  cv[cv.MomentumXIndex()] = rho * this->U();
  cv[cv.MomentumYIndex()] = rho * this->V();
  cv[cv.MomentumZIndex()] = rho * this->W();
  cv[cv.EnergyIndex()] = rho * this->Energy(eqnState, thermo);
  for (auto ii = 0; ii < cv.NumTurbulence(); ++ii) {
    cv[cv.TurbulenceIndex() + ii] = rho * this->TurbN(ii);
  }
  return cv;
}

// function to return the state of the appropriate ghost cell
primitive GetGhostState(const primitive &, const string &,
                        const vector3d<double> &, const double &, const int &,
                        const input &, const int &, const unique_ptr<eos> &,
                        const unique_ptr<thermodynamic> &,
                        const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &, wallVars &, const int &,
                        const double & = 0.0, const primitive & = {0, 0},
                        const vector3d<double> & = {},
                        const tensor<double> & = {}, const double & = 0.0,
                        const double & = 0.0);

ostream &operator<<(ostream &os, const primitive &);

primitive RoeAveragedState(const primitive&, const primitive&);

#endif

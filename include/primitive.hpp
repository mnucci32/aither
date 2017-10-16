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
#include "arrayView.hpp"

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

class primitive : public varArray {
 public:
  // constructors
  primitive() : varArray() {}
  primitive(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}
  primitive(const int &numEqns, const int &numSpecies, const double &val)
      : varArray(numEqns, numSpecies, val) {}
  template <typename T,
            typename = std::enable_if_t<std::is_base_of<varArray, T>::value ||
                                        std::is_same<conservedView, T>::value>>
  primitive(const T &, const unique_ptr<eos> &,
            const unique_ptr<thermodynamic> &, const unique_ptr<turbModel> &);
  primitive(const vector<double>::const_iterator &b,
            const vector<double>::const_iterator &e, const int &numSpecies)
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
  primitive Squared() const {
    auto sq = (*this);
    return sq *= sq;
  }

  arrayView<primitive, double> GetView() const {
    return {this->begin(), this->end(), this->NumSpecies()};
  }

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
                               const varArrayView &,
                               const unique_ptr<turbModel> &) const;

  void ApplyFarfieldTurbBC(const vector3d<double> &, const double &,
                           const double &, const unique_ptr<transport> &,
                           const unique_ptr<eos> &,
                           const unique_ptr<turbModel> &);
  void LimitTurb(const unique_ptr<turbModel> &);

  // destructor
  ~primitive() noexcept {}
};

// ---------------------------------------------------------------------------
// constructors
template <typename T, typename TT>
primitive::primitive(const T &cons, const unique_ptr<eos> &eqnState,
                     const unique_ptr<thermodynamic> &thermo,
                     const unique_ptr<turbModel> &turb) {
  // cons -- array of conserved variables
  // eqnState -- equation of state
  // thermo -- thermodynamic model
  // turb -- turbulence model

  *this = primitive(cons.Size(), cons.NumSpecies());

  for (auto ii = 0; ii < this->NumSpecies(); ++ii) {
    (*this)[ii] = cons.SpeciesN(ii);
  }
  
  const auto rho = cons.SpeciesSum();
  (*this)[this->MomentumXIndex()] = cons.MomentumX() / rho;
  (*this)[this->MomentumYIndex()] = cons.MomentumY() / rho;
  (*this)[this->MomentumZIndex()] = cons.MomentumZ() / rho;
  
  const auto energy = cons.Energy() / rho;
  (*this)[this->EnergyIndex()] =
      eqnState->PressFromEnergy(thermo, rho, energy, this->Velocity().Mag());

  for (auto ii = 0; ii < this->NumTurbulence(); ++ii) {
    (*this)[this->TurbulenceIndex() + ii] = cons.TurbulenceN(ii) / rho;
  }

  // Adjust turbulence variables to be above minimum if necessary
  this->LimitTurb(turb);
}

// ---------------------------------------------------------------------------
// non member functions
// function to calculate conserved variables from primitive variables
template <typename T>
conserved PrimToCons(const T &state, const unique_ptr<eos> &eqnState,
                     const unique_ptr<thermodynamic> &thermo) {
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");
  conserved cv(state.Size(), state.NumSpecies());
  for (auto ii = 0; ii < cv.NumSpecies(); ++ii) {
    cv[ii] = state[ii];
  }
  const auto rho = state.Rho();
  cv[cv.MomentumXIndex()] = rho * state.U();
  cv[cv.MomentumYIndex()] = rho * state.V();
  cv[cv.MomentumZIndex()] = rho * state.W();
  cv[cv.EnergyIndex()] = rho * state.Energy(eqnState, thermo);
  for (auto ii = 0; ii < cv.NumTurbulence(); ++ii) {
    cv[cv.TurbulenceIndex() + ii] = rho * state.TurbN(ii);
  }
  return cv;
}

// function to take in a genArray of updates to the conservative
// variables, and update the primitive variables with it.
// this is used in the implicit solver
template <typename T>
primitive UpdatePrimWithCons(const T &state, const unique_ptr<eos> &eqnState,
                             const unique_ptr<thermodynamic> &thermo,
                             const varArrayView &du,
                             const unique_ptr<turbModel> &turb) {
  // eqnState -- equation of state
  // du -- updates to conservative variables
  // turb -- turbulence model
  static_assert(std::is_same<primitive, T>::value ||
                    std::is_same<primitiveView, T>::value,
                "T requires primitive or primativeView type");

  // convert primitive to conservative and update
  const auto consUpdate = state.ConsVars(eqnState, thermo) + du;
  return primitive(consUpdate, eqnState, thermo, turb);
}
                               


// ---------------------------------------------------------------------------
// member function to calculate conserved variables from primitive variables
conserved primitive::ConsVars(const unique_ptr<eos> &eqnState,
                            const unique_ptr<thermodynamic> &thermo) const {
  return PrimToCons((*this), eqnState, thermo);
}

ostream &operator<<(ostream &os, const primitive &);

// function to calculate the Roe averaged state
template <typename T1, typename T2>
primitive RoeAveragedState(const T1 &left, const T2 &right) {
  static_assert(std::is_same<primitive, T1>::value ||
                    std::is_same<primitiveView, T1>::value,
                "T1 requires primitive or primativeView type");
  static_assert(std::is_same<primitive, T2>::value ||
                    std::is_same<primitiveView, T2>::value,
                "T2 requires primitive or primativeView type");
                
  // compute Rho averaged quantities
  primitive rhoState(left.Size(), left.NumSpecies());
  // density ratio
  const auto denRatio = sqrt(right.Rho() / left.Rho());
  // Roe averaged density
  for (auto ii = 0; ii < rhoState.NumSpecies(); ++ii) {
    rhoState[ii] = left.RhoN(ii) * denRatio;
  }
  // Roe averaged velocities - u, v, w
  rhoState[rhoState.MomentumXIndex()] =
      (left.U() + denRatio * right.U()) / (1.0 + denRatio);
  rhoState[rhoState.MomentumYIndex()] =
      (left.V() + denRatio * right.V()) / (1.0 + denRatio);
  rhoState[rhoState.MomentumZIndex()] =
      (left.W() + denRatio * right.W()) / (1.0 + denRatio);

  // Roe averaged pressure
  rhoState[rhoState.EnergyIndex()] =
      (left.P() + denRatio * right.P()) / (1.0 + denRatio);

  // Roe averaged turbulence variables
  for (auto ii = 0; ii < rhoState.NumTurbulence(); ++ii) {
    rhoState[rhoState.TurbulenceIndex() + ii] =
        (left.TurbN(ii) + denRatio * right.TurbN(ii)) / (1.0 + denRatio);
  }

  return rhoState;
}


#endif

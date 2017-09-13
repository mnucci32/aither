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

#ifndef PRIMATIVEHEADERDEF
#define PRIMATIVEHEADERDEF

/* This header contains the primative class.

The primative class stores the primative variables for the Euler and
Navier-Stokes equations [rho, u, v, w, P]. It contains several member functions
to manipulate these primative varibles. It also contains member functions to
extrapolate the primative variables from the cell center to the cell
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

class primative : public varArray {
 public:
  // constructors
  primative(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}

  primative(const conserved &, const unique_ptr<eos> &,
           const unique_ptr<thermodynamic> &, const unique_ptr<turbModel> &);

  // move constructor and assignment operator
  primative(primative&&) noexcept = default;
  primative& operator=(primative&&) noexcept = default;

  // copy constructor and assignment operator
  primative(const primative&) = default;
  primative& operator=(const primative&) = default;

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

  primative Abs() const;

  inline vector3d<double> Velocity() const;

  inline double Energy(const unique_ptr<eos> &,
                       const unique_ptr<thermodynamic> &) const;
  inline double Enthalpy(const unique_ptr<eos> &,
                         const unique_ptr<thermodynamic> &) const;
  inline double Temperature(const unique_ptr<eos> &) const;
  inline double SoS(const unique_ptr<thermodynamic> &,
                    const unique_ptr<eos> &) const;

  inline conserved ConsVars(const unique_ptr<eos> &,
                            const unique_ptr<thermodynamic> &) const;
  primative UpdateWithConsVars(const unique_ptr<eos> &,
                               const unique_ptr<thermodynamic> &,
                               const conserved &,
                               const unique_ptr<turbModel> &) const;

  void ApplyFarfieldTurbBC(const vector3d<double> &, const double &,
                           const double &, const unique_ptr<transport> &,
                           const unique_ptr<eos> &,
                           const unique_ptr<turbModel> &);
  void LimitTurb(const unique_ptr<turbModel> &);

  double InvCellSpectralRadius(const unitVec3dMag<double> &,
                               const unitVec3dMag<double> &,
                               const unique_ptr<thermodynamic> &,
                               const unique_ptr<eos> &) const;
  double InvFaceSpectralRadius(const unitVec3dMag<double> &,
                               const unique_ptr<thermodynamic> &,
                               const unique_ptr<eos> &) const;

  double ViscCellSpectralRadius(const unitVec3dMag<double> &,
                                const unitVec3dMag<double> &,
                                const unique_ptr<thermodynamic> &,
                                const unique_ptr<eos> &,
                                const unique_ptr<transport> &, const double &,
                                const double &, const double &,
                                const unique_ptr<turbModel> &) const;
  double ViscFaceSpectralRadius(const unitVec3dMag<double> &,
                                const unique_ptr<thermodynamic> &,
                                const unique_ptr<eos> &,
                                const unique_ptr<transport> &, const double &,
                                const double &, const double &,
                                const unique_ptr<turbModel> &) const;

  double CellSpectralRadius(const unitVec3dMag<double> &,
                            const unitVec3dMag<double> &,
                            const unique_ptr<thermodynamic> &,
                            const unique_ptr<eos> &,
                            const unique_ptr<transport> &, const double &,
                            const double &, const double &,
                            const unique_ptr<turbModel> &, const bool &) const;
  double FaceSpectralRadius(const unitVec3dMag<double> &,
                            const unique_ptr<thermodynamic> &,
                            const unique_ptr<eos> &,
                            const unique_ptr<transport> &, const double &,
                            const double &, const double &,
                            const unique_ptr<turbModel> &, const bool &) const;

  // member function to calculate reconstruction of state variables from cell
  // center to cell face assuming value at cell center is constant over cell
  // volume; zeroth order reconstruction results in first order accuracy
  const primative &FaceReconConst() const { return *this; }

  // member function to calculate reconstruction of state variables from cell
  // center to cell face this function uses muscle extrapolation resulting in
  // higher order accuracy
  primative FaceReconMUSCL(const primative &, const primative &, const double &,
                          const string &, const double &,
                          const double &, const double &) const;

  // calculate face reconstruction using 5th order weno scheme
  primative FaceReconWENO(const primative &, const primative &, const primative &,
                         const primative &, const double &, const double &,
                         const double &, const double &, const double &,
                         const bool &) const;

  // member function to calculate Van Albada limiter function
  primative LimiterVanAlbada(const primative &) const;
  primative LimiterMinmod(const primative &, const primative &,
                         const double &) const;
  primative LimiterNone() const;

  // member function to return the state of the appropriate ghost cell
  primative GetGhostState(const string &, const vector3d<double> &,
                         const double &, const int &, const input &,
                         const int &, const unique_ptr<eos> &,
                         const unique_ptr<thermodynamic> &,
                         const unique_ptr<transport> &,
                         const unique_ptr<turbModel> &, wallVars &, const int &,
                         const double & = 0.0, const primative & = {0, 0},
                         const vector3d<double> & = {},
                         const tensor<double> & = {}, const double & = 0.0,
                         const double & = 0.0) const;

  // destructor
  ~primative() noexcept {}
};

// function definitions
// member function to calculate temperature from primative variables and
// equation of state
double primative::Temperature(const unique_ptr<eos> &eqnState) const {
  return eqnState->Temperature(this->P(), this->Rho());
}

// member function to calculate velocity from primative variables
vector3d<double> primative::Velocity() const {
  vector3d<double> vel(this->U(), this->V(), this->W());
  return vel;
}

// member function to calculate total energy from primative variables
double primative::Energy(const unique_ptr<eos> &eqnState,
                        const unique_ptr<thermodynamic> &thermo) const {
  const auto t = this->Temperature(eqnState);
  return eqnState->Energy(eqnState->SpecEnergy(thermo, t),
                          this->Velocity().Mag());
}

// member function to calculate speed of sound from primative variables
double primative::SoS(const unique_ptr<thermodynamic> &thermo,
                     const unique_ptr<eos> &eqnState) const {
  return sqrt(thermo->Gamma(this->Temperature(eqnState)) * this->P() /
              this->Rho());
}

// member function to calculate enthalpy from conserved variables and equation
// of state
double primative::Enthalpy(const unique_ptr<eos> &eqnState,
                          const unique_ptr<thermodynamic> &thermo) const {
  const auto t = this->Temperature(eqnState);
  return eqnState->Enthalpy(thermo, t, this->Velocity().Mag());
}

// member function to calculate conserved variables from primative variables
conserved primative::ConsVars(const unique_ptr<eos> &eqnState,
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

ostream &operator<<(ostream &os, const primative &);

primative RoeAveragedState(const primative&, const primative&);

#endif

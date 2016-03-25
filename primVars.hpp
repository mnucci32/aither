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

#ifndef PRIMVARSHEADERDEF  // only if the macro PRIMVARSHEADERDEF is not defined
                           // execute these lines of code
#define PRIMVARSHEADERDEF  // define the macro

/* This header contains the primVars class.

The primVars class stores the primative variables for the Euler and
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
#include "eos.hpp"                 // idealGas, sutherland
#include "multiArray3d.hpp"        // multiArray3d
#include "genArray.hpp"            // genArray
#include "macros.hpp"

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

class primVars {
  double data_[NUMVARS];  // primative variables at cell center

 public:
  // constructors
  primVars(const double &a, const double &b, const double &c, const double &d,
           const double &e, const double &f, const double &g)
      : data_{a, b, c, d, e, f, g} {}
  primVars(const double &a, const double &b, const double &c, const double &d,
           const double &e)
      : primVars(a, b, c, d, e, 0.0, 0.0) {}
  primVars() : primVars(0.0, 0.0, 0.0, 0.0, 0.0) {}
  explicit primVars(const double &a) : primVars(a, a, a, a, a, a, a) {}
  primVars(const double &r, const double &p, const vector3d<double> &v)
      : primVars(r, v.X(), v.Y(), v.Z(), p) {}
  primVars(const genArray &, const bool &, const idealGas &,
           const unique_ptr<turbModel> &);

  // move constructor and assignment operator
  primVars(primVars&&) noexcept = default;
  primVars& operator=(primVars&&) noexcept = default;

  // copy constructor and assignment operator
  primVars(const primVars&) = default;
  primVars& operator=(const primVars&) = default;

  // member functions
  double Rho() const { return data_[0]; }
  double U() const { return data_[1]; }
  double V() const { return data_[2]; }
  double W() const { return data_[3]; }
  double P() const { return data_[4]; }
  double Tke() const { return data_[5]; }
  double Omega() const { return data_[6]; }

  void NondimensionalInitialize(const idealGas&, const double&, const input&,
                                const sutherland&);
  bool IsZero() const;

  inline vector3d<double> Velocity() const;

  inline double Energy(const idealGas &) const;
  inline double Enthalpy(const idealGas &) const;
  inline double Temperature(const idealGas &) const;
  inline double SoS(const idealGas &) const;

  inline genArray ConsVars(const idealGas &) const;
  primVars UpdateWithConsVars(const idealGas &, const genArray &,
                              const unique_ptr<turbModel> &) const;

  void ApplyFarfieldTurbBC(const vector3d<double> &, const double &,
                           const double &, const sutherland &,
                           const idealGas &);
  void LimitTurb(const unique_ptr<turbModel> &);

  double CellSpectralRadius(const unitVec3dMag<double> &,
                            const unitVec3dMag<double> &,
                            const idealGas &) const;
  double ViscCellSpectralRadius(const unitVec3dMag<double> &,
                                const unitVec3dMag<double> &, const idealGas &,
                                const sutherland &, const double &,
                                const unique_ptr<turbModel> &) const;

  // operator overloads for addition and subtraction of states
  inline primVars & operator+=(const primVars &);
  inline primVars & operator-=(const primVars &);
  inline primVars & operator*=(const primVars &);
  inline primVars & operator/=(const primVars &);

  inline primVars & operator+=(const double &);
  inline primVars & operator-=(const double &);
  inline primVars & operator*=(const double &);
  inline primVars & operator/=(const double &);

  inline primVars operator+(const double &s) const {
    auto lhs = *this;
    return lhs += s;
  }
  inline primVars operator-(const double &s) const {
    auto lhs = *this;
    return lhs -= s;
  }
  inline primVars operator*(const double &s) const {
    auto lhs = *this;
    return lhs *= s;
  }
  inline primVars operator/(const double &s) const {
    auto lhs = *this;
    return lhs /= s;
  }

  friend inline const primVars operator-(const double &lhs, primVars rhs);
  friend inline const primVars operator/(const double &lhs, primVars rhs);

  // member function to calculate reconstruction of state variables from cell
  // center to cell face assuming value at cell center is constant over cell
  // volume; zeroth order reconstruction results in first order accuracy
  primVars FaceReconConst() const { return *this; }

  // member function to calculate reconstruction of state variables from cell
  // center to cell face this function uses muscle extrapolation resulting in
  // higher order accuracy
  primVars FaceReconMUSCL(const primVars &, const primVars &, const double &,
                          const string &, const double &,
                          const double &, const double &) const;

  // member function to calculate Van Albada limiter function
  primVars LimiterVanAlbada(const primVars &) const;
  primVars LimiterMinmod(const primVars &, const primVars &,
                         const double &) const;
  primVars LimiterNone() const;

  // member function to return the state of the appropriate ghost cell
  primVars GetGhostState(const string &, const vector3d<double> &,
                         const double &, const string &, const input &,
                         const idealGas &, const sutherland &,
                         const unique_ptr<turbModel> &, const int = 1) const;

  // destructor
  ~primVars() noexcept {}
};

// function definitions
multiArray3d<primVars> GetGhostStates(
    const multiArray3d<primVars> &, const string &,
    const multiArray3d<unitVec3dMag<double>> &, const multiArray3d<double> &,
    const string &, const input &, const idealGas &, const sutherland &,
    const unique_ptr<turbModel> &, const int = 1);

// member function to calculate temperature from conserved variables and
// equation of state
double primVars::Temperature(const idealGas &eqnState) const {
  return eqnState.Temperature(data_[4], data_[0]);
}

// member function to calculate velocity from conserved variables
vector3d<double> primVars::Velocity() const {
  vector3d<double> vel(data_[1], data_[2], data_[3]);
  return vel;
}

// member function to calculate total enthalpy from conserved variables
double primVars::Energy(const idealGas &eqnState) const {
  return eqnState.Energy(eqnState.SpecEnergy(data_[4], data_[0]),
                            (*this).Velocity().Mag());
}

// member function to calculate speed of sound from primative varialbes
double primVars::SoS(const idealGas &eqnState) const {
  return sqrt(eqnState.Gamma() * data_[4] / data_[0]);
}

// member function to calculate enthalpy from conserved variables and equation
// of state
double primVars::Enthalpy(const idealGas &eqnState) const {
  return eqnState.Enthalpy((*this).Energy(eqnState), data_[4], data_[0]);
}

// member function to calculate conserved variables from primative variables
genArray primVars::ConsVars(const idealGas &eqnState) const {
  genArray cv(data_[0],
              data_[0] * data_[1],
              data_[0] * data_[2],
              data_[0] * data_[3],
              data_[0] * (*this).Energy(eqnState),
              data_[0] * data_[5],
              data_[0] * data_[6]);
  return cv;
}

// operator overload for addition
primVars & primVars::operator+=(const primVars &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] += arr.data_[rr];
  }
  return *this;
}

// operator overload for subtraction with a scalar
primVars & primVars::operator-=(const primVars &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] -= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise multiplication
primVars & primVars::operator*=(const primVars &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] *= arr.data_[rr];
  }
  return *this;
}

// operator overload for elementwise division
primVars & primVars::operator/=(const primVars &arr) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    data_[rr] /= arr.data_[rr];
  }
  return *this;
}

inline const primVars operator+(primVars lhs, const primVars &rhs) {
  return lhs += rhs;
}

inline const primVars operator-(primVars lhs, const primVars &rhs) {
  return lhs -= rhs;
}

inline const primVars operator*(primVars lhs, const primVars &rhs) {
  return lhs *= rhs;
}

inline const primVars operator/(primVars lhs, const primVars &rhs) {
  return lhs /= rhs;
}

// operator overloads for double -------------------------------------
// operator overload for addition
primVars & primVars::operator+=(const double &scalar) {
  for (auto &val : data_) {
    val += scalar;
  }
  return *this;
}

// operator overload for subtraction with a scalar
primVars & primVars::operator-=(const double &scalar) {
  for (auto &val : data_) {
    val -= scalar;
  }
  return *this;
}

// operator overload for elementwise multiplication
primVars & primVars::operator*=(const double &scalar) {
  for (auto &val : data_) {
    val *= scalar;
  }
  return *this;
}

// operator overload for elementwise division
primVars & primVars::operator/=(const double &scalar) {
  for (auto &val : data_) {
    val /= scalar;
  }
  return *this;
}

inline const primVars operator+(const double &lhs, primVars rhs) {
  return rhs += lhs;
}

inline const primVars operator-(const double &lhs, primVars rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs - rhs.data_[rr];
  }
  return rhs;
}

inline const primVars operator*(const double &lhs, primVars rhs) {
  return rhs *= lhs;
}

inline const primVars operator/(const double &lhs, primVars rhs) {
  for (auto rr = 0; rr < NUMVARS; rr++) {
    rhs.data_[rr] = lhs / rhs.data_[rr];
  }
  return rhs;
}

ostream &operator<<(ostream &os, const primVars &);

primVars RoeAveragedState(const primVars&, const primVars&, const idealGas&);

#endif

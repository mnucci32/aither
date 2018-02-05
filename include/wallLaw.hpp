/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef WALLLAWHEADERDEF
#define WALLLAWHEADERDEF

/* This defines the classes for the turbulent wall law implementation */
#include <memory>

#include "vector3d.hpp"
#include "tensor.hpp"
#include "primitive.hpp"
#include "arrayView.hpp"

using std::unique_ptr;

// forward class declaration
class physics;
struct wallVars;

class wallLaw {
  const bool isRANS_;
  const double vonKarmen_;
  const double wallConst_;
  const double wallDist_;
  const primitive state_;

  double yplus0_;
  double beta_;
  double gamma_;
  double q_;
  double phi_;
  double yplusWhite_;
  double uStar_;
  double uplus_;
  double tW_;
  double rhoW_;
  double muW_;
  double mutW_;
  double kW_;
  double recoveryFactor_;

  // private member functions
  void UpdateConstants(const double &);
  void UpdateGamma(const unique_ptr<thermodynamic> &, const double &);
  void CalcYplusWhite();
  double CalcHeatFlux(const unique_ptr<eos> &) const;
  void SetWallVars(const double &, const physics &);
  void EddyVisc(const physics &);
  void CalcVelocities(const double &, const double &);
  void CalcTurbVars(const physics &, double &, double &);
  double CalcYplusRoot(const double &) const;
  double ShearStressMag() const {return uStar_ * uStar_ * rhoW_;};
  void CalcRecoveryFactor(const unique_ptr<thermodynamic> &, const double &);
  double CalcWallTemperature(const physics &, const double &) const;

 public:
  // constructor
  template <typename T>
  wallLaw(const double &k, const double &c, const T &s, const double &d,
          const bool &isRANS)
      : isRANS_(isRANS),
        vonKarmen_(k),
        wallConst_(c),
        wallDist_(d),
        state_(s.begin(), s.end(), s.NumSpecies()),
        yplus0_(std::exp(-k * c)),
        beta_(0.0),
        gamma_(0.0),
        q_(0.0),
        phi_(0.0),
        yplusWhite_(0.0),
        uStar_(0.0),
        uplus_(0.0),
        tW_(0.0),
        rhoW_(0.0),
        muW_(0.0),
        mutW_(0.0),
        kW_(0.0),
        recoveryFactor_(0.0) {
    static_assert(std::is_same<primitive, T>::value ||
                      std::is_same<primitiveView, T>::value,
                  "T requires primitive or primativeView type");
  }

  // move constructor and assignment operator
  wallLaw(wallLaw&&) = default;
  wallLaw& operator=(wallLaw&&) = default;

  // copy constructor and assignment operator
  wallLaw(const wallLaw&) = default;
  wallLaw& operator=(const wallLaw&) = default;

  // member functions
  double VonKarmen() const { return vonKarmen_; }
  double WallConstant() const { return wallConst_; }
  wallVars AdiabaticBCs(const vector3d<double> &, const vector3d<double> &,
                        const vector<double> &, const physics &, const bool &);
  wallVars HeatFluxBCs(const vector3d<double> &, const vector3d<double> &,
                       const vector<double> &, const physics &, const double &,
                       const bool &);
  wallVars IsothermalBCs(const vector3d<double> &, const vector3d<double> &,
                         const vector<double> &, const physics &,
                         const double &, const bool &);

  // destructor
  ~wallLaw() noexcept {}
};


// function declarations



#endif

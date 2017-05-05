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

#ifndef WALLLAWHEADERDEF
#define WALLLAWHEADERDEF

/* This defines the classes for the turbulent wall law implementation */
#include <memory>

#include "vector3d.hpp"
#include "tensor.hpp"
#include "primVars.hpp"

using std::unique_ptr;

// forward class declaration
class eos;
class transport;
class turbModel;
struct wallVars;

class wallLaw {
  const bool isRANS_;
  const double vonKarmen_;
  const double wallConst_;
  const double wallDist_;
  const primVars state_;

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
  void UpdateGamma(const unique_ptr<eos> &);
  void CalcYplusWhite();
  double CalcHeatFlux(const unique_ptr<eos> &) const;
  void SetWallVars(const double &, const unique_ptr<eos> &,
                   const unique_ptr<transport> &);
  void EddyVisc(const unique_ptr<eos> &, const unique_ptr<transport> &);
  void CalcVelocities(const double &, const double &);
  void CalcTurbVars(const unique_ptr<turbModel> &, const unique_ptr<eos> &,
                    const unique_ptr<transport> &, double &, double &);
  double CalcYplusRoot(const double &) const;
  double ShearStressMag() const {return uStar_ * uStar_ * rhoW_;};
  void CalcRecoveryFactor(const unique_ptr<eos> &);
  double CalcWallTemperature(const unique_ptr<eos> &, const double &) const;

 public:
  // constructor
  wallLaw(const double &k, const double &c, const primVars &s, const double &d,
          const bool &isRANS)
      : isRANS_(isRANS),
        vonKarmen_(k),
        wallConst_(c),
        wallDist_(d),
        state_(s),
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
        recoveryFactor_(0.0) {}

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
                        const unique_ptr<eos> &, const unique_ptr<transport> &,
                        const unique_ptr<turbModel> &, const bool &);
  wallVars HeatFluxBCs(const vector3d<double> &, const vector3d<double> &,
                       const unique_ptr<eos> &, const unique_ptr<transport> &,
                       const unique_ptr<turbModel> &, const double &,
                       const bool &);
  wallVars IsothermalBCs(const vector3d<double> &, const vector3d<double> &,
                         const unique_ptr<eos> &, const unique_ptr<transport> &,
                         const unique_ptr<turbModel> &, const double &,
                         const bool &);

  // destructor
  ~wallLaw() noexcept {}
};


// function declarations



#endif

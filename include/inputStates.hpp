/*  This file is part of aither.
    Copyright (C) 2015-16  Michael Nucci (michael.nucci@gmail.com)

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

#ifndef INPUTSTATESHEADERDEF
#define INPUTSTATESHEADERDEF

/* This header file contains the bcState and icState classes and their children

These states are used to store data from the input file and apply boundary
conditions and initial conditions.
*/

#include <string>  // string
#include "vector3d.hpp"

using std::string;

class icState {
  vector3d<double> velocity_;
  double density_;
  double pressure_;
  double turbIntensity_;
  double eddyViscRatio_;
  int tag_;

 public:
  // Constructor
  icState(const vector3d<double> &v, const double &r, const double &p,
          const double &ti, const double &evr, const int &t) :
      velocity_(v), density_(r), pressure_(p), turbIntensity_(ti),
      eddyViscRatio_(evr), tag_(t) {}
  icState(const vector3d<double> &v, const double &r, const double &p,
          const double &ti, const double &evr) :
      velocity_(v), density_(r), pressure_(p), turbIntensity_(ti),
      eddyViscRatio_(evr), tag_(-1) {}
  icState(const vector3d<double> &v, const double &r, const double &p,
          const int &t) :
      velocity_(v), density_(r), pressure_(p), turbIntensity_(0.0),
      eddyViscRatio_(0.0), tag_(t) {}
  icState(const vector3d<double> &v, const double &r, const double &p) :
      velocity_(v), density_(r), pressure_(p), turbIntensity_(0.0),
      eddyViscRatio_(0.0), tag_(-1) {}

  // move constructor and assignment operator
  icState(icState&&) noexcept = default;
  icState& operator=(icState&&) noexcept = default;

  // copy constructor and assignment operator
  icState(const icState&) = default;
  icState& operator=(const icState&) = default;

  // Member functions
  const vector3d<double> & Velocity() const {return velocity_;}
  const double & Density() const {return density_;}
  const double & Pressure() const {return pressure_;}
  const double & TurbulenceIntensity() const {return turbIntensity_;}
  const double & EddyViscosityRatio() const {return eddyViscRatio_;}
  const int & Tag() const {return tag_;}

  // Destructor
  ~icState() noexcept {}
};

// function declarations
ostream &operator<<(ostream &, const icState &);

#endif

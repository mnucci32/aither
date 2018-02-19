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

#ifndef REACTIONHEADERDEF
#define REACTIONHEADERDEF

#include <vector>
#include <cmath>

using std::vector;

// class to hold reaction data
class reaction {
  vector<double> stoichReactants_;
  vector<double> stoichProducts_;
  double arrheniusC_;
  double arrheniusEta_;
  double arrheniusTheta_;

 public:
  // Constructor
  reaction(const int& ns)
      : stoichReactants_(ns, 0.0), stoichProducts_(ns, 0.0) {}

  // move constructor and assignment operator
  reaction(reaction&&) noexcept = default;
  reaction& operator=(reaction&&) noexcept = default;

  // copy constructor and assignment operator
  reaction(const reaction&) = default;
  reaction& operator=(const reaction&) = default;

  // Member functions for abstract base class
  const double &StoichReactant(const int &ss) const {
    return stoichReactants_[ss];
  }
  const double &StoichProduct(const int &ss) const {
    return stoichProducts_[ss];
  }
  double ForwardRate(const double &t) const {
    return arrheniusC_ * pow(t, arrheniusEta_) * std::exp(-arrheniusTheta_ / t);
  }
  double BackwardRate(const double &t) const {
    return this->ForwardRate(t) / this->EquilibriumRate(t);
  }
  double EquilibriumRate(const double &t) const { return 0.0; }

  // Destructor
  ~reaction() noexcept {}
};

// --------------------------------------------------------------------------
// function declarations



#endif

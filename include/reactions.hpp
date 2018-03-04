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
#include <string>
#include <iostream>
#include <memory>

using std::vector;
using std::string;
using std::unique_ptr;

// forward class declarations
class input;
class thermodynamic;

// class to hold reaction data
class reaction {
  vector<double> stoichReactants_;
  vector<double> stoichProducts_;
  vector<double> modifyReactants_;
  vector<string> species_;
  double arrheniusC_;
  double arrheniusEta_;
  double arrheniusTheta_;
  bool isForwardOnly_;
  bool isNondimensional_ = false;

 public:
  // Constructor
  reaction(const string &str, const input &inp);

  // move constructor and assignment operator
  reaction(reaction&&) noexcept = default;
  reaction& operator=(reaction&&) noexcept = default;

  // copy constructor and assignment operator
  reaction(const reaction&) = default;
  reaction& operator=(const reaction&) = default;

  // Member functions
  void Print(std::ostream &os) const;
  bool IsForwardOnly() const { return isForwardOnly_; }
  const double &StoichReactant(const int &ss) const {
    return stoichReactants_[ss];
  }
  const double &StoichProduct(const int &ss) const {
    return stoichProducts_[ss];
  }
  double ForwardRate(const double &t) const {
    return arrheniusC_ * pow(t, arrheniusEta_) * std::exp(-arrheniusTheta_ / t);
  }
  double BackwardRate(const double &t, const double &R, const double &refP,
                      const unique_ptr<thermodynamic> &thermo) const {
    return isForwardOnly_ ? 0.0
                          : this->ForwardRate(t) /
                                this->EquilibriumRate(t, R, refP, thermo);
  }
  double EquilibriumRate(const double &, const double &, const double &,
                         const unique_ptr<thermodynamic> &) const;
  void Nondimensionalize(const double &tref, const double &lref,
                         const double &aref) {
    if (!isNondimensional_) {
      const auto tauRef = lref / aref;
      arrheniusTheta_ /= tref;
      arrheniusC_ *= tauRef;
    }
    isNondimensional_ = true;
  }
  bool IsNondimensional() const { return isNondimensional_; }

  // Destructor
  ~reaction() noexcept {}
};


// --------------------------------------------------------------------------
// function declarations

std::ostream &operator<<(std::ostream &, const reaction &);

#endif

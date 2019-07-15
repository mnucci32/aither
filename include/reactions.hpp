/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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
#include <numeric>
#include <memory>

using std::vector;
using std::string;
using std::unique_ptr;

// forward class declarations
class input;

// class to hold reaction data
class reaction {
  vector<double> stoichReactants_;
  vector<double> stoichProducts_;
  vector<double> modifyReactants_;
  vector<string> species_;
  double arrheniusC_;
  double arrheniusEta_;
  double arrheniusTheta_;
  double universalGasConst_;
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
  double BackwardRate(const double &t, const double &refP,
                      const vector<double> &omega) const {
    return isForwardOnly_
               ? 0.0
               : this->ForwardRate(t) / this->EquilibriumRate(t, refP, omega);
  }
  double EquilibriumRate(const double &, const double &,
                         const vector<double> &) const;
  void Nondimensionalize(const double &tref, const double &lref,
                         const double &aref) {
    if (!isNondimensional_) {
      // forward rate units are (mol / m^3)^(1 - nu_reac_sum) / s
      // backward rate units are (mol / m^3)^(1 - nu_prod_sum) / s
      arrheniusTheta_ /= tref;
      const auto tauRef = lref / aref;
      const auto nuReacSum = std::accumulate(std::begin(stoichReactants_),
                                             std::end(stoichReactants_), 0.0);
      const auto conRef = pow(1.0 / pow(lref, 3.0), 1.0 - nuReacSum);
      arrheniusC_ *= tauRef * pow(tref, arrheniusEta_) / conRef;
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

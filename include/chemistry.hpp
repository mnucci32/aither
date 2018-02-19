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

#ifndef CHEMISTRYHEADERDEF
#define CHEMISTRYHEADERDEF

// This header file contains the classes for chemistry models

#include <vector>
#include "reactions.hpp"

using std::vector;

// abstract base class for chemistry models
class chemistry {
  int numReactions_;

 public:
  // Constructor
  chemistry(const int &nr) : numReactions_(nr) {}
  chemistry() : chemistry(0) {}

  // move constructor and assignment operator
  chemistry(chemistry&&) noexcept = default;
  chemistry& operator=(chemistry&&) noexcept = default;

  // copy constructor and assignment operator
  chemistry(const chemistry&) = default;
  chemistry& operator=(const chemistry&) = default;

  // Member functions for abstract base class
  const int& NumReactions() const { return numReactions_; }
  virtual vector<double> SourceTerms(const vector<double>& rho,
                                     const double& t) const = 0;
  virtual void AddReaction(const reaction &r) {}
  virtual bool IsReacting() const { return false; }
  virtual double SrcSpecRad(const vector<double>& rho, const double& t,
                            const double& vol) const {
    return 0.0;
  }

  // Destructor
  virtual ~chemistry() noexcept {}
};

// this class models no reacting chemistry
class frozen : public chemistry {

 public:
  // Constructors
  frozen() : chemistry(0) {}

  // move constructor and assignment operator
  frozen(frozen&&) noexcept = default;
  frozen& operator=(frozen&&) noexcept = default;

  // copy constructor and assignment operator
  frozen(const frozen&) = default;
  frozen& operator=(const frozen&) = default;

  // Member functions
  vector<double> SourceTerms(const vector<double>& rho,
                             const double& t) const override {
    return vector<double>(rho.size(), 0.0);
  }

  // Destructor
  ~frozen() noexcept {}
};


// this class models reacting chemistry
class reacting : public chemistry {
  vector<reaction> reactions_;
  vector<double> molarMass_;

 public:
  // Constructors
  reacting(const int& nr, const vector<double>& m)
      : chemistry(nr), molarMass_(m) {
    reactions_.reserve(nr);
  }

  // move constructor and assignment operator
  reacting(reacting&&) noexcept = default;
  reacting& operator=(reacting&&) noexcept = default;

  // copy constructor and assignment operator
  reacting(const reacting&) = default;
  reacting& operator=(const reacting&) = default;

  // Member functions
  vector<double> SourceTerms(const vector<double>& rho,
                             const double& t) const override;
  void AddReaction(const reaction& r) override { reactions_.push_back(r); }
  bool IsReacting() const override { return true; }
  double SrcSpecRad(const vector<double>& rho, const double& t,
                    const double& vol) const override {
    return 0.0;
  }

  // Destructor
  ~reacting() noexcept {}
};

// --------------------------------------------------------------------------
// function declarations



#endif

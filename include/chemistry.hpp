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

  // Destructor
  virtual ~chemistry() noexcept {}
};

// this class models no chemistry
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

  // Destructor
  ~frozen() noexcept {}
};

// --------------------------------------------------------------------------
// function declarations



#endif

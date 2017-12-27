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

#ifndef DIFFUSIONHEADERDEF
#define DIFFUSIONHEADERDEF

// This header file contains the classes for all diffusion models

#include <vector>

using std::vector;

// forward class declarations
class fluid;

// abstract base class for diffusion models
class diffusion {

 public:
  // Constructor
  diffusion() {}

  // move constructor and assignment operator
  diffusion(diffusion&&) noexcept = default;
  diffusion& operator=(diffusion&&) noexcept = default;

  // copy constructor and assignment operator
  diffusion(const diffusion&) = default;
  diffusion& operator=(const diffusion&) = default;

  // Member functions for abstract base class
  virtual double SpeciesDiffusion(const int &ii, const int &jj) const = 0;

  // Destructor
  virtual ~diffusion() noexcept {}
};

// this class models no diffusion
class diffNone : public diffusion {

 public:
  // Constructors
  diffNone() {}

  // move constructor and assignment operator
  diffNone(diffNone&&) noexcept = default;
  diffNone& operator=(diffNone&&) noexcept = default;

  // copy constructor and assignment operator
  diffNone(const diffNone&) = default;
  diffNone& operator=(const diffNone&) = default;

  // Member functions
  double SpeciesDiffusion(const int &ii, const int &jj) const override {
    return 0.0;
  }

  // Destructor
  ~diffNone() noexcept {}
};



// this class models diffusion using Schmidt number
class schmidt : public diffusion {
  vector<double> diffCoeff_;

  // private member functions
  double WilkesDiff(const vector<double> &, const vector<double> &) const;

 public:
  // Constructors
  schmidt(const vector<fluid> &);

  // move constructor and assignment operator
  schmidt(schmidt&&) noexcept = default;
  schmidt& operator=(schmidt&&) noexcept = default;

  // copy constructor and assignment operator
  schmidt(const schmidt&) = default;
  schmidt& operator=(const schmidt&) = default;

  // Member functions
  int NumSpecies() const { return diffCoeff_.size(); }
  double SpeciesDiffusion(const int &, const int &) const override;

  // Destructor
  ~schmidt() noexcept {}
};


// --------------------------------------------------------------------------
// function declarations



#endif

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

#ifndef PHYSICSMODELSHEADERDEF
#define PHYSICSMODELSHEADERDEF

// This header file contains the classes for all physics models

#include <memory>
#include "eos.hpp"
#include "transport.hpp"
#include "thermodynamic.hpp"
#include "diffusion.hpp"
#include "turbulence.hpp"

using std::unique_ptr;

// forward class declarations
class fluid;

// class to store all physics models
class physics {
  unique_ptr<eos> eos_;
  unique_ptr<transport> transport_;
  unique_ptr<thermodynamic> thermodynamic_;
  unique_ptr<diffusion> diffusion_;
  unique_ptr<turbModel> turbulence_;

 public:
  // Constructor
  physics(unique_ptr<eos> &eqnState, unique_ptr<transport> &trans,
          unique_ptr<thermodynamic> &thermo, unique_ptr<diffusion> &diff,
          unique_ptr<turbModel> &turb)
      : eos_(std::move(eqnState)),
        transport_(std::move(trans)),
        thermodynamic_(std::move(thermo)),
        diffusion_(std::move(diff)),
        turbulence_(std::move(turb)) {}

  // move constructor and assignment operator
  physics(physics&&) noexcept = default;
  physics& operator=(physics&&) noexcept = default;

  // copy constructor and assignment operator
  // class is non-copyable b/c it holds unique_ptrs
  physics(const physics&) = delete;
  physics& operator=(const physics&) = delete;

  // Member functions
  const unique_ptr<eos> &EoS() const { return eos_; }
  const unique_ptr<transport> &Transport() const { return transport_; }
  const unique_ptr<thermodynamic> &Thermodynamic() const {
    return thermodynamic_;
  }
  const unique_ptr<diffusion> &Diffusion() const { return diffusion_; }
  const unique_ptr<turbModel> &Turbulence() const { return turbulence_; }

  // Destructor
  virtual ~physics() noexcept {}
};


// --------------------------------------------------------------------------
// function declarations



#endif

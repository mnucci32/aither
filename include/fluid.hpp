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

#ifndef FLUIDHEADERDEF
#define FLUIDHEADERDEF

// This header file contains the properties for a given fluid

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using std::vector;
using std::ifstream;
using std::string;
using std::ostream;

// class for fluid properties
class fluid {
  double n_ = 2.5;
  double molarMass_ = 0.02897;  // kg/mol
  double vibTemp_ = 3056.0;    // K
  double universalGasConst_ = 8.3144598;  // J / mol-K
  string name_ = "air";
  bool nondimensional_ = false;

  // private member functions
  void SetNondimensional(const bool &nd) { nondimensional_ = nd; }

 public:
  // Constructor
  fluid() {}
  fluid(string& str, const string = "fluid");

  // move constructor and assignment operator
  fluid(fluid&&) noexcept = default;
  fluid& operator=(fluid&&) noexcept = default;

  // copy constructor and assignment operator
  fluid(const fluid&) = default;
  fluid& operator=(const fluid&) = default;

  // Member functions for abstract base class
  double N() const { return n_; }
  double MolarMass() const { return molarMass_; }
  double VibrationalTemperature() const { return vibTemp_; }
  double GasConstant() const { return universalGasConst_ / molarMass_; }
  double UniversalGasConstant() const { return universalGasConst_; }
  string Name() const { return name_; }
  bool IsNondimensional() const { return nondimensional_; }

  void Nondimensionalize(const double&);

  // Destructor
  ~fluid() noexcept {}
};

// function declarations
ostream &operator<<(ostream &, const fluid &);
vector<fluid> ReadFluidList(ifstream &, string &);


#endif

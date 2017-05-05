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

#ifndef THERMOHEADERDEF
#define THERMOHEADERDEF

// This header file contains the thermodynamic model classes

// abstract base class for thermodynamic model
class thermodynamic {

 public:
  // Constructor
  thermodynamic() {}

  // move constructor and assignment operator
  thermodynamic(thermodynamic&&) noexcept = default;
  thermodynamic& operator=(thermodynamic&&) noexcept = default;

  // copy constructor and assignment operator
  thermodynamic(const thermodynamic&) = default;
  thermodynamic& operator=(const thermodynamic&) = default;

  // Member functions for abstract base class
  virtual double Gamma() const = 0;
  virtual double Prandtl() const = 0;
  virtual double SpecificHeat() const = 0;

  // Destructor
  virtual ~thermodynamic() noexcept {}
};


class caloricallyPerfect : public thermodynamic {
  const double gamma_;

 public:
  // Constructor
  explicit caloricallyPerfect(const double &a) : gamma_(a) {}
  caloricallyPerfect() : caloricallyPerfect(1.4) {}

  // move constructor and assignment operator
  caloricallyPerfect(caloricallyPerfect&&) noexcept = default;
  caloricallyPerfect& operator=(caloricallyPerfect&&) noexcept = default;

  // copy constructor and assignment operator
  caloricallyPerfect(const caloricallyPerfect&) = default;
  caloricallyPerfect& operator=(const caloricallyPerfect&) = default;

  // Member functions
  double Gamma() const override {return gamma_;}
  double Prandtl() const override {return (4.0 * gamma_) / (9.0 * gamma_ - 5.0);}
  double SpecificHeat() const override {return 1.0 / (gamma_ - 1.0);}

  // Destructor
  ~caloricallyPerfect() noexcept {}
};


#endif

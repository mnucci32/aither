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

#ifndef CONSERVEDHEADERDEF
#define CONSERVEDHEADERDEF

#include <iostream>
#include "varArray.hpp"

using std::ostream;

class conserved : public varArray {
 public:
  // constructor
  conserved(const int &numEqns, const int &numSpecies)
      : varArray(numEqns, numSpecies) {}
  conserved(const vector<double>::const_iterator &b,
            const vector<double>::const_iterator &e, const int &numSpecies)
      : varArray(b, e, numSpecies) {}

  // member functions
  const double & RhoN(const int &ii) const { return this->SpeciesN(ii); }
  double Rho() const { return this->SpeciesSum(); }
  const double & RhoU() const { return this->MomentumX(); }
  const double & RhoV() const { return this->MomentumY(); }
  const double & RhoW() const { return this->MomentumZ(); }
  const double & RhoE() const { return this->Energy(); }
  const double & RhoTke() const { return this->TurbulenceN(0); }
  const double & RhoOmega() const { return this->TurbulenceN(1); }
  const double & RhoTurbN(const int &ii) const { return this->TurbulenceN(ii); }

  // move constructor and assignment operator
  conserved(conserved&&) noexcept = default;
  conserved& operator=(conserved&&) noexcept = default;

  // copy constructor and assignment operator
  conserved(const conserved&) = default;
  conserved& operator=(const conserved&) = default;

  // destructor
  ~conserved() noexcept {}
};

ostream &operator<<(ostream &os, const conserved &);



#endif

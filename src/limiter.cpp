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

#include <algorithm>
#include <cmath>
#include "limiter.hpp"
#include "primitive.hpp"


// function to calculate minmod limiter
primitive LimiterMinmod(const primitive &r) {
  // r -- ratio of divided differences

  primitive limiter(r.Size(), r.NumSpecies());
  // calculate minmod limiter
  for (auto ii = 0; ii < limiter.Size(); ++ii) {
    limiter[ii] = std::max(0.0, std::min(1.0, r[ii]));
  }
  return limiter;
}

// function to calculate Van Albada limiter
primitive LimiterVanAlbada(const primitive &r) {
  // r -- ratio of divided differences

  const auto r2 = r * r;
  auto limiter = (r + r2) / (1.0 + r2);
  // if value is negative, return zero
  for (auto ii = 0; ii < limiter.Size(); ++ii) {
    limiter[ii] = std::max(0.0, limiter[ii]);
  }
  return limiter;
}

// function to return no limiter
primitive LimiterNone(const int &numEq, const int &numSpec) {
  // for no limiter return all 1s
  primitive limiter(numEq, numSpec, 1.0);
  return limiter;
}


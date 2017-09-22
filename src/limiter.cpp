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

#include <algorithm>
#include <cmath>
#include "limiter.hpp"
#include "primitive.hpp"
#include "utility.hpp"    // Sign


// function to calculate minmod limiter
primitive LimiterMinmod(const primitive &upwind, const primitive &downwind,
                        const double &kap) {
  // upwind -- upwind state (primitive)
  // downwind -- downwind state (primitive)
  // kap -- MUSCL parameter kappa

  primitive limiter(upwind.Size(), upwind.NumSpecies());

  // calculate minmod parameter beta
  const auto beta = (3.0 - kap) / (1.0 - kap);

  // calculate minmod limiter
  for (auto ii = 0; ii < limiter.Size(); ++ii) {
    auto sign = Sign(upwind[ii]);
    limiter[ii] = sign * std::max(0.0, std::min(fabs(upwind[ii]),
                                                sign * downwind[ii] * beta));
  }
  return limiter;
}

// function to calculate Van Albada limiter
primitive LimiterVanAlbada(const primitive &r) {
  // r -- ratio of divided differences

  auto limiter = (r + r * r) / (1.0 + r * r);
  // if value is negative, return zero
  for (auto ii = 0; ii < limiter.Size(); ++ii) {
    limiter[ii] = std::max(0.0, limiter[ii]);
  }
  return limiter;
}

// function to return no limiter
primitive LimiterNone(const primitive &state) {
  // for no limiter return all 1s
  primitive limiter(state.Size(), state.NumSpecies(), 1.0);
  return limiter;
}


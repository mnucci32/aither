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

#include <cstdlib>  // exit()
#include <iostream>  // cout
#include <cmath>
#include "varArray.hpp"
#include "mpi.h"

using std::cout;
using std::endl;
using std::cerr;

// ------------------------------------------------------------------
// functions for varArray class

varArray Squared() const {
  auto sq = (*this);
  sq *= sq;
  return sq;
}

// member function to sum the residuals from all processors
void residual::GlobalReduceMPI(const int &rank) {
  // Get residuals from all processors
  if (rank == ROOTP) {
    MPI_Reduce(MPI_IN_PLACE, &(*this)[0], this->Size(), MPI_DOUBLE, MPI_SUM,
               ROOTP, MPI_COMM_WORLD);
  } else {
    MPI_Reduce((*this)[0], &(*this)[0], this->Size(), MPI_DOUBLE, MPI_SUM, ROOTP,
               MPI_COMM_WORLD);
  }
}

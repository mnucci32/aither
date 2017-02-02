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

#include "resid.hpp"

// member function to get the addresses of a resid to create an MPI_Datatype
void resid::GetAddressesMPI(MPI_Aint (&displacement)[2]) const {
  // get addresses of each field
  MPI_Get_address(&(*this).linf_, &displacement[0]);
  MPI_Get_address(&(*this).blk_, &displacement[1]);
}

// member function to update class when a new maximum residual is found
void resid::UpdateMax(const double &a, const int &b, const int &c, const int &d,
                      const int &e, const int &f) {
  linf_ = a;
  blk_ = b;
  i_ = c;
  j_ = d;
  k_ = e;
  eqn_ = f;
}

// Member function to calculate the maximum residual from all processors
void resid::GlobalReduceMPI(const int &rank,
                            const MPI_Datatype &MPI_DOUBLE_5INT,
                            const MPI_Op &MPI_MAX_LINF) {
  // Get residuals from all processors
  if (rank == ROOTP) {
    MPI_Reduce(MPI_IN_PLACE, &(*this), 1, MPI_DOUBLE_5INT, MPI_MAX_LINF, ROOTP,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&(*this), &(*this), 1, MPI_DOUBLE_5INT, MPI_MAX_LINF, ROOTP,
               MPI_COMM_WORLD);
  }
}

/* Function to calculate the maximum of two resid instances and allow access to
all the data_ in the resid instance. This is used to create an operation
for MPI_Reduce.
*/
void MaxLinf(resid *in, resid *inout, int *len, MPI_Datatype *MPI_DOUBLE_5INT) {
  // *in -- pointer to all input residuals (from all procs)
  // *inout -- pointer to input and output residuals. The answer is stored here
  // *len -- pointer to array size of *in and *inout
  // *MPI_DOUBLE_5INT -- pointer to MPI_Datatype of double followed by 5 Ints,
  // which represents the resid class

  resid resLinf;  // intialize a resid

  for (auto ii = 0; ii < *len; ii++) {  // loop over array of resids
    if (in->Linf() >= inout->Linf()) {  // if linf from input is greater than or
                                        // equal to linf from output, then new
                                        // max has been found
      resLinf = *in;
    } else {  // no new max
      resLinf = *inout;
    }

    *inout = resLinf;  // assign max to output

    // increment to next entry in array
    in++;
    inout++;
  }
}

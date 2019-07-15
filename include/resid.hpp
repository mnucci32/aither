/*  This file is part of aither.
    Copyright (C) 2015-19  Michael Nucci (mnucci@pm.me)

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

#ifndef RESIDHEADERDEF  // only if the macro RESIDHEADERDEF is not defined
                         // execute these lines of code
#define RESIDHEADERDEF  // define the macro

#include "mpi.h"       // parallelism
#include "macros.hpp"

class resid {
  double linf_;
  int blk_;
  int i_;
  int j_;
  int k_;
  int eqn_;

 public:
  // constructor
  resid() : linf_(0.0), blk_(0), i_(0), j_(0), k_(0), eqn_(0) {}
  resid(const double &a, const int &b, const int &c, const int &d, const int &e,
        const int &f)
      : linf_(a), blk_(b), i_(c), j_(d), k_(e), eqn_(f) {}

  // move constructor and assignment operator
  resid(resid&&) noexcept = default;
  resid& operator=(resid&&) noexcept = default;

  // copy constructor and assignment operator
  resid(const resid&) = default;
  resid& operator=(const resid&) = default;

  // member functions
  double Linf() const { return linf_; }
  int Block() const { return blk_; }
  int ILoc() const { return i_; }
  int JLoc() const { return j_; }
  int KLoc() const { return k_; }
  int Eqn() const { return eqn_; }

  void UpdateMax(const double &, const int &, const int &, const int &,
                 const int &, const int &);
  void GetAddressesMPI(MPI_Aint (&)[2]) const;
  void GlobalReduceMPI(const int &, const MPI_Datatype &, const MPI_Op &);

  void MaxLinf(resid *, resid *, int *, MPI_Datatype *);

  void Zero() {
    linf_ = 0.0;
    blk_ = 0;
    i_ = 0;
    j_ = 0;
    k_ = 0;
    eqn_ = 0;
  }

  // destructor
  ~resid() noexcept {}
};

// function declarations

#endif

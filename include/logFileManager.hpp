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

#ifndef LOG_FILE_MANAGER_HEADERDEF
#define LOG_FILE_MANAGER_HEADERDEF

#include <iostream>
#include <fstream>
#include <chrono>
#include "varArray.hpp"
#include "macros.hpp"

// forward class declaration
class input;
class resid;

class logFileManager {
  std::ofstream resid_;
  std::ofstream time_;
  int rank_;
  residual l2First_;
  std::chrono::high_resolution_clock::time_point simStartTime_;
  std::chrono::high_resolution_clock::time_point iterStartTime_;

 public:
  // Constructor
  logFileManager(const input& inp, const int &rank);

  // Member functions
  residual &L2First() { return l2First_; }
  const residual &L2First() const { return l2First_; }
  void WriteResiduals(const input &, const residual &, const resid &,
                      const double &, const int &, const int &);
  void GetIterStart();
  void WriteTime(const int &nn);

  // Destructor
  ~logFileManager();
};



#endif

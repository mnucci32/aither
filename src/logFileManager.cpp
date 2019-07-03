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

#include <iostream>     // cout
#include <iomanip>
#include <cstdlib>      // exit()
#include <fstream>
#include <chrono>
#include "logFileManager.hpp"
#include "input.hpp"
#include "output.hpp"
#include "resid.hpp"
#include "varArray.hpp"
#include "macros.hpp"

using std::cout;
using std::endl;
using std::cerr;

// constructor
logFileManager::logFileManager(const input &inp, const int &rank)
    : l2First_(inp.NumEquations(), inp.NumSpecies()) {
  rank_ = rank;
  if (rank_ != ROOTP) {
    return;
  }

  simStartTime_ = std::chrono::high_resolution_clock::now();

  // open residual file
  if (inp.IsRestart()) {
    resid_.open(inp.SimNameRoot() + ".resid", std::ios::app);
  } else {
    resid_.open(inp.SimNameRoot() + ".resid", std::ios::out);
  }
  if (resid_.fail()) {
    cerr << "ERROR: Could not open residual file!" << endl;
    exit(EXIT_FAILURE);
  }
  PrintHeaders(inp, resid_);

  // open timing file
  time_.open(inp.SimNameRoot() + ".tme", std::ios::out);
  if (time_.fail()) {
    cerr << "ERROR: Could not open residual file!" << endl;
    exit(EXIT_FAILURE);
  }
  time_ << std::left << std::setw(7) << "Step" << std::setw(16) << "Iter-Time"
        << std::setw(16) << "Sim-Time" << endl;
}

// destructor
logFileManager::~logFileManager() { 
  if (rank_ != ROOTP) {
    return;
  }
  resid_.close();
  time_.close();
}

// function to write out residual information
void logFileManager::WriteResiduals(const input &inp, const residual &residL2,
                                    const resid &residLinf,
                                    const double &matrixResid, const int &nn,
                                    const int &mm) {
  if (rank_ != ROOTP) {
    return;
  }

  // write out column headers every 100 iterations to standard out
  if (nn % 100 == 0 && mm == 0) {
    PrintHeaders(inp, cout);
  }

  // print residuals to standard out
  PrintResiduals(inp, l2First_, residL2, residLinf, matrixResid, nn, mm, cout);
  // print residuals to residual file
  PrintResiduals(inp, l2First_, residL2, residLinf, matrixResid, nn, mm,
                 resid_);
}

void logFileManager::GetIterStart() {
  if (rank_ != ROOTP) {
    return;
  }
  iterStartTime_ = std::chrono::high_resolution_clock::now();
}

void logFileManager::WriteTime(const int &nn) {
  if (rank_ != ROOTP) {
    return;
  }
  const auto currTime = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> duration = currTime - simStartTime_;
  const std::chrono::duration<double> iterDuration = currTime - iterStartTime_;
  time_ << std::left << std::setw(7) << nn << std::setw(16)
        << std::setprecision(6) << std::scientific << iterDuration.count()
        << std::setw(16) << duration.count() << endl;
  time_.unsetf(std::ios::fixed | std::ios::scientific);
}
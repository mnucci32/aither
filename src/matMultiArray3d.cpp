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
#include <vector>
#include "matMultiArray3d.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;

// operation overload for << - allows use of cout, cerr, etc.
ostream &operator<<(ostream &os, const matMultiArray3d &arr) {
  os << "Size: " << arr.NumI() << ", " << arr.NumJ() << ", " << arr.NumK()
     << endl;
  os << "Number of ghost layers: " << arr.GhostLayers() << endl;

  auto GetVal = [](const auto &mat, const int &r, const int &c,
                   const int &size) -> decltype(auto) {
    return *(mat + r * size + c);
  };

  for (auto kk = arr.StartK(); kk < arr.EndK(); ++kk) {
    for (auto jj = arr.StartJ(); jj < arr.EndJ(); ++jj) {
      for (auto ii = arr.StartI(); ii < arr.EndI(); ++ii) {
        os << ii << ", " << jj << ", " << kk << endl;

        // print flow jacobian
        for (auto rf = 0; rf < arr.FlowSize(); ++rf) {
          for (auto cf = 0; cf < arr.FlowSize(); ++cf) {
            os << GetVal(arr.begin() + arr.GetLoc1D(ii, jj, kk), rf, cf,
                         arr.FlowSize())
               << " ";
            if (cf == arr.FlowSize() - 1) {
              os << endl;
            }
          }
        }

        // print turbulence jacobian
        for (auto rt = 0; rt < arr.TurbSize(); ++rt) {
          for (auto ct = 0; ct < arr.TurbSize(); ++ct) {
            os << GetVal(arr.begin() + arr.GetLoc1D(ii, jj, kk) +
                             arr.FlowSize() * arr.FlowSize(),
                         rt, ct, arr.FlowSize())
               << " ";
            if (ct == arr.TurbSize() - 1) {
              os << endl;
            }
          }
        }
      }
    }
  }
  return os;
}

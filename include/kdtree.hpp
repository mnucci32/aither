/*  This file is part of aither.
    Copyright (C) 2015-18  Michael Nucci (mnucci@pm.me)

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

// Only if the macro KDTREEHEADERDEF is not defined execute these lines of code
#ifndef KDTREEHEADERDEF
#define KDTREEHEADERDEF  // define the macro

#include <vector>        // vector
#include <string>        // string
#include <utility>       // pair
#include "vector3d.hpp"

using std::vector;
using std::string;
using std::pair;

/* This class implements a k-d tree for the purpose of efficiently performing
a nearest neighbor search. The search is used to find the minimum distance of
all cells to a viscous wall. The coordinates of cell faces contacting a 
viscous wall are stored in the data structure in the nodes_ variable. When
the tree is built the nodes are reordered such that a tree branch to the left
is the next indice down, and the index of a tree branch to the right is 
stored in the right_ variable. A simple 2D example will be shown to explain
how the data is ordered.

Given a set of points: (7,1), (0,2), (5,4), (3,9), (2,7), (8,6), (1,3), 
(6,5), (4,0), (9,8)

Using a bin size of 2, the k-d tree would look like the following:

Root                         (4,0)
                            /     \
Split - X                (1,3)   (6,5)
                        /    \    /   \
Split - Y           (0,2) (2,7) (7,1) (8,6)
                          (3,9) (5,4) (9,8)

This tree would be stored in the nodes_ and right_ vectors in the following
way:
nodes_ = [(4,0), (1,3), (0,2), (2,7), (3,9), (6,5), (7,1), (5,4), (8,6), (9.8)]
right_ = [  5,     3,    -1,    -1     -1,     8,    -1,    -1,    -1,    -1  ]

To traverse the tree, the branch to the left is given by the next indice in 
the nodes_ vector. For example point (1,3) is at index 1, and its left branch
is at index 2. The branch to the right is given by the right_ vector. For
example at point (1,3) is at index 1, its right branch is given by the
corresponding value at index 1 in the right_ vector. In this case its right
branch is at index 3.
*/
class kdtree {
  // all points to search through in k-d tree order
  vector<pair<vector3d<double>, int>> nodes_;
  vector<int> right_;           // right branch indices
  const int dim_ = 3;           // dimension of space to search
  const double binSize_ = 32;   // max number of points in a leaf node

  // private member functions
  int FindMedian(const int &, const int &, const int &);
  void BuildKdtree(const int &, const int &, const int &);
  void NearestNeighbor(const int &, const int &, const int &,
                       const vector3d<double> &, pair<vector3d<double>, int> &,
                       double &) const;
  bool SphereInside(const vector3d<double>&, const double &,
                    const vector3d<double>&, const int &) const;

 public:
  // constructor
  explicit kdtree(const vector<vector3d<double>> &);

  // move constructor and assignment operator
  kdtree(kdtree&&) noexcept = default;
  kdtree& operator=(kdtree&&) noexcept = default;

  // copy constructor and assignment operator
  kdtree(const kdtree&) = default;
  kdtree& operator=(const kdtree&) = default;

  // member functions
  double NearestNeighbor(const vector3d<double> &, vector3d<double> &,
                         int &) const;
  int Size() const { return nodes_.size(); }

  // destructor
  ~kdtree() noexcept {}
};

#endif

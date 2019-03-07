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

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "plot3d.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ios;

// plot 3d block member functions

// member function to calculate the centroid of a given cell
vector3d<double> plot3dBlock::Centroid(const int &ii, const int &jj,
                                       const int &kk) const {
  // centroid is average of all coordinates forming cell
  return 0.125 *
         (coords_(ii, jj, kk) + coords_(ii + 1, jj, kk) +
          coords_(ii, jj + 1, kk) + coords_(ii + 1, jj + 1, kk) +
          coords_(ii, jj, kk + 1) + coords_(ii + 1, jj, kk + 1) +
          coords_(ii, jj + 1, kk + 1) + coords_(ii + 1, jj + 1, kk + 1));
}

// plot3dBlock member function that calcualtes the volume of each cell
/*
All cells are assumed to be hexahedra. The 8 points that make up each hexahedron
are used to split the cell into 6 pyramids. The area of each pyramid is then
calculated and the volume of the 6 pyramids are summed to get the volume of the
hexahedra. The point of each pyramid is the centroid, and the six sides of the
cell make up the six bases.

Vp = Havg * (D1 (cross) D2))

The equation above shows how the volume of a pyramid (Vp) is calculated. Havg is
the average distance from the four base points to the top of the pyramid. D1 and
D2 are the diagonal distances of the four base points.
*/
multiArray3d<double> plot3dBlock::Volume() const {
  // Allocate multiArray3d to store cell volumes in
  multiArray3d<double> vol(this->NumI() - 1, this->NumJ() - 1,
                           this->NumK() - 1, 0);

  // Loop over all cells
  for (auto kk = 0; kk < vol.NumK(); kk++) {
    for (auto jj = 0; jj < vol.NumJ(); jj++) {
      for (auto ii = 0; ii < vol.NumI(); ii++) {
        // get centroid
        const auto centroid = this->Centroid(ii, jj, kk);

        // calculate area of i-lower pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii, jj, kk), coords_(ii, jj, kk + 1),
            coords_(ii, jj + 1, kk + 1), coords_(ii, jj + 1, kk));

        // calculate area of i-upper pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii + 1, jj, kk), coords_(ii + 1, jj + 1, kk),
            coords_(ii + 1, jj + 1, kk + 1), coords_(ii + 1, jj, kk + 1));

        // calculate area of j-lower pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii, jj, kk), coords_(ii + 1, jj, kk),
            coords_(ii + 1, jj, kk + 1), coords_(ii, jj, kk + 1));

        // calculate area of j-upper pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii, jj + 1, kk), coords_(ii, jj + 1, kk + 1),
            coords_(ii + 1, jj + 1, kk + 1), coords_(ii + 1, jj + 1, kk));

        // calculate area of k-lower pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii, jj, kk), coords_(ii, jj + 1, kk),
            coords_(ii + 1, jj + 1, kk), coords_(ii + 1, jj, kk));

        // calculate area of k-upper pyramid
        vol(ii, jj, kk) += PyramidVolume(
            centroid, coords_(ii, jj, kk + 1), coords_(ii + 1, jj, kk + 1),
            coords_(ii + 1, jj + 1, kk + 1), coords_(ii, jj + 1, kk + 1));

        // Check for negative volumes
        if (vol(ii, jj, kk) <= 0) {
          cerr << "ERROR: Negative volume in PLOT3D block!!!" << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj
               << ", k-dim = " << kk << ", vol = " << vol(ii, jj, kk) << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  return vol;
}

// plot3dBlock member function that calcualtes the centroid of each cell
// the centroid of the hexahedron is the average of the 8 points that define it
multiArray3d<vector3d<double>> plot3dBlock::Centroid() const {
  // Allocate multiArray3d to store cell centroids in
  multiArray3d<vector3d<double>> centroid(this->NumI() - 1, this->NumJ() - 1,
                                          this->NumK() - 1, 0);

  // loop over all cells
  for (auto kk = 0; kk < centroid.NumK(); kk++) {
    for (auto jj = 0; jj < centroid.NumJ(); jj++) {
      for (auto ii = 0; ii < centroid.NumI(); ii++) {
        // Calculate the centroid of the cell
        centroid(ii, jj, kk) = this->Centroid(ii, jj, kk);
      }
    }
  }

  return centroid;
}

/* plot3dBlock member function that calcualtes the area of each face normal to
   the i-direction the area is calculated as half the cross product of the
   diagonal vectors of the 4-sided face

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B
to C. The normal vector points in the direction of increasing i.
*/
multiArray3d<unitVec3dMag<double>> plot3dBlock::FaceAreaI() const {
  // Allocate multiArray3d to store cell face areas
  multiArray3d<unitVec3dMag<double>> fArea(this->NumI(), this->NumJ() - 1,
                                           this->NumK() - 1, 0);

  // loop over all i-faces
  for (auto kk = 0; kk < fArea.NumK(); kk++) {
    for (auto jj = 0; jj < fArea.NumJ(); jj++) {
      for (auto ii = 0; ii < fArea.NumI(); ii++) {
        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        auto xac = coords_(ii, jj + 1, kk + 1) - coords_(ii, jj, kk);
        auto xbd = coords_(ii, jj + 1, kk) - coords_(ii, jj, kk + 1);

        // area vector is calculated so that normal points nominally in
        // direction of increasing i-coordinate
        fArea(ii, jj, kk) = unitVec3dMag<double>(0.5 * xbd.CrossProd(xac));

        if (fArea(ii, jj, kk).Mag() <= 0) {
          cerr << "ERROR: Negative i-face area in PLOT3D block!!!" << endl;
          cerr << "Face area = " << fArea(ii, jj, kk).Mag() << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj << ", k-dim = " << kk
               << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and "
               << xbd << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  return fArea;
}

// plot3dBlock member function that calcualtes the center of each face normal to
// the i-direction the face center is calculated as the average of the 4 points
// that comprise it
multiArray3d<vector3d<double>> plot3dBlock::FaceCenterI() const {
  // Allocate multiArray3d to store cell face centers
  multiArray3d<vector3d<double>> fCenter(this->NumI(), this->NumJ() - 1,
                                         this->NumK() - 1, 0);
  // loop over all i-faces
  for (auto kk = 0; kk < fCenter.NumK(); kk++) {
    for (auto jj = 0; jj < fCenter.NumJ(); jj++) {
      for (auto ii = 0; ii < fCenter.NumI(); ii++) {
        // Calculate face center by averaging four points that make up the face
        fCenter(ii, jj, kk) =
            0.25 * (coords_(ii, jj, kk) + coords_(ii, jj + 1, kk) +
                    coords_(ii, jj, kk + 1) + coords_(ii, jj + 1, kk + 1));
      }
    }
  }
  return fCenter;
}

/* plot3dBlock member function that calcualtes the area of each face normal to
// the j-direction the area is calculated as half the cross product of the
// diagonal vectors of the 4-sided face

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B
to C. The normal points in the direction of increasing j.
*/
multiArray3d<unitVec3dMag<double>> plot3dBlock::FaceAreaJ() const {
  // Allocate multiArray3d to store cell face areas
  multiArray3d<unitVec3dMag<double>> fArea(this->NumI() - 1, this->NumJ(),
                                           this->NumK() - 1, 0);

  // loop over all j-faces
  for (auto kk = 0; kk < fArea.NumK(); kk++) {
    for (auto jj = 0; jj < fArea.NumJ(); jj++) {
      for (auto ii = 0; ii < fArea.NumI(); ii++) {
        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        auto xac = coords_(ii, jj, kk + 1) - coords_(ii + 1, jj, kk);
        auto xbd = coords_(ii, jj, kk) - coords_(ii + 1, jj, kk + 1);

        // area vector is calculated so that normal nominally points in
        // direction of increasing j-coordinate
        fArea(ii, jj, kk) = unitVec3dMag<double>(0.5 * xbd.CrossProd(xac));

        if (fArea(ii, jj, kk).Mag() <= 0) {
          cerr << "ERROR: Negative j-face area in PLOT3D block!!!" << endl;
          cerr << "Face area = " << fArea(ii, jj, kk).Mag() << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj
               << ", k-dim = " << kk << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and "
               << xbd << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  return fArea;
}

// plot3dBlock member function that calcualtes the area of each face normal to
// the j-direction
// the face center is calculated as the average of the 4 points that comprise it
multiArray3d<vector3d<double>> plot3dBlock::FaceCenterJ() const {
  // Allocate multiArray3d to store cell face centers
  multiArray3d<vector3d<double>> fCenter(this->NumI() - 1, this->NumJ(),
                                         this->NumK() - 1, 0);

  // loop over all j-faces
  for (auto kk = 0; kk < fCenter.NumK(); kk++) {
    for (auto jj = 0; jj < fCenter.NumJ(); jj++) {
      for (auto ii = 0; ii < fCenter.NumI(); ii++) {
        // Calculate face center by averaging the four points that make up the
        // face
        fCenter(ii, jj, kk) =
            0.25 * (coords_(ii, jj, kk) + coords_(ii + 1, jj, kk) +
                    coords_(ii, jj, kk + 1) + coords_(ii + 1, jj, kk + 1));
      }
    }
  }
  return fCenter;
}

/* plot3dBlock member function that calcualtes the area of each face normal to
// the k-direction the area is calculated as half the cross product of the
// diagonal vectors of the 4-sided face

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B
to C. The normal vector points in the direction of increasing k.
*/
multiArray3d<unitVec3dMag<double>> plot3dBlock::FaceAreaK() const {
  // Allocate multiArray3d to store cell face areas
  multiArray3d<unitVec3dMag<double>> fArea(this->NumI() - 1, this->NumJ() - 1,
                                           this->NumK(), 0);

  // loop over all k-faces
  for (auto kk = 0; kk < fArea.NumK(); kk++) {
    for (auto jj = 0; jj < fArea.NumJ(); jj++) {
      for (auto ii = 0; ii < fArea.NumI(); ii++) {
        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        auto xac = coords_(ii, jj + 1, kk) - coords_(ii + 1, jj, kk);
        auto xbd = coords_(ii + 1, jj + 1, kk) - coords_(ii, jj, kk);

        // area vector is calculated so that normal nominally points in
        // direction of increasing k-coordinate
        fArea(ii, jj, kk) = unitVec3dMag<double>(0.5 * xbd.CrossProd(xac));

        if (fArea(ii, jj, kk).Mag() <= 0) {
          cerr << "ERROR: Negative k-face area in PLOT3D block!!!" << endl;
          cerr << "Face area = " << fArea(ii, jj, kk).Mag() << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj
               << ", k-dim = " << kk << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and "
               << xbd << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  return fArea;
}

// plot3dBlock member function that calcualtes the area of each face normal to
// the k-direction
// the face center is calculated as the average of the 4 points that comprise it
multiArray3d<vector3d<double>> plot3dBlock::FaceCenterK() const {
  // Allocate multiArray3d to store cell face centers
  multiArray3d<vector3d<double>> fCenter(this->NumI() - 1, this->NumJ() - 1,
                                         this->NumK(), 0);

  // loop over all k-faces
  for (auto kk = 0; kk < fCenter.NumK(); kk++) {
    for (auto jj = 0; jj < fCenter.NumJ(); jj++) {
      for (auto ii = 0; ii < fCenter.NumI(); ii++) {
        // baseline location
        auto botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        auto botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the j-direction
        auto botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the i and j directions
        auto botPortFore = coords_(ii + 1, jj + 1, kk);

        // Calculate face center by averaging four points that make up cell face
        fCenter(ii, jj, kk) =
            0.25 * (coords_(ii, jj, kk) + coords_(ii + 1, jj, kk) +
                    coords_(ii, jj + 1, kk) + coords_(ii + 1, jj + 1, kk));
      }
    }
  }
  return fCenter;
}

//------------------------------------------------------------------------------
// function to read in a plot3d grid and assign it to a plot3dMesh data type
vector<plot3dBlock> ReadP3dGrid(const string &gridName, const double &LRef,
                                double &numCells) {
  // open binary plot3d grid file
  ifstream fName;
  string fPostfix = ".xyz";
  auto readName = gridName + fPostfix;
  fName.open(readName, ios::in | ios::binary);

  // check to see if file opened correctly
  if (fName.fail()) {
    cerr << "ERROR: Error in plot3d.cpp:ReadP3dGrid(). Grid file " << readName
         << " did not open correctly!!!" << endl;
    exit(EXIT_FAILURE);
  }

  // read the number of plot3d blocks in the file
  cout << "Reading grid file..." << endl << endl;
  auto numBlks = 1;
  fName.read(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));
  cout << "Number of blocks: " << numBlks << endl << endl;

  // read the number of i, j, k coordinates in each plot3d block
  cout << "Size of each block is..." << endl;
  vector<vector3d<int>> blkSize(numBlks);
  auto tempInt = 0;
  numCells = 0;

  // loop over all blocks and fill i, j, k vectors with block sizes
  for (auto ii = 0; ii < numBlks; ii++) {
    cout << "Block Number: " << ii << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    blkSize[ii][0] = tempInt;
    cout << "I-DIM: " << tempInt << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    blkSize[ii][1] = tempInt;
    cout << "J-DIM: " << tempInt << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    blkSize[ii][2] = tempInt;
    cout << "K-DIM: " << tempInt << endl;

    // calculate total number of cells (subtract 1 because number of cells is 1
    // less than number of points)
    numCells +=
        (blkSize[ii][0] - 1) * (blkSize[ii][1] - 1) * (blkSize[ii][2] - 1);
  }
  cout << endl;

  // read each block and add it to the vector of plot3dBlocks
  auto tempDouble = 0.0;
  vector<plot3dBlock> mesh;
  mesh.reserve(numBlks);

  for (auto ii = 0; ii < numBlks; ii++) {
    multiArray3d<vector3d<double>> coordinates(blkSize[ii][0], blkSize[ii][1],
                                               blkSize[ii][2], 0);

    for (auto jj = 0; jj < coordinates.Size(); jj++) {
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      coordinates(jj)[0] = tempDouble / LRef;
    }
    for (auto jj = 0; jj < coordinates.Size(); jj++) {
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      coordinates(jj)[1] = tempDouble / LRef;
    }
    for (auto jj = 0; jj < coordinates.Size(); jj++) {
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      coordinates(jj)[2] = tempDouble / LRef;
    }

    // create single plot3dBlock and assign it appropriate location in vector
    mesh.push_back(plot3dBlock(coordinates));

    cout << "Block " << ii << " read" << endl;
  }

  cout << endl << "Grid file read" << endl;
  cout << "Total number of cells is " << numCells << endl;

  // close plot3d grid file
  fName.close();

  return mesh;
}


/* Member function to split a plot3dBlock along a plane defined by a direction
and an index.
*/
plot3dBlock plot3dBlock::Split(const string &dir, const int &ind) {
  // dir -- plane to split along, either i, j, or k
  // ind -- index (face) to split at (w/o counting ghost cells)

  const auto blk2 = plot3dBlock(coords_.Slice(dir, {ind, coords_.End(dir)}));
  (*this) = plot3dBlock(coords_.Slice(dir, {coords_.Start(dir), ind + 1}));
  return blk2;
}

/* Member function to join a plot3dBlock along a plane defined by a direction.
The calling instance will be the lower portion of the joined block,
and the input instance will be the upper portion of the joined block.
*/
void plot3dBlock::Join(const plot3dBlock &blk, const string &dir) {
  auto iTot = this->NumI();
  auto jTot = this->NumJ();
  auto kTot = this->NumK();

  if (dir == "i") {
    iTot += blk.NumI() - 1;
  } else if (dir == "j") {
    jTot += blk.NumJ() - 1;
  } else if (dir == "k") {
    kTot += blk.NumK() - 1;
  } else {
    cerr << "ERROR: Error in plot3dBlock::Join(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(EXIT_FAILURE);
  }

  plot3dBlock newBlk(iTot, jTot, kTot);

  newBlk.coords_.Insert(dir, {coords_.Start(dir), coords_.End(dir)}, coords_);
  // lower and upper splits share one plane of nodes
  newBlk.coords_.Insert(dir, {coords_.End(dir) - 1, newBlk.coords_.End(dir)},
                        blk.coords_);
  (*this) = newBlk;
}

double PyramidVolume(const vector3d<double> &p, const vector3d<double> &a,
                     const vector3d<double> &b, const vector3d<double> &c,
                     const vector3d<double> &d) {
  auto xp = 0.25 * ((a - p) + (b - p) + (c - p) + (d - p));
  // vectors along diagonal of base
  auto xac = c - a;
  auto xbd = d - b;
  return 1.0 / 6.0 * xp.DotProd(xac.CrossProd(xbd));
}
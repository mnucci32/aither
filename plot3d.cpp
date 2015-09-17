/*  An open source Navier-Stokes CFD solver.
    Copyright (C) 2015  Michael Nucci (michael.nucci@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
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

// plot3dBlock member function that calcualtes the volume of each cell
/*
All cells are assumed to be hexahedra. The 8 points that make up each hexahedron
are used to split the cell into 3 pyramids. The area of each pyramid is then
calculated and the volume of the 3 pyramids are summed to get the volume of the
hexahedra. This method is outlined in Hirsch.

Vp = Havg * (D1 (cross) D2))

The equation above shows how the volume of a pyramid (Vp) is calculated. Havg is
the average distance from the four base points to the top of the pyramid. D1 and
D2 are the diagonal distances of the four base points.
*/
multiArray3d<double> plot3dBlock::Volume() const {
  // Allocate multiArray3d to store cell volumes in
  multiArray3d<double> vol(this->NumI() - 1, this->NumJ() - 1,
                           this->NumK() - 1);

  // Loop over all cells
  for (int kk = 0; kk < vol.NumK(); kk++) {
    for (int jj = 0; jj < vol.NumJ(); jj++) {
      for (int ii = 0; ii < vol.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the i and j directions
        vector3d<double> botPortFore = coords_(ii + 1, jj + 1, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the i and k directions
        vector3d<double> topStarFore = coords_(ii + 1, jj, kk + 1);
        // up 1 in the j and k directions
        vector3d<double> topPortAft = coords_(ii, jj + 1, kk + 1);
        // up 1 in the i, j, and k directions
        vector3d<double> topPortFore = coords_(ii + 1, jj + 1, kk + 1);

        // Point of all three pyramids is located at the top, starboard, aft
        // corner of the cell
        // Calculate volume for pyramid 1 - quad face is bottom side
        // xp is average vector from 4 base points to peak of pyramid
        vector3d<double> xp = 0.25 * ((botStarAft - topStarAft) +
                                      (botPortAft - topStarAft) +
                                      (botStarFore - topStarAft) +
                                      (botPortFore - topStarAft));
        // vectors along diagonal of base
        vector3d<double> xac = botPortFore - botStarAft;
        vector3d<double> xbd = botStarFore - botPortAft;
        double pyramidVol = 1.0 / 6.0 * xp.DotProd(xac.CrossProd(xbd));

        // Calculate volume for pyramid2 - quad face is fore side
        xp = 0.25 * ((botStarFore - topStarAft) + (botPortFore - topStarAft) +
                     (topStarFore - topStarAft) + (topPortFore - topStarAft));
        xac = topPortFore - botStarFore;
        xbd = topStarFore - botPortFore;
        pyramidVol += 1.0 / 6.0 * xp.DotProd(xac.CrossProd(xbd));

        // Calculate volume for pyramid3 - quad face is port side
        xp = 0.25 * ((botPortFore - topStarAft) + (botPortAft - topStarAft) +
                     (topPortFore - topStarAft) + (topPortAft - topStarAft));
        xac = topPortFore - botPortAft;
        xbd = botPortFore - topPortAft;
        pyramidVol += 1.0 / 6.0 * xp.DotProd(xac.CrossProd(xbd));

        // Assign volume to appropriate location
        vol(ii, jj, kk) = pyramidVol;

        // Check for negative volumes
        if (pyramidVol <= 0) {
          cerr << "ERROR: Negative volume in PLOT3D block!!!" << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj
               << ", k-dim = " << kk << endl;
          exit(0);
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
  multiArray3d<vector3d<double>> centroid(this->NumI() - 1,
                                           this->NumJ() - 1,
                                           this->NumK() - 1);

  // loop over all cells
  for (int kk = 0; kk < centroid.NumK(); kk++) {
    for (int jj = 0; jj < centroid.NumJ(); jj++) {
      for (int ii = 0; ii < centroid.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the i and j directions
        vector3d<double> botPortFore = coords_(ii + 1, jj + 1, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the i and k directions
        vector3d<double> topStarFore = coords_(ii + 1, jj, kk + 1);
        // up 1 in the j and k directions
        vector3d<double> topPortAft = coords_(ii, jj + 1, kk + 1);
        // up 1 in the i, j, and k directions
        vector3d<double> topPortFore = coords_(ii + 1, jj + 1, kk + 1);

        // Calculate the centroid of the cell
        centroid(ii, jj, kk) = (0.125 * (botStarAft + botStarFore + botPortAft +
                                         botPortFore + topStarAft + topStarFore
                                         + topPortAft + topPortFore));
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
                                            this->NumK() - 1);

  // loop over all i-faces
  for (int kk = 0; kk < fArea.NumK(); kk++) {
    for (int jj = 0; jj < fArea.NumJ(); jj++) {
      for (int ii = 0; ii < fArea.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the j and k directions
        vector3d<double> topPortAft = coords_(ii, jj + 1, kk + 1);

        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        vector3d<double> xac = topPortAft - botStarAft;
        vector3d<double> xbd = botPortAft - topStarAft;

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
          exit(0);
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
  multiArray3d<vector3d<double>> fCenter(this->NumI(),
                                          this->NumJ() - 1,
                                          this->NumK() - 1);
  // loop over all i-faces
  for (int kk = 0; kk < fCenter.NumK(); kk++) {
    for (int jj = 0; jj < fCenter.NumJ(); jj++) {
      for (int ii = 0; ii < fCenter.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the j and k directions
        vector3d<double> topPortAft = coords_(ii, jj + 1, kk + 1);

        // Calculate face center by averaging four points that make up the face
        fCenter(ii, jj, kk) = 0.25 *
            (botStarAft + botPortAft + topStarAft + topPortAft);
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
                                            this->NumK() - 1);

  // loop over all j-faces
  for (int kk = 0; kk < fArea.NumK(); kk++) {
    for (int jj = 0; jj < fArea.NumJ(); jj++) {
      for (int ii = 0; ii < fArea.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the i and k directions
        vector3d<double> topStarFore = coords_(ii + 1, jj, kk + 1);

        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        vector3d<double> xac = topStarAft - botStarFore;
        vector3d<double> xbd = botStarAft - topStarFore;

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
          exit(0);
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
  multiArray3d<vector3d<double>> fCenter(this->NumI() - 1,
                                          this->NumJ(),
                                          this->NumK() - 1);
  // loop over all j-faces
  for (int kk = 0; kk < fCenter.NumK(); kk++) {
    for (int jj = 0; jj < fCenter.NumJ(); jj++) {
      for (int ii = 0; ii < fCenter.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the k direction
        vector3d<double> topStarAft = coords_(ii, jj, kk + 1);
        // up 1 in the i and k directions
        vector3d<double> topStarFore = coords_(ii + 1, jj, kk + 1);

        // Calculate face center by averaging the four points that make up the
        // face
        fCenter(ii, jj, kk) = 0.25 *
            (botStarAft + botStarFore + topStarAft + topStarFore);
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
  multiArray3d<unitVec3dMag<double>> fArea(this->NumI() - 1,
                                            this->NumJ() - 1,
                                            this->NumK());
  // loop over all k-faces
  for (int kk = 0; kk < fArea.NumK(); kk++) {
    for (int jj = 0; jj < fArea.NumJ(); jj++) {
      for (int ii = 0; ii < fArea.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the i and j directions
        vector3d<double> botPortFore = coords_(ii + 1, jj + 1, kk);

        // Calculate area for face by taking 1/2 of the cross product between
        // opposite diagonals
        // vectors from opposite corners of face
        vector3d<double> xac = botPortAft - botStarFore;
        vector3d<double> xbd = botPortFore - botStarAft;

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
          exit(0);
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
  multiArray3d<vector3d<double>> fCenter(this->NumI() - 1,
                                          this->NumJ() - 1,
                                          this->NumK());
  // loop over all k-faces
  for (int kk = 0; kk < fCenter.NumK(); kk++) {
    for (int jj = 0; jj < fCenter.NumJ(); jj++) {
      for (int ii = 0; ii < fCenter.NumI(); ii++) {
        // baseline location
        vector3d<double> botStarAft = coords_(ii, jj, kk);
        // up 1 in the i-direction
        vector3d<double> botStarFore = coords_(ii + 1, jj, kk);
        // up 1 in the j-direction
        vector3d<double> botPortAft = coords_(ii, jj + 1, kk);
        // up 1 in the i and j directions
        vector3d<double> botPortFore = coords_(ii + 1, jj + 1, kk);

        // Calculate face center by averaging four points that make up cell face
        fCenter(ii, jj, kk) = 0.25 *
            (botStarAft + botStarFore + botPortAft + botPortFore);
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
  string readName = gridName + fPostfix;
  fName.open(readName.c_str(), ios::in | ios::binary);

  // check to see if file opened correctly
  if (fName.fail()) {
    cerr << "ERROR: Error in plot3d.cpp:ReadP3dGrid(). Grid file " << readName
         << " did not open correctly!!!" << endl;
    exit(0);
  }

  // read the number of plot3d blocks in the file
  cout << "Reading grid file..." << endl << endl;
  int numBlks;
  fName.read(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));
  cout << "Number of blocks: " << numBlks << endl << endl;

  // read the number of i, j, k coordinates in each plot3d block
  cout << "Size of each block is..." << endl;
  vector<int> vecI(numBlks, 0);
  vector<int> vecJ(numBlks, 0);
  vector<int> vecK(numBlks, 0);
  int tempInt = 0;
  numCells = 0;

  // loop over all blocks and fill i, j, k vectors with block sizes
  for (int ii = 0; ii < numBlks; ii++) {
    cout << "Block Number: " << ii << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    vecI[ii] = tempInt;
    cout << "I-DIM: " << tempInt << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    vecJ[ii] = tempInt;
    cout << "J-DIM: " << tempInt << "     ";
    fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
    vecK[ii] = tempInt;
    cout << "K-DIM: " << tempInt << endl;

    // calculate total number of cells (subtract 1 because number of cells is 1
    // less than number of points)
    numCells += (vecI[ii] - 1) * (vecJ[ii] - 1) * (vecK[ii] - 1);
  }
  cout << endl;

  // read each block and add it to the vector of plot3dBlocks
  double tempDouble = 0;
  vector<plot3dBlock> mesh;
  mesh.reserve(numBlks);

  for (int ii = 0; ii < numBlks; ii++) {
    multiArray3d<vector3d<double>> coordinates(vecI[ii], vecJ[ii], vecK[ii]);

    for (int jj = 0; jj < coordinates.Size(); jj++) {
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      coordinates(jj)[0] = tempDouble / LRef;
    }
    for (int jj = 0; jj < coordinates.Size(); jj++) {
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      coordinates(jj)[1] = tempDouble / LRef;
    }
    for (int jj = 0; jj < coordinates.Size(); jj++) {
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

// function to reorder the block by hyperplanes
/*A hyperplane is a plane of i+j+k=constant within an individual block. The
LUSGS solver must sweep along these hyperplanes to avoid
calculating a flux jacobian. Ex. The solver must visit all points on hyperplane
1 before visiting any points on hyperplane 2.
*/
vector<vector3d<int>> HyperplaneReorder(const int &imax, const int &jmax,
                                         const int &kmax) {
  int numPlanes =
      imax + jmax + kmax - 2;  // total number of hyperplanes in a given block
  vector<vector3d<int>> reorder;
  reorder.reserve(imax * jmax * kmax);

  int count = 0;
  for (int pp = 0; pp < numPlanes; pp++) {
    for (int kk = 0; kk < kmax; kk++) {
      for (int jj = 0; jj < jmax; jj++) {
        for (int ii = 0; ii < imax; ii++) {
          if (ii + jj + kk == pp) {  // if sum of ii, jj, and kk equals pp than
                                     // point is on hyperplane pp
            reorder[count] = vector3d<int>(ii, jj, kk);
            count++;
          }
        }
      }
    }
  }

  return reorder;
}

/* Member function to split a plot3dBlock along a plane defined by a direction
and an index.
*/
void plot3dBlock::Split(const string &dir, const int &ind, plot3dBlock &blk1,
                        plot3dBlock &blk2) const {
  // dir -- plane to split along, either i, j, or k
  // ind -- index (face) to split at (w/o counting ghost cells)

  if (dir == "i") {  // split along i-plane
    blk1 = plot3dBlock(coords_.Slice(0, ind,
                                     0, coords_.NumJ() - 1,
                                     0, coords_.NumK() - 1));
    blk2 = plot3dBlock(coords_.Slice(ind, coords_.NumI() - 1,
                                     0, coords_.NumJ() - 1,
                                     0, coords_.NumK() - 1));
  } else if (dir == "j") {  // split along j-plane
    blk1 = plot3dBlock(coords_.Slice(0, coords_.NumI() - 1,
                                     0, ind,
                                     0, coords_.NumK() - 1));
    blk2 = plot3dBlock(coords_.Slice(0, coords_.NumI() - 1,
                                     ind, coords_.NumJ() - 1,
                                     0, coords_.NumK() - 1));
  } else if (dir == "k") {  // split along k-plane
    blk1 = plot3dBlock(coords_.Slice(0, coords_.NumI() - 1,
                                     0, coords_.NumJ() - 1,
                                     0, ind));
    blk2 = plot3dBlock(coords_.Slice(0, coords_.NumI() - 1,
                                     0, coords_.NumJ() - 1,
                                     ind, coords_.NumK() - 1));
  } else {
    cerr << "ERROR: Error in plot3dBlock::Split(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }
}

/* Member function to join a plot3dBlock along a plane defined by a direction.
The calling instance will be the lower portion of the joined block,
and the input instance will be the upper portion of the joined block.
*/
void plot3dBlock::Join(const plot3dBlock &blk, const string &dir) {
  if (dir == "i") {
    int newNumI = this->NumI() + blk.NumI() - 1;

    plot3dBlock newBlk(newNumI, this->NumJ(), this->NumK());

    newBlk.coords_.Insert(0, this->NumI(), 0, this->NumJ(), 0,
                          this->NumK(), coords_);

    newBlk.coords_.Insert(this->NumI(), newNumI, 0, this->NumJ(), 0,
                          this->NumK(), blk.coords_);
    (*this) = newBlk;
  } else if (dir == "j") {
    int newNumJ = this->NumJ() + blk.NumJ() - 1;

    plot3dBlock newBlk(this->NumI(), newNumJ, this->NumK());

    newBlk.coords_.Insert(0, this->NumI(), 0, this->NumJ(), 0,
                          this->NumK(), coords_);

    newBlk.coords_.Insert(0, this->NumI(), this->NumJ(), newNumJ, 0,
                          this->NumK(), blk.coords_);
    (*this) = newBlk;
  } else if (dir == "k") {
    int newNumK = this->NumK() + blk.NumK() - 1;

    plot3dBlock newBlk(this->NumI(), this->NumJ(), newNumK);

    newBlk.coords_.Insert(0, this->NumI(), 0, this->NumJ(), 0,
                          this->NumK(), coords_);

    newBlk.coords_.Insert(0, this->NumI(), 0, this->NumJ(), this->NumK(),
                          newNumK, blk.coords_);
    (*this) = newBlk;
  } else {
    cerr << "ERROR: Error in plot3dBlock::Join(). Direction " << dir
         << " is not recognized! Choose either i, j, or k." << endl;
    exit(0);
  }
}

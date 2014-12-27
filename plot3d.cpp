#include "plot3d.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ios;

//----------------------------------------------------------------------------------------------------------------//
//plot 3d block constructor, member functions
//constructor -- assign passed variables to create plot3d block
plot3dBlock::plot3dBlock( int i, int j, int k, vector<double> &xCoord, vector<double> &yCoord, vector<double> &zCoord){
  numi = i;
  numj = j;
  numk = k;
  x = xCoord;
  y = yCoord;
  z = zCoord;
}
//constructor with no arguments
plot3dBlock::plot3dBlock(){
  numi = 0;
  numj = 0;
  numk = 0;
  vector<double> xCoord(1,0);
  vector<double> yCoord(1,0);
  vector<double> zCoord(1,0);
  x = xCoord;
  y = yCoord;
  z = zCoord;
}

//plot3dBlock member function that calcualtes the volume of each cell
/*
All cells are assumed to be hexahedra. The 8 points that make up each hexahedron are used to 
split the cell into 3 pyramids. The area of each pyramid is then calculated and the volume of the 
3 pyramids are summed to get the volume of the hexahedra. This method is outlined in Hirsch.

Vp = Havg * (D1 (cross) D2))

The equation above shows how the volume of a pyramid (Vp) is calculated. Havg is the average distance from 
the four base points to the top of the pyramid. D1 and D2 are the diagonal distances of the four base points.
*/
const vector<double> plot3dBlock::Volume() const {

  int len = (numi-1)*(numj-1)*(numk-1);  //number of cells in the block

  vector<double> vol(len,0);             //initially assign a volume of 0

  //each hex cell is broken up into 3 pyramids and the volume of each pyramid is calculated
  double pyramid1 = 0;
  double pyramid2 = 0;
  double pyramid3 = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botStarFore(0,0,0);
  vector3d<double> botPortFore(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topStarFore(0,0,0);
  vector3d<double> topPortFore(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  //vectors of the components that go into calculating the volume of each pyramid
  vector3d<double> xp(0,0,0);            //average vector from 4 base points to peak point
  vector3d<double> xac(0,0,0);           //vector along diagonal of base
  vector3d<double> xbd(0,0,0);           //vector along diagonal of base

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in the i-direction
        botStarFore.SetX(x[loc+1]);
        botStarFore.SetY(y[loc+1]);
        botStarFore.SetZ(z[loc+1]);

	//up 1 in the j-direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);

	//up 1 in the i and j directions
        botPortFore.SetX(x[loc+numi+1]);
        botPortFore.SetY(y[loc+numi+1]);
        botPortFore.SetZ(z[loc+numi+1]);
              
	//up 1 in the k direction
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

	//up 1 in the i and k directions
        topStarFore.SetX(x[loc+numi*numj+1]);
        topStarFore.SetY(y[loc+numi*numj+1]);
        topStarFore.SetZ(z[loc+numi*numj+1]);

	//up 1 in the j and k directions
        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

	//up 1 in the i, j, and k directions
        topPortFore.SetX(x[loc+numi+numi*numj+1]);
        topPortFore.SetY(y[loc+numi+numi*numj+1]);
        topPortFore.SetZ(z[loc+numi+numi*numj+1]);

        //Point of all three pyramids is located at the top, starboard, aft corner of the cell
        //Calculate volume for pyramid 1 - quad face is bottom side
        xp = 0.25*((botStarAft-topStarAft) + (botPortAft-topStarAft) + (botStarFore-topStarAft) + (botPortFore-topStarAft));
        xac = botPortFore - botStarAft;
        xbd = botStarFore - botPortAft;
        pyramid1 = 1.0/6.0 * xp.DotProd(xac.CrossProd(xbd));

        //Calculate volume for pyramid2 - quad face is fore side
        xp = 0.25*((botStarFore-topStarAft) + (botPortFore-topStarAft) + (topStarFore-topStarAft) + (topPortFore-topStarAft));
        xac = topPortFore - botStarFore;
        xbd = topStarFore - botPortFore;
        pyramid2 = 1.0/6.0 * xp.DotProd(xac.CrossProd(xbd));

        //Calculate volume for pyramid3 - quad face is port side
        xp = 0.25*((botPortFore-topStarAft) + (botPortAft-topStarAft) + (topPortFore-topStarAft) + (topPortAft-topStarAft));
        xac = topPortFore - botPortAft;
        xbd = botPortFore - topPortAft;
        pyramid3 = 1.0/6.0 * xp.DotProd(xac.CrossProd(xbd));

        vol[index] = pyramid1 + pyramid2 + pyramid3;

        if (vol[index] == 0){
          cerr << "ERROR: Negative volume in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          exit(1);
	}
        index++;

      }
    }
  }

  return vol;
}

//plot3dBlock member function that calcualtes the centroid of each cell
//the centroid of the hexahedron is the average of the 8 points that define it
const vector<vector3d<double> > plot3dBlock::Centroid() const {

  int len = (numi-1)*(numj-1)*(numk-1);                 //number of cells in the block
  vector<vector3d<double> > centroid(len);             //preallocate centroid vector

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botStarFore(0,0,0);
  vector3d<double> botPortFore(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topStarFore(0,0,0);
  vector3d<double> topPortFore(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in the i direction
        botStarFore.SetX(x[loc+1]);
        botStarFore.SetY(y[loc+1]);
        botStarFore.SetZ(z[loc+1]);

	//up 1 in the j direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);

	//up 1 in the i and j directions
        botPortFore.SetX(x[loc+numi+1]);
        botPortFore.SetY(y[loc+numi+1]);
        botPortFore.SetZ(z[loc+numi+1]);
              
	//up 1 in the k direction
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

	//up 1 in the i and k directions
        topStarFore.SetX(x[loc+numi*numj+1]);
        topStarFore.SetY(y[loc+numi*numj+1]);
        topStarFore.SetZ(z[loc+numi*numj+1]);

	//up 1 in the j and k directions
        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

	//up 1 in the i, j, and k directions
        topPortFore.SetX(x[loc+numi+numi*numj+1]);
        topPortFore.SetY(y[loc+numi+numi*numj+1]);
        topPortFore.SetZ(z[loc+numi+numi*numj+1]);

	//Calculate the centroid of the cell
        centroid[index] = 0.125 * (botStarAft + botStarFore + botPortAft + botPortFore + topStarAft + topStarFore + topPortAft + topPortFore);
	index++;
      }
    }
  }

  return centroid;
}

//plot3dBlock member function that calcualtes the area of each face normal to the i-direction
//the area is calculated as half the cross product of the diagonal vectors of the 4-sided face
/*

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B to C. The normal vector
points in the direction of increasing i.

*/
const vector<vector3d<double> > plot3dBlock::FaceAreaI() const {

  int len = numi * (numj - 1) * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi; ii++){

        //assign coordinates to 4 points that make up the face being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in j direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
	//up 1 in k direction
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

	//up 1 in j and k directions
        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        vector3d<double> xac = topPortAft - botStarAft; //vector from opposite corners of face
        vector3d<double> xbd = botPortAft - topStarAft; //vector from opposite corners of face

        fArea[index] = 0.5 * xbd.CrossProd(xac);     //area vector is calculated so that normal points nominally in direction of increasing i-coordinate

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(0);
	}
	index++;

      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the center of each face normal to the i-direction
//the face center is calculated as the average of the 4 points that comprise it
const vector<vector3d<double> > plot3dBlock::FaceCenterI() const {

  int len = numi * (numj - 1) * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi; ii++){

        //assign coordinates to 4 points that make up the face being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in j direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
	//up 1 in k direction
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

	//up 1 in j and k directions
        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

        //Calculate face center by averaging four points that make up the face
        fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
	index++;
      }
    }
  }

  return fCenter;
}

//plot3dBlock member function that calcualtes the area of each face normal to the j-direction
//the area is calculated as half the cross product of the diagonal vectors of the 4-sided face
/*

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B to C. The normal
points in the direction of increasing j.

*/
const vector<vector3d<double> > plot3dBlock::FaceAreaJ() const {

  int len = (numi -1) * numj * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 4 points that make up the face being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//up 1 in i direction
        botStarAft.SetX(x[loc+1]);
        botStarAft.SetY(y[loc+1]);
        botStarAft.SetZ(z[loc+1]);

	//baseline location
        botPortAft.SetX(x[loc]);
        botPortAft.SetY(y[loc]);
        botPortAft.SetZ(z[loc]);
              
	//up 1 in i and k directions
        topStarAft.SetX(x[loc+numi*numj+1]);
        topStarAft.SetY(y[loc+numi*numj+1]);
        topStarAft.SetZ(z[loc+numi*numj+1]);

	//up 1 in k direction
        topPortAft.SetX(x[loc+numi*numj]);
        topPortAft.SetY(y[loc+numi*numj]);
        topPortAft.SetZ(z[loc+numi*numj]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        vector3d<double> xac = topPortAft - botStarAft; //vector from opposite corners of face
        vector3d<double> xbd = botPortAft - topStarAft; //vector from opposite corners of face

        fArea[index] = 0.5 * xbd.CrossProd(xac);       //area vector is calculated so that normal nominally points in direction of increasing j-coordinate

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(0);
	}
	index++;
      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the area of each face normal to the j-direction
//the face center is calculated as the average of the 4 points that comprise it
const vector<vector3d<double> > plot3dBlock::FaceCenterJ() const {

  int len = (numi -1) * numj * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk-1; kk++){
    for ( int jj = 0; jj < numj; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 4 points that make up the face being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//up 1 in i direction
        botStarAft.SetX(x[loc+1]);
        botStarAft.SetY(y[loc+1]);
        botStarAft.SetZ(z[loc+1]);

	//baseline location
        botPortAft.SetX(x[loc]);
        botPortAft.SetY(y[loc]);
        botPortAft.SetZ(z[loc]);
              
	//up 1 in i and k directions
        topStarAft.SetX(x[loc+numi*numj+1]);
        topStarAft.SetY(y[loc+numi*numj+1]);
        topStarAft.SetZ(z[loc+numi*numj+1]);

	//up 1 in k direction
        topPortAft.SetX(x[loc+numi*numj]);
        topPortAft.SetY(y[loc+numi*numj]);
        topPortAft.SetZ(z[loc+numi*numj]);

        //Calculate face center by averaging the four points that make up the face
	fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
	index++;
      }
    }
  }

  return fCenter;
}

//plot3dBlock member function that calcualtes the area of each face normal to the k-direction
//the area is calculated as half the cross product of the diagonal vectors of the 4-sided face
/*

  A---------------B
  |               |
  |               |
  |               |
  |               |
  C---------------D

A = 0.5 * rAD (cross) rCB

In the equation above rAD is the vector from D to A and rCD is the vector from B to C. The normal vector
points in the direction of increasing k.

*/
const vector<vector3d<double> > plot3dBlock::FaceAreaK() const {

  int len = (numi - 1) * (numj - 1) * numk;         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 4 points that make up the face being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in j direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
	//up 1 in i direction
        topStarAft.SetX(x[loc+1]);
        topStarAft.SetY(y[loc+1]);
        topStarAft.SetZ(z[loc+1]);

	//up 1 in i and j directions
        topPortAft.SetX(x[loc+numi+1]);
        topPortAft.SetY(y[loc+numi+1]);
        topPortAft.SetZ(z[loc+numi+1]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        vector3d<double> xac = topPortAft - botStarAft; //vector from opposite corners of face
        vector3d<double> xbd = botPortAft - topStarAft; //vector from opposite corners of face

        fArea[index] = 0.5 * xac.CrossProd(xbd);   //area vector is calculated so that normal nominally points in direction of increasing k-coordinate

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(0);
	}
	index++;
      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the area of each face normal to the k-direction
//the face center is calculated as the average of the 4 points that comprise it
const vector<vector3d<double> > plot3dBlock::FaceCenterK() const {

  int len = (numi - 1) * (numj - 1) * numk;         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter

  for ( int kk = 0; kk < numk; kk++){
    for ( int jj = 0; jj < numj-1; jj++){
      for ( int ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        int loc = GetLoc1D(ii, jj, kk, numi, numj); 

	//baseline location
        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

	//up 1 in j direction
        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
	//up 1 in i direction
        topStarAft.SetX(x[loc+1]);
        topStarAft.SetY(y[loc+1]);
        topStarAft.SetZ(z[loc+1]);

	//up 1 in i and j directions
        topPortAft.SetX(x[loc+numi+1]);
        topPortAft.SetY(y[loc+numi+1]);
        topPortAft.SetZ(z[loc+numi+1]);

        //Calculate face center by averaging four points that make up cell face
        fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
	index++;
      }
    }
  }

  return fCenter;
}

//---------------------------------------------------------------------------------------------------------------//
//plot 3d solution block constructor
//constructor -- create a plot3d solution block when passed a plot3dblock
plot3dQBlock::plot3dQBlock(const  plot3dBlock &singleBlock ):q1(singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK()),
                                                             q2(singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK()),
                                                             q3(singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK()),
                                                             q4(singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK()),
                                                             q5(singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK()){
             
  numi = singleBlock.NumI();
  numj = singleBlock.NumJ();
  numk = singleBlock.NumK();
  mach = 0.0;
  alpha = 0.0;
  Re = 0.0;
  ti = 0.0;
  int ii = 0;
  int length = singleBlock.NumI()*singleBlock.NumJ()*singleBlock.NumK();
  for(ii = 0; ii < length; ii++){
    q1[ii] = 1.0;
    q2[ii] = 2.0;
    q3[ii] = 3.0;
    q4[ii] = 4.0;
    q5[ii] = 5.0;
  }
}

//constructor for plot3dQBlock with no arguments
plot3dQBlock::plot3dQBlock( ) {
  numi = 1;
  numj = 1;
  numk = 1;
  mach = 0.0;
  alpha = 0.0;
  Re = 0.0;
  ti = 0.0;

  vector<double> cons1(1,0);
  vector<double> cons2(1,0);
  vector<double> cons3(1,0);
  vector<double> cons4(1,0);
  vector<double> cons5(1,0);

  q1 = cons1;
  q2 = cons2;
  q3 = cons3;
  q4 = cons4;
  q5 = cons5;
}

//---------------------------------------------------------------------------------------------------------------//
//plot 3d mesh constructor, member functions
//constructor -- assign the passed single plot3dBlock into the first position in the plot3dMesh
plot3dMesh::plot3dMesh(const plot3dBlock &singleBlock ){
  blocks.assign(1,singleBlock);
}
//constructor with no arguements
plot3dMesh::plot3dMesh(){
  plot3dBlock singleBlock;
  blocks.assign(1,singleBlock);
}
//constructor with arugument for number of blocks
plot3dMesh::plot3dMesh(int num){
  plot3dBlock singleBlock;
  blocks.assign(num,singleBlock);
}

//plot3dMesh member function that adds a single plot3d block to the end of the plot3dMesh
void plot3dMesh::AddP3dBlock(const plot3dBlock &singleBlock) {
  return blocks.push_back(singleBlock);
}

//plot3dMesh member function that replaces a single plot3d block with the block that is passed
void plot3dMesh::ReplaceBlock(int index, const plot3dBlock &singleBlock) {
  blocks[index]=singleBlock;
  return;
}

//---------------------------------------------------------------------------------------------------------------//
//function to read in a plot3d grid and assign it to a plot3dMesh data type
plot3dMesh ReadP3dGrid(const string &gridName, double &numCells) {

  //open binary plot3d grid file 
  ifstream fName;
  string fPostfix = ".xyz";
  string readName = gridName + fPostfix;
  fName.open(readName.c_str(), ios::in | ios::binary);

  //check to see if file opened correctly
  if (fName.fail()) {
    cerr << "ERROR: Grid file " << readName << " did not open correctly!!!" << endl;
    exit(0);
  }

  //read the number of plot3d blocks in the file
  cout << "Reading grid file..." << endl << endl;
  int numBlks;
  fName.read(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));
  cout << "Number of blocks: " << numBlks << endl << endl;

  //read the number of i, j, k coordinates in each plot3d block
  cout << "Size of each block is..." << endl;
  vector<int> vecI, vecJ, vecK;
  int tempInt = 0;
  int ii = 0;
  numCells = 0;
  for ( ii=0; ii < numBlks; ii++){
  
    if (ii==0){
      cout << "Block Number: " << ii << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecI.assign(1,tempInt);
      cout << "I-DIM: " << tempInt << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecJ.assign(1,tempInt);
      cout << "J-DIM: " << tempInt << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecK.assign(1,tempInt);
      cout << "K-DIM: " << tempInt << endl;
    }
    else {    
      cout << "Block Number: " << ii << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecI.push_back(tempInt);
      cout << "I-DIM: " << tempInt << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecJ.push_back(tempInt);
      cout << "J-DIM: " << tempInt << "     ";
      fName.read(reinterpret_cast<char *>(&tempInt), sizeof(tempInt));
      vecK.push_back(tempInt);
      cout << "K-DIM: " << tempInt << endl;
    }

    //calculate total number of cells
    numCells += (vecI[ii] - 1) * (vecJ[ii] - 1) * (vecK[ii] - 1);

  }
  cout << endl;

  //read each block and add it to the plot3dMesh
  double tempDouble=0;
  vector<double> tempX, tempY, tempZ;
  plot3dMesh mesh(numBlks);

  for ( ii=0; ii < numBlks; ii++){
  
    int size = vecI[ii]*vecJ[ii]*vecK[ii];
    tempX.assign(size,0);
    tempY.assign(size,0);
    tempZ.assign(size,0);

    int jj = 0;

    for (jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempX[jj] = tempDouble;
    }
    for (jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempY[jj] = tempDouble;
    }
    for (jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempZ[jj] = tempDouble;
    }

    plot3dBlock singleBlock(vecI[ii],vecJ[ii],vecK[ii],tempX,tempY,tempZ);

    mesh.ReplaceBlock(ii,singleBlock);

    cout << "Block " << ii << " read" << endl;

  }

  cout << endl << "Grid file read" << endl;
  cout << "Total number of cells is " << numCells << endl;

  //close plot3d grid file
  fName.close();

  return mesh;
}


//---------------------------------------------------------------------------------------------------------------//
//plot 3d solution mesh constructor, member functions
//constructor -- create a plot3d solution mesh from a passed plot3d mesh
plot3dQMesh::plot3dQMesh(const plot3dMesh &mesh ):qblocks(mesh.NumBlocks()){
  int ii = 0;
  for (ii = 0; ii< mesh.NumBlocks(); ii++){
    cout << "writing solution block " << ii << endl;
    qblocks[ii] = plot3dQBlock(mesh.Blocks(ii));
  }
}

//member function to write data to a file
void plot3dQMesh::WriteData() const {

  string solName = "sol.q";

  ofstream outFile;
  outFile.precision(16);
  outFile.open(solName.c_str(), ios::out | ios::binary);

  //write number of blocks to file
  int numBlks = (*this).QBlocks().size();
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  //write i, j, k dimension for each block
  int ll = 0;
  int dumInt = 0;
  for ( ll=0; ll < numBlks; ll++ ){
    dumInt = (*this).QBlocks()[ll].NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = (*this).QBlocks()[ll].NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = (*this).QBlocks()[ll].NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
  }


  int q = 0;
  int ii = 0;
  int jj = 0;
  int kk = 0;
  double mach = 0.0;
  double alpha = 0.0;
  double Re = 0.0;
  double ti = 0.0;
  int maxi = 0;
  int maxj = 0;
  int maxk = 0;
  double dumDouble = 0.0;


  for ( ll=0; ll < numBlks; ll++ ){
      mach = (*this).QBlocks()[ll].Mach();
      alpha = (*this).QBlocks()[ll].Alpha();
      Re = (*this).QBlocks()[ll].Reynolds();
      ti = (*this).QBlocks()[ll].Ti();

      outFile.write(reinterpret_cast<char *>(&mach), sizeof(mach));
      outFile.write(reinterpret_cast<char *>(&alpha), sizeof(alpha));
      outFile.write(reinterpret_cast<char *>(&Re), sizeof(Re));
      outFile.write(reinterpret_cast<char *>(&ti), sizeof(ti));

      maxi = (*this).QBlocks()[ll].NumI();
      maxj = (*this).QBlocks()[ll].NumJ();
      maxk = (*this).QBlocks()[ll].NumK();
  for ( q=0; q < 5; q++ ) {  
    for ( kk=0; kk < maxk; kk++ ){
        for ( jj = 0; jj < maxj; jj++ ){
          for ( ii = 0; ii < maxi; ii++ ){
            //write out q-data
            dumDouble = q + 1.0;
            outFile.write(reinterpret_cast<char *>(&dumDouble), sizeof(dumDouble));
          }
        }
      }
    }
  }

  outFile.close();

}


//function to take in imax, jmax, kmax, and vector index and return the i, j, k indices of that cell
vector3d<int> GetIJK(const int &imax, const int &jmax, const int &kmax, const int &index){

  vector3d<int> ijk;
  int ii = 0;
  int jj = 0;
  int kk = 0;
  int loc = 0;
  int breakFlag = 0;
  for ( kk = 0; kk < kmax; kk++ ){
    for ( jj = 0; jj < jmax; jj++ ){
      for ( ii = 0; ii < imax; ii++ ){
	loc = ii + jj*imax + kk*imax*jmax;
	if ( loc == index ){
	  ijk.SetX(ii);
	  ijk.SetY(jj);
	  ijk.SetZ(kk);
	  breakFlag = 1;
	  break;
	}
      }
      if (breakFlag) break;
    }
    if (breakFlag) break;
  }


  return ijk;
}

//------------------------------------------------------------------------------------------------------------------
//functions to take in i, j, k indexes of a cell and return the 1D index of the face belonging to that cell
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the upper i-direction
int GetUpperFaceI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i+1+(num-1) + j*(imax+1) + k*(imax+1)*jmax; 
}
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the lower i-direction
int GetLowerFaceI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i+1-num + j*(imax+1) + k*(imax+1)*jmax; 
}
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the upper j-direction
int GetUpperFaceJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j+1+(num-1))*imax + k*imax*(jmax+1);
}
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the lower j-direction
int GetLowerFaceJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j+1-num)*imax + k*imax*(jmax+1);
}
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the upper k-direction
int GetUpperFaceK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k+1+(num-1))*imax*jmax;
}
//function to take in i, j, k indexes of a cell and return the 1D index of the face in the lower k-direction
int GetLowerFaceK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k+1-num)*imax*jmax;
}

//-------------------------------------------------------------------------------------------------------------------
//functions to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the lower i face
int GetNeighborLowI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i-num + j*imax + k*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the upper i face
int GetNeighborUpI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i+num + j*imax + k*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the lower j face
int GetNeighborLowJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j-num)*imax + k*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the upper j face
int GetNeighborUpJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j+num)*imax + k*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the lower k face
int GetNeighborLowK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k-num)*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the upper k face
int GetNeighborUpK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k+num)*imax*jmax;;
}

//------------------------------------------------------------------------------------------------------------------------------
//functions to take in i, j, k, imax, jmax of a face and return the 1D index of the neighbor cell
//function to take in i, j, k indexes of a k-face and return the 1D index of the neighbor cell in the lower k direction
int GetCellFromFaceLowerK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k-1-(num-1))*imax*jmax;
}
//function to take in i, j, k indexes of a k-face and return the 1D index of the neighbor cell in the upper k direction
int GetCellFromFaceUpperK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + j*imax + (k-1+num)*imax*jmax;
}
//function to take in i, j, k indexes of a j-face and return the 1D index of the neighbor cell in the lower j direction
int GetCellFromFaceLowerJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j-1-(num-1))*imax + k*imax*(jmax-1);
}
//function to take in i, j, k indexes of a j-face and return the 1D index of the neighbor cell in the upper j direction
int GetCellFromFaceUpperJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j-1+num)*imax + k*imax*(jmax-1);
}
//function to take in i, j, k indexes of a i-face and return the 1D index of the neighbor cell in the lower i direction
int GetCellFromFaceLowerI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return (i-1-(num-1)) + j*(imax-1) + k*(imax-1)*jmax;
}
//function to take in i, j, k indexes of a i-face and return the 1D index of the neighbor cell in the upper i direction
int GetCellFromFaceUpperI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return (i-1+num) + j*(imax-1) + k*(imax-1)*jmax;
}


//function to take in i, j, k, imax, jmax and return corresponding location inside 1D array
int GetLoc1D(const int &i, const int &j, const int &k, const int &imax, const int &jmax){
  return i + j*imax + k*imax*jmax;
}

//function to take in the 1D index of a cell and return the 1D index of the upper i diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagUpperFromMainI(const int &main){
  return main - 1;
}
//function to take in the 1D index of a cell and return the 1D index of the lower i diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagLowerFromMainI(const int &main){
  return main + 1;
}
//function to take in the 1D index of a cell and return the 1D index of the upper j diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagUpperFromMainJ(const int &main, const int &imax){
  return main - imax;
}
//function to take in the 1D index of a cell and return the 1D index of the lower j diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagLowerFromMainJ(const int &main, const int &imax){
  return main + imax;
}
//function to take in the 1D index of a cell and return the 1D index of the upper k diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagUpperFromMainK(const int &main, const int &imax, const int &jmax){
  return main - imax*jmax;
}
//function to take in the 1D index of a cell and return the 1D index of the lower k diagonal to place the flux jacobian for the implicit matrix
int GetMatrixDiagLowerFromMainK(const int &main, const int &imax, const int &jmax){
  return main + imax*jmax;
}


//function to take in 1D index of a cell and return the 1D index of the upper i diagonal (in the same equation row) for the implicit matrix
int GetDiagPosUpperI(const int &main){
  return main + 1;
}
//function to take in 1D index of a cell and return the 1D index of the lower i diagonal (in the same equation row) for the implicit matrix
int GetDiagPosLowerI(const int &main){
  return main - 1;
}
//function to take in 1D index of a cell and return the 1D index of the upper j diagonal (in the same equation row) for the implicit matrix
int GetDiagPosUpperJ(const int &main, const int &imax){
  return main + imax;
}
//function to take in 1D index of a cell and return the 1D index of the lower j diagonal (in the same equation row) for the implicit matrix
int GetDiagPosLowerJ(const int &main, const int &imax){
  return main - imax;
}
//function to take in 1D index of a cell and return the 1D index of the upper k diagonal (in the same equation row) for the implicit matrix
int GetDiagPosUpperK(const int &main, const int &imax, const int &jmax){
  return main + imax*jmax;
}
//function to take in 1D index of a cell and return the 1D index of the lower k diagonal (in the same equation row) for the implicit matrix
int GetDiagPosLowerK(const int &main, const int &imax, const int &jmax){
  return main - imax*jmax;
}


//function to take in j, j, k indexes of a cell and the diagonal label, and determine if the matrix should have data at the cell
bool IsMatrixData(const int &i, const int &j, const int &k, const int &imax, const int &jmax, const int &kmax, const string &label){

  bool isData;

  if (label == "il"){
    if (i == 0){
      isData = false;  //at lower i boundary, there is no lower i diagonal
    }
    else{
      isData = true;
    }
  }
  else if (label == "iu"){
    if (i == imax-1){
      isData = false;  //at upper i boundary, there is no upper i diagonal
    }
    else{
      isData = true;
    }
  }
  else if (label == "jl"){
    if (j == 0){
      isData = false;  //at lower j boundary, there is no lower j diagonal
    }
    else{
      isData = true;
    }
  }
  else if (label == "ju"){
    if (j == jmax-1){
      isData = false;  //at upper j boundary, there is no upper j diagonal
    }
    else{
      isData = true;
    }
  }
  else if (label == "kl"){
    if (k == 0){
      isData = false;  //at lower k boundary, there is no lower k diagonal
    }
    else{
      isData = true;
    }
  }
  else if (label == "ku"){
    if (k == kmax-1){
      isData = false;  //at upper k boundary, there is no upper k diagonal
    }
    else{
      isData = true;
    }
  }
  else{
    cerr << "ERROR: Cannot determine position in matrix because block surface " << label << " is not recognized!" << endl;
    isData = false;
    exit(0);
  }
  return isData;
}

//function to reorder the block by hyperplanes
vector<vector3d<int> > HyperplaneReorder( const int &imax, const int &jmax, const int &kmax){

  int numPlanes = imax + jmax + kmax - 2;
  vector<vector3d<int> > reorder(imax*jmax*kmax);

  int count = 0;
  for ( int pp = 0; pp < numPlanes; pp++ ){
    for ( int kk = 0; kk < kmax; kk++ ){
      for ( int jj = 0; jj < jmax; jj++ ){
	for ( int ii = 0; ii < imax; ii++ ){
	  if ( ii+jj+kk == pp){
	    vector3d<int> loc(ii,jj,kk);
	    reorder[count] = loc;
	    count++;
	  }
	}
      }
    }
  }

  return reorder;
}

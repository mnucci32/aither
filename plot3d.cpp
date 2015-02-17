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
plot3dBlock::plot3dBlock( const int &i, const int &j, const int &k, const vector<double> &xCoord, const vector<double> &yCoord, const vector<double> &zCoord){
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

  vector<double> vol;             //declare vector and reserve necessary space
  vol.reserve(len);

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

        //vol[index] = pyramid1 + pyramid2 + pyramid3;
        vol.push_back( pyramid1 + pyramid2 + pyramid3 );

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
  vector<vector3d<double> > centroid;                   //declare and preallocate centroid vector
  centroid.reserve(len);

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botStarFore(0,0,0);
  vector3d<double> botPortFore(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topStarFore(0,0,0);
  vector3d<double> topPortFore(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  //int index = 0;                         //cell index counter

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
        //centroid[index] = 0.125 * (botStarAft + botStarFore + botPortAft + botPortFore + topStarAft + topStarFore + topPortAft + topPortFore);
        centroid.push_back( 0.125 * (botStarAft + botStarFore + botPortAft + botPortFore + topStarAft + topStarFore + topPortAft + topPortFore) );
	//index++;
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
  vector<vector3d<double> > fArea;                  //declare and reserve space for vector
  fArea.reserve(len);

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

        //fArea[index] = 0.5 * xbd.CrossProd(xac);     //area vector is calculated so that normal points nominally in direction of increasing i-coordinate
	fArea.push_back( 0.5 * xbd.CrossProd(xac) );     //area vector is calculated so that normal points nominally in direction of increasing i-coordinate

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii << ", j-dim = " << jj << ", k-dim = " << kk << endl;
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
  vector<vector3d<double> > fCenter;                //declare and reserve space for vector
  fCenter.reserve(len);

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  //int index = 0;                         //cell index counter

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
        //fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
        fCenter.push_back( 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft) );
	//index++;
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
  vector<vector3d<double> > fArea;                 //declare and reserve space for vector
  fArea.reserve(len);

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

        //fArea[index] = 0.5 * xbd.CrossProd(xac);       //area vector is calculated so that normal nominally points in direction of increasing j-coordinate
        fArea.push_back( 0.5 * xbd.CrossProd(xac) );       //area vector is calculated so that normal nominally points in direction of increasing j-coordinate

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
  vector<vector3d<double> > fCenter;               //declare and reserve space for vector
  fCenter.reserve(len);

  //vectors of coordinates of the 4 points that make up each face
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  //int index = 0;                         //cell index counter

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
	//fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
	fCenter.push_back( 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft) );
	//index++;
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
  vector<vector3d<double> > fArea;                  //declare and reserve space for vector
  fArea.reserve(len);

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

        //fArea[index] = 0.5 * xac.CrossProd(xbd);   //area vector is calculated so that normal nominally points in direction of increasing k-coordinate
        fArea.push_back( 0.5 * xac.CrossProd(xbd) );   //area vector is calculated so that normal nominally points in direction of increasing k-coordinate

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
  vector<vector3d<double> > fCenter;                //declare and reserve space for vector
  fCenter.reserve(len);

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  //int index = 0;                         //cell index counter

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
        //fCenter[index] = 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft);
        fCenter.push_back( 0.25 * (botStarAft + botPortAft + topStarAft + topPortAft) );
	//index++;
      }
    }
  }

  return fCenter;
}

//---------------------------------------------------------------------------------------------------------------//
//function to read in a plot3d grid and assign it to a plot3dMesh data type
vector<plot3dBlock> ReadP3dGrid(const string &gridName, double &numCells) {

  //open binary plot3d grid file 
  ifstream fName;
  string fPostfix = ".xyz";
  string readName = gridName + fPostfix;
  fName.open(readName.c_str(), ios::in | ios::binary);

  //check to see if file opened correctly
  if (fName.fail()) {
    cerr << "ERROR: Error in plot3d.cpp:ReadP3dGrid(). Grid file " << readName << " did not open correctly!!!" << endl;
    exit(0);
  }

  //read the number of plot3d blocks in the file
  cout << "Reading grid file..." << endl << endl;
  int numBlks;
  fName.read(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));
  cout << "Number of blocks: " << numBlks << endl << endl;

  //read the number of i, j, k coordinates in each plot3d block
  cout << "Size of each block is..." << endl;
  vector<int> vecI(numBlks,0);
  vector<int> vecJ(numBlks,0);
  vector<int> vecK(numBlks,0);
  int tempInt = 0;
  numCells = 0;

  //loop over all blocks and fill i, j, k vectors with block sizes
  for ( int ii=0; ii < numBlks; ii++){
  
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

    //calculate total number of cells (subtract 1 because number of cells is 1 less than number of points)
    numCells += (vecI[ii] - 1) * (vecJ[ii] - 1) * (vecK[ii] - 1);

  }
  cout << endl;

  //read each block and add it to the vector of plot3dBlocks
  double tempDouble=0;
  vector<plot3dBlock> mesh;
  mesh.reserve(numBlks);

  for ( int ii=0; ii < numBlks; ii++){
  
    int size = vecI[ii]*vecJ[ii]*vecK[ii];
    vector<double> tempX(size,0.0);
    vector<double> tempY(size,0.0);
    vector<double> tempZ(size,0.0);

    for (int jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempX[jj] = tempDouble;
    }
    for (int jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempY[jj] = tempDouble;
    }
    for (int jj=0; jj < size; jj++){
      fName.read(reinterpret_cast<char *>(&tempDouble), sizeof(tempDouble));
      tempZ[jj] = tempDouble;
    }

    //create single plot3dBlock and assign it appropriate location in vector
    plot3dBlock singleBlock(vecI[ii],vecJ[ii],vecK[ii],tempX,tempY,tempZ);
    mesh.push_back(singleBlock);

    cout << "Block " << ii << " read" << endl;
  }

  cout << endl << "Grid file read" << endl;
  cout << "Total number of cells is " << numCells << endl;

  //close plot3d grid file
  fName.close();

  return mesh;
}

//function to take in imax, jmax, kmax, and vector index and return the i, j, k indices of that cell
vector3d<int> GetIJK(const int &imax, const int &jmax, const int &kmax, const int &index){

  vector3d<int> ijk;
  bool breakFlag = false;

  //Loop over all possible i, j, k combinations. When calculated index matches given index, i, j, k location
  //is found. Once found, break out of loops and return
  for ( int kk = 0; kk < kmax; kk++ ){
    for ( int jj = 0; jj < jmax; jj++ ){
      for ( int ii = 0; ii < imax; ii++ ){
        int loc = GetLoc1D(ii, jj, kk, imax, jmax); 
	if ( loc == index ){
	  ijk.SetX(ii);
	  ijk.SetY(jj);
	  ijk.SetZ(kk);
	  breakFlag = true;
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
/*
    _________________________________________
    |         |         |          |         |
    |         |         |          |         |
    |   Ui-1  |   Ui    |   Ui+1   |   Ui+2  |  
    |         |         |          |         |
    |         |         |          |         |
    |_________|_________|__________|_________|
  Ui-1       Ui       Ui+1       Ui+2      Ui+3            

If the given i, j, k indices correspond to the location of cell Ui, the GetUpper functions will return the index of
the face at face location Ui+1 and the GetLower functions will return the index of the face at face location Ui.

The default is the return the index of the adjacent face, but this can be changed by supplying the last argument (num) to
the function.

*/
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
/*
    _________________________________________
    |         |         |          |         |
    |         |         |          |         |
    |   Ui-1  |   Ui    |   Ui+1   |   Ui+2  |  
    |         |         |          |         |
    |         |         |          |         |
    |_________|_________|__________|_________|
  Ui-1       Ui       Ui+1       Ui+2      Ui+3            

If the given i, j, k indices correspond to the location of cell Ui, the GetUpper functions will return the index of
the cell at location Ui+1 and the GetLower functions will return the index of the cell at location Ui-1.

If the given i, j, k indices correspond to the location of face Ui, the GetUpper functions will return the index of
the face at location Ui+1 and the GetLower functions will return the index of the face at location Ui-1.

The default is to return the index of the adjacent neighbor, but this can be changed by supplying the last argument (num) to
the function.

*/
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
/*
    _________________________________________
    |         |         |          |         |
    |         |         |          |         |
    |   Ui-1  |   Ui    |   Ui+1   |   Ui+2  |  
    |         |         |          |         |
    |         |         |          |         |
    |_________|_________|__________|_________|
  Ui-1       Ui       Ui+1       Ui+2      Ui+3            

If the given i, j, k indices correspond to the location of face Ui, the GetUpper functions will return the index of
the cell at location Ui and the GetLower functions will return the index of the cell at location Ui-1.

The default is to return the index of the adjacent cell, but this can be changed by supplying the last argument (num) to
the function.

*/
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

//function to take in j, j, k indexes of a cell and the diagonal label, and determine if the matrix should have data at the cell
//NOT currently USED, but useful for debugging
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
/*A hyperplane is a plane of i+j+k=constant within an individual block. The LUSGS solver must sweep along these hyperplanes to avoid
calculating a flux jacobian. Ex. The solver must visit all points on hyperplane 1 before visiting any points on hyperplane 2.
*/
vector<vector3d<int> > HyperplaneReorder( const int &imax, const int &jmax, const int &kmax){

  int numPlanes = imax + jmax + kmax - 2; //total number of hyperplanes in a given block
  vector<vector3d<int> > reorder;
  reorder.reserve(imax*jmax*kmax);

  int count = 0;
  for ( int pp = 0; pp < numPlanes; pp++ ){
    for ( int kk = 0; kk < kmax; kk++ ){
      for ( int jj = 0; jj < jmax; jj++ ){
	for ( int ii = 0; ii < imax; ii++ ){
	  if ( ii+jj+kk == pp){ //if sum of ii, jj, and kk equals pp than point is on hyperplane pp
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

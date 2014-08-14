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
const vector<double> plot3dBlock::Volume() const {

  int len = (numi-1)*(numj-1)*(numk-1);  //number of cells in the block
  vector<double> vol(len,0);             //initially assign a volume of 0

  //each hex cell is broken up into 3 pyramids and the volume of each pyramid is calculated
  double pyramid1 = 0;
  double pyramid2 = 0;
  double pyramid3 = 0;

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

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
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botStarFore.SetX(x[loc+1]);
        botStarFore.SetY(y[loc+1]);
        botStarFore.SetZ(z[loc+1]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);

        botPortFore.SetX(x[loc+numi+1]);
        botPortFore.SetY(y[loc+numi+1]);
        botPortFore.SetZ(z[loc+numi+1]);
              
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

        topStarFore.SetX(x[loc+numi*numj+1]);
        topStarFore.SetY(y[loc+numi*numj+1]);
        topStarFore.SetZ(z[loc+numi*numj+1]);

        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

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
const vector<vector3d<double> > plot3dBlock::Centroid() const {

  int len = (numi-1)*(numj-1)*(numk-1);                 //number of cells in the block
  vector<vector3d<double> > centroid(len);             //preallocate centroid vector
  vector3d<double> curCentroid(0,0,0);

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

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
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botStarFore.SetX(x[loc+1]);
        botStarFore.SetY(y[loc+1]);
        botStarFore.SetZ(z[loc+1]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);

        botPortFore.SetX(x[loc+numi+1]);
        botPortFore.SetY(y[loc+numi+1]);
        botPortFore.SetZ(z[loc+numi+1]);
              
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

        topStarFore.SetX(x[loc+numi*numj+1]);
        topStarFore.SetY(y[loc+numi*numj+1]);
        topStarFore.SetZ(z[loc+numi*numj+1]);

        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

        topPortFore.SetX(x[loc+numi+numi*numj+1]);
        topPortFore.SetY(y[loc+numi+numi*numj+1]);
        topPortFore.SetZ(z[loc+numi+numi*numj+1]);

	//Calculate the centroid of the cell
        curCentroid = (1.0/8.0) * (botStarAft + botStarFore + botPortAft + botPortFore + topStarAft + topStarFore + topPortAft + topPortFore);
        centroid[index] = curCentroid;

        index++;

      }
    }
  }

  return centroid;
}

//plot3dBlock member function that calcualtes the area of each face normal to the i-direction
const vector<vector3d<double> > plot3dBlock::FaceAreaI() const {

  int len = numi * (numj - 1) * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  vector3d<double> tempArea(0,0,0);

  //vectors of the components that go into calculating the volume of each pyramid
  vector3d<double> xac(0,0,0);           //vector from opposite corners of base
  vector3d<double> xbd(0,0,0);           //vector from opposite corners of base

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        xac = topPortAft - botStarAft;
        xbd = botPortAft - topStarAft;

        tempArea = 0.5 * xbd.CrossProd(xac);     //area vector is calculated so that normal points nominally in direction of increasing i-coordinate
        fArea[index].SetX(tempArea.X());
        fArea[index].SetY(tempArea.Y());
        fArea[index].SetZ(tempArea.Z());

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(1);
	}
        index++;

      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the center of each face normal to the i-direction
const vector<vector3d<double> > plot3dBlock::FaceCenterI() const {

  int len = numi * (numj - 1) * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
        topStarAft.SetX(x[loc+numi*numj]);
        topStarAft.SetY(y[loc+numi*numj]);
        topStarAft.SetZ(z[loc+numi*numj]);

        topPortAft.SetX(x[loc+numi+numi*numj]);
        topPortAft.SetY(y[loc+numi+numi*numj]);
        topPortAft.SetZ(z[loc+numi+numi*numj]);

        //Calculate face center by averaging four points that make up the face
        fCenter[index] = (1.0/4.0) * (botStarAft + botPortAft + topStarAft + topPortAft);

        index++;

      }
    }
  }

  return fCenter;
}

//plot3dBlock member function that calcualtes the area of each face normal to the j-direction
const vector<vector3d<double> > plot3dBlock::FaceAreaJ() const {

  int len = (numi -1) * numj * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  vector3d<double> tempArea(0,0,0);

  //vectors of the components that go into calculating the volume of each pyramid
  vector3d<double> xac(0,0,0);           //vector from opposite corners of base
  vector3d<double> xbd(0,0,0);           //vector from opposite corners of base

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc+1]);
        botStarAft.SetY(y[loc+1]);
        botStarAft.SetZ(z[loc+1]);

        botPortAft.SetX(x[loc]);
        botPortAft.SetY(y[loc]);
        botPortAft.SetZ(z[loc]);
              
        topStarAft.SetX(x[loc+numi*numj+1]);
        topStarAft.SetY(y[loc+numi*numj+1]);
        topStarAft.SetZ(z[loc+numi*numj+1]);

        topPortAft.SetX(x[loc+numi*numj]);
        topPortAft.SetY(y[loc+numi*numj]);
        topPortAft.SetZ(z[loc+numi*numj]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        xac = topPortAft - botStarAft;
        xbd = botPortAft - topStarAft;

        tempArea = 0.5 * xbd.CrossProd(xac);       //area vector is calculated so that normal nominally points in direction of increasing j-coordinate
        fArea[index].SetX(tempArea.X());
        fArea[index].SetY(tempArea.Y());
        fArea[index].SetZ(tempArea.Z());

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(1);
	}
        index++;

      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the area of each face normal to the j-direction
const vector<vector3d<double> > plot3dBlock::FaceCenterJ() const {

  int len = (numi -1) * numj * (numk - 1);         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk-1; kk++){
    for ( jj = 0; jj < numj; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc+1]);
        botStarAft.SetY(y[loc+1]);
        botStarAft.SetZ(z[loc+1]);

        botPortAft.SetX(x[loc]);
        botPortAft.SetY(y[loc]);
        botPortAft.SetZ(z[loc]);
              
        topStarAft.SetX(x[loc+numi*numj+1]);
        topStarAft.SetY(y[loc+numi*numj+1]);
        topStarAft.SetZ(z[loc+numi*numj+1]);

        topPortAft.SetX(x[loc+numi*numj]);
        topPortAft.SetY(y[loc+numi*numj]);
        topPortAft.SetZ(z[loc+numi*numj]);

        //Calculate face center by averaging the four points that make up the face
	fCenter[index] = (1.0/4.0) * (botStarAft + botPortAft + topStarAft + topPortAft);
        index++;

      }
    }
  }

  return fCenter;
}

//plot3dBlock member function that calcualtes the area of each face normal to the k-direction
const vector<vector3d<double> > plot3dBlock::FaceAreaK() const {

  int len = (numi - 1) * (numj - 1) * numk;         //number of i-faces in the block
  vector<vector3d<double> > fArea(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  vector3d<double> tempArea(0,0,0);

  //vectors of the components that go into calculating the volume of each pyramid
  vector3d<double> xac(0,0,0);           //vector from opposite corners of base
  vector3d<double> xbd(0,0,0);           //vector from opposite corners of base

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
        topStarAft.SetX(x[loc+1]);
        topStarAft.SetY(y[loc+1]);
        topStarAft.SetZ(z[loc+1]);

        topPortAft.SetX(x[loc+numi+1]);
        topPortAft.SetY(y[loc+numi+1]);
        topPortAft.SetZ(z[loc+numi+1]);

        //Calculate area for face by taking 1/2 of the cross product between opposite diagonals
        xac = topPortAft - botStarAft;
        xbd = botPortAft - topStarAft;

        tempArea = 0.5 * xac.CrossProd(xbd);   //area vector is calculated so that normal nominally points in direction of increasing k-coordinate
        fArea[index].SetX(tempArea.X());
        fArea[index].SetY(tempArea.Y());
        fArea[index].SetZ(tempArea.Z());

        if (fArea[index].Mag() <= 0){
          cerr << "ERROR: Negative face area in PLOT3D block at index " << index << "!!!" << endl;
          cerr << "Face area = " << fArea[index].Mag() << endl;
          cerr << "i-dim = " << ii+1 << ", j-dim = " << jj+1 << ", k-dim = " << kk+1 << endl;
          cerr << "Vectors to opposite diagonals are : " << xac << " and " << xbd << endl;
          exit(1);
	}
        index++;

      }
    }
  }

  return fArea;
}

//plot3dBlock member function that calcualtes the area of each face normal to the k-direction
const vector<vector3d<double> > plot3dBlock::FaceCenterK() const {

  int len = (numi - 1) * (numj - 1) * numk;         //number of i-faces in the block
  vector<vector3d<double> > fCenter(len);             //initially assign an area of 0

  //iteration counters
  int ii = 0;
  int jj = 0;
  int kk = 0;

  //vectors of coordinates of the 8 points that make up each hex cell
  vector3d<double> botStarAft(0,0,0);
  vector3d<double> botPortAft(0,0,0);
  vector3d<double> topStarAft(0,0,0);
  vector3d<double> topPortAft(0,0,0);

  int index = 0;                         //cell index counter
  int loc = 0;                           //location within vector of nodal coordinates

  for ( kk = 0; kk < numk; kk++){
    for ( jj = 0; jj < numj-1; jj++){
      for ( ii = 0; ii < numi-1; ii++){

        //assign coordinates to 8 points that make up the cell being analyzed
        loc = ii + jj*numi + kk*numi*numj;

        botStarAft.SetX(x[loc]);
        botStarAft.SetY(y[loc]);
        botStarAft.SetZ(z[loc]);

        botPortAft.SetX(x[loc+numi]);
        botPortAft.SetY(y[loc+numi]);
        botPortAft.SetZ(z[loc+numi]);
              
        topStarAft.SetX(x[loc+1]);
        topStarAft.SetY(y[loc+1]);
        topStarAft.SetZ(z[loc+1]);

        topPortAft.SetX(x[loc+numi+1]);
        topPortAft.SetY(y[loc+numi+1]);
        topPortAft.SetZ(z[loc+numi+1]);

        //Calculate face center by averaging four points that make up cell face
        fCenter[index] = (1.0/4.0) * (botStarAft + botPortAft + topStarAft + topPortAft);
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
plot3dMesh ReadP3dGrid(const string &gridName) {

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
  return i + (j-1)*imax + k*imax*jmax;
}
//function to take in i, j, k indexes of a cell/face and return the 1D index of the neighbor cell/face to the upper j face
int GetNeighborUpJ(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  return i + (j+1)*imax + k*imax*jmax;
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
  if( (k-1-(num-1)) < 0 ){
    cerr << "ERROR: Tried to access cell outside of range!" << endl;
    cerr << "Tried to access cell " << i << ", " << j << ", " << k-1-(num-1) << endl;
    exit(0);
  }

  return i + j*imax + (k-1-(num-1))*imax*jmax;
}
//function to take in i, j, k indexes of a k-face and return the 1D index of the neighbor cell in the upper k direction
int GetCellFromFaceUpperK(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){
  if( (k-1+num) < 0 ){
    cerr << "ERROR: Tried to access cell outside of range!" << endl;
    cerr << "Tried to access cell " << i << ", " << j << ", " << k-1+num << endl;
    exit(0);
  }

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

  if( (i-1-(num-1)) < 0 || (i-1-(num-1)) > (imax-2) ){
    cerr << "ERROR: Tried to access cell outside of range!" << endl;
    cerr << "Tried to access cell " << i-1-(num-1) << ", " << j << ", " << k << endl;
    exit(0);
  }

  return (i-1-(num-1)) + j*(imax-1) + k*(imax-1)*jmax;
}
//function to take in i, j, k indexes of a i-face and return the 1D index of the neighbor cell in the upper i direction
int GetCellFromFaceUpperI(const int &i, const int &j, const int &k, const int &imax, const int &jmax, int num){

  if( (i-1+num) < 0 || (i-1+num) > (imax-2) ){
    cerr << "ERROR: Tried to access cell outside of range!" << endl;
    cerr << "Tried to access cell " << i-1+num << ", " << j << ", " << k << endl;
    exit(0);
  }

  return (i-1+num) + j*(imax-1) + k*(imax-1)*jmax;
}


//function to take in i, j, k, imax, jmax and return corresponding location inside 1D array
int GetLoc1D(const int &i, const int &j, const int &k, const int &imax, const int &jmax){
  return i + j*imax + k*imax*jmax;
}

//function to take in i, j, k indexes of a i-face and return the 1D index of the upper i diagonal for the implicit matrix
int GetMatrixDiagUpperFromMainI(const int &main){
  return main - 1;
}
//function to take in i, j, k indexes of a i-face and return the 1D index of the lower i diagonal for the implicit matrix
int GetMatrixDiagLowerFromMainI(const int &main){
  return main + 1;
}
//function to take in i, j, k indexes of a j-face and return the 1D index of the upper j diagonal for the implicit matrix
int GetMatrixDiagUpperFromMainJ(const int &main, const int &imax){
  return main - imax;
}
//function to take in i, j, k indexes of a j-face and return the 1D index of the lower j diagonal for the implicit matrix
int GetMatrixDiagLowerFromMainJ(const int &main, const int &imax){
  return main + imax;
}
//function to take in i, j, k indexes of a k-face and return the 1D index of the upper k diagonal for the implicit matrix
int GetMatrixDiagUpperFromMainK(const int &main, const int &imax, const int &jmax){
  return main - imax*jmax;
}
//function to take in i, j, k indexes of a k-face and return the 1D index of the lower k diagonal for the implicit matrix
int GetMatrixDiagLowerFromMainK(const int &main, const int &imax, const int &jmax){
  return main + imax*jmax;
}

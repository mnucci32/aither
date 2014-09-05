#include "output.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::to_string;
using std::max;

//---------------------------------------------------------------------------------------------------------------//
//function declarations
//function to write out cell centers of grid in plot3d format
void WriteCellCenter(const string &gridName, const vector<blockVars > &vars) {

  //open binary plot3d grid file 
  ofstream outFile;
  string fEnd = "_center";
  string fPostfix = ".xyz";
  string writeName = gridName + fEnd + fPostfix;
  outFile.open(writeName.c_str(), ios::out | ios::binary);

  //check to see if file opened correctly
  if (outFile.fail()) {
    cerr << "ERROR: Grid file " << writeName << " did not open correctly!!!" << endl;
    exit(0);
  }

  //write number of blocks to file
  int numBlks = vars.size();
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  //write i, j, k dimension for each block
  int ll = 0;
  int dumInt = 0;
  int jj = 0;
  int kk = 0;
  vector3d<double> dumVec;
  double dumDouble=0.0;

  for ( ll=0; ll < numBlks; ll++ ){ //loop over all blocks
    //subtract 1 from max values because blockVars maxes are in terms of nodes, not cells
    dumInt = vars[ll].NumI()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumJ()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumK()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
  }

  //write out x, y, z coordinates of cell centers
  for ( ll = 0; ll < numBlks; ll++ ){  //loop over all blocks
    int maxi = vars[ll].NumI()-1;
    int maxj = vars[ll].NumJ()-1;
    int maxk = vars[ll].NumK()-1;
    int blkLen = maxi * maxj * maxk;

    for ( kk = 0; kk < 3; kk++ ){  //loop over dimensions (3)
      for ( jj = 0; jj < blkLen; jj++ ){ //loop over length of block

	dumVec = vars[ll].Center(jj);  //at a given cell, get the cell center coordinates

	//for a given block, first write out all x coordinates, then all y coordinates, then all z coordinates
	if (kk == 0 ) {
	  dumDouble = dumVec.X();
	}
	else if (kk == 1){
	  dumDouble = dumVec.Y();
	}
	else {
	  dumDouble = dumVec.Z();
	}
	//write to file
	outFile.write(reinterpret_cast<char *>(&dumDouble), sizeof(dumDouble));
      }
    }
  }

  //close plot3d grid file
  outFile.close();

}

//---------------------------------------------------------------------------------------------------------------//
//function to write out variables in function file format
void WriteFun(const string &gridName, const vector<blockVars> &vars, const vector<viscBlockVars> &vVars, const idealGas &eqnState, const double &solTime, const double &refRho, const double &refSoS, const double &refT) {

  //open binary plot3d function file 
  ofstream outFile;
  string fEnd = "_center";
  string fPostfix = ".fun";
  string writeName = gridName + to_string(static_cast<long long int>(solTime)) + fEnd + fPostfix;
  outFile.open(writeName.c_str(), ios::out | ios::binary);

  //check to see if file opened correctly
  if (outFile.fail()) {
    cerr << "ERROR: Function file " << writeName << " did not open correctly!!!" << endl;
    exit(0);
  }

  //write number of blocks to file
  int numBlks = vars.size();
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  //write i, j, k, vars dimension for each block
  int numVars = 23;            //number of variables to write out
  int ll = 0;
  int dumInt = 0;
  int jj = 0;
  int kk = 0;
  int ii = 0;

  vector3d<double> vel;
  double dumDouble=0.0;

  for ( ll=0; ll < numBlks; ll++ ){ //loop over all blocks and write out imax, jmax, kmax, numVars
    //subtract 1 from maxes because blockVars maxes are in terms of nodes, not cells
    dumInt = vars[ll].NumI()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumJ()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumK()-1;
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));

    outFile.write(reinterpret_cast<char *>(&numVars), sizeof(numVars));
  }

  //write out variables
  for ( ll = 0; ll < numBlks; ll++ ){ //loop over all blocks
    int maxi = vars[ll].NumI()-1;
    int maxj = vars[ll].NumJ()-1;
    int maxk = vars[ll].NumK()-1;
    int blkLen = maxi * maxj * maxk;
    vector<double> dumVec(blkLen);

    for ( kk = 0; kk < numVars; kk++ ){ //loop over the number of variables to write out
      //store nondimensional variable in dumVec for a given block in order. i.e. var1 var2 var3 etc
      if (kk == 0) {
	for ( ii = 0; ii < blkLen; ii++){         //density
	  dumVec[ii] = vars[ll].State(ii).Rho();
	}
      }
      else if (kk == 1) {
	for ( ii = 0; ii < blkLen; ii++){         //vel-x
	  dumVec[ii] = vars[ll].State(ii).U();
	}
      }
      else if (kk == 2) {
	for ( ii = 0; ii < blkLen; ii++){         //vel-y
	  dumVec[ii] = vars[ll].State(ii).V();
	}
      }
      else if (kk == 3) {                         //vel-z
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].State(ii).W();
	}
      }
      else if (kk == 4) {                        //pressure
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].State(ii).P();
	}
      }
      else if (kk == 5) {                      //mach
	for ( ii = 0; ii < blkLen; ii++){
	  vel = vars[ll].State(ii).Velocity();
	  dumVec[ii] = vel.Mag() / eqnState.GetSoS( vars[ll].State(ii).P(), vars[ll].State(ii).Rho() );
	}
      }
      else if (kk == 6) {                     //speed of sound
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = eqnState.GetSoS( vars[ll].State(ii).P(), vars[ll].State(ii).Rho() );
	}
      }
      else if (kk == 7) {                         
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Dt(ii);    //time step
	}
      }
      else if (kk == 8) {                         
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].State(ii).Temperature(eqnState);    //temperature
	}
      }
      else if (kk == 9) {
	for ( ii = 0; ii < blkLen; ii++){       //du/dx
	  dumVec[ii] = vVars[ll].VelGrad(ii).XX();
	}
      }
      else if (kk == 10) {
	for ( ii = 0; ii < blkLen; ii++){       //dv/dx
	  dumVec[ii] = vVars[ll].VelGrad(ii).XY();
	}
      }
      else if (kk == 11) {
	for ( ii = 0; ii < blkLen; ii++){       //dw/dx
	  dumVec[ii] = vVars[ll].VelGrad(ii).XZ();
	}
      }
      else if (kk == 12) {
	for ( ii = 0; ii < blkLen; ii++){       //du/dy
	  dumVec[ii] = vVars[ll].VelGrad(ii).YX();
	}
      }
      else if (kk == 13) {
	for ( ii = 0; ii < blkLen; ii++){       //dv/dy
	  dumVec[ii] = vVars[ll].VelGrad(ii).YY();
	}
      }
      else if (kk == 14) {
	for ( ii = 0; ii < blkLen; ii++){       //dw/dy
	  dumVec[ii] = vVars[ll].VelGrad(ii).YZ();
	}
      }
      else if (kk == 15) {
	for ( ii = 0; ii < blkLen; ii++){       //du/dz
	  dumVec[ii] = vVars[ll].VelGrad(ii).ZX();
	}
      }
      else if (kk == 16) {
	for ( ii = 0; ii < blkLen; ii++){       //dv/dz
	  dumVec[ii] = vVars[ll].VelGrad(ii).ZY();
	}
      }
      else if (kk == 17) {
	for ( ii = 0; ii < blkLen; ii++){       //dw/dz
	  dumVec[ii] = vVars[ll].VelGrad(ii).ZZ();
	}
      }
      else if (kk == 18) {                        //mass residual
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,0);
	}
      }
      else if (kk == 19) {                        //momentum-x residual
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,1);
	}
      }
      else if (kk == 20) {                        //momentum-y residual
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,2);
	}
      }
      else if (kk == 21) {                        //momentum-z residual
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,3);
	}
      }
      else if (kk == 22) {                        //energy residual
	for ( ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,4);
	}
      }
      else {
	cerr << "ERROR: Variable to write to function file is not defined!" << endl;
        exit(0);
      }

      for ( jj = 0; jj < blkLen; jj++ ){                              //write out dimensional variables -- loop over block length

        dumDouble = dumVec[jj];

	if (kk == 0){                                        //density
	  dumDouble = dumDouble * refRho;
	}
	else if (kk == 1){                                  //velocity x
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 2){                                 //velocity y
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 3){                                 //velocity z
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 4){                                 // pressure
	  dumDouble = dumDouble * refRho * refSoS * refSoS ; 
	}
	else if (kk == 5){                                 //mach is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (kk == 6){                                 // speed of sound
	  dumDouble = dumDouble * refSoS ; 
	}
	else if (kk == 7){                                 //time step
	  dumDouble = dumDouble / refSoS;                        
	}
	else if (kk == 8){                                 //temperature
	  dumDouble = dumDouble * refT;                        
	}
	else if (kk == 9){                                  //du/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 10){                                 //dv/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 11){                                 //dw/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 12){                                  //du/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 13){                                 //dv/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 14){                                 //dw/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 15){                                  //du/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 16){                                 //dv/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 17){                                 //dw/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (kk == 18){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (kk == 19){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (kk == 20){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (kk == 21){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (kk == 22){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	outFile.write(reinterpret_cast<char *>(&dumDouble), sizeof(dumDouble));
      }
    }
  }


  //close plot3d function file
  outFile.close();

}

//function to write out results file for ensight
void WriteRes(const string &gridName, const int &iter, const int &outFreq) {

  //open results file
  ofstream resFile;
  string fResPostfix = ".res";
  string fPostfix = ".fun";
  string fEnd = "_center";
  string resName = gridName + fEnd + fResPostfix;
  resFile.open(resName.c_str(), ios::out);

  string writeName = gridName + "*" + fEnd + fPostfix;

  //check to see if file opened correctly
  if (resFile.fail()) {
    cerr << "ERROR: Results file " << resName << " did not open correctly!!!" << endl;
    exit(0);
  }

  //write number of scalars and number of vectors
  int numScalar = 18;
  int numVector = 1;
  resFile << numScalar << "     " << numVector << "     " << 0 << endl;

  //write number of time points that there is solution data at
  int numTime = iter/outFreq;
  resFile << numTime << endl;

  //Write solution times or iteration numbers
  int solTime = 0;
  int ii = 0;
  int count = 1;
  for (ii = 0; ii < numTime; ii++){
    solTime += outFreq;
    if ( count%10 == 0 ){
      resFile << endl;
    }
    resFile << solTime << "   ";
    count++;
  }
  resFile << endl;

  //Write starting iteration and iteration increment
  resFile << outFreq << "  " << outFreq << endl;


  //Write out variables
  resFile << writeName << " F 0001 density" << endl;
  resFile << writeName << " F 0002 Vx" << endl;
  resFile << writeName << " F 0003 Vy" << endl;
  resFile << writeName << " F 0004 Vz" << endl;
  resFile << writeName << " F 0005 pressure" << endl;
  resFile << writeName << " F 0006 mach" << endl;
  resFile << writeName << " F 0007 sos" << endl;
  resFile << writeName << " F 0008 dt" << endl;
  resFile << writeName << " F 0009 temperature" << endl;
  resFile << writeName << " F 0010 du/dx" << endl;
  resFile << writeName << " F 0011 dv/dx" << endl;
  resFile << writeName << " F 0012 dw/dx" << endl;
  resFile << writeName << " F 0013 du/dy" << endl;
  resFile << writeName << " F 0014 dv/dy" << endl;
  resFile << writeName << " F 0015 dw/dy" << endl;
  resFile << writeName << " F 0016 du/dz" << endl;
  resFile << writeName << " F 0017 dv/dz" << endl;
  resFile << writeName << " F 0018 dw/dz" << endl;
  resFile << writeName << " F 0019 massRes" << endl;
  resFile << writeName << " F 0020 momxRes" << endl;
  resFile << writeName << " F 0021 momyRes" << endl;
  resFile << writeName << " F 0022 momzRes" << endl;
  resFile << writeName << " F 0023 engyRes" << endl;
  resFile << writeName << " F 0002 0003 0004 velocity" << endl;

  //Close results file
  resFile.close();

}


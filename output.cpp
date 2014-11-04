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
void WriteCellCenter(const string &gridName, const vector<procBlock > &vars) {

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
  int dumInt = 0;
  vector3d<double> dumVec;
  double dumDouble=0.0;

  for ( int ll=0; ll < numBlks; ll++ ){ //loop over all blocks
    //subtract 1 from max values because procBlock maxes are in terms of nodes, not cells
    dumInt = vars[ll].NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
  }

  //write out x, y, z coordinates of cell centers
  for ( int ll = 0; ll < numBlks; ll++ ){  //loop over all blocks
    int maxi = vars[ll].NumI();
    int maxj = vars[ll].NumJ();
    int maxk = vars[ll].NumK();
    int maxiG = vars[ll].NumI() + 2 * vars[ll].NumGhosts();
    int maxjG = vars[ll].NumJ() + 2 * vars[ll].NumGhosts();

    for ( int nn = 0; nn < 3; nn++ ){  //loop over dimensions (3)

      for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	  for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){ 

	    int loc = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	    dumVec = vars[ll].Center(loc);  //at a given cell, get the cell center coordinates

	    //for a given block, first write out all x coordinates, then all y coordinates, then all z coordinates
	    if (nn == 0 ) {
	      dumDouble = dumVec.X();
	    }
	    else if (nn == 1){
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

    }
  }

  //close plot3d grid file
  outFile.close();

}

//---------------------------------------------------------------------------------------------------------------//
//function to write out variables in function file format
void WriteFun(const string &gridName, const vector<procBlock> &vars, const idealGas &eqnState, const double &solTime, const double &refRho, const double &refSoS, const double &refT) {

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
  int dumInt = 0;

  vector3d<double> vel;
  double dumDouble=0.0;

  for ( int ll=0; ll < numBlks; ll++ ){ //loop over all blocks and write out imax, jmax, kmax, numVars
    dumInt = vars[ll].NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = vars[ll].NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));

    outFile.write(reinterpret_cast<char *>(&numVars), sizeof(numVars));
  }

  //write out variables
  for ( int ll = 0; ll < numBlks; ll++ ){ //loop over all blocks
    int maxi = vars[ll].NumI();
    int maxj = vars[ll].NumJ();
    int maxk = vars[ll].NumK();
    int blkLen = vars[ll].NumCells();
    vector<double> dumVec(blkLen);
    int maxiG = vars[ll].NumI() + 2 * vars[ll].NumGhosts();
    int maxjG = vars[ll].NumJ() + 2 * vars[ll].NumGhosts();

    for ( int vv = 0; vv < numVars; vv++ ){ //loop over the number of variables to write out
      //store nondimensional variable in dumVec for a given block in order. i.e. var1 var2 var3 etc
      if (vv == 0) {   //density
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){ 
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).Rho(); 
	    }
	  }
	}
      }
      else if (vv == 1) {  //vel-x
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).U();  
	    }
	  }
	}
      }
      else if (vv == 2) {  //vel-y
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).V();
	    }
	  }
	}
      }
      else if (vv == 3) {   //vel-z
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).W(); 
	    }
	  }
	}
      }
      else if (vv == 4) {                        //pressure
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).P();
	    }
	  }
	}
      }
      else if (vv == 5) {                      //mach
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      vel = vars[ll].State(locG).Velocity();
	      dumVec[loc] = vel.Mag() / eqnState.GetSoS( vars[ll].State(locG).P(), vars[ll].State(locG).Rho() );
	    }
	  }
	}
      }
      else if (vv == 6) {                     //speed of sound
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = eqnState.GetSoS( vars[ll].State(locG).P(), vars[ll].State(locG).Rho() );
	    }
	  }
	}
      }
      else if (vv == 7) {                    //time step - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Dt(ii);  
	}
      }
      else if (vv == 8) {                     //temperature
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].State(locG).Temperature(eqnState);   
	    }
	  }
	}
      }
      else if (vv == 9) {    // du/dx
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).XX();
	    }
	  }
	}
      }
      else if (vv == 10) {  // dv/dx
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).XY();
	    }
	  }
	}
      }
      else if (vv == 11) {  // dw/dx
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).XZ();
	    }
	  }
	}
      }
      else if (vv == 12) {   // du/dy
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).YX();
	    }
	  }
	}
      }
      else if (vv == 13) {  // dv/dy
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).YY();
	    }
	  }
	}
      }
      else if (vv == 14) {   // dw/dy
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).YZ();
	    }
	  }
	}
      }
      else if (vv == 15) {   // du/dz
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).ZX();
	    }
	  }
	}
      }
      else if (vv == 16) {   // dv/dz
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).ZY();
	    }
	  }
	}
      }
      else if (vv == 17) {  // dw/dz
	for ( int kk = vars[ll].NumGhosts(); kk < maxk + vars[ll].NumGhosts(); kk++ ){
	  for ( int jj = vars[ll].NumGhosts(); jj < maxj + vars[ll].NumGhosts(); jj++ ){
	    for ( int ii = vars[ll].NumGhosts(); ii < maxi + vars[ll].NumGhosts(); ii++){         
	      int loc = GetLoc1D(ii - vars[ll].NumGhosts(), jj - vars[ll].NumGhosts(), kk - vars[ll].NumGhosts(), maxi, maxj);
	      int locG = GetLoc1D(ii, jj, kk, maxiG, maxjG);
	      dumVec[loc] = vars[ll].VelGrad(locG).ZZ();
	    }
	  }
	}
      }
      else if (vv == 18) {                        //mass residual - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,0);
	}
      }
      else if (vv == 19) {                        //momentum-x residual - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,1);
	}
      }
      else if (vv == 20) {                        //momentum-y residual - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,2);
	}
      }
      else if (vv == 21) {                        //momentum-z residual - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,3);
	}
      }
      else if (vv == 22) {                        //energy residual - no ghost cells
	for ( int ii = 0; ii < blkLen; ii++){
	  dumVec[ii] = vars[ll].Residual(ii,4);
	}
      }
      else {
	cerr << "ERROR: Variable to write to function file is not defined!" << endl;
        exit(0);
      }

      for ( int nn = 0; nn < blkLen; nn++ ){                              //write out dimensional variables -- loop over block length

        dumDouble = dumVec[nn];

	if (vv == 0){                                        //density
	  dumDouble = dumDouble * refRho;
	}
	else if (vv == 1){                                  //velocity x
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 2){                                 //velocity y
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 3){                                 //velocity z
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 4){                                 // pressure
	  dumDouble = dumDouble * refRho * refSoS * refSoS ; 
	}
	else if (vv == 5){                                 //mach is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (vv == 6){                                 // speed of sound
	  dumDouble = dumDouble * refSoS ; 
	}
	else if (vv == 7){                                 //time step
	  dumDouble = dumDouble / refSoS;                        
	}
	else if (vv == 8){                                 //temperature
	  dumDouble = dumDouble * refT;                        
	}
	else if (vv == 9){                                  //du/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 10){                                 //dv/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 11){                                 //dw/dx
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 12){                                  //du/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 13){                                 //dv/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 14){                                 //dw/dy
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 15){                                  //du/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 16){                                 //dv/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 17){                                 //dw/dz
	  dumDouble = dumDouble * refSoS;
	}
	else if (vv == 18){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (vv == 19){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (vv == 20){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (vv == 21){                                 //residual is already nondimensional
	  dumDouble = dumDouble ;                        
	}
	else if (vv == 22){                                 //residual is already nondimensional
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


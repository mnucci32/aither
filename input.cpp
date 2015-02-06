#include <iostream>     //cout
#include <fstream>      //ifstream
#include <cstdlib>      //exit()
#include <sstream>      //istringstream
#include <iterator>     //istring_iterator
#include <time.h>       //strftime
#include "input.h"

#define ROOT 0

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using std::cerr;
using std::istringstream;
using std::istream_iterator;

//constructor for input class
//initialize vector to have length of number of acceptable inputs to the code
input::input(): vars(32){
  vector<double> stagProps(6,0.0);
  vector<double> pOutlet(2,0.0);

  //default values for each variable
  gName = "";
  dt = -1.0;
  iterations = 1;
  pRef = -1.0;
  rRef = -1.0;
  lRef = 1.0;
  gamma = 1.4;
  gasConst = 287.058;
  velRef.SetX(1.0);
  velRef.SetY(0.0);
  velRef.SetZ(0.0);
  vector<boundaryConditions> dummy(1);
  bc = dummy;
  timeIntegration = "explicitEuler";
  cfl = -1.0;
  kappa = -2.0;         //default to value outside of range to tell if higher order or constant method should be used
  limiter = "none";
  outputFrequency = 1;
  equationSet = "euler";
  tRef = -1.0;
  matrixSolver = "lusgs";
  matrixSweeps = 1;
  matrixRelaxation = 1.0;  //default is symmetric Gauss-Seidel with no overrelaxation
  timeIntTheta = 1.0;      //default results in implicit euler
  timeIntZeta = 0.0;       //default results in implicit euler
  nonlinearIterations = 1;    //default is 1 (steady)
  cflMax = 1.0;
  cflStep = 0.0;
  cflStart = 1.0;
  invFluxJac = "laxFriedrichs";  //default is approximate roe flux
  dualTimeCFL = -1.0;          //default value of -1; negative value means dual time stepping is not used
  inviscidFlux = "roe";        //default value is roe flux
  stagInletProps = stagProps;
  pressureOutlet = pOutlet;
  decompMethod = "manual";     //default is manual decomposition

  //keywords in the input file that the parser is looking for to define variables
  vars[0] = "gridName:";
  vars[1] = "timeStep:";
  vars[2] = "iterations:";
  vars[3] = "pressureRef:";
  vars[4] = "densityRef:";
  vars[5] = "lengthRef:";
  vars[6] = "gamma:";
  vars[7] = "gasConstant:";
  vars[8] = "velocity:";
  vars[9] = "timeIntegration:";
  vars[10] = "cfl:";
  vars[11] = "faceReconstruction:";
  vars[12] = "limiter:";
  vars[13] = "outputFrequency:";
  vars[14] = "equationSet:";
  vars[15] = "temperatureRef:";
  vars[16] = "matrixSolver:";
  vars[17] = "matrixSweeps:";
  vars[18] = "matrixRelaxation:";
  vars[19] = "timeIntTheta:";
  vars[20] = "timeIntZeta:";
  vars[21] = "nonlinearIterations:";
  vars[22] = "cflMax:";
  vars[23] = "cflStep:";
  vars[24] = "cflStart:";
  vars[25] = "inviscidFluxJacobian:";
  vars[26] = "dualTimeCFL:";
  vars[27] = "inviscidFlux:";
  vars[28] = "stagnationInlet:";
  vars[29] = "pressureOutlet:";
  vars[30] = "decompositionMethod:";

  vars[31] = "boundaryConditions:";  //boundary conditions should be listed last
}

//member function to set vector holding boundary conditions for each block
void input::SetBCVec(const int &a){
  bc.resize(a);
}

//function to trim leading and trailing whitespace from a string, and also remove data after a comment
string trim(const string &s, const string &whitespace = " \t"){

  const string comment = "#";                                  //# is comment character for input file

  if (s.empty()){
    return ""; //string is empty
  }
  else{
    const unsigned int sBegin = s.find_first_not_of(whitespace);          //find index of first non whitespace character
    const int sEnd = s.find_last_not_of(whitespace);            //find index of last non whitespace character
    const int sRange = sEnd - sBegin +1;                        //range to trim string to
    string temp = s.substr(sBegin,sRange);

    const int tempComment = temp.find(comment);                 //find index of first comment character
    const int tempRange = tempComment - 0;                      //find range of string to return

    return temp.substr(0,tempRange);
  }
}

//function to print the time
void PrintTime(){

  time_t rawtime;
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);

  int bufLen = 100; //hard coded max buffer length
  char buffer [bufLen];
  strftime(buffer, bufLen, "%c", timeinfo);
  cout << buffer << endl;
}

//function to read the input file and return the data as a member of the input class
//read the input file
input ReadInput(const string &inputName, const int &rank){
  // inputName -- name of input file
  // rank -- rank of processor

  if (rank == ROOT){
    cout << "###########################################################################################################################" << endl;
    PrintTime();
    cout << endl;
    cout << "Parsing input file " << inputName << endl << endl;
    cout << "Solver Inputs" << endl;
  }

  input inputVars;

  //open input file
  ifstream inFile;
  inFile.open(inputName.c_str(), ios::in);
  if (inFile.fail()) {
    cerr << "ERROR: Error in input.cpp:ReadInput(). Input file " << inputName << " did not open correctly!!!" << endl;
    exit(0);
  }

  string line;
  vector3d<double> temp;
  int readingBCs = 0;
  int blk = 0;
  int surfCounter = 0;
  int numBCBlks = 0;
  vector<boundaryConditions> tempBC(1);
  int lEnd = 0;
  int numSurf = 0;

  while(getline(inFile,line)){   //while there are still lines in the input file, execute loop

    line = trim(line);  //remove leading and trailing whitespace and ignore comments

    //split line into words
    istringstream buf(line);
    istream_iterator<string> beg(buf), end;
    vector<string> tokens(beg,end);
   
    //search to see if first token corresponds to any keywords
    if (tokens.size() >= 2) { //line must contain at least 2 tokens (keyword and value)

      for( int ii=0; ii < inputVars.NumVars(); ii++){ //loop over all input variables and see if they match keywords

        if (tokens[0] == inputVars.Vars(ii) || readingBCs > 0){ //if first token matches a keyword or reading boundary conditions
	  
	  //if not yet reading BCs (readingBCs == 0), set variable in input class to corresponding value and print assignment to std out
          if (ii==0 && readingBCs == 0){
            inputVars.SetGridName(tokens[1]);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.GridName() << endl;}
            continue;
          }
          else if (ii==1 && readingBCs == 0){
            inputVars.SetDt(atof(tokens[1].c_str()));                          //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.Dt() << endl;}
            continue;
	  }
          else if (ii==2 && readingBCs == 0){
            inputVars.SetIterations(atoi(tokens[1].c_str()));
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.Iterations() << endl;}
            continue;
	  }
          else if (ii==3 && readingBCs == 0){
            inputVars.SetPRef(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.PRef() << endl;}
            continue;
	  }
          else if (ii==4 && readingBCs == 0){
            inputVars.SetRRef(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.RRef() << endl;}
            continue;
	  }
          else if (ii==5 && readingBCs == 0){
            inputVars.SetLRef(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.LRef() << endl;}
            continue;
	  }
          else if (ii==6 && readingBCs == 0){
            inputVars.SetGamma(atof(tokens[1].c_str()));                      //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.Gamma() << endl;}
            continue;
	  }
          else if (ii==7 && readingBCs == 0){
            inputVars.SetR(atof(tokens[1].c_str()));                         //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.R() << endl;}
            continue;
	  }
          else if (ii==8 && readingBCs == 0){
            temp.SetX(atof(tokens[1].c_str()));                              //double variable (atof)
            temp.SetY(atof(tokens[2].c_str()));
            temp.SetZ(atof(tokens[3].c_str()));
            inputVars.SetVelRef(temp);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.VelRef() << endl;}
            continue;
	  }
          else if (ii==9 && readingBCs == 0){
            inputVars.SetTimeIntegration(tokens[1]);
	    if (inputVars.TimeIntegration() == "implicitEuler"){
	      inputVars.SetTheta(1.0);
	      inputVars.SetZeta(0.0);
	    }
	    else if (inputVars.TimeIntegration() == "crankNicholson"){
	      inputVars.SetTheta(0.5);
	      inputVars.SetZeta(0.0);
	    }
	    else if (inputVars.TimeIntegration() == "bdf2"){
	      inputVars.SetTheta(1.0);
	      inputVars.SetZeta(0.5);
	    }

            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.TimeIntegration() << endl;}
            continue;
          }
          else if (ii==11 && readingBCs == 0){
	    if (tokens[1] == "upwind"){
	      inputVars.SetKappa(-1.0);
	    }
	    else if (tokens[1] == "fromm"){
	      inputVars.SetKappa(0.0);
	    }
	    else if (tokens[1] == "quick"){
	      inputVars.SetKappa(0.5);
	    }
	    else if (tokens[1] == "central"){
	      inputVars.SetKappa(1.0);
	    }
	    else if (tokens[1] == "thirdOrder"){
	      inputVars.SetKappa(1.0/3.0);
	    }
	    else if (tokens[1] == "constant"){
	      inputVars.SetKappa(-2.0);
	    }
	    else{
	      inputVars.SetKappa(atof(tokens[1].c_str()));         //if string is not recognized, set kappa to number input
	    }

            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << tokens[1] << " kappa = " << inputVars.Kappa() << endl;}
	    if ( (inputVars.Kappa() < -1.0) || (inputVars.Kappa() > 1.0) ){
	      cerr << "ERROR: Error in input.cpp:ReadInput(). Kappa value of " << inputVars.Kappa() 
		   << " is not valid! Choose a value between -1.0 and 1.0." << endl;
	      exit(0);
	    }
            continue;
          }
          else if (ii==12 && readingBCs == 0){
            inputVars.SetLimiter(tokens[1]);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.Limiter() << endl;}
            continue;
          }
          else if (ii==13 && readingBCs == 0){
            inputVars.SetOutputFrequency(atoi(tokens[1].c_str()));
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.OutputFrequency() << endl;}
            continue;
	  }
          else if (ii==14 && readingBCs == 0){
            inputVars.SetEquationSet(tokens[1]);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.EquationSet() << endl;}
            continue;
          }
          else if (ii==15 && readingBCs == 0){
            inputVars.SetTRef(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.TRef() << endl;}
            continue;
	  }
          else if (ii==16 && readingBCs == 0){
            inputVars.SetMatrixSolver(tokens[1]);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.MatrixSolver() << endl;}
            continue;
	  }
          else if (ii==17 && readingBCs == 0){
            inputVars.SetMatrixSweeps(atoi(tokens[1].c_str()));
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.MatrixSweeps() << endl;}
            continue;
	  }
          else if (ii==18 && readingBCs == 0){
            inputVars.SetMatrixRelaxation(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.MatrixRelaxation() << endl;}
            continue;
	  }
          else if (ii==21 && readingBCs == 0){
            inputVars.SetNonlinearIterations(atoi(tokens[1].c_str()));
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.NonlinearIterations() << endl;}
            continue;
	  }
          else if (ii==22 && readingBCs == 0){
            inputVars.SetCFLMax(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.CFLMax() << endl;}
            continue;
          }
          else if (ii==23 && readingBCs == 0){
            inputVars.SetCFLStep(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.CFLStep() << endl;}
            continue;
          }
          else if (ii==24 && readingBCs == 0){
            inputVars.SetCFLStart(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.CFLStart() << endl;}
            continue;
          }
          else if (ii==25 && readingBCs == 0){
            inputVars.SetInvFluxJac(tokens[1]);
	    if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.InvFluxJac() << endl;}
	  }
          else if (ii==26 && readingBCs == 0){
            inputVars.SetDualTimeCFL(atof(tokens[1].c_str()));                       //double variable (atof)
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.DualTimeCFL() << endl;}
            continue;
          }
          else if (ii==27 && readingBCs == 0){
            inputVars.SetInviscidFlux(tokens[1]);
	    if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.InviscidFlux() << endl;}
	  }
          else if (ii==28 && readingBCs == 0){
            inputVars.SetStagInletTag(atoi(tokens[1].c_str()));
            inputVars.SetStagInletP0(atof(tokens[2].c_str()));
            inputVars.SetStagInletT0(atof(tokens[3].c_str()));
            inputVars.SetStagInletDx(atof(tokens[4].c_str()));
            inputVars.SetStagInletDy(atof(tokens[5].c_str()));
            inputVars.SetStagInletDz(atof(tokens[6].c_str()));
	    if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.StagInletTag() << " " << inputVars.StagInletP0() <<
		" " << inputVars.StagInletT0() << " " << inputVars.StagInletDx() << " " << inputVars.StagInletDy() <<
		" " << inputVars.StagInletDz() << endl;}
	  }
          else if (ii==29 && readingBCs == 0){
            inputVars.SetPressureOutletTag(atoi(tokens[1].c_str()));
            inputVars.SetPressureOutletP(atof(tokens[2].c_str()));
	    if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.PressureOutletTag() << " " << inputVars.PressureOutletP() << endl;}
	  }
          else if (ii==30 && readingBCs == 0){
            inputVars.SetDecompMethod(tokens[1]);
            if (rank == ROOT) {cout << inputVars.Vars(ii) << " " << inputVars.DecompMethod() << endl;}
            continue;
          }

	  //reading BCs -------------------------------------------------------------------------------------------------------
          else if (ii==inputVars.NumVars()-1 || readingBCs > 0){
	    //read in boundary conditions and assign to boundaryConditions class
	    if (readingBCs == 0) {  //variable read must be number of blocks if first line of BCs
	      numBCBlks = atoi(tokens[1].c_str());
	      tempBC.resize(numBCBlks); //resize vector holding BCs to correct number of blocks
	      readingBCs++;
	      lEnd++;
	    }	    
	    else if (readingBCs == lEnd){ //variables must be number of i, j, k surfaces for each block 
	      //set number of i, j, k surfaces
	      tempBC[blk].SetNumSurfI(atoi(tokens[0].c_str()));
	      tempBC[blk].SetNumSurfJ(atoi(tokens[1].c_str()));
	      tempBC[blk].SetNumSurfK(atoi(tokens[2].c_str()));

	      //set total number of surfaces
	      numSurf = atoi(tokens[0].c_str()) + atoi(tokens[1].c_str()) + atoi(tokens[2].c_str());
	      tempBC[blk].ResizeVecs(numSurf);

	      lEnd += (numSurf + 1); //set new target for when block data is finished
	      readingBCs++;

	    }
	    else {  //assign BC block variables
	      tempBC[blk].SetBCTypes(tokens[0], surfCounter);
	      tempBC[blk].SetIMin(atoi(tokens[1].c_str()), surfCounter);
	      tempBC[blk].SetIMax(atoi(tokens[2].c_str()), surfCounter);
	      tempBC[blk].SetJMin(atoi(tokens[3].c_str()), surfCounter);
	      tempBC[blk].SetJMax(atoi(tokens[4].c_str()), surfCounter);
	      tempBC[blk].SetKMin(atoi(tokens[5].c_str()), surfCounter);
	      tempBC[blk].SetKMax(atoi(tokens[6].c_str()), surfCounter);
	      tempBC[blk].SetTag(atoi(tokens[7].c_str()), surfCounter);

	      surfCounter++;
	      if (surfCounter == numSurf){ //at end of block data, increment block index, reset surface counter
		blk++;
		surfCounter = 0;
	      }

	      readingBCs++;
	    }

	    //if block counter reaches number of blocks, BCs are finished (b/c counter starts at 0), so assign BCs and write them out
	    if (blk == numBCBlks){ 
	      inputVars.SetBC(tempBC);
	      if (rank == ROOT) {
		cout << inputVars.Vars(inputVars.NumVars()-1) << " " << inputVars.NumBC() << endl << endl;
		for ( int ll = 0; ll < inputVars.NumBC(); ll++ ){
		  cout << "Block: " << ll << endl;
		  cout << inputVars.BC(ll) << endl;
		}
	      }
	    }

	    break;
	  }

	}
      }
    }
    else{ //if there aren't 2 or more tokens just skip line
      continue;
    }
  }

  if (rank == ROOT) {
    cout << endl;
    cout << "Input file parse complete" << endl;
    cout << "###########################################################################################################################" << endl << endl;
  }

  return inputVars;

}


//member function to calculate the cfl value for the step from the starting, ending, and step values
void input::CalcCFL(const int &i){

  if ( i == 0 ){  //first time step
    cfl = cflStart;
  }
  else if ( cflStart + i*cflStep > cflMax ){ //if calculated value is higher than max, set to max
    cfl = cflMax;
  }
  else{
    cfl = cflStart + i * cflStep;
  }

}

#ifndef OUTPUTHEADERDEF             //only if the macro OUTPUTHEADERDEF is not defined execute these lines of code
#define OUTPUTHEADERDEF             //define the macro

/* This header contains the function declarations for file output.

It contains function headers to write out the grid at the cell centers in Plot3D format, as well as the Plot3D function files.
It also writes out a master file in Ensight format to name the Plot3D functions.  */

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "plot3d.h" //plot3d
#include "eos.h"
#include "primVars.h" //primVars
#include "procBlock.h" //procBlock
#include "inviscidFlux.h" //inviscidFlux
#include "input.h" //inputVars
#include "boundaryConditions.h" //decomposition
#include <fstream>
#include <iostream>

using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::cout;
using std::endl;
using std::cerr;

//function definitions
void WriteCellCenter(const string&, const vector<procBlock> &, const decomposition&);
void WriteCellCenterGhost(const string&, const vector<procBlock> &);
void WriteFun(const string&, const vector<procBlock> &, const idealGas&, const int&, const double&, const double&, const double&, const decomposition&);
void WriteRes(const string&, const int&, const int&);

void WriteResiduals(const input&, genArray&, genArray&, const resid&, const double&, const int&, const int&);

vector<procBlock> Recombine( const vector<procBlock>&, const decomposition& );

#endif

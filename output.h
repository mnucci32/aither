#ifndef OUTPUTHEADERDEF             //only if the macro OUTPUTHEADERDEF is not defined execute these lines of code
#define OUTPUTHEADERDEF             //define the macro

#include <vector>  //vector
#include <string>  //string
#include "vector3d.h" //vector3d
#include "tensor.h" //tensor
#include "plot3d.h" //plot3d
#include "eos.h"
#include "primVars.h" //primVars
#include "blockVars.h" //blockVars
#include "viscBlockVars.h" //viscBlockVars
#include "inviscidFlux.h" //inviscidFlux
#include "input.h" //inputVars
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
void WriteCellCenter(const string&, const vector<blockVars> &);
void WriteFun(const string&, const vector<blockVars> &, const vector<viscBlockVars> &, const idealGas&, const double&, const double&, const double&, const double&);
void WriteRes(const string&, const int&, const int&);


#endif

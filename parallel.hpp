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

#ifndef PARALLELHEADERDEF             //only if the macro PARALLELHEADERDEF is not defined execute these lines of code
#define PARALLELHEADERDEF             //define the macro

/* This header contains the function declarations for many of the parallel functions in the code */

#include "mpi.h" //parallelism
#include <vector>  //vector
#include <string>  //string
#include "vector3d.hpp" //vector3d
#include "plot3d.hpp" //plot3d
#include "primVars.hpp" //primVars
#include "procBlock.hpp" //procBlock
#include "boundaryConditions.hpp" //interblock
#include <iostream>
#include "macros.hpp"

using std::vector;
using std::string;
using std::ios;
using std::cout;
using std::endl;
using std::cerr;
using std::ostream;


//function definitions
decomposition ManualDecomposition(vector<plot3dBlock>&, vector<boundaryConditions>&, const int& );
decomposition CubicDecomposition(vector<plot3dBlock>&, vector<boundaryConditions>&, const int& );

void SendNumProcBlocks(const vector<int>&, const int&, int&);
void SendConnections(vector<interblock>&, const MPI_Datatype& );

void SetDataTypesMPI(const int&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype&, MPI_Datatype& );

vector<procBlock> SendProcBlocks( const vector<procBlock>&, const int&, const int&, const MPI_Datatype&, const MPI_Datatype& );
void GetProcBlocks( vector<procBlock>&, const vector<procBlock>&, const int&, const MPI_Datatype& );

void MaxLinf( resid*, resid*, int*, MPI_Datatype* );

void BroadcastString( string &str );
#endif

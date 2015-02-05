#include "parallel.h"

/* Function to return processor list for manual decomposition. Manual decomposition assumes that each block will reside on it's own processor.
*/
void ManualDecomposition(vector<procBlock> &grid, const int &numProc){
  // grid -- vector of procBlocks (no need to split procBlocks or combine them with manual decomposition)
  // numProc -- number of processors in run

  if ( (int)grid.size() != numProc ){
    cerr << "ERROR: Error in parallel.cpp:ManualDecomposition(). Manual decomposition assumes that the number of processors used is equal to the " <<
      "number of blocks in the grid. This grid has " << grid.size() << " blocks and the simulation is using " << numProc << " processors." << endl;
    exit(0);
  }

  cout << "Using manual grid decomposition." << endl;

  //assign processor rank for each procBlock
  for ( int ii = 0; ii < numProc; ii++ ){
    grid[ii].SetRank(ii);
  }

}


void SetDataTypesMPI(const int &numEqn, MPI_Datatype &MPI_vec3d, MPI_Datatype &MPI_cellData, MPI_Datatype &MPI_procBlockInts){

  //create vector3d<double> MPI datatype
  MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_vec3d);
  MPI_Type_commit(&MPI_vec3d);

  //create MPI datatype for states (primVars), residuals (colMatrix), etc
  MPI_Type_contiguous(numEqn, MPI_DOUBLE, &MPI_cellData);
  MPI_Type_commit(&MPI_cellData);

  //create MPI datatype for all the integers in the procBlock class
  MPI_Type_contiguous(14, MPI_INT, &MPI_procBlockInts);
  MPI_Type_commit(&MPI_procBlockInts);

}


//function to send procBlocks to their appropriate processor
vector<procBlock> SendProcBlocks( const vector<procBlock> &blocks, const int &rank, const MPI_Datatype &MPI_procBlockInts, const MPI_Datatype &MPI_cellData, const MPI_Datatype &MPI_vec3d ){

  //PROBLEM -- function is passed blocks, but blocks is only available on root!!!


  for ( unsigned int ii = 0; ii < blocks.size(); ii++ ){
    int numRank = 0;
    if ( blocks[ii].Rank() != ROOT ){ //send/receive data
      if ( rank == ROOT ){ //send data
	int tempS[14];
	tempS[0]  = blocks[ii].NumCells();
	tempS[1]  = blocks[ii].NumVars();
	tempS[2]  = blocks[ii].NumI();
	tempS[3]  = blocks[ii].NumJ();
	tempS[4]  = blocks[ii].NumK();
	tempS[5]  = blocks[ii].NumGhosts();
	tempS[6]  = blocks[ii].ParentBlock();
	tempS[7]  = blocks[ii].ParentBlockStartI();
	tempS[8]  = blocks[ii].ParentBlockEndI();
	tempS[9]  = blocks[ii].ParentBlockStartJ();
	tempS[10] = blocks[ii].ParentBlockEndJ();
	tempS[11] = blocks[ii].ParentBlockStartK();
	tempS[12]  = blocks[ii].ParentBlockEndK();
	tempS[13]  = blocks[ii].Rank();

	MPI_Request request;
	MPI_Isend(&tempS, 1, MPI_procBlockInts, blocks[ii].Rank(), 1, MPI_COMM_WORLD, &request); //send integers

	vector<double> doubleVecS = blocks[ii].VolVec();
	MPI_Isend(&doubleVecS, doubleVecS.size(), MPI_DOUBLE, blocks[ii].Rank(), 2, MPI_COMM_WORLD, &request); //send volumes

	vector<primVars> primVecS = blocks[ii].StateVec(); 
	MPI_Isend(&primVecS, primVecS.size(), MPI_cellData, blocks[ii].Rank(), 3, MPI_COMM_WORLD, &request); //send states

	vector<vector3d<double> > vec3dVecS = blocks[ii].CenterVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 4, MPI_COMM_WORLD, &request); //send cell centers

	vec3dVecS = blocks[ii].FAreaIVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 5, MPI_COMM_WORLD, &request); //send face area Is

	vec3dVecS = blocks[ii].FAreaJVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 6, MPI_COMM_WORLD, &request); //send face area Js

	vec3dVecS = blocks[ii].FAreaKVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 7, MPI_COMM_WORLD, &request); //send face area Ks

	vec3dVecS = blocks[ii].FCenterIVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 8, MPI_COMM_WORLD, &request); //send face center Is

	vec3dVecS = blocks[ii].FCenterJVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 9, MPI_COMM_WORLD, &request); //send face center Js

	vec3dVecS = blocks[ii].FCenterKVec(); 
	MPI_Isend(&vec3dVecS, vec3dVecS.size(), MPI_vec3d, blocks[ii].Rank(), 10, MPI_COMM_WORLD, &request); //send face center Ks

	//send boundary conditions
	boundaryConditions bound = blocks[ii].BC();
	vector<int> intVecS = bound.GetIMinVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 11, MPI_COMM_WORLD, &request); //send i-min coordinates
	intVecS = bound.GetIMaxVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 12, MPI_COMM_WORLD, &request); //send i-max coordinates
	intVecS = bound.GetJMinVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 13, MPI_COMM_WORLD, &request); //send j-min coordinates
	intVecS = bound.GetJMaxVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 14, MPI_COMM_WORLD, &request); //send j-max coordinates
	intVecS = bound.GetKMinVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 15, MPI_COMM_WORLD, &request); //send k-min coordinates
	intVecS = bound.GetKMaxVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 16, MPI_COMM_WORLD, &request); //send k-max coordinates
	intVecS = bound.GetTagVec();
	MPI_Isend(&intVecS, intVecS.size(), MPI_INT, blocks[ii].Rank(), 17, MPI_COMM_WORLD, &request); //send tags
	int dims[3] = {bound.NumSurfI(), bound.NumSurfJ(), bound.NumSurfK()};
	MPI_Isend(&dims, 3, MPI_INT, blocks[ii].Rank(), 1, MPI_COMM_WORLD, &request); //send number of surfaces
	vector<string> strVecS = bound.GetBCTypesVec();
	MPI_Isend(strVecS[0].c_str(), strVecS.size(), MPI_CHAR, blocks[ii].Rank(), 18, MPI_COMM_WORLD, &request); //send bc types
      }
      else{ //receive data
	int tempR[14];
	MPI_Status status;
	MPI_Recv(&tempR, 1, MPI_procBlockInts, ROOT, 1, MPI_COMM_WORLD, &status); //receive integers

      }

    }
  }



}

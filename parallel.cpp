#include "parallel.h"

/* Function to return processor list for manual decomposition. Manual decomposition assumes that each block will reside on it's own processor.
*/
vector<int> ManualDecomposition(vector<procBlock> &grid, const int &numProc){
  // grid -- vector of procBlocks (no need to split procBlocks or combine them with manual decomposition)
  // numProc -- number of processors in run

  if ( (int)grid.size() != numProc ){
    cerr << "ERROR: Error in parallel.cpp:ManualDecomposition(). Manual decomposition assumes that the number of processors used is equal to the " <<
      "number of blocks in the grid. This grid has " << grid.size() << " blocks and the simulation is using " << numProc << " processors." << endl;
    exit(0);
  }

  cout << "Using manual grid decomposition." << endl;

  vector<int> loadBal(numProc,1); //vector containing number of procBlocks for each processor

  //assign processor rank for each procBlock
  for ( int ii = 0; ii < numProc; ii++ ){
    grid[ii].SetRank(ii);
  }

  return loadBal;
}

//function to send each processor the number of procBlocks that it should contain
void SendNumProcBlocks(const vector<int> &loadBal, const int &rank, int &numProcBlock){
  MPI_Scatter(&loadBal[0], 1, MPI_INT, &numProcBlock, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
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
vector<procBlock> SendProcBlocks( const vector<procBlock> &blocks, const int &rank, const int &numProcBlock, const MPI_Datatype &MPI_procBlockInts, 
				  const MPI_Datatype &MPI_cellData, const MPI_Datatype &MPI_vec3d ){

  vector<procBlock> localBlocks;
  localBlocks.reserve(numProcBlock);

  if ( rank == ROOT ){ // may have to pack and send data
    for ( unsigned int ii = 0; ii < blocks.size(); ii++ ){

      if (blocks[ii].Rank() == ROOT ){ //no need to send data because it is already on ROOT processor
	localBlocks.push_back(blocks[ii]);
      }
      else{ //send data
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

	vector<primVars> primVecS = blocks[ii].StateVec(); 
	vector<vector3d<double> > centerVecS = blocks[ii].CenterVec(); 
	vector<vector3d<double> > fAreaIVecS = blocks[ii].FAreaIVec(); 
	vector<vector3d<double> > fAreaJVecS = blocks[ii].FAreaJVec(); 
	vector<vector3d<double> > fAreaKVecS = blocks[ii].FAreaKVec(); 
	vector<vector3d<double> > fCenterIVecS = blocks[ii].FCenterIVec(); 
	vector<vector3d<double> > fCenterJVecS = blocks[ii].FCenterJVec(); 
	vector<vector3d<double> > fCenterKVecS = blocks[ii].FCenterKVec(); 
	vector<double> volVecS = blocks[ii].VolVec(); 
	boundaryConditions bound = blocks[ii].BC();
	int surfs [3] = {bound.NumSurfI(), bound.NumSurfJ(), bound.NumSurfK()};
	vector<string> names = bound.GetBCTypesVec();
	vector<int> iMin = bound.GetIMinVec();
	vector<int> iMax = bound.GetIMaxVec();
	vector<int> jMin = bound.GetJMinVec();
	vector<int> jMax = bound.GetJMaxVec();
	vector<int> kMin = bound.GetKMinVec();
	vector<int> kMax = bound.GetKMaxVec();
	vector<int> tags = bound.GetTagVec();

	//string lengths
	vector<int> strLengthS(bound.NumSurfI() + bound.NumSurfJ() + bound.NumSurfK() );
	for ( unsigned int jj = 0; jj < strLengthS.size(); jj++ ){
	  strLengthS[jj] = bound.GetBCTypes(jj).size();
	}

	int sendBufSize = 0;
	int tempSize = 0;
	MPI_Pack_size(14, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for ints in class
	sendBufSize += tempSize;
	MPI_Pack_size(primVecS.size(), MPI_cellData, MPI_COMM_WORLD, &tempSize); //add size for states
	sendBufSize += tempSize;
	MPI_Pack_size(centerVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for cell centers
	sendBufSize += tempSize;
	MPI_Pack_size(fAreaIVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face area I
	sendBufSize += tempSize;
	MPI_Pack_size(fAreaJVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face area J
	sendBufSize += tempSize;
	MPI_Pack_size(fAreaKVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face area K
	sendBufSize += tempSize;
	MPI_Pack_size(fCenterIVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face center I
	sendBufSize += tempSize;
	MPI_Pack_size(fCenterJVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face center J
	sendBufSize += tempSize;
	MPI_Pack_size(fCenterKVecS.size(), MPI_vec3d, MPI_COMM_WORLD, &tempSize); //add size for face center K
	sendBufSize += tempSize;
	MPI_Pack_size(volVecS.size(), MPI_DOUBLE, MPI_COMM_WORLD, &tempSize); //add size for volumes
	sendBufSize += tempSize;
	MPI_Pack_size(3, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for number of surfaces
	sendBufSize += tempSize;
	//8x because iMin, iMax, jMin, jMax, kMin, kMax, tags, string sizes
	MPI_Pack_size( (bound.NumSurfI() + bound.NumSurfJ() + bound.NumSurfK()) * 8, MPI_INT, MPI_COMM_WORLD, &tempSize); //add size for states
	sendBufSize += tempSize;
	int stringSize = 0;
	for ( int jj = 0; jj < (bound.NumSurfI() + bound.NumSurfJ() + bound.NumSurfK()); jj++ ){
	  MPI_Pack_size(bound.GetBCTypes(jj).size(), MPI_CHAR, MPI_COMM_WORLD, &tempSize); //add size for states
	  stringSize += tempSize;
	}
	sendBufSize += stringSize;

	char *sendBuffer = new char[sendBufSize]{};

	int position = 0;
	MPI_Pack(&tempS[0], 14, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&primVecS[0], primVecS.size(), MPI_cellData, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&centerVecS[0], centerVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fAreaIVecS[0], fAreaIVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fAreaJVecS[0], fAreaJVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fAreaKVecS[0], fAreaKVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fCenterIVecS[0], fCenterIVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fCenterJVecS[0], fCenterJVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&fCenterKVecS[0], fCenterKVecS.size(), MPI_vec3d, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&volVecS[0], volVecS.size(), MPI_DOUBLE, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);

	MPI_Pack(&surfs[0], 3, MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&iMin[0], iMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&iMax[0], iMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&jMin[0], jMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&jMax[0], jMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&kMin[0], kMin.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&kMax[0], kMax.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&tags[0], tags.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	MPI_Pack(&strLengthS[0], strLengthS.size(), MPI_INT, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	for ( int jj = 0; jj < (bound.NumSurfI() + bound.NumSurfJ() + bound.NumSurfK()); jj++ ){
	  MPI_Pack(names[jj].c_str(), names[jj].size(), MPI_CHAR, sendBuffer, sendBufSize, &position, MPI_COMM_WORLD);
	  cout << "sent char is " << names[jj].c_str() << ", size of " << strLengthS[jj] << endl;
	}

	MPI_Send(sendBuffer, sendBufSize, MPI_PACKED, blocks[ii].Rank(), 2, MPI_COMM_WORLD);

	delete [] sendBuffer;

      }
    }

  }
  

  if (rank != ROOT) { // receive and unpack data
    for ( int ii = 0; ii < numProcBlock; ii++ ){

      MPI_Status status;
      int tempR[14];

      int bufRecvData = 0;
      MPI_Probe(ROOT, 2, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_CHAR, &bufRecvData);

      char *recvBuffer = new char[bufRecvData]{};

      cout << "processor " << rank << " has a buffer of size " << bufRecvData << endl;

      MPI_Recv(recvBuffer, bufRecvData, MPI_PACKED, ROOT, 2, MPI_COMM_WORLD, &status);

      int position = 0;
      MPI_Unpack(recvBuffer, bufRecvData, &position, &tempR[0], 14, MPI_INT, MPI_COMM_WORLD); //unpack ints

      procBlock tempBlock(tempR[2], tempR[3], tempR[4], tempR[5]);
      tempBlock.SetParentBlock(tempR[6]);
      tempBlock.SetParentBlockStartI(tempR[7]);
      tempBlock.SetParentBlockEndI(tempR[8]);
      tempBlock.SetParentBlockStartJ(tempR[9]);
      tempBlock.SetParentBlockEndJ(tempR[10]);
      tempBlock.SetParentBlockStartK(tempR[11]);
      tempBlock.SetParentBlockEndK(tempR[12]);
      tempBlock.SetRank(rank);

      int numCells = (tempR[2] + 2 * tempR[5]) * (tempR[3] + 2 * tempR[5]) * (tempR[4] + 2 * tempR[5]);
      int numFaceI = (tempR[2] + 2 * tempR[5] + 1) * (tempR[3] + 2 * tempR[5]) * (tempR[4] + 2 * tempR[5]);
      int numFaceJ = (tempR[2] + 2 * tempR[5]) * (tempR[3] + 2 * tempR[5] + 1) * (tempR[4] + 2 * tempR[5]);
      int numFaceK = (tempR[2] + 2 * tempR[5]) * (tempR[3] + 2 * tempR[5]) * (tempR[4] + 2 * tempR[5] + 1);

      vector<primVars> primVecR(numCells);
      vector<vector3d<double> > centerVecR(numCells);
      vector<vector3d<double> > fAreaIVecR(numFaceI);
      vector<vector3d<double> > fAreaJVecR(numFaceJ);
      vector<vector3d<double> > fAreaKVecR(numFaceK);
      vector<vector3d<double> > fCenterIVecR(numFaceI);
      vector<vector3d<double> > fCenterJVecR(numFaceJ);
      vector<vector3d<double> > fCenterKVecR(numFaceK);
      vector<double> volVecR(numCells);

      MPI_Unpack(recvBuffer, bufRecvData, &position, &primVecR[0], primVecR.size(), MPI_cellData, MPI_COMM_WORLD); //unpack states
      MPI_Unpack(recvBuffer, bufRecvData, &position, &centerVecR[0], centerVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack cell centers
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fAreaIVecR[0], fAreaIVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face area I
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fAreaJVecR[0], fAreaJVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face area J
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fAreaKVecR[0], fAreaKVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face area K
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fCenterIVecR[0], fCenterIVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face center I
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fCenterJVecR[0], fCenterJVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face center J
      MPI_Unpack(recvBuffer, bufRecvData, &position, &fCenterKVecR[0], fCenterKVecR.size(), MPI_vec3d, MPI_COMM_WORLD); //unpack face center K
      MPI_Unpack(recvBuffer, bufRecvData, &position, &volVecR[0], volVecR.size(), MPI_DOUBLE, MPI_COMM_WORLD); //unpack volumes

      tempBlock.SetStateVec(primVecR); //assign states
      tempBlock.SetCenterVec(centerVecR); //assign cell centers
      tempBlock.SetFAreaIVec(fAreaIVecR); //assign face area I
      tempBlock.SetFAreaJVec(fAreaJVecR); //assign face area J
      tempBlock.SetFAreaKVec(fAreaKVecR); //assign face area K
      tempBlock.SetFCenterIVec(fCenterIVecR); //assign face center I
      tempBlock.SetFCenterJVec(fCenterJVecR); //assign face center J
      tempBlock.SetFCenterKVec(fCenterKVecR); //assign face center K
      tempBlock.SetVolVec(volVecR); //assign volumes

      int surfs[3];
      MPI_Unpack(recvBuffer, bufRecvData, &position, &surfs[0], 3, MPI_INT, MPI_COMM_WORLD); //unpack number of surfaces

      boundaryConditions bound(surfs[0], surfs[1], surfs[2]);
      vector<int> iMin(surfs[0] + surfs[1] + surfs[2]);
      vector<int> iMax(surfs[0] + surfs[1] + surfs[2]);
      vector<int> jMin(surfs[0] + surfs[1] + surfs[2]);
      vector<int> jMax(surfs[0] + surfs[1] + surfs[2]);
      vector<int> kMin(surfs[0] + surfs[1] + surfs[2]);
      vector<int> kMax(surfs[0] + surfs[1] + surfs[2]);
      vector<int> tags(surfs[0] + surfs[1] + surfs[2]);
      vector<int> strLength(surfs[0] + surfs[1] + surfs[2]);
      vector<string> names(surfs[0] + surfs[1] + surfs[2]);

      MPI_Unpack(recvBuffer, bufRecvData, &position, &iMin[0], iMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack i min coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &iMax[0], iMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack i max coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &jMin[0], jMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack j min coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &jMax[0], jMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack j max coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &kMin[0], kMin.size(), MPI_INT, MPI_COMM_WORLD); //unpack k min coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &kMax[0], kMax.size(), MPI_INT, MPI_COMM_WORLD); //unpack k max coordinates
      MPI_Unpack(recvBuffer, bufRecvData, &position, &tags[0], tags.size(), MPI_INT, MPI_COMM_WORLD); //unpack tags
      MPI_Unpack(recvBuffer, bufRecvData, &position, &strLength[0], strLength.size(), MPI_INT, MPI_COMM_WORLD); //unpack string sizes
      for ( unsigned int jj = 0; jj < strLength.size(); jj++ ){
	char nameBuf[100]{};
      	cout << iMin[jj] << ", " << iMax[jj] << ", " << jMin[jj] << ", " << jMax[jj] << ", " << kMin[jj] << ", " << kMax[jj] << ", " << tags[jj] << ", " << strLength[jj] << endl;
      	MPI_Unpack(recvBuffer, bufRecvData, &position, &nameBuf[0], strLength[jj], MPI_CHAR, MPI_COMM_WORLD); //unpack bc types
      	string bcName(nameBuf);
      	names[jj] = bcName;
      	cout << "I am processor " << rank << ". I received boundary condition " << names[jj] << endl;
      }

      delete [] recvBuffer;

      //set up boundary conditions
      bound.SetBCTypesVec(names);
      bound.SetIMinVec(iMin);
      bound.SetIMaxVec(iMax);
      bound.SetJMinVec(jMin);
      bound.SetJMaxVec(jMax);
      bound.SetKMinVec(kMin);
      bound.SetKMaxVec(kMax);
      bound.SetTagVec(tags);
      tempBlock.SetBCs(bound);

      localBlocks.push_back(tempBlock);

    }
  }

  return localBlocks;

}

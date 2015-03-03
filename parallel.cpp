#include "parallel.h"
#include "output.h"

/* Function to return processor list for manual decomposition. Manual decomposition assumes that each block will reside on it's own processor.
The processor list tells how many procBlocks a processor will have.
*/
vector<int> ManualDecomposition(vector<plot3dBlock> &grid, vector<vector3d<int> > &rankParPos, const int &numProc, const int &totalCells ){
  // grid -- vector of procBlocks (no need to split procBlocks or combine them with manual decomposition)
  // blkRank -- rank of processor that plot3dBlock will go on
  // blkPar -- parent block of plot3dBlock
  // numProc -- number of processors in run
  // totalCells -- total number of cells in the grid; used for load balancing metrics

  if ( (int)grid.size() != numProc ){
    cerr << "ERROR: Error in parallel.cpp:ManualDecomposition(). Manual decomposition assumes that the number of processors used is equal to the " <<
      "number of blocks in the grid. This grid has " << grid.size() << " blocks and the simulation is using " << numProc << " processors." << endl;
    exit(0);
  }

  cout << "--------------------------------------------------------------------------------" << endl;
  cout << "Using manual grid decomposition." << endl;

  double idealLoad = (double)totalCells / (double)numProc; //average number of cells per processor
  int maxLoad = 0;

  //vector containing number of procBlocks for each processor
  //in manual decomp, each proc gets 1 block
  vector<int> loadBal(numProc, 1);

  //assign processor rank for each plot3dBlock; global position is index, ii
  rankParPos.resize(grid.size());
  for ( unsigned int ii = 0; ii < rankParPos.size(); ii++ ){
    rankParPos[ii][0] = ii;   //rank
    rankParPos[ii][1] = ii;   //parent block
    rankParPos[ii][2] = 0;    //local block
  }

  //find maximum number of cells on a processor
  for ( unsigned int ii = 0; ii < grid.size(); ii++ ){
    maxLoad = std::max(maxLoad, grid[ii].NumCells());
  }

  cout << "Ratio of most loaded processor to average processor is : " << (double)maxLoad / idealLoad << endl;
  cout << "--------------------------------------------------------------------------------" << endl << endl;;

  return loadBal;
}

/* Function to return processor list for cubic decomposition. 
The processor list tells how many procBlocks a processor will have.
*/
vector<int> CubicDecomposition(vector<plot3dBlock> &grid, vector<vector3d<int> >&rankParPos, vector<boundaryConditions> &bcs, const int &numProc, const int &totalCells ){
  // grid -- vector of procBlocks (no need to split procBlocks or combine them with manual decomposition)
  // blkRank -- rank of processor that plot3dBlock will go on
  // blkPar -- parent block of plot3dBlock
  // bcs -- vector of boundary conditions for all blocks
  // numProc -- number of processors in run
  // totalCells -- total number of cells in the grid; used for load balancing metrics

  cout << "--------------------------------------------------------------------------------" << endl;
  cout << "Using cubic grid decomposition." << endl;

  double idealLoad = (double)totalCells / (double)numProc; //average number of cells per processor
  int maxLoad = 0;

  ///////////////////////////////////////////////////

  //vector containing number of procBlocks for each processor
  //in cubic decomp, each proc gets 1 block
  //vector<int> loadBal(numProc, 1);
  vector<int> loadBal(1, 5);

  rankParPos.resize(grid.size());
  for ( unsigned int ii = 0; ii < rankParPos.size(); ii++ ){
    //blkRank[ii] = ii;
    rankParPos[ii][0] = 0;      //rank
    //blkPar[ii] = ii;
    rankParPos[ii][1] = 0;      //parent block
    rankParPos[ii][2] = ii;     //local position
  }

  string dir = "j";
  int ind = 50;
  int blkNum = 3;
  int newBlk = 4;
  vector<boundarySurface> altSurf;
  plot3dBlock lBlk, uBlk; 
  grid[blkNum].Split(dir,ind, lBlk, uBlk);
  grid.push_back(uBlk);
  vector3d<int> temp(0, blkNum, newBlk);
  rankParPos.push_back(temp);
  boundaryConditions newBcs = bcs[blkNum].Split(dir,ind, blkNum, newBlk, altSurf);
  bcs.push_back(newBcs);

  for ( unsigned int ii = 0; ii < altSurf.size(); ii++ ){
    bcs[altSurf[ii].PartnerBlock()].DependentSplit(altSurf[ii], grid[blkNum], grid[altSurf[ii].PartnerBlock()], altSurf[ii].PartnerBlock(), dir, ind, blkNum, newBlk);
  }
  //reassign split grid
  grid[blkNum] = lBlk;



  cout << "updated BCs are:" << endl;
  for ( unsigned int ii = 0; ii < bcs.size(); ii++ ){
    cout << bcs[ii] << endl;
  }

  cout << "updated block sizes are:" << endl;
  for ( unsigned int ii = 0; ii < grid.size(); ii++ ){
    cout << "Block: " << ii << "  " << grid[ii].NumI() << ", " << grid[ii].NumJ() << ", " << grid[ii].NumK() << endl;
  }




  // boundaryConditions myBC = bcs[2];
  // boundarySurface mySurf1 = myBC.GetSurface(3);
  // cout << "original surface: " << mySurf1 << endl;
  // bool split;
  // boundarySurface mySurf2 = mySurf1.Split("j", 20, 1, 7, split);
  // cout << "split surfaces:" << endl;
  // cout << "was split: " << split << endl;
  // cout << mySurf1 << endl;
  // cout << mySurf2 << endl;

  //exit(0);


  //find maximum number of cells on a processor
  for ( unsigned int ii = 0; ii < grid.size(); ii++ ){
    maxLoad = std::max(maxLoad, grid[ii].NumCells());
  }

  ////////////////////////////////////////////////////////////

  cout << "Ratio of most loaded processor to average processor is : " << (double)maxLoad / idealLoad << endl;
  cout << "--------------------------------------------------------------------------------" << endl << endl;;

  return loadBal;
}

//function to send each processor the number of procBlocks that it should contain
void SendNumProcBlocks(const vector<int> &loadBal, const int &rank, int &numProcBlock){
  MPI_Scatter(&loadBal[0], 1, MPI_INT, &numProcBlock, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
}

//function to send each processor the vector of interblocks it needs to compute its boundary conditions
void SendConnections(vector<interblock> &connections, const MPI_Datatype &MPI_interblock) {

  //first determine the number of interblocks and send that to all processors
  int numCon = connections.size();
  MPI_Bcast(&numCon, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  connections.resize(numCon); //allocate space to receive the interblocks

  //broadcast all interblocks to all processors
  MPI_Bcast(&connections[0], connections.size(), MPI_interblock, ROOT, MPI_COMM_WORLD);
}

/* Function to set custom MPI datatypes to allow for easier data transmission */
void SetDataTypesMPI(const int &numEqn, MPI_Datatype &MPI_vec3d, MPI_Datatype &MPI_cellData, MPI_Datatype &MPI_procBlockInts, MPI_Datatype &MPI_interblock,
		     MPI_Datatype & MPI_DOUBLE_5INT ){
  // numEqn -- number of equations being solved
  // MPI_vec3d -- output MPI_Datatype for a vector3d<double>
  // MPI_cellData -- output MPI_Datatype for primVars or genArray
  // MPI_procBlockInts -- output MPI_Datatype for 14 INTs (14 INTs in procBlock class)
  // MPI_interblock -- output MPI_Datatype to send interblock class
  // MPI_DOUBLE_5INT -- output MPI_Datatype for a double followed by 5 ints

  //create vector3d<double> MPI datatype
  MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_vec3d);
  MPI_Type_commit(&MPI_vec3d);

  //create MPI datatype for states (primVars), residuals (genArray), etc
  MPI_Type_contiguous(numEqn, MPI_DOUBLE, &MPI_cellData);
  MPI_Type_commit(&MPI_cellData);

  //create MPI datatype for all the integers in the procBlock class
  MPI_Type_contiguous(15, MPI_INT, &MPI_procBlockInts);
  MPI_Type_commit(&MPI_procBlockInts);


  //create MPI datatype for a double followed by 5 ints
  int fieldCounts[2] = {1,5}; //number of entries per field
  MPI_Datatype fieldTypes[2] = {MPI_DOUBLE, MPI_INT}; //field types
  MPI_Aint displacement[2], lBound, ext;
  resid res; //dummy resid to get layout of class
  res.GetAddressesMPI(displacement);

  //make addresses relative to first field
  for ( int ii = 1; ii >= 0; ii-- ){
    displacement[ii] -= displacement[0];
  }
  MPI_Type_create_struct(2, fieldCounts, displacement, fieldTypes, &MPI_DOUBLE_5INT);

  //check that datatype has the correct extent, if it doesn't change the extent
  //this is necessary to portably send an array of this type
  MPI_Type_get_extent(MPI_DOUBLE_5INT, &lBound, &ext);
  if ( ext != sizeof(res) ){
    MPI_Datatype temp = MPI_DOUBLE_5INT;
    MPI_Type_create_resized(temp, 0, sizeof(res), &MPI_DOUBLE_5INT);
    MPI_Type_free(&temp);
  }

  MPI_Type_commit(&MPI_DOUBLE_5INT);


  //create MPI datatype for interblock class
  int counts[11] = {2,2,2,2,2,2,2,2,2,8,1}; //number of entries per field
  MPI_Datatype types[11] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_C_BOOL, MPI_INT}; //field types
  MPI_Aint disp[11], lowerBound, extent;
  interblock inter; //dummy interblock to get layout of class
  inter.GetAddressesMPI(disp);

  //make addresses relative to first field
  for ( int ii = 10; ii >= 0; ii-- ){
    disp[ii] -= disp[0];
  }
  MPI_Type_create_struct(11, counts, disp, types, &MPI_interblock);

  //check that datatype has the correct extent, if it doesn't change the extent
  //this is necessary to portably send an array of this type
  MPI_Type_get_extent(MPI_interblock, &lowerBound, &extent);
  if ( extent != sizeof(inter) ){
    MPI_Datatype temp = MPI_interblock;
    MPI_Type_create_resized(temp, 0, sizeof(inter), &MPI_interblock);
    MPI_Type_free(&temp);
  }

  MPI_Type_commit(&MPI_interblock);

}

/* Function to send procBlocks to their appropriate processor. This function is called after the decomposition has been run. The procBlock data all 
resides on the ROOT processor. In this function, the ROOT processor packs the procBlocks and sends them to the appropriate processor. All the non-ROOT
processors receive and unpack the data from ROOT. This is used to send the geometric block data from ROOT to all the processors at the beginning of the
simulation.
*/
vector<procBlock> SendProcBlocks( const vector<procBlock> &blocks, const int &rank, const int &numProcBlock, 
				  const MPI_Datatype &MPI_cellData, const MPI_Datatype &MPI_vec3d ){

  // blocks -- full vector of all procBlocks. This is only used on ROOT processor, all other processors just need a dummy variable to call the function
  // rank -- processor rank. Used to determine if process should send or receive
  // numProcBlock -- number of procBlocks that the processor should have. (All processors may give different values).
  // MPI_cellData -- MPI_Datatype used for primVars and genArray transmission
  // MPI_vec3d -- MPI_Datatype used for vector3d<double>  transmission

  vector<procBlock> localBlocks; //vector of procBlocks for each processor
  localBlocks.reserve(numProcBlock); //each processor may allocate for a different size

  //------------------------------------------------------------------------------------------------------------------------------------------------
  //                                                                ROOT
  //------------------------------------------------------------------------------------------------------------------------------------------------
  if ( rank == ROOT ){ // may have to pack and send data
    for ( unsigned int ii = 0; ii < blocks.size(); ii++ ){ //loop over ALL blocks

      if (blocks[ii].Rank() == ROOT ){ //no need to send data because it is already on ROOT processor
	localBlocks.push_back(blocks[ii]);
      }
      else{ //send data to receiving processors
	//pack and send procBlock
	blocks[ii].PackSendGeomMPI( MPI_cellData, MPI_vec3d );
      }
    }

  }
  //------------------------------------------------------------------------------------------------------------------------------------------------
  //                                                                NON - ROOT
  //------------------------------------------------------------------------------------------------------------------------------------------------
  else { // receive and unpack data (non-root)
    for ( int ii = 0; ii < numProcBlock; ii++ ){
      //recv and unpack procBlock
      procBlock tempBlock;
      tempBlock.RecvUnpackGeomMPI( MPI_cellData, MPI_vec3d );

      localBlocks.push_back(tempBlock); //add procBlock to output vector
    }
  }

  return localBlocks;

}


/* Function to send procBlocks to the root processor. In this function, the non-ROOT processors pack the procBlocks and send them to the ROOT processor. 
The ROOT processor receives and unpacks the data from the non-ROOT processors. This is used to get all the data on the ROOT processor to write out results.
*/
void GetProcBlocks( vector<procBlock> &blocks, const vector<procBlock> &localBlocks, const int &rank, const MPI_Datatype &MPI_cellData ){

  // blocks -- full vector of all procBlocks. This is only used on ROOT processor, all other processors just need a dummy variable to call the function
  // localBlocks -- procBlocks local to each processor. These are sent to ROOT
  // rank -- processor rank. Used to determine if process should send or receive
  // MPI_cellData -- MPI_Datatype used for primVars and genArray transmission

  //------------------------------------------------------------------------------------------------------------------------------------------------
  //                                                                ROOT
  //------------------------------------------------------------------------------------------------------------------------------------------------
  if ( rank == ROOT ){ // may have to recv and unpack data
    int locNum = 0;

    for ( unsigned int ii = 0; ii < blocks.size(); ii++ ){ //loop over ALL blocks

      if (blocks[ii].Rank() == ROOT ){ //no need to recv data because it is already on ROOT processor
	//assign local state block to global state block in order of local state vector
	blocks[ii] = localBlocks[locNum];
	locNum++;
      }
      else{ //recv data from sending processors
	blocks[ii].RecvUnpackSolMPI(MPI_cellData);
      }
    }

  }
  
  //------------------------------------------------------------------------------------------------------------------------------------------------
  //                                                                NON - ROOT
  //------------------------------------------------------------------------------------------------------------------------------------------------
  else { // pack and send data (non-root)
    for ( unsigned int ii = 0; ii < localBlocks.size(); ii++ ){
      localBlocks[ii].PackSendSolMPI(MPI_cellData);
    }
  }

}

/*function to broadcast a string from ROOT to all processors. This is needed because it is not garunteed in the MPI standard that the commmand
line arguments will be on any processor but ROOT.  */
void BroadcastString( string &str ){
  // str -- string to broadcast to all processors

  int strSize = str.size() + 1;  //get size of string (+1 for c_str end character)
  MPI_Bcast(&strSize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);  //broadcast string size

  char *buf = new char[strSize]; //allcate a char buffer of string size
  strcpy(buf, str.c_str()); //copy string into buffer
  MPI_Bcast(&buf[0], strSize, MPI_CHAR, ROOT, MPI_COMM_WORLD);  //broadcast string as char

  //create new string and assign to old string
  string newStr(buf, strSize - 1); //-1 to not include c_str end character
  str = newStr;

  delete [] buf; //deallocate buffer
}



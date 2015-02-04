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

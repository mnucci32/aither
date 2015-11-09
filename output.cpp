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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <utility>  // pair
#include <cmath>
#include "output.hpp"
#include "turbulence.hpp"
#include "vector3d.hpp"  // vector3d
#include "multiArray3d.hpp"  // multiArray3d
#include "tensor.hpp"    // tensor
#include "plot3d.hpp"    // plot3d
#include "eos.hpp"
#include "primVars.hpp"            // primVars
#include "procBlock.hpp"           // procBlock
#include "inviscidFlux.hpp"        // inviscidFlux
#include "input.hpp"               // inputVars
#include "boundaryConditions.hpp"  // decomposition
#include "resid.hpp"               // resid

#define VAROUT 15

using std::cout;
using std::endl;
using std::cerr;
using std::vector;
using std::string;
using std::ios;
using std::ofstream;
using std::to_string;
using std::max;
using std::pair;
using std::setw;
using std::setprecision;
using std::unique_ptr;

//-----------------------------------------------------------------------
// function declarations
// function to write out cell centers of grid in plot3d format
void WriteCellCenter(const string &gridName, const vector<procBlock> &vars,
                     const decomposition &decomp, const double &LRef) {
  // recombine procblocks into original configuration
  auto recombVars = Recombine(vars, decomp);

  // open binary output file
  ofstream outFile;
  string fEnd = "_center";
  string fPostfix = ".xyz";
  auto writeName = gridName + fEnd + fPostfix;
  outFile.open(writeName, ios::out | ios::binary);

  // check to see if file opened correctly
  if (outFile.fail()) {
    cerr << "ERROR: Grid file " << writeName << " did not open correctly!!!"
         << endl;
    exit(0);
  }

  // write number of blocks to file
  auto numBlks = static_cast<int>(recombVars.size());
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  // write i, j, k dimension for each block
  for (auto ll = 0; ll < numBlks; ll++) {  // loop over all blocks
    auto dumInt = recombVars[ll].NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = recombVars[ll].NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = recombVars[ll].NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
  }

  // write out x, y, z coordinates of cell centers
  for (auto ll = 0; ll < numBlks; ll++) {  // loop over all blocks
    auto maxi = recombVars[ll].NumI();
    auto maxj = recombVars[ll].NumJ();
    auto maxk = recombVars[ll].NumK();

    for (auto nn = 0; nn < 3; nn++) {  // loop over dimensions (3)
      for (auto kk = recombVars[ll].NumGhosts();
           kk < maxk + recombVars[ll].NumGhosts(); kk++) {
        for (auto jj = recombVars[ll].NumGhosts();
             jj < maxj + recombVars[ll].NumGhosts(); jj++) {
          for (auto ii = recombVars[ll].NumGhosts();
               ii < maxi + recombVars[ll].NumGhosts(); ii++) {
            // get the cell center coordinates (dimensionalized)
            auto dumVec = recombVars[ll].Center(ii, jj, kk) * LRef;

            // for a given block, first write out all x coordinates, then all y
            // coordinates, then all z coordinates
            auto dumDouble = dumVec[nn];
            // write to file
            outFile.write(reinterpret_cast<char *>(&dumDouble),
                          sizeof(dumDouble));
          }
        }
      }
    }
  }

  // close output file
  outFile.close();
}

//----------------------------------------------------------------------
// function to write out variables in function file format
void WriteFun(const vector<procBlock> &vars, const idealGas &eqnState,
              const sutherland &suth, const int &solIter,
              const decomposition &decomp, const input &inp,
              const unique_ptr<turbModel> &turb) {
  // define reference speed of sound
  auto refSoS = eqnState.SoS(inp.PRef(), inp.RRef());

  auto recombVars = Recombine(vars, decomp);

  // open binary plot3d function file
  ofstream outFile;
  string fEnd = "_center";
  string fPostfix = ".fun";
  auto writeName = inp.SimNameRoot() + "_" + to_string(solIter) + fEnd +
      fPostfix;
  outFile.open(writeName, ios::out | ios::binary);

  // check to see if file opened correctly
  if (outFile.fail()) {
    cerr << "ERROR: Function file " << writeName << " did not open correctly!!!"
         << endl;
    exit(0);
  }

  // write number of blocks to file
  auto numBlks = static_cast<int>(recombVars.size());
  outFile.write(reinterpret_cast<char *>(&numBlks), sizeof(numBlks));

  // write i, j, k, recombVars dimension for each block
  auto numVars = VAROUT;  // number of variables to write out

  // loop over all blocks and write out imax, jmax, kmax, numVars
  for (auto ll = 0; ll < numBlks; ll++) {
    auto dumInt = recombVars[ll].NumI();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = recombVars[ll].NumJ();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));
    dumInt = recombVars[ll].NumK();
    outFile.write(reinterpret_cast<char *>(&dumInt), sizeof(dumInt));

    outFile.write(reinterpret_cast<char *>(&numVars), sizeof(numVars));
  }

  // write out variables
  for (auto ll = 0; ll < numBlks; ll++) {  // loop over all blocks
    multiArray3d<double> dumArr(recombVars[ll].NumI(), recombVars[ll].NumJ(),
                                recombVars[ll].NumK());

    // loop over the number of variables to write out
    for (auto vv = 0; vv < numVars; vv++) {
      // store nondimensional variable in dumArr for a given block in order.
      // i.e. var1 var2 var3 etc
      if (vv == 0) {  // density
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).Rho();
            }
          }
        }
      } else if (vv == 1) {  // vel-x
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).U();
            }
          }
        }
      } else if (vv == 2) {  // vel-y
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).V();
            }
          }
        }
      } else if (vv == 3) {  // vel-z
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).W();
            }
          }
        }
      } else if (vv == 4) {  // pressure
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).P();
            }
          }
        }
      } else if (vv == 5) {  // mach
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              vector3d<double> vel =
                  recombVars[ll].State(ii.g, jj.g, kk.g).Velocity();
              dumArr(ii.p, jj.p, kk.p) = vel.Mag() /
                  eqnState.SoS(recombVars[ll].State(ii.g, jj.g, kk.g).P(),
                               recombVars[ll].State(ii.g, jj.g, kk.g).Rho());
            }
          }
        }
      } else if (vv == 6) {  // speed of sound
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  eqnState.SoS(recombVars[ll].State(ii.g, jj.g, kk.g).P(),
                               recombVars[ll].State(ii.g, jj.g, kk.g).Rho());
            }
          }
        }
      } else if (vv == 7) {  // time step - no ghost cells
        for (int kk = 0; kk < recombVars[ll].NumK(); kk++) {
          for (int jj = 0; jj < recombVars[ll].NumJ(); jj++) {
            for (int ii = 0; ii < recombVars[ll].NumI(); ii++) {
              dumArr(ii, jj, kk) = recombVars[ll].Dt(ii, jj, kk);
            }
          }
        }
      } else if (vv == 8) {  // temperature
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).Temperature(eqnState);
            }
          }
        }
      } else if (vv == 9) {  // processor rank
        for (auto kk = 0; kk < recombVars[ll].NumK(); kk++) {
          for (auto jj = 0; jj < recombVars[ll].NumJ(); jj++) {
            for (auto ii = 0; ii < recombVars[ll].NumI(); ii++) {
              dumArr(ii, jj, kk) = vars[SplitBlockNumber(recombVars, decomp, ll,
                                                         ii, jj, kk)].Rank();
            }
          }
        }
      } else if (vv == 10) {  // global position
        for (auto kk = 0; kk < recombVars[ll].NumK(); kk++) {
          for (auto jj = 0; jj < recombVars[ll].NumJ(); jj++) {
            for (auto ii = 0; ii < recombVars[ll].NumI(); ii++) {
              dumArr(ii, jj, kk) =
                  vars[SplitBlockNumber(recombVars, decomp, ll,
                                        ii, jj, kk)].GlobalPos();
            }
          }
        }
      } else if (vv == 11) {  // viscosity ratio
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  turb->EddyViscNoLim(recombVars[ll].State(ii.g, jj.g, kk.g)) /
                  suth.Viscosity(recombVars[ll].State(ii.g, jj.g, kk.g).
                                 Temperature(eqnState));
            }
          }
        }
      } else if (vv == 12) {  // tke
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).Tke();
            }
          }
        }
      } else if (vv == 13) {  // omega
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].State(ii.g, jj.g, kk.g).Omega();
            }
          }
        }
      } else if (vv == 14) {  // wall distance
        for (struct {int p; int g;} kk = {0, recombVars[ll].NumGhosts()};
             kk.p < recombVars[ll].NumK(); kk.g++, kk.p++) {
          for (struct {int p; int g;} jj = {0, recombVars[ll].NumGhosts()};
               jj.p < recombVars[ll].NumJ(); jj.g++, jj.p++) {
            for (struct {int p; int g;} ii = {0, recombVars[ll].NumGhosts()};
                 ii.p < recombVars[ll].NumI(); ii.g++, ii.p++) {
              dumArr(ii.p, jj.p, kk.p) =
                  recombVars[ll].WallDist(ii.g, jj.g, kk.g);
            }
          }
        }
      } else {
        cerr << "ERROR: Variable to write to function file is not defined!"
             << endl;
        exit(0);
      }

      // write out dimensional variables -- loop over block length
      for (auto kk = 0; kk < recombVars[ll].NumK(); kk++) {
        for (auto jj = 0; jj < recombVars[ll].NumJ(); jj++) {
          for (auto ii = 0; ii < recombVars[ll].NumI(); ii++) {
            auto dumDouble = dumArr(ii, jj, kk);

            if (vv == 0) {  // density
              dumDouble = dumDouble * inp.RRef();
            } else if (vv == 1) {  // velocity x
              dumDouble = dumDouble * refSoS;
            } else if (vv == 2) {  // velocity y
              dumDouble = dumDouble * refSoS;
            } else if (vv == 3) {  // velocity z
              dumDouble = dumDouble * refSoS;
            } else if (vv == 4) {  // pressure
              dumDouble = dumDouble * inp.RRef() * refSoS * refSoS;
            } else if (vv == 5) {  // mach is already nondimensional
              dumDouble = dumDouble;
            } else if (vv == 6) {  // speed of sound
              dumDouble = dumDouble * refSoS;
            } else if (vv == 7) {  // time step
              dumDouble = dumDouble / refSoS * inp.LRef();
            } else if (vv == 8) {  // temperature
              dumDouble = dumDouble * inp.TRef();
            } else if (vv == 12) {  // tke
              dumDouble = dumDouble * refSoS * refSoS;
            } else if (vv == 13) {  // omega
              dumDouble = dumDouble * refSoS * refSoS * inp.RRef() /
                  suth.MuRef();
            } else if (vv == 14) {  // wall distance
              dumDouble = dumDouble * inp.LRef();
            }
            outFile.write(reinterpret_cast<char *>(&dumDouble),
                          sizeof(dumDouble));
          }
        }
      }
    }
  }

  // close plot3d function file
  outFile.close();
}

// function to write out results file for ensight
void WriteRes(const string &gridName, const int &iter, const int &outFreq) {
  // open results file
  ofstream resFile;
  string fResPostfix = ".res";
  string fPostfix = ".fun";
  string fEnd = "_center";
  auto resName = gridName + fEnd + fResPostfix;
  resFile.open(resName, ios::out);

  auto writeName = gridName + "_*" + fEnd + fPostfix;

  // check to see if file opened correctly
  if (resFile.fail()) {
    cerr << "ERROR: Results file " << resName << " did not open correctly!!!"
         << endl;
    exit(0);
  }

  // write number of scalars and number of vectors
  constexpr auto numScalar = VAROUT;
  auto numVector = 1;
  resFile << numScalar << "     " << numVector << "     " << 0 << endl;

  // write number of time points that there is solution data at
  auto numTime = iter / outFreq;
  resFile << numTime << endl;

  // Write solution times or iteration numbers
  auto solTime = 0;
  auto count = 1;
  for (auto ii = 0; ii < numTime; ii++) {
    solTime += outFreq;
    if (count % 10 == 0) {
      resFile << endl;
    }
    resFile << solTime << "   ";
    count++;
  }
  resFile << endl;

  // Write starting iteration and iteration increment
  resFile << outFreq << "  " << outFreq << endl;

  // Write out variables
  resFile << writeName << " F 0001 density" << endl;
  resFile << writeName << " F 0002 Vx" << endl;
  resFile << writeName << " F 0003 Vy" << endl;
  resFile << writeName << " F 0004 Vz" << endl;
  resFile << writeName << " F 0005 pressure" << endl;
  resFile << writeName << " F 0006 mach" << endl;
  resFile << writeName << " F 0007 sos" << endl;
  resFile << writeName << " F 0008 dt" << endl;
  resFile << writeName << " F 0009 temperature" << endl;
  resFile << writeName << " F 0010 procRank" << endl;
  resFile << writeName << " F 0011 procBlockID" << endl;
  resFile << writeName << " F 0012 viscRatio" << endl;
  resFile << writeName << " F 0013 tke" << endl;
  resFile << writeName << " F 0014 omega" << endl;
  resFile << writeName << " F 0015 wallDistance" << endl;
  resFile << writeName << " F 0002 0003 0004 velocity" << endl;

  // Close results file
  resFile.close();
}

// function to write out residual information
void WriteResiduals(const input &inp, genArray &residL2First, genArray &residL2,
                    const resid &residLinf, const double &matrixResid,
                    const int &nn, const int &mm) {
  // determine normalization
  if (nn == 0 && mm == 0) {  // if at first iteration, normalize by itself
    residL2First = residL2;
  // if within first 5 iterations reset normalization
  } else if ((nn < 5) && mm == 0) {
    for (auto cc = 0; cc < NUMVARS; cc++) {
      if (residL2[cc] > residL2First[cc]) {
        residL2First[cc] = residL2[cc];
      }
    }
  }

  // normalize residuals
  residL2 = (residL2 + EPS) / (residL2First + EPS);

  // write out column headers every 100 iterations
  if (nn % 100 == 0 && mm == 0) {
    cout << std::left << setw(7) << "Step" << setw(8) << "NL-Iter";
    if (inp.Dt() > 0.0) {
      cout << std::left << setw(12) << "Time-Step";
    } else if (inp.CFL() > 0.0) {
      cout << std::left << setw(12) << "CFL";
    }
    cout << std::left << setw(12) << "Res-Mass" << setw(12)
         << "Res-Mom-X" << setw(12) << "Res-Mom-Y" << setw(12) << "Res-Mom-Z"
         << setw(12) << "Res-Energy";
    if (inp.IsTurbulent()) {
      cout << std::left << setw(12) << "Res-Tke" << setw(12) << "Res-Omega";
    }
    cout << std::left << setw(8) << "Max-Eqn" << setw(8)
         << "Max-Blk" << setw(8) << "Max-I" << setw(8)
         << "Max-J" << setw(8) << "Max-K" << setw(12) << "Max-Res"
         << setw(12) << "Res-Matrix" << endl;
  }

  cout << std::left << setw(7) << nn << setw(8) << mm;
  if (inp.Dt() > 0.0) {
    cout << std::left << setw(12) << setprecision(4) << std::scientific
         << inp.Dt();
  } else if (inp.CFL() > 0.0) {
    cout << std::left << setw(12) << setprecision(4) << std::scientific
         << inp.CFL();
  }
  cout << std::left << setw(12) << residL2[0] << setw(12) << residL2[1]
       << setw(12) << residL2[2] << setw(12) << residL2[3] << setw(12)
       << residL2[4];
  if (inp.IsTurbulent()) {
    cout << std::left << setw(12) << residL2[5] << setw(12) << residL2[6];
  }
  cout.unsetf(std::ios::fixed | std::ios::scientific);
  cout << std::left << setw(8) << residLinf.Eqn() << setw(8)
       << residLinf.Block() << setw(8) << residLinf.ILoc() << setw(8)
       << residLinf.JLoc() << setw(8) << residLinf.KLoc() << setw(12)
       << setprecision(4) << std::scientific << residLinf.Linf() << setw(12)
       << matrixResid << endl;

  cout.unsetf(std::ios::fixed | std::ios::scientific);
}

/*Function to take in a vector of split procBlocks and return a vector of joined
 * procblocks (in their original configuration before grid decomposition).*/
vector<procBlock> Recombine(const vector<procBlock> &vars,
                            const decomposition &decomp) {
  // vars -- vector of split procBlocks
  // decomp -- decomposition

  auto recombVars = vars;
  vector<boundarySurface> dumSurf;
  for (auto ii = decomp.NumSplits() - 1; ii >= 0; ii--) {
    // recombine blocks and resize vector
    recombVars[decomp.SplitHistBlkLower(ii)]
        .Join(recombVars[decomp.SplitHistBlkUpper(ii)], decomp.SplitHistDir(ii),
              dumSurf);
    recombVars.resize(recombVars.size() - 1);
  }

  return recombVars;
}

/*Function to take in indices from the recombined procBlocks and determine which
 * split procBlock index the cell is associated with.*/
int SplitBlockNumber(const vector<procBlock> &vars, const decomposition &decomp,
                     const int &blk, const int &ii, const int &jj,
                     const int &kk) {
  // vars -- vector of recombined procblocks
  // decomp -- decomposition data structure
  // blk -- block number
  // ii -- i index of cell in recombined block to find split block number
  // jj -- j index of cell in recombined block to find split block number
  // kk -- k index of cell in recombined block to find split block number

  // Get block dimensions (both lower and upper extents)
  vector<pair<vector3d<int>, vector3d<int>>> blkDims(vars.size());
  vector3d<int> initialLower(0, 0, 0);
  for (auto bb = 0; bb < static_cast<int>(blkDims.size()); bb++) {
    vector3d<int> dims(vars[bb].NumI(), vars[bb].NumJ(), vars[bb].NumK());
    blkDims[bb].first = initialLower;
    blkDims[bb].second = dims;
  }

  auto ind = blk;

  // no splits, cell must be in parent block already
  if (decomp.NumSplits() == 0) {
    return ind;
  } else {  // cell is in lower split already
    for (auto ss = 0; ss < decomp.NumSplits(); ss++) {  // loop over all splits
      // wrong parent block - split won't effect search so use dummy value
      if (blk != decomp.ParentBlock(ss + vars.size())) {
        pair<vector3d<int>, vector3d<int>> dumBlk(initialLower, initialLower);
        blkDims.push_back(dumBlk);
      } else {
        // "split" blocks - change lower limits of block
        if (decomp.SplitHistDir(ss) == "i") {
          pair<vector3d<int>, vector3d<int>> splitBlk =
              blkDims[decomp.SplitHistBlkLower(ss)];
          splitBlk.first[0] += decomp.SplitHistIndex(ss);
          blkDims.push_back(splitBlk);
        } else if (decomp.SplitHistDir(ss) == "j") {
          pair<vector3d<int>, vector3d<int>> splitBlk =
              blkDims[decomp.SplitHistBlkLower(ss)];
          splitBlk.first[1] += decomp.SplitHistIndex(ss);
          blkDims.push_back(splitBlk);
        } else {  // direction is k
          pair<vector3d<int>, vector3d<int>> splitBlk =
              blkDims[decomp.SplitHistBlkLower(ss)];
          splitBlk.first[2] += decomp.SplitHistIndex(ss);
          blkDims.push_back(splitBlk);
        }

        // test to see if block is in upper split
        if (!(ii <= blkDims[decomp.SplitHistBlkUpper(ss)].second.X() &&
              jj <= blkDims[decomp.SplitHistBlkUpper(ss)].second.Y() &&
              kk <= blkDims[decomp.SplitHistBlkUpper(ss)].second.Z() &&
              ii >= blkDims[decomp.SplitHistBlkUpper(ss)].first.X() &&
              jj >= blkDims[decomp.SplitHistBlkUpper(ss)].first.Y() &&
              kk >= blkDims[decomp.SplitHistBlkUpper(ss)].first.Z())) {
          // cell not in upper split, but in lower split - found block index
          return decomp.SplitHistBlkLower(ss);
        } else {  // cell in upper split (and lower split)
          ind = decomp.SplitHistBlkUpper(ss);
        }
      }
    }
  }

  return ind;  // cell was in uppermost split for given parent block
}

#   This file is part of aither.
#   Copyright (C) 2015-17  Michael Nucci (michael.nucci@gmail.com)
#
#   Aither is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Aither is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
#   This script runs regression tests to test builds on linux and osx for
#   travis ci.

import os
import optparse
import shutil
import sys
import datetime
import subprocess

class regressionTest:
    def __init__(self):
        self.caseName = "none"
        self.iterations = 100
        self.procs = 1
        self.residuals = [1.0, 1.0, 1.0, 1.0, 1.0]
        self.ignoreIndices = []
        self.location = os.getcwd()
        self.runDirectory = "."
        self.aitherPath = "aither"
        self.mpirunPath = "mpirun"
        self.percentTolerance = 0.01
        self.isRestart = False
        self.restartFile = "none"
        
    def SetRegressionCase(self, name):
        self.caseName = name
        
    def SetNumberOfIterations(self, num):
        self.iterations = num
        
    def SetNumberOfProcessors(self, num):
        self.procs = num
        
    def SetResiduals(self, resid):
        self.residuals = resid
        
    def SetRunDirectory(self, path):
        self.runDirectory = path

    def SetAitherPath(self, path):
        self.aitherPath = path

    def SetMpirunPath(self, path):
        self.mpirunPath = path
        
    def SetIgnoreIndices(self, ind):
        self.ignoreIndices.append(ind)
        
    def SetPercentTolerance(self, per):
        self.percentTolerance = per

    def GoToRunDirectory(self):
        os.chdir(self.runDirectory)

    def SetRestart(self, resFlag):
        self.isRestart = resFlag

    def SetRestartFile(self, resFile):
        self.restartFile = resFile
        
    def ReturnToHomeDirectory(self):
        os.chdir(self.location)
        
    def GetTestCaseResiduals(self):
        fname = self.caseName + ".resid"
        file = open(fname, "r")
        lastLine = file.readlines()[-1]
        file.close()
        tokens = lastLine.split()
        resids = [float(ii) for ii in tokens[3:3+len(self.residuals)]]
        return resids

    def CompareResiduals(self, returnCode):
        testResids = self.GetTestCaseResiduals()
        resids = []
        truthResids = []
        for ii in range(0, len(testResids)):
            if ii not in self.ignoreIndices:
                resids.append(testResids[ii])
                truthResids.append(self.residuals[ii])
        if (returnCode == 0):
            passing = [abs(resid - truthResids[ii]) <= self.percentTolerance * truthResids[ii]
                       for ii, resid in enumerate(resids)]
        else:
            passing = [False for ii in resids]
        return passing, resids
        
    def GetResiduals(self):
        return self.residuals
        
    # change input file to have number of iterations specified for test
    def ModifyInputFile(self):
        fname = self.caseName + ".inp"
        fnameBackup = fname + ".old"
        shutil.move(fname, fnameBackup)
        with open(fname, "w") as fout:
            with open(fnameBackup, "r") as fin:
                for line in fin:
                    if "iterations:" in line:
                        fout.write("iterations: " + str(self.iterations) + "\n")
                    elif "outputFrequency:" in line:
                        fout.write("outputFrequency: " + str(self.iterations) + "\n")
                    else:
                        fout.write(line)
                
    # modify the input file and run the test
    def RunCase(self):
        self.GoToRunDirectory()
        print("---------- Starting Test:", self.caseName, "----------")
        print("Current directory:", os.getcwd())
        print("Modifying input file...")
        self.ModifyInputFile()
        if self.isRestart:
            cmd = self.mpirunPath + " -np " + str(self.procs) + " " + self.aitherPath \
                  + " " + self.caseName + ".inp " + self.restartFile + " > " + self.caseName \
                  + ".out"
        else:
            cmd = self.mpirunPath + " -np " + str(self.procs) + " " + self.aitherPath \
                  + " " + self.caseName + ".inp > " + self.caseName + ".out"
        print(cmd)
        start = datetime.datetime.now()
        process = subprocess.Popen(cmd, shell=True)
        returnCode = process.wait()
        if (returnCode == 0):
            print("Simulation completed with no errors")
        else:
            print("ERROR: Simulation terminated with errors")
        duration = datetime.datetime.now() - start
        
        # test residuals for pass/fail
        passed, resids = self.CompareResiduals(returnCode)
        if (all(passed)):
            print("All tests for", self.caseName, "passed!")
        else:
            print("Tests for", self.caseName, "failed!")
            print("Residuals should be:", self.GetResiduals())
            print("Residuals are:", resids)        
            
        print("Test Duration:",duration)
        print("---------- End Test:", self.caseName, "----------")
        self.ReturnToHomeDirectory()
        return passed
        
        
def main():
    # Set up options
    parser = optparse.OptionParser()
    parser.add_option("-a", "--aitherPath", action="store", dest="aitherPath",
                      default="aither", help="Path to aither executable.")
    parser.add_option("-o", "--operatingSystem", action="store",
                      dest="operatingSystem", default="linux",
                      help="Operating system that tests will run on [linux/osx]")
    parser.add_option("-m", "--mpirunPath", action="store",
                      dest="mpirunPath", default="",
                      help="Path to mpirun")
                      
    options, remainder = parser.parse_args()

    # travis osx images have 1 proc, ubuntu have 2
    if (options.operatingSystem == "linux"):
        maxProcs = 2
    else:
        maxProcs = 1
        
    numIterations = 100
    numIterationsRestart = 50
    totalPass = True
    
    # ------------------------------------------------------------------
    # Regression tests
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # subsonic cylinder
    # laminar, inviscid, lu-sgs
    subCyl = regressionTest()
    subCyl.SetRegressionCase("subsonicCylinder")
    subCyl.SetAitherPath(options.aitherPath)
    subCyl.SetRunDirectory("subsonicCylinder")
    subCyl.SetNumberOfProcessors(1)
    subCyl.SetNumberOfIterations(numIterations)
    subCyl.SetResiduals([1.5394e-1, 1.4989e-1, 1.5909e-1, 8.1415e-1, 1.5295e-1])
    subCyl.SetIgnoreIndices(3)
    subCyl.SetMpirunPath(options.mpirunPath)
    
    # run regression case
    passed = subCyl.RunCase()   
    totalPass = totalPass and all(passed)
        
    # ------------------------------------------------------------------
    # multi-block subsonic cylinder
    # laminar, inviscid, lusgs, multi-block
    multiCyl = regressionTest()
    multiCyl.SetRegressionCase("multiblockCylinder")
    multiCyl.SetAitherPath(options.aitherPath)
    multiCyl.SetRunDirectory("multiblockCylinder")
    multiCyl.SetNumberOfProcessors(maxProcs)
    multiCyl.SetNumberOfIterations(numIterations)
    if (options.operatingSystem == "linux"):
        multiCyl.SetResiduals([2.3188e-1, 2.9621e-1, 4.5868e-1, 1.2813, 2.3009e-1])
    else:
        multiCyl.SetResiduals([2.3188e-1, 2.9621e-1, 4.5868e-1, 1.2813, 2.3009e-1])
    multiCyl.SetIgnoreIndices(3)
    multiCyl.SetMpirunPath(options.mpirunPath)
    
    # run regression case
    passed = multiCyl.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # sod shock tube
    # laminar, inviscid, bdf2, weno
    shockTube = regressionTest()
    shockTube.SetRegressionCase("shockTube")
    shockTube.SetAitherPath(options.aitherPath)
    shockTube.SetRunDirectory("shockTube")
    shockTube.SetNumberOfProcessors(1)
    shockTube.SetNumberOfIterations(numIterations)
    shockTube.SetResiduals([5.0503e-1, 4.4569e-1, 1.0e0, 1.0e0, 2.6181e-1])
    shockTube.SetIgnoreIndices(2)
    shockTube.SetIgnoreIndices(3)
    shockTube.SetMpirunPath(options.mpirunPath)
    
    # run regression case
    passed = shockTube.RunCase()   
    totalPass = totalPass and all(passed)
        
    # ------------------------------------------------------------------
    # sod shock tube restart
    # laminar, inviscid, bdf2, weno
    shockTubeRestart = shockTube
    shockTubeRestart.SetNumberOfIterations(numIterationsRestart)
    shockTubeRestart.SetRestart(True)
    shockTubeRestart.SetRestartFile("shockTube_50.rst")

    # run regression case
    passed = shockTubeRestart.RunCase()   
    totalPass = totalPass and all(passed)
        
    # ------------------------------------------------------------------
    # supersonic wedge
    # laminar, inviscid, explicit euler
    supWedge = regressionTest()
    supWedge.SetRegressionCase("supersonicWedge")
    supWedge.SetAitherPath(options.aitherPath)
    supWedge.SetRunDirectory("supersonicWedge")
    supWedge.SetNumberOfProcessors(1)
    supWedge.SetNumberOfIterations(numIterations)
    supWedge.SetResiduals([4.1813e-1, 4.2549e-1, 3.6525e-1, 3.8013e-1, 4.0998e-1])
    supWedge.SetIgnoreIndices(3)
    supWedge.SetMpirunPath(options.mpirunPath)
        
    # run regression case
    passed = supWedge.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # transonic bump in channel
    # laminar, inviscid, dplur
    transBump = regressionTest()
    transBump.SetRegressionCase("transonicBump")
    transBump.SetAitherPath(options.aitherPath)
    transBump.SetRunDirectory("transonicBump")
    transBump.SetNumberOfProcessors(1)
    transBump.SetNumberOfIterations(numIterations)
    transBump.SetResiduals([1.1839e-1, 6.8615e-2, 8.4925e-2, 1.0398, 9.9669e-2])        
    transBump.SetIgnoreIndices(3)
    transBump.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = transBump.RunCase()
    totalPass = totalPass and all(passed)
        
    # ------------------------------------------------------------------
    # viscous flat plate
    # laminar, viscous, lu-sgs
    viscPlate = regressionTest()
    viscPlate.SetRegressionCase("viscousFlatPlate")
    viscPlate.SetAitherPath(options.aitherPath)
    viscPlate.SetRunDirectory("viscousFlatPlate")
    viscPlate.SetNumberOfProcessors(maxProcs)
    viscPlate.SetNumberOfIterations(numIterations)
    if (options.operatingSystem == "linux"):
        viscPlate.SetResiduals([7.7265e-2, 2.4712e-1, 5.6413e-2, 1.0228, 7.9363e-2])
    else:
        viscPlate.SetResiduals([7.6468e-2, 2.4713e-1, 4.0109e-2, 9.8730e-1, 7.9237e-2])
    viscPlate.SetIgnoreIndices(3)
    viscPlate.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = viscPlate.RunCase()
    totalPass = totalPass and all(passed)        

    # ------------------------------------------------------------------
    # turbulent flat plate
    # viscous, lu-sgs, k-w wilcox
    turbPlate = regressionTest()
    turbPlate.SetRegressionCase("turbFlatPlate")
    turbPlate.SetAitherPath(options.aitherPath)
    turbPlate.SetRunDirectory("turbFlatPlate")
    turbPlate.SetNumberOfProcessors(maxProcs)
    turbPlate.SetNumberOfIterations(numIterations)
    if (options.operatingSystem == "linux"):
        turbPlate.SetResiduals([4.1174e-2, 4.2731e-2, 1.0641, 8.3686e-2, 3.9585e-2,
                                4.5098e-8, 1.1416e-5])
    else:
        turbPlate.SetResiduals([3.9338e-2, 4.2745e-2, 1.0167, 7.4604e-2, 3.8146e-2,
                                4.7610e-8, 1.1583e-5])
    turbPlate.SetIgnoreIndices(2)
    turbPlate.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = turbPlate.RunCase()
    totalPass = totalPass and all(passed)        

    # ------------------------------------------------------------------
    # rae2822
    # turbulent, k-w sst, c-grid
    rae2822 = regressionTest()
    rae2822.SetRegressionCase("rae2822")
    rae2822.SetAitherPath(options.aitherPath)
    rae2822.SetRunDirectory("rae2822")
    rae2822.SetNumberOfProcessors(maxProcs)
    rae2822.SetNumberOfIterations(numIterations)
    if (options.operatingSystem == "linux"):
        rae2822.SetResiduals([6.3720e-1, 1.0432, 6.1385e-1, 4.8859e-1, 5.8621e-1,
                              2.5389e-5, 4.6321e-5])
    else:
        rae2822.SetResiduals([6.3424e-1, 1.0518, 6.1904e-1, 6.0576e-1, 5.8718e-1,
                              2.5387e-5, 4.6450e-5])
    rae2822.SetIgnoreIndices(3)
    rae2822.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = rae2822.RunCase()
    totalPass = totalPass and all(passed)        
    
    # ------------------------------------------------------------------
    # regression test overall pass/fail
    # ------------------------------------------------------------------
    if (totalPass):
        print("All tests passed!")
        sys.exit(0)
    else:
        print("ERROR: Some tests failed")
        sys.exit(1)
        
    
if __name__ == "__main__":
    main()

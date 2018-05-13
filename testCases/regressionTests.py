#   This file is part of aither.
#   Copyright (C) 2015-18  Michael Nucci (michael.nucci@gmail.com)
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
#   This script runs regression tests to test builds on linux and macOS for
#   travis ci, and windows for appveyor

import os
import optparse
import shutil
import sys
import datetime
import subprocess
import time

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
        self.passedStatus = "none"
        self.isProfile = False

    def SetRegressionCase(self, name):
        self.caseName = name

    def SetNumberOfIterations(self, num):
        self.iterations = num

    def SetNumberOfProcessors(self, num):
        self.procs = num

    def Processors(self):
        return self.procs

    def PassedStatus(self):
        return self.passedStatus

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

    def SetProfile(self, profFlag):
        self.isProfile = profFlag

    def SetRestartFile(self, resFile):
        self.restartFile = resFile

    def ReturnToHomeDirectory(self):
        os.chdir(self.location)

    def GetTestCaseResiduals(self):
        fname = self.caseName + ".resid"
        rfile = open(fname, "r")
        lastLine = rfile.readlines()[-1]
        rfile.close()
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
        return passing, resids, truthResids

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
                    elif "restartFrequency:" in line and self.isProfile:
                        fout.write("restartFrequency: " + str(self.iterations) + "\n")
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
        interval = start
        process = subprocess.Popen(cmd, shell=True)
        while process.poll() is None:
            current = datetime.datetime.now()
            if (current - interval).total_seconds() > 60.:
                print("----- Run Time: %s -----" % (current - start))
                interval = current
            time.sleep(0.5)
        returnCode = process.poll()

        if (returnCode == 0):
            print("Simulation completed with no errors")
            # test residuals for pass/fail
            if not self.isProfile:
                passed, resids, truth = self.CompareResiduals(returnCode)
                if all(passed):
                    print("All tests for", self.caseName, "PASSED!")
                    self.passedStatus = "PASSED"
                else:
                    print("Tests for", self.caseName, "FAILED!")
                    print("Residuals should be:", truth)
                    print("Residuals are:", resids)
                    self.passedStatus = "MISMATCH"
            else:
              passed = [True]
              self.passedStatus = "PROFILE"
        else:
            print("ERROR: Simulation terminated with errors")
            self.passedStatus = "ERRORS"
        duration = datetime.datetime.now() - start

        print("Test Duration:", duration)
        print("---------- End Test:", self.caseName, "----------")
        print("")
        print("")
        self.ReturnToHomeDirectory()
        return passed


def main():
    # Set up options
    parser = optparse.OptionParser()
    parser.add_option("-a", "--aitherPath", action="store", dest="aitherPath",
                      default="aither",
                      help="Path to aither executable. Default = aither")
    parser.add_option("-o", "--operatingSystem", action="store",
                      dest="operatingSystem", default="linux",
                      help="Operating system that tests will run on [linux/macOS/windows]. Default = linux")
    parser.add_option("-m", "--mpirunPath", action="store",
                      dest="mpirunPath", default="mpirun",
                      help="Path to mpirun. Default = mpirun")
    parser.add_option("-b", "--build", action="store",
                      dest="build", default="release",
                      help="build type used in compilation. Default = release")

    options, remainder = parser.parse_args()

    # travis macOS images have 1 proc, ubuntu have 2
    # appveyor windows images have 2 procs
    maxProcs = 2
    if (options.operatingSystem == "macOS"):
        maxProcs = 1

    isProfile = options.build == "debug"
    numIterations = 100
    numIterationsShort = 20
    numIterationsRestart = 50
    if isProfile:
      numIterations = 1
      numIterationsShort = 1
      numIterationsRestart = 1

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
    subCyl.SetProfile(isProfile)
    subCyl.SetNumberOfProcessors(1)
    subCyl.SetNumberOfIterations(numIterations)
    subCyl.SetResiduals(
        [1.8751e-01, 2.6727e-01, 3.1217e-01, 7.9662e-01, 1.8639e-01])
    subCyl.SetIgnoreIndices(3)
    subCyl.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = subCyl.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # multi-block subsonic cylinder
    # laminar, inviscid, lusgs, multi-block, ausmpw+
    multiCyl = regressionTest()
    multiCyl.SetRegressionCase("multiblockCylinder")
    multiCyl.SetAitherPath(options.aitherPath)
    multiCyl.SetRunDirectory("multiblockCylinder")
    multiCyl.SetProfile(isProfile)
    multiCyl.SetNumberOfProcessors(maxProcs)
    multiCyl.SetNumberOfIterations(numIterations)
    multiCyl.SetResiduals(
        [2.0529e-01, 3.4540e-01, 5.0153e-01, 1.0180e+00, 1.9997e-01])
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
    shockTube.SetProfile(isProfile)
    shockTube.SetNumberOfProcessors(1)
    shockTube.SetNumberOfIterations(numIterations)
    shockTube.SetResiduals(
        [4.8537e-01, 4.5855e-01, 1.0000e+00, 1.0000e+00, 2.6434e-01])
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
    shockTubeRestart.SetRestartFile("shockTube_" + str(numIterationsRestart) + ".rst")

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
    supWedge.SetProfile(isProfile)
    supWedge.SetNumberOfProcessors(1)
    supWedge.SetNumberOfIterations(numIterations)
    supWedge.SetResiduals([4.1813e-1, 4.2549e-1, 3.6525e-1, 3.9971e-1, 4.0998e-1])
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
    transBump.SetProfile(isProfile)
    transBump.SetNumberOfProcessors(1)
    transBump.SetNumberOfIterations(numIterations)
    transBump.SetResiduals([1.1901e-01, 7.0606e-02, 8.4288e-02, 1.0000e+00, 1.0032e-01])
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
    viscPlate.SetProfile(isProfile)
    viscPlate.SetNumberOfProcessors(maxProcs)
    viscPlate.SetNumberOfIterations(numIterations)
    if viscPlate.Processors() == 2:
        viscPlate.SetResiduals(
            [7.6770e-02, 2.4712e-01, 5.2446e-02, 1.0000e+00, 7.9490e-02])
    else:
        viscPlate.SetResiduals(
            [7.4673e-02, 2.4711e-01, 3.8960e-02, 1.0000e+00, 7.7683e-02])
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
    turbPlate.SetProfile(isProfile)
    turbPlate.SetNumberOfProcessors(maxProcs)
    turbPlate.SetNumberOfIterations(numIterationsShort)
    if turbPlate.Processors() == 2:
        turbPlate.SetResiduals([2.2801e-01, 2.9863e-01, 1.0000e+00, 3.2381e-01,
                                2.2326e-01, 2.5206e-07, 3.3015e-06])
    else:
        turbPlate.SetResiduals([2.2309e-01, 2.9862e-01, 1.0000e+00, 3.2376e-01,
                                2.1910e-01, 2.5208e-07, 3.3009e-06])
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
    rae2822.SetProfile(isProfile)
    rae2822.SetNumberOfProcessors(maxProcs)
    rae2822.SetNumberOfIterations(numIterationsShort)
    if rae2822.Processors() == 2:
        rae2822.SetResiduals([5.5892e-01, 6.7268e-01, 5.3250e-01, 1.0000e+00,
                              5.0058e-01, 2.5771e-09, 3.4059e-10])
    else:
        rae2822.SetResiduals([5.5618e-01, 6.6813e-01, 5.3620e-01, 1.0000e+00,
                              4.9726e-01, 2.5769e-09, 3.4032e-10])
    rae2822.SetIgnoreIndices(3)
    rae2822.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = rae2822.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # couette flow
    # laminar, viscous, periodic bcs, moving wall, isothermal wall
    couette = regressionTest()
    couette.SetRegressionCase("couette")
    couette.SetAitherPath(options.aitherPath)
    couette.SetRunDirectory("couette")
    couette.SetProfile(isProfile)
    couette.SetNumberOfProcessors(1)
    couette.SetNumberOfIterations(numIterations)
    couette.SetResiduals([1.1816e-01, 5.0725e-01, 6.9807e-02, 5.5916e-01,
                          2.3024e-01])
    couette.SetIgnoreIndices(3)
    couette.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = couette.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # wall law
    # wall law bc, turbulent, blusgs
    wallLaw = regressionTest()
    wallLaw.SetRegressionCase("wallLaw")
    wallLaw.SetAitherPath(options.aitherPath)
    wallLaw.SetRunDirectory("wallLaw")
    wallLaw.SetProfile(isProfile)
    wallLaw.SetNumberOfProcessors(maxProcs)
    wallLaw.SetNumberOfIterations(numIterationsShort)
    if wallLaw.Processors() == 2:
        wallLaw.SetResiduals([7.3745e-01, 1.5345e-01, 3.1677e-01, 9.2831e-01,
                              7.1928e-01, 2.6861e-02, 2.6255e-07])
    else:
        wallLaw.SetResiduals([7.4098e-01, 1.4914e-01, 3.1463e-01, 9.2837e-01,
                              7.2133e-01, 2.6860e-02, 2.6250e-07])
    wallLaw.SetIgnoreIndices(1)
    wallLaw.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = wallLaw.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # thermally perfect gas
    # turbulent, thermally perfect, supersonic
    thermallyPerfect = regressionTest()
    thermallyPerfect.SetRegressionCase("thermallyPerfect")
    thermallyPerfect.SetAitherPath(options.aitherPath)
    thermallyPerfect.SetRunDirectory("thermallyPerfect")
    thermallyPerfect.SetProfile(isProfile)
    thermallyPerfect.SetNumberOfProcessors(maxProcs)
    thermallyPerfect.SetNumberOfIterations(numIterationsShort)
    if thermallyPerfect.Processors() == 2:
        thermallyPerfect.SetResiduals([5.8177e-01, 3.8066e-01, 4.8670e-01,
                                       1.0000e+00, 5.9931e-01, 1.2830e-06,
                                       3.5031e-04])
    else:
        thermallyPerfect.SetResiduals([5.8177e-01, 3.8066e-01, 4.8670e-01,
                                       1.0000e+00, 5.9931e-01, 1.2830e-06,
                                       3.5031e-04])
    thermallyPerfect.SetIgnoreIndices(3)
    thermallyPerfect.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = thermallyPerfect.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # uniform flow
    # turbulent, all 8 block-to-block orientations
    uniform = regressionTest()
    uniform.SetRegressionCase("uniformFlow")
    uniform.SetAitherPath(options.aitherPath)
    uniform.SetRunDirectory("uniformFlow")
    uniform.SetProfile(isProfile)
    uniform.SetNumberOfProcessors(1)
    uniform.SetNumberOfIterations(numIterationsShort)
    uniform.SetResiduals([1.0342e+00, 9.0115e-01, 2.7211e-01, 8.8273e-01,
                          9.0131e-01, 1.4756e-07, 1.8748e-07])
    uniform.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = uniform.RunCase()
    # only care if this case ran, since nothing is changing, residuals have
    # some variation with compiler, os, optimiziation level, etc

    # ------------------------------------------------------------------
    # convecting vortex
    # initialization from file, nonreflecting boundary
    vortex = regressionTest()
    vortex.SetRegressionCase("convectingVortex")
    vortex.SetAitherPath(options.aitherPath)
    vortex.SetRunDirectory("convectingVortex")
    vortex.SetProfile(isProfile)
    vortex.SetNumberOfProcessors(1)
    vortex.SetNumberOfIterations(numIterations)
    vortex.SetResiduals([5.2772e+00, 6.3732e-01, 7.0928e-01, 1.0000e+00, 
                         7.9563e-01])
    vortex.SetIgnoreIndices(3)
    vortex.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = vortex.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # supersonic mixing with diffusion
    # turbulent, schmidt diffusion, 4th central, bdplur
    supersonicMixing = regressionTest()
    supersonicMixing.SetRegressionCase("supersonicMixing")
    supersonicMixing.SetAitherPath(options.aitherPath)
    supersonicMixing.SetRunDirectory("supersonicMixing")
    supersonicMixing.SetProfile(isProfile)
    supersonicMixing.SetNumberOfProcessors(maxProcs)
    supersonicMixing.SetNumberOfIterations(numIterationsShort)
    if supersonicMixing.Processors() == 2:
        supersonicMixing.SetResiduals([2.1642e-01, 1.5503e-01, 1.3670e+00,
                                       8.2043e-02, 3.3908e-01, 3.6563e-04,
                                       1.2388e-05])
    else:
        supersonicMixing.SetResiduals([2.1360e-01, 1.5278e-01, 1.3632e+00,
                                       7.8807e-02, 3.3470e-01, 3.6610e-04,
                                       1.2393e-05])
    supersonicMixing.SetIgnoreIndices(3)
    supersonicMixing.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = supersonicMixing.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # dissociation
    # reacting chemistry
    dissociation = regressionTest()
    dissociation.SetRegressionCase("dissociation")
    dissociation.SetAitherPath(options.aitherPath)
    dissociation.SetRunDirectory("dissociation")
    dissociation.SetProfile(isProfile)
    dissociation.SetNumberOfProcessors(1)
    dissociation.SetNumberOfIterations(numIterations)
    dissociation.SetResiduals([4.6027e-01, 4.6120e-01, 2.6470e+07, 1.0000e+00,
                               2.3423e-01])
    dissociation.SetIgnoreIndices(2)
    dissociation.SetIgnoreIndices(3)
    dissociation.SetMpirunPath(options.mpirunPath)

    # run regression case
    passed = dissociation.RunCase()
    totalPass = totalPass and all(passed)

    # ------------------------------------------------------------------
    # regression test overall pass/fail
    # ------------------------------------------------------------------
    errorCode = 0
    if totalPass and uniform.PassedStatus() != "ERRORS":
        print("All tests passed!")
    else:
        print("ERROR: Some tests failed")
        errorCode = 1
    print("--------------------------------------------------")
    print("subsonicCylinder:", subCyl.PassedStatus())
    print("multiblockCylinder:", multiCyl.PassedStatus())
    print("shockTube:", shockTube.PassedStatus())
    print("shockTubeRestart:", shockTubeRestart.PassedStatus())
    print("supersonicWedge:", supWedge.PassedStatus())
    print("transonicBump:", transBump.PassedStatus())
    print("viscousFlatPlate:", viscPlate.PassedStatus())
    print("turbulentFlatPlate:", turbPlate.PassedStatus())
    print("rae2822:", rae2822.PassedStatus())
    print("couette:", couette.PassedStatus())
    print("wallLaw:", wallLaw.PassedStatus())
    print("thermallyPerfect:", thermallyPerfect.PassedStatus())
    print("uniform:", uniform.PassedStatus())
    print("convectingVortex:", vortex.PassedStatus())
    print("supersonicMixing:", supersonicMixing.PassedStatus())
    print("dissociation:", dissociation.PassedStatus())
    sys.exit(errorCode)

if __name__ == "__main__":
    main()

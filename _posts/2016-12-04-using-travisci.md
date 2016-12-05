---
layout: post
title: "Using Travis CI For Regression Tests"
date: 2016-12-03 12:00
tags: [CFD, Aither, C++, v0.4.0, travis, continuous integration, travisci, regression, C++11, C++14, Cmake]
comments: true
---
## Why Use A Continuous Integration Service?
Continuous integration services allow code updates to be built and tested on a variety of platforms. This saves the developer
a lot of time by not having to manually test the code. As Aither grows larger it becomes more and more beneficial to use a
continuous integration service. For example, say a more efficient way to calculate the inviscid flux was found, and a new
branch was created to refactor the invisicid flux code to use this new method. This change should result in the same solution,
but should take less time to complete. To be thorough, before merging the code back into the **develop** branch, unit tests
covering all of the code's various functionality should be completed. These tests should still show that the solution is the
same as it was prior to the refactor. It can be tedious and time consuming to manually run these tests, not to mention
the tests should be run on different operating systems, and with different compilers as well. This is where continuous
integration saves the day! A continuous integration service will automatically build the most updated code on a variety of
operating systems with a variety of compilers, and can be made to run regression tests. This way it can easily be determined
if the refactor introduced any bugs.


## Aither's Requirements For Continuous Integration
Ok, so it is clear that continuous integration is a good thing, but which service should be used? Ideally, a continuous integration
service would provide the following:

* Free (Aither is not a money making venture after all)
* Testing on multiple operating systems (Aither is cross platform)
* Testing with multiple compilers
* Ability to use modern C++ (Aither uses C++14)
* Support for required dependencies (Aither requires an MPI implementation and Cmake)
* Ability to run regression tests in parallel
* Easy to use within [Github](https://github.com)

After a brief survery of available options, Aither recently started using [Travis CI](https://travis-ci.org) for continuous
integration. Travis CI meets all of the above requirements. It is free for open source projects, widely used in the Github 
community (i.e. [SU2](https://github.com/su2code/SU2)), and supports builds on Ubuntu and macOS. 


## Using Travis CI
Once an account has been created with Travis CI it is easy to integrate with Github. All that is required is to add a **.travis.yml** file
to the repository. This file instructs Travis CI on how to build the code and run any regression tests. For Aither, a matrix of five
builds is setup (Ubuntu/gcc-5, Ubuntu/gcc-6, Ubuntu/clang, macOS/gcc-6, macOS/clang). These builds are setup under the `matrix` data
field of the **.travis.yml** file. An abbreviated build matrix is shown below; each of the builds is marked by the `- os:` line. 


```yaml
# set up build matrix
matrix:
  include:
    # build for Ubuntu/gcc-6
    - os: linux
      dist: trusty
      sudo: required
      compiler: gcc
      # add toolchains for newer, C++14 supporting gcc-6
      addons:
        apt:
          sources:
            - ubuntu-toolchain-r-test
          packages:
            - g++-6 gcc-6 libstdc++-6-dev
      # change default compiler to newer gcc-6
      env:
        - CXX_COMPILER=g++-6
        - C_COMPILER=gcc-6
    # build for macOS/clang
    - os: osx
      osx_image: xcode8
      compiler: clang
      # change defualt and homebrew compilers to clang
      env:
        - CXX_COMPILER=clang++
        - C_COMPILER=clang
        - HOMEBREW_CC=clang
        - HOMEBREW_CXX=clang++
```

### Installing MPI
During the build matrix setup, environment variables for the C/C++ compilers are changed to reflect the newer C++14 supporting
compiler to be used in the build. Travis CI only offers Ubuntu 14.04 as its newest linux offering. Since this version of the
operating systems is a few years old, updated compilers are needed for the latest C++ standard. However, when the Ubuntu
package manager is used, it installs binaries that were created with the system C/C++ compilers. For compatability purposes,
it would be best if Aither used a version of MPI that was compiled with the same compiler that will be used to compile
Aither itself. For this reason OpenMPI is compiled from source using the updated compilers.

For macOS, things are a little different. The macOS virtual machines from Travis CI come preinstalled with the 
[homebrew](http://brew.sh) package manager. With homebrew the `HOMEBREW_CC` and `HOMEBREW_CXX` environment variables
control the compiler that new packages are built with. This means that installing MPI is easier because the package
manager can do it automatically.

This means that the **.travis.yml** script has to tell Travis CI to install MPI in a different way depending on which
operating system the build is happening on. This can easily be done with a simple bash script as shown below.

```bash
#!/bin/bash

# for macOS builds use OpenMPI from homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    cd openmpi
    # check to see if OpenMPI is cached from previous build
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
        echo "Installing OpenMPI with homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/openmpi
        brew install open-mpi
    fi
else
    # for Ubuntu builds install OpenMPI from source
    # check to see if OpenMPI is cached from previous build
    if [ -f "openmpi/bin/mpirun" ] && [ -f "openmpi-2.0.1/config.log" ]; then
	echo "Using cached OpenMPI"
	echo "Configuring OpenMPI"
	cd openmpi-2.0.1
	./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
    else
        # install OpenMPI from source
	echo "Downloading OpenMPI Source"
	wget https://www.open-mpi.org/software/ompi/v2.0/downloads/openmpi-2.0.1.tar.gz
	tar zxf openmpi-2.0.1.tar.gz
	echo "Configuring and building OpenMPI"
	cd openmpi-2.0.1
	./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=$C_COMPILER CXX=$CXX_COMPILER &> openmpi.configure
	make -j4 &> openmpi.make
	make install &> openmpi.install
	cd ..
    fi
    # recommended by Travis CI documentation to unset these for MPI builds
    test -n $CC && unset CC
    test -n $CXX && unset CXX
fi
```


### Caching Dependencies
Builds can be sped up on Travis CI by caching dependencies. Aither depends on MPI which can take a while to build from
source. However, this really only needs to be done once if the MPI installation can be cached and retrieved from build
to build. Fortunately, Travis CI allows this capability even for their free tier of services. To cache the MPI install
directory is simple. Only the following few lines need to be added to the **.travis.yml** file. This caches the OpenMPI
source code directory, as well as the installation directory.


```yaml
cache:
  directories:
    - openmpi
    - openmpi-2.0.1
```


## Regression Tests
Once the build completes Travis CI will run the Aither regression tests. The regression tests are located in the **testCases**
directory of the repository. Travis CI will run each case for 100 iterations and compare the residuals to some "truth" values.
If the residuals differ by less than a given ammount (1% for Aither), the test passes. The idea is that the regression tests
cover most or all of the code's functionality. On Ubuntu builds there are two processors available, so the tests are run in
parallel. On macOS there is only one processor available, so the tests are run in serial. The Aither repository includes a
python script to automate the running of these regression tests. After Travis CI builds the code, the script is invoked to
run the tests.

## Conclusion
Travis CI is now used by Aither to test builds on Ubuntu and macOS using gcc-5, gcc-6, and clang. Regression tests are run for
all of Aither's test cases to ensure that no existing functionality is broken with changes to the code. For more information on
how the whole thing is set up, visit the [repository](https://github.com/mnucci32/aither) and check out the **.travis.yml** and
**travis/installMPI** files.
---
layout: post
title: "Version 0.1.0 Released"
date: 2016-01-20 23:13
tags: [CFD, Aither, C++, v0.1.0, Release]
comments: true
---
# Release Notes
The first official release of Aither is now available. Check it out on [GitHub](https://github.com/mnucci32/aither/releases). This first release contains enough basic features to be able to run most standard CFD problems. Aither has the following functionality:

* **Time Integration**
  * Backwards Euler - First order, implicit
  * BDF2 - Second order, implicit
  * Forward Euler - First order, explicit
  * Low Storage RK4 - Second order, explicit
* **State Reconstruction**
  * MUSCL reconstruction for inviscid fluxes
  * Central reconstruction for viscous fluxes
* **Green-Gauss Gradient Formulation**
* **LU-SGS Implicit Solver**
* **Turbulence Models**
  * Wilcox K-Omega 2006
  * SST 2003
* **Input/Output**
  * Binary PLOT3D grid
  * Binary PLOT3D function files for output
* **K-D Tree Wall Distance Solver**

Have a look and see if you can make any improvements! For feature requests comment here, on Github, or on Twitter.


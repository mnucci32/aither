---
layout: post
title: "Introducing Aither"
date: 2015-11-24 10:24:00
tags: [CFD, Aither]
comments: true
---
# Introduction
Aither (**A**ccurate **I**mplicit **T**hree-dimensional **H**igh-order **E**ffcient **R**ANS), pronounced [*ee-thir*](https://en.wikipedia.org/wiki/Aether_(mythology)), is an open source CFD code written in C++11/14. Aither is a structured multi-block Reynolds-Averaged Navier-Stokes code. The Aither project is hosted on Github.

# Motivation
Aither began as side project to learn more about the inner workings of a CFD code. It is very thoroughly commented in hopes that this would allow someone else to learn from the code and/or contribute to the project. 

# Features
- **Parallel**: Domain decomposition via MPI
- **Implicit**: LU-SGS method used for implicit solver; Euler and Bdf2 time integration methods available
- **Explicit**: Explicit Euler and 4-stage low-storage Runge-Kutta methods available
- **Inviscid Flux**: Roe approximate riemann solver used for inviscid flux
- **Viscous Flux**: Central difference method used for viscous flux
- **Spatial Accuracy**: MUSCL extrapolation used for 2nd order accuracy in space
- **Temporal Accuracy**: Bdf2 (implicit) and RK4 (explicit) methods for 2nd order accuracy in time
- **Turbulence Modeling**: Wilcox K-Omega 2006 and Menter's K-Omega SST 2003 models
- **Equation of State**: Ideal gas equation
- **Input**: Input file is a text file and input grid is PLOT3D format
- **Output**: Output is a PLOT3D function file which is readable from Paraview

# Direction Forward
Aither will be continuously improved. The next couple of items on the to-do list are to implement an implicit solver with more accurate jacobians to facilitate faster convergence. A block LU-SGS method will most likely be used. The SST turbulence model will be extended to work with the DES and IDDES variants. A higher order face-state reconstruction will also be implemented to allow the code to accurately run Hybrid RANS/LES and LES simulations.
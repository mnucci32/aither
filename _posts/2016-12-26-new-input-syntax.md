---
layout: post
title: "New Input File Syntax: Vectors, States, & Lists"
date: 2016-12-26 02:00
tags: [CFD, Aither, C++, develop, vectors, states, lists, input, syntax]
comments: true
---
## New Input File Syntax
There is a new input file syntax for Aither now in use in the **develop** branch of the code. This
syntax makes it easy to specify initial conditions by grid block, and boundary conditions by
boundary condtion tag. This is a huge upgrade in usability as it now allows for problems such as
Sod's shock tube to be simulated. It also allows for easy implementation of various **viscousWall**
boundary conditions such as *adiabatic*, *isothermal*, and *constant heat flux*. The new input file
syntax is based off of three new objects (**vectors**, **lists**, & **states**) which will be 
discussed in detail below.


### Vectors
Vectors are now input in a comma separated list enclosed in brackets like below. Vectors must be
defined entirely on one line in the input file.

```
velocityRef: [1.0, 0.0, 0.0]
```

Valid vector inputs have three components. If three components are not specified, Aither will throw
and error. Vector inputs are now used wherever vector quantities are needed such as for velocity
(above), or specifiying a direction as is done with the **stagnationInlet** boundary condition.

### States
States are a group of properties that apply to an initial condition state or a boundary condtion
state. States are identified by name, enclosed in parenthesis, and individual properties within a
state are assigned with the equals operator and separated by semicolons. States must be defined
entirely on one line in the input file. Below is an example of *icState* which is used to specify
a flow state for initial condtions.

```
icState(tag=0; pressure=101325; density=1.225; velocity=[100, 0, 0])
```

The supported properties depend on the type of state. There are optional turbulence properities
*tubulenceIntensity* and *eddyViscosityRatio* that may be specified for states used for inflow
boundary conditions or *icState*. If either of the optional turbulence properties are specified,
both must be specified. For *icState* the *tag* property is special. It refers to the block
number in which the *icState* will be applied. A value of -1 functions as the default state in
the event that there is not an *icState* with a tag pointing to a given block. An explicity
specified tag takes precedence over the default state. For example for a four block grid with
two *icState*s defined, one with a tag of -1, and another with a tag of 0, blocks 1-3 will use 
the default *icState* with tag -1, and block 0 will use the *icState* with the tag of 0.

In addition to initial conditions, states are used for boundary conditions that may require
additional information. An example of each such boundary condition is shown below. For boundary
conditions, the tag property in each state refers to the boundary surface tag that is specified
in the boundary condition definition. 

Inflow boundary conditions. These may optionally specify the turbulence properties.

```
characteristic(tag=0; pressure=101325; density=1.225; velocity=[100, 0, 0])
```

```
stagnationInlet(tag=0; p0=101325; t0=300; direction=[1, 0, 0])
```

```
supersonicInflow(tag=0; pressure=101325; density=1.225; velocity=[100, 0, 0])
```

```
subsonicInflow(tag=0; density=1.225; velocity=[100, 0, 0])
```

Outflow boundary conditions.

```
pressureOutlet(tag=0; pressure=101325)
```

```
subsonicOutflow(tag=0; pressure=101325)
```

Wall boundary condtions. One of *heatFlux* or *temperature* may be specified. The default behavior is zero
velocity and zero heat flux which corresponds to a stationary adiabatic wall.

```
viscousWall(tag=0; heatFlux=100)
```

```
viscousWall(tag=0; temperature=400)
```

```
viscousWall(tag=0; velocity=[10, 0, 0])
```

### Lists
Lists are a comma separated group of properties that are enclosed in angle brackets. Lists may be specified
across multiple lines. Lists are most commonly used to specify the variables to output, the initial condition
states, and the boundary condition states. Examples are shown below.

```
outputVariables: <density, vel_x, vel_y, vel_z, pressure, temperature, mach>
```

```
initialConditions: <icState(tag=0; pressure=101325; density=1.225; velocity=[0, 0, 0]),
                    icState(tag=1; pressure=10132.5; density=0.153125; velocity=[0, 0, 0])>
```

```
boundaryStates: <characteristic(tag=0; pressure=101325; density=1.225; velocity=[100, 0, 0]),
                 viscousWall(tag=1; velocity=[10, 0, 0])>
```

## Summary
The new input syntax in Aither is more intuitive and now allows for a wider variety of problems to
easily be simulated. All of the test cases in the **develop** branch have been updated to support this
new syntax. Grab the **develop** branch from Github and try it out today. This will be merging into the
**master** branch shortly.


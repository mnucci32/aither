---
layout: post
title: "Couette Flow & New Boundary Conditions"
date: 2017-03-04 17:00
tags: [CFD, Aither, C++, couette, periodic, isothermal]
comments: true
---
## Couette Flow
[Couette flow](https://en.wikipedia.org/wiki/Couette_flow) is viscous laminar flow
between two parallel plates, one of which is moving relative to the other. Due to its
simple nature and the existance of an analytical solution, it is a common validation
case for CFD codes. An example validation case is shown in Hirsch [1]. Couette flow
results in a constant shear stress which has a linear velocity profile and a parabolic
temperature profile as shown below.

$$ v(y) = \frac{y}{L} v_{wall} $$

$$ P_r E_c = \frac{\mu v^2_{wall}}{k \Delta T} $$

$$ T(y) = T_{low} + \Delta T \frac{y}{L} \left[1 + \frac{1}{2} P_r E_c \left(1 - \frac{y}{L} \right) \right] $$


Three of the newer features in Aither are periodic boundary conditions, moving walls,
and isothermal walls. A couette flow simulation can make use of all of these features,
so it makes a great addition to the test cases suite.

## Problem Setup
The parallel plates are placed a distance of 0.001 meters apart. The bottom plate is
held at a temperature of 288 K and is stationary. The top plate is held at a temperature
of 289 K and is moving at 75.4 m/s. For this setup, the product of the Prandtl and Eckert
numbers is 4. This means that for the temperature profile, the maximum temperature will
not be at the plate, but in the flow instead. The exact solution dictates that the
maximum temperature should be three fourths of the way between the cold and hot plates.

The CFD domain is a rectangular prism with the top and bottom modeled as viscous walls,
the sides as slip walls, and the front / back as periodic. Isothermal walls can be
specified in Aither by adding the **temperature** parameter to the boundary state list.
Similarly moving walls can be specified by adding the **velocity** parameter. Periodic
boundary conditions are specified by indicating which boundary condition tags should be
paired as periodic. This is done through the **startTag** and **endTag** parameters. For
each periodic boundary condition a transformation must be done to get from one periodic
face to the other. Currently a translation can be specified by adding the **translation**
parameter which is a vector specifying how the boundary at the **startTag** should be
translated to get to the boundary at the **endTag**. Alternatively a rotation can be
specified by using the **axis**, **point**, and **rotation** parameters. The **axis**
parameter is a vector defining the axis of rotation. The **point** parameter is a vector
defining a point about which to rotate. The **rotation** parameter is a scalar defining
the rotation angle in radians. An example of these new boundary condition options is
shown below.

```
boundaryStates: <periodic(startTag=4; endTag=5; translation=[0.01, 0, 0]),
                 viscousWall(tag=1; temperature=288),
		 viscousWall(tag=2; temperature=289; velocity=[75.4, 0, 0])>
```

## Results and Summary
The results from Aither show a linear velocity profile and a parabolic temperature
profile as expected. The results agree very well with the exact solution. These results
can be reproduced by running the [Aither](https://github.com/mnucci32/aither) code. The
grid and input file for the couette flow case can be found in the **testCases** directory
of the repository.

![Couette]({{site.baseurl}}/downloads/couette.png){: .center-image}

<center>Velocity and temperature profiles for Couette flow.</center>

## References
[[1]](https://www.amazon.com/Numerical-Computation-Internal-External-Flows/dp/0750665947)
Hirsch, Charles. "Numerical Computation of Internal and External Flows". 2nd Edition.
Butterworth-Heinemann. 2006.



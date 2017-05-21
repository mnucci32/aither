---
layout: post
title: "Thermally Perfect Thermodynamic Model"
date: 2017-05-20 17:00
tags: [CFD, Aither, C++, thermally perfect, calorically perfect, ideal gas, thermodynamics]
comments: true
---
## Calorically Perfect Gas
The default thermodynamic model in Aither is the calorically perfect gas model. 
Calorically perfect gases have constant specific heats 
$$\left( c_p, c_v \right)$$, and therefore a constant $$\gamma$$. In general 
the calorically perfect gas model is a good assumption for air at lower 
temperatures. The molecules of a calorically perfect gas are assumed to be rigid
so there is no vibrational component of internal energy. The equations below
show the implementation of the calorically perfect gas thermodynamic model in
Aither.


$$ e = e_{translational} + e_{rotational} + e_{vibrational} = c_v T$$

$$ e_{translational} =  \left( n - 1 \right)  R T 
\,\,\,\,\,\,\,
 e_{rotational} = R T 
\,\,\,\,\,\,\,
 e_{vibrational} = 0 $$

$$ h = c_p T $$

$$ c_v = n R
\,\,\,\,\,\,\,
c_p = \left( n + 1 \right) R
\,\,\,\,\,\,\,
\gamma = \frac{c_p}{c_v} = \frac{1}{n} + 1 $$


## Thermally Perfect Gas
As the temperature of a gas increases the vibrational modes of its molecules 
are activated [1,2]. This means that some of the energy of the gas will move 
into these vibrational modes instead of raising the temperature of the gas. 
Therefore at higher temperatures where the vibrational modes are activated 
a thermally perfect gas will have a cooler temperature than a calorically 
perfect gas. The temperature at which the vibrational modes become significant 
is determined by the vibrational temperature $$T_v$$ of the gas. As the 
temperature of the gas approaches this value the vibrational component of 
internal energy becomes more and more significant.

The activiation of vibrational modes also means the the specific heats are no
longer constant. Instead they are assumed to be functions of temperature only.
The equations below show the implementation of the thermally perfect 
thermodynamic model in Aither.

$$ e = e_{translational} + e_{rotational} + e_{vibrational} = 
\int_0^T c_v \left( t \right) dt 
\,\,\,\,\,\,\,
e_{vibrational} = \frac{R T_v}{e^{\frac{T_v}{T}} - 1} $$

$$ h = \int_0^T c_p \left( t \right) dt $$

$$ c_v = \frac{\partial e}{\partial t} = 
R \left( n + \left[ \frac{\theta_v}{sinh \left( \theta_v \right)} \right] 
^ 2 \right) 
\,\,\,\,\,\,\,
c_p = \frac{\partial h}{\partial t} = 
R \left( n + 1 + \left[ \frac{\theta_v}{sinh \left( \theta_v \right)} \right]
 ^ 2 \right) 
\,\,\,\,\,\,\,
\theta_v = \frac{T_v}{2 T}$$

$$ \gamma = \frac{c_p \left( t \right)}{c_v \left( t \right)} $$

In Aither the thermally perfect thermodynamic model can be activated as shown
below.

```
fluids: <fluid(name=air; n=2.5; molarMass=0.02897; vibrationalTemperature=3056)>
thermodynamicModel: thermallyPerfect
```

## Example Problem
An example problem of when the thermally perfect model is needed is now a part
of the test cases that come with the Aither repository. It is currently only 
available on the **develop** branch, but will be available on **master** after 
the next release. The test case involves hot supersonic flow over a 20 degree 
ramp. The freestream conditions of the flow are Mach 3, static temperature of 
2000 K, and static pressure of 229,600 Pa.

![ThermallyPerfect]({{site.baseurl}}/downloads/cpg_tpg.png){: .center-image}

<center>Comparison of calorically perfect and thermally perfect thermodynamic
models.</center>

The results show that behind the shocks the thermally perfect gas model predicts
a cooler flow than the calorically perfect gas model. This is expected as with
the thermally perfect model some of the energy goes into the vibrational modes
of the gas molecules.

## References
[[1]](http://www.donnerflug.de/thesis/Lampe_MS_Thesis.pdf)
Lampe, Dietrich Rudolf. "Thermally Perfect, Calorically Imperfect 
Taylor-Maccoll Flow". 1994.

[[2]](https://www.amazon.com/Hypersonic-High-Temperature-Dynamics-Second-
Education/dp/1563477807/ref=asap_bc?ie=UTF8) 
Anderson, John. "Hypersonic and High Temperature Gas Dynamics". 2nd Edition. 
AIAA. 2006.


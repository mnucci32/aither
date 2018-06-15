---
layout: post
title: "Reacting Flow"
date: 2018-05-21 08:00
tags: [CFD, Aither, C++, finite-rate, chemsitry, reactions]
comments: true
---
## Finite Rate Chemistry in Aither
As of [v0.9.0](https://github.com/mnucci32/aither/releases/tag/v0.9.0) Aither 
has a finite rate reacting flow capability. This post will detail how the 
chemistry source terms are calculated in Aither. The governing equations for 
multispecies reacting flow are given below.

$$\frac{\partial}{\partial t}\int_V W\,dV\ +
\int_{A} \left[ \frac{\partial F_i(W)}{\partial x_j} - \frac{\partial F_v(W)}{\partial x_j} \right]\,dA\
= \int_V S\,dV\$$

In order to calculate the source term $$ S $$ due to chemical reactions, some
information about each species in the simulation is needed. In Aither the
fluid database specifies the enthalpy of formation $$h_f$$, molar mass $$M$$,
vibrational temperatures $$T_v$$, and the linear component of internal energy
$$n$$. The fluid database also supplies reference properties for pressure
$$p_{ref}$$, temperature $$T_{ref}$$, and entropy $$s_{ref}$$. It is possible
to add to Aither's fluid database by specifying these properties for a new 
species. This information can be found in a variety of sources, including the
[NIST Chemistry WebBook](https://webbook.nist.gov).

## Chemistry Source Terms
Reacting flow simulations differ from nonreacting simulations in that the 
species mass equations have a source term, $$ S_s $$. This source term is
calculated via the law of mass action. In the equation below 
$$ \nu_{s,r}^{''} $$ represents the stoichiometric coefficient on the product 
side of the reaction for species $$ s $$ in reaction $$ r $$. Similarly, 
$$ \nu_{s,r}^{'} $$ represents the stoichiometric coefficient on the reactant
side of the reaction for species $$ s $$ in reaction $$ r $$.

$$ S_s = M_s \sum_{r=1}^{NR} \left( \nu_{s,r}^{''} - \nu_{s,r}^{'} \right)
   \left[ k_{f,r} \prod_{i=1}^{NS} \left(\frac{\rho_i}{M_i} \right)^{\nu_{i,r}^{'}} - 
   k_{b,r} \prod_{i=1}^{NS} \left(\frac{\rho_i}{M_i} \right)^{\nu_{i,r}^{''}} \right]
$$

The source term depends on the forward and backward reaction rates. These rates
are functions of temparature only. In Aither the forward reaction rate 
$$ k_f $$ is calculated using an Arrhenius curve fit. The backward reaction 
rate $$ k_b $$ is calculated from thermodynamic properties using the 
equilibrium rate $$ k_e $$. The units of the forward rate are 
$$ \frac{1}{s} \left( \frac{mol}{m^3} \right)^{1 - \sum_{s=1}^{NS} \nu_s^{'}} $$
and the units of the backward rate are
$$ \frac{1}{s} \left( \frac{mol}{m^3} \right)^{1 - \sum_{s=1}^{NS} \nu_s^{''}} $$.

$$ k_{f,r} = C T^{\eta} e^{\frac{T_a}{T}} \qquad
 k_{b,r} = \frac{k_{f,r}}{k_{e,r}} $$

## Calculation of Equilibrium Rate
The equilibrium rate of the reaction is calculated via the minimization of the
Gibbs free energy [1]. For a thermally perfect gas using the vibrational 
equilibrium model in Aither, the Gibbs free energy for a given species is shown 
below. In the equations below $$NV$$ represents the number of vibrational modes
for a given species.

$$ G_s = H_s - T S_s = h_{f,s} + \int_0^T c_{p,s} \left( t \right) dt - 
T \int_0^T \frac{c_{p,s} \left( t \right)}{t} dt - T s_{0,s} $$

$$ G_s = h_{f,s} + \int_0^T R_s \left( n + 1 + 
\sum_{v=1}^{NV} \left[ \frac{T_{v,s}}{2 t \cdot sinh \left( \frac{T_{v,s}}{2t} \right)} \right] ^ 2 \right) dt - 
T \int_0^T \frac{R_s}{t} \left( n + 1 + 
\sum_{v=1}^{NV} \left[ \frac{T_{v,s}}{2 t \cdot sinh \left( \frac{T_{v,s}}{2t} \right)} \right] ^ 2 \right) dt - T s_{0,s} $$

$$ G_s = h_{f,s} + R_s T \left( n + 1 \right) + R_s
\sum_{v=1}^{NV} \frac{T_{v,s}}{e^{\frac{T_{v,s}}{T}} - 1} - 
R_s T \left( n + 1 \right) ln \left( T \right) - T R_s
\int_0^T \sum_{v=1}^{NV} \left[ \frac{T_{v,s}}{2 t \cdot sinh \left( \frac{T_{v,s}}{2t} \right)} \right] ^ 2 dt - T s_{0,s}
$$

The remaining integral represents the vibrational contribution to the entropy
term. It can be solved with integration by parts.

$$ \int_0^T \left[ \frac{T_v}{2 t \cdot sinh \left( \frac{T_v}{2t} \right)} \right] ^ 2 dt =
\frac{1}{T} \frac{T_v}{e^{\frac{T_v}{T}} - 1} -
\int_0^T \frac{T_v}{e^{\frac{T_v}{t}} - 1} \frac{-1}{t^2} dt =
\frac{T_v}{T \cdot \left( e^{\frac{T_v}{T}} - 1 \right)} - 
ln \left( e^{\frac{T_v}{T}} - 1 \right) +
\frac{T_v}{T}$$

$$ \frac{T_v}{T \cdot \left( e^{\frac{T_v}{T}} - 1 \right)} - 
ln \left( e^{\frac{T_v}{T}} - 1 \right) +
\frac{T_v}{T} =
\frac{T_v}{T \cdot \left( e^{\frac{T_v}{T}} - 1 \right)} - 
ln \left( 1 - e^{\frac{-T_v}{T}} \right) $$

Substituting this expression back into Gibbs free energy equation yields the 
following.

$$ G_s = h_{f,s} + R_s T \left( n + 1 \right) + R_s
\sum_{v=1}^{NV} \frac{T_{v,s}}{e^{\frac{T_{v,s}}{T}} - 1} - 
R_s T \left( n + 1 \right) ln \left( T \right) - R_s
\sum_{v=1}^{NV} \left[ \frac{T_{v,s}}{e^{\frac{T_{v,s}}{T}} - 1} - T
ln \left( 1 - e^{\frac{-T_{v,s}}{T}} \right) \right] - T s_{0,s} $$

$$ G_s = h_{f,s} + R_s T \left( n + 1 \right) \left[ 1 - ln \left( T \right) \right] -
R_s T \sum_{v=1}^{NV} ln \left( 1 - e^{\frac{-T_{v,s}}{T}} \right) - T s_{0,s} $$

Where $$ s_{0,s} $$ is a constant to account for any difference between the
given reference entropy of the species and the calculated entropy at the 
species' reference temperature.

$$ s_{0,s} = s_{ref,s} - R_s \left( n + 1 \right) ln \left( T_{ref} \right) -
\sum_{v=1}^{NV} \frac{T_{v,s}}{T_{ref} \cdot \left( e^{\frac{T_{v,s}}{T_{ref}}} - 1 \right)} - ln \left( 1 - e^{\frac{-T_{v,s}}{T_{ref}}} \right) $$

From this the equilibrium rate can be calculated as shown below.

$$ k_{e,r} = \left( \frac{p_{ref}}{\tilde{R} T} \right) 
^{\sum_{s=1}^{NS} \nu_{s,r}^{''} - \nu_{s,r}^{'}}
e^{-\sum_{s=1}^{NS} \frac{G_s}{R_s T} \left( \nu_{s,r}^{''} - \nu_{s,r}^{'} \right)}
$$

## Implicit Treatment
For the implicit solvers in Aither the chemistry jacobian and/or spectral 
radius of the jacobian are needed. The analytical chemistry jacobian is very 
complicated and somewhat expensive to calculate. For this reason Aither uses
a numerical jacobian. The derivative of the chemistry source terms is needed 
with respect to $$ \rho_s $$ and $$ \rho E $$. This is done by perturbing the 
conserved variable of interest by a small amount $$ \epsilon $$, while keeping
the other conserved variables constant, and recalucuating the chemistry source 
terms.

$$ \frac{\partial S}{\partial \rho_s} = 
\frac{S\left( \rho_s + \epsilon \right) - S \left( \rho_s \right)}{\epsilon}
$$

For the scalar matrix implicit solvers in Aither (LU-SGS and DPLUR) only the 
spectral radius of the chemistry jacobian is needed. This is approximated using
the procedure in [2] where only the negative terms that would increase the
diagonal dominance of the implicit matrix are treated implicitly.

$$ \lambda_{chem} = min \left( \frac{S_s^{-}}{Y_s} \right) $$

$$ S_s^{-} = M_s \sum_{r=1}^{NR} \left( \nu_{s,r}^{''} - \nu_{s,r}^{'} \right)
\left[ - k_{b,r} \prod_{i=1}^{NS} \left(\frac{\rho_i}{M_i} \right)^{\nu_{i,r}^{''}} \right]
$$


## References
[1] Luke, E. A. "A Rule-Based Specification System for Computational Fluid
Dynamics". Ph. D. Thesis. 1999.

[2] Savard, B. et al. "A Computationally-Efficient Semi-Implicit Iterative
Method for the Time Integration of Reacting Flows with Stiff Chemistry". 
Journal of Computational Physics. 2015.



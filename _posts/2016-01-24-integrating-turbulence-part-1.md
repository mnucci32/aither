---
layout: post
title: "Integrating Turbulence Models Into Aither - Part 1"
date: 2016-01-24 18:48
tags: [CFD, Aither, C++, RANS, SST, k-omega, turbulence, Menter, Roe, Wilcox, explicit, implicit]
comments: true
---
## Part 1 - Explicit Formulation
This post will focus on integrating turbulence models into Aither's explicit solvers. This includes calculation of the
residual which is necessary for the implicit solvers. The next post in this series will cover integrating turbulence
models into the LU-SGS implicit solver.

## Reynolds-Averaged Navier-Stokes Equations
To solve the Navier-Stokes equations with turbulence models we use the Favre and Reynolds-averaged equations as shown below.
The mean flow equations are identical to the Navier-Stokes equations with the exception of the viscosity and thermal
conductivity terms. The viscosity $$\mu$$ is replaced by its sum with the turbulent eddy viscosity $$\mu_t$$. The thermal
conductivity term $$\frac{\mu}{Pr(\gamma - 1)}$$ is replaced by its sum with the turbulent thermal conductivity
$$\frac{\mu_t}{Pr_t(\gamma - 1)}$$. The turbulence equations themselves only couple to the mean flow equations through the
turbulent eddy viscosity.

$$\frac{\partial}{\partial t}\int_V W\,dV\ +
\int_{A} \left[ \frac{\partial F_i(W)}{\partial x_j} - \frac{\partial F_v(W)}{\partial x_j} \right]\,dA\
= \int_V S(W)\,dV\$$

$$W = \left[ \begin{array}{c}
\rho \\
\rho v_x \\
\rho v_y \\
\rho v_z \\
\rho E \\
\rho k \\
\rho \omega
\end{array} \right] \;

F_i = \left[ \begin{array}{c}
\rho v_i n_i \\
\rho v_i n_i v_x + p n_x \\
\rho v_i n_i v_y + p n_y \\
\rho v_i n_i v_z + p n_z \\
\rho v_i n_i H \\
\rho v_i n_i k \\
\rho v_i n_i \omega
\end{array} \right] \;

F_v = \left[ \begin{array}{c}
0 \\
\tau_{xx} n_x + \tau_{xy} n_y + \tau_{xz} n_z \\
\tau_{yx} n_x + \tau_{yy} n_y + \tau_{yz} n_z \\
\tau_{zx} n_x + \tau_{zy} n_y + \tau_{zz} n_z \\
\tau_{ij} n_i v_j - (\frac{\mu}{Pr(\gamma - 1)} + \frac{\mu_t}{Pr_t(\gamma - 1)}) n_i \frac{\partial T}{\partial x_i}\\
(\mu + \sigma_k \mu_t) \frac{\partial k}{\partial x_i} n_i\\
(\mu + \sigma_{\omega} \mu_t) \frac{\partial \omega}{\partial x_i} n_i
\end{array} \right] \;

S = \left[ \begin{array}{c}
0 \\
0 \\
0 \\
0 \\
0 \\
P_k - D_k \\
P_{\omega} - D_{\omega} + CD_{k\omega}
\end{array} \right]$$

$$\tau_{ij} = \lambda \frac{\partial v_i}{\partial x_i} \delta_{jk} +
(\mu + \mu_t) \left[ \frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i} \right]$$

In the above equations the turbulence model shown is a two equation k-$$\omega$$ model. The turbulent source terms consist
of production, destruction, and cross diffusion terms. These terms are model dependent and their values for the SST model
can be seen [here](http://turbmodels.larc.nasa.gov/sst.html).

## Numerical Solution Strategy
We solve the turbulence equations separately from the mean flow equations instead of solving both sets of equations
simultaneously. This is done because the coupling between the equation sets is relatively weak, and solving the equations
simultaneously requires more work. To solve the equations simultaneously using an implicit method requires a new flux
jacobian to be calculated for each turbulence model in the code. Solving the equations separately allows the same flux
jacobian to be used for all turbulence models.

To fully implement the turbulence models into the Aither code the inviscid and viscous flux calculations must be extended
for the turbulence equations. Also the code must now calculate the source terms of the turbulence models. The source terms
can provide a lot of issues numerically because they are characteristically *stiff*. These terms can severely limit the stable
time increment for a given time integration scheme. For this reason it is usually only practical to solve these equations
with an implicit method. Aither uses the Lower-Upper Symmetric Gauss Seidel (LU-SGS) implicit solver, so this must be extended
for the turbulence equations (covered in next post).

## Flux Calculations With Turbulence Equations
Aither uses the Roe flux difference splitting method to calculate the inviscid fluxes, and a central difference method to
calculate the viscous fluxes. These flux calculations must be extended for use with the turbulence equations.

# Roe Flux With Turbulence Equations
To calculate the inviscid fluxes the primative variables are reconstructed at the cell faces. A given face will have two
adjacent cells, so there will be two separated reconstructed states. These states may not be equal to each other and therefore
form a Riemann problem. Roe's approximate Riemann solver is used to determine the inviscid flux at the cell face. The flux is
calculated using the convective fluxes from the left and right reconstructed states as well as a dissipation matrix as shown
below. In the equations below variables marked with a ~ indicated Roe averaged quantities, and the $$\Delta$$ refers to the
right state minus the left state. The dissipation matrix ($$D$$) is calculated from the left eigenvalues of the Roe matrix ($$T$$),
the eigenvalues of the Roe matrix ($$\Lambda$$), and the wave strengths ($$\Delta C$$).

$$F_i = \frac{1}{2} (F_{c_l} + F_{c_r} - D)$$

$$D = \frac{\partial \tilde{F}}{\partial W} (W_r - W_l) = \tilde{A} \Delta W = T \Lambda \Delta C$$

$$T = \left[ \begin{array}{cccccc}
1 & 1 & 1 & 0 & 0 & 0 \\
\tilde{v}_x - \tilde{a} n_x & \tilde{v}_x & \tilde{v}_x + \tilde{a} n_x & \Delta v_x - \Delta v_i n_i n_x & 0 & 0  \\
\tilde{v}_y - \tilde{a} n_y & \tilde{v}_y & \tilde{v}_y + \tilde{a} n_y & \Delta v_y - \Delta v_i n_i n_y & 0 & 0  \\
\tilde{v}_z - \tilde{a} n_z & \tilde{v}_z & \tilde{v}_z + \tilde{a} n_z & \Delta v_z - \Delta v_i n_i n_z & 0 & 0  \\
\tilde{h} - \tilde{a} \tilde{v}_i n_i & \frac{1}{2} \tilde{v}_i \tilde{v}_i & \tilde{h} + \tilde{a} \tilde{v}_i n_i &
\tilde{v}_i \Delta v_i - \tilde{v}_i n_i \Delta v_j n_j & 0 & 0  \\
\tilde{k} & 0 & \tilde{k} & 0 & 1 & 0 \\
 \tilde{\omega} & 0 & \tilde{\omega} & 0 & 0 & 1
\end{array} \right]

\Lambda = \left[ \begin{array}{c}
\tilde{v}_i n_i - \tilde{a} \\
\tilde{v}_i n_i \\
\tilde{v}_i n_i + \tilde{a} \\
\tilde{v}_i n_i \\
\tilde{v}_i n_i \\
\tilde{v}_i n_i 
\end{array} \right]

\Delta C = \left[ \begin{array}{c}
\frac{\Delta P - \tilde{\rho} \tilde{a} \Delta v_i n_i}{2 \tilde{a}^2} \\
\Delta \rho - \frac{\Delta P}{\tilde{a}^2} \\
\frac{\Delta P + \tilde{\rho} \tilde{a} \Delta v_i n_i}{2 \tilde{a}^2} \\
\tilde{\rho} \\
\tilde{\rho} \Delta k + \tilde{k} \Delta \rho - \frac{\Delta P \tilde{k}}{\tilde{a}^2} \\
\tilde{\rho} \Delta \omega + \tilde{\omega} \Delta \rho - \frac{\Delta \tilde{\omega}}{\tilde{a}^2} 
\end{array} \right]$$

As can be seen from the above equations there is no coupling from the turbulence equations to the mean flow equations in
the inviscid flux calculation. However, the mean flow equations have some coupling to the turbulence equations. The above
formulation is written in a way that is independent of the face tangent vectors. More information on this formulation and
its derivation can be found [here](http://ossanworld.com/cfdbooks//cfdcodes/threed_euler_fluxes_v3.f90).

# Viscous Flux With Turbulence Equations
The viscous fluxes are calculated in the same way whether turbulence equations are present or not. A central difference is
used to reconstruct the state at the cell face. The viscous flux itself is calculated as in $$F_v$$ above. In order to calculate
the viscous flux for the turbulence equations gradients of $$k$$ and $$\omega$$ are needed. Therefore the existing gradient
calculation methods in the code are extended to calculate these additional gradients. All gradients are calculated using the
Green-Gauss method.

## Residual Calculation
The residual calculation with turbulence equations is identical with to that without it with the exception of the source terms.
There are no source terms in the mean flow equations, but they are present in the turbulence equations. The source terms are
multiplied by the cell volume, not the cell face area. It is important to note that the source terms start out on the right hand
side of the equation, opposite of the inviscid and viscous fluxes. The equation below shows how the source terms contribute to
the residual calculation.

$$W^{n+1} = W^n - \frac{\Delta t}{V} R = W^n - \frac{\Delta t}{V} \left[ (F_i - F_v) A  - S V \right]$$

With the changes above and the addition of appropriate boundary conditions, the code has been extended to solve the turbulence
equations in an explicit manner.



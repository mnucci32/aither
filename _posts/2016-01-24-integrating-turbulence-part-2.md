---
layout: post
title: "Integrating Turbulence Models Into Aither - Part 2"
date: 2016-01-24 19:00
tags: [CFD, Aither, C++, RANS, SST, k-omega, turbulence, Menter, Roe, Wilcox, explicit, implicit]
comments: true
---
## Part 2 - Implicit Formulation
This post will focus on integrating turbulence models into Aither's implicit solver. This includes calculation of the
flux jacobians which are necessary for the implicit solver. A previous post covered integration with the explicit solver,
which is a prerequisite to this post.

## LU-SGS Method With Turbulence Equations
An excellent description of the LU-SGS method for the mean flow equations can be found in
[Blazek](http://www.amazon.com/Computational-Fluid-Dynamics-Principles-Applications/dp/0080445063). 
The mean flow procedure remains unchanged because the turbulence equations are solved separately. The procedure for solving the
turbulence equations with the LU-SGS is identical to that of the mean flow equations with the exception that the flux jacobians
are different. As a review, the solution update is obtained with two sweeps, one forward and one backward over hyperplanes through
the domain. Sweeping across hyperplanes allows the off-diagonal terms to be computed on-the-fly. The two sweeps along with the
diagonal and off-diagonal matrices are shown below.

$$D \Delta W^{n^*} = - R^n - L \Delta W^{n^*}$$

$$D \Delta W^n = D \Delta W^{n^*} - U \Delta W^n$$

$$D = \frac{V}{\Delta t} + \omega \left( \frac{F_{i_i}}{\partial W} A_i + \frac{F_{i_j}}{\partial W} A_j + \frac{F_{i_k}}{\partial W} A_k \right)+
2 \left( \frac{\partial F_{v_i}}{\partial W} A_i+ \frac{\partial F_{v_j}}{\partial W} A_j+ \frac{\partial F_{v_k}}{\partial W} A_k \right)
- \frac{\partial S}{\partial W} V$$

$$L \Delta W^n = \frac{A}{2} \left( \Delta F_i + \omega \frac{\partial F_i}{\partial W} \Delta W^n + 2 \frac{\partial F_v}{\partial W} \right)$$

$$U \Delta W^n = \frac{A}{2} \left( \Delta F_i - \omega \frac{\partial F_i}{\partial W} \Delta W^n - 2 \frac{\partial F_v}{\partial W} \right)$$

The equations above show that the source terms only contribute to the main diagonal. The LU-SGS method requires inversion of the main diagonal
matrix. However a great computational savings can be incurred if the flux jacobians on the main diagonal are approximated by their spectral
radii. Using this approximation makes the elements of the main diagonal scalars instead of matrices, so their inversion is trivial. These flux
jacobian approximations are shown below. 

$$\frac{\partial F_i}{\partial W} = \left[ \begin{array}{cc}
0.5 \left(\left| v_i n_i \right| + v_i n_i \right) & 0 \\
0 & 0.5 \left(\left| v_i n_i \right| + v_i n_i \right)
\end{array} \right] \approx v_i n_i $$

$$\frac{\partial F_v}{\partial W} = \left[ \begin{array}{cc}
\frac{A}{V} \frac{\mu + \sigma_k \mu_t}{\rho} & 0 \\
0 & \frac{A}{V} \frac{\mu + \sigma_{\omega} \mu_t}{\rho}
\end{array} \right] \approx \frac{A}{V} \frac{\mu + \sigma_k \mu_t}{\rho} $$

$$\frac{\partial S}{\partial W} = \left[ \begin{array}{cc}
-2 \beta^* \omega & 0 \\
0 & -2 \beta \omega
\end{array} \right] \approx -2 \beta^* \omega $$

As can be seen above, only the destruction terms are included in the source term flux jacobian. This treatment is detailed in 
[Wilcox](http://www.amazon.com/Turbulence-Modeling-CFD-Third-Edition/dp/1928729088).
It is beneficial because the source terms will only contribute to the diagonal dominance of the system, helping with stability.

---
layout: post
title: "Evaluation of Eigen"
date: 2016-10-15 12:00
tags: [CFD, Aither, C++, v0.4.0, Eigen, matrix, linear algebra, blas, linpack, lapack, OpenBLAS, ATLAS, intel MKL]
comments: true
---
## Should Aither Use Eigen?
Aither, like many CFD codes requires matrix and vector operations when calculating the flow solution. When using block implicit
methods, Aither uses 5x5 matrices for the flow equations and 2x2 matrices if a turbulence model is selected. These matrices
are mutiplied with other matrices, multiplied with vectors, scaled by scalar values, added, subtracted, and inverted. If these
operations could be performed more efficiently, it would result in a great performance improvement. There are many third party
linear algebra libraries available that could be used by Aither such as [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page),
[Armadillo](http://arma.sourceforge.net/), and [PETSc](https://www.mcs.anl.gov/petsc/). PETSc is widely used and has the ability
to run in parallel. It contains many robust matrix solvers as well. However it is written in C, must be linked to, and is probably
overkill for small matrix/vector operations as described above. Eigen is also widely used and writen in C++. It has the added
advantage of being entirely header-only, so there is no need to link to anything. It could therefore be distributed with the Aither
source code eliminating the issue of findind a dependency on various computer systems. Armadillo is another linear algebra library
written in C++. It links with libraries such as LAPACK, OpenBLAS, MKL, or ATLAS. Since Eigen claims comprable or better performance
than many linear algebra libraries, and it has the best ease of use, it was choosen for evalution. Armadillo may be evaluated at
a later point.

## Test Cases
The following five test cases were used to evaluate Eigen versus the linear algebra code already in Aither.

* Matrix-matrix multiplication $$\left( A B = C \right)$$
* Matrix-vector multiplication $$\left( A \vec{x} = \vec{b} \right)$$
* Matrix multiplication with a scalar and addition $$\left( A s + B = C \right)$$
* Vector multiplication with a scalar and addition $$\left( \vec{x} s + \vec{y} = \vec{z} \right)$$
* Matrix inverse $$\left( A^{-1} = B \right)$$

The above tests were repeated 10 million times each using 5x5 and 2x2 matrices. Eigen has predefined classes for small (< 4)
matrices and vectors which statically allocate their memory on the stack. These can be quite a bit faster than the more general n-dimensional
matrices and vectors which dynamically allocate their memory on the heap. Aither requires general sized matrices because the matrix size is
determined at run time. For scalar implicit methods like LU-SGS and DPLUR the matrix size is 1. For block implicit methods like
BLU-SGS and BDPLUR the matrix size is 5. For the tests using the 2x2 matrices both the static `Eigen::Matrix2d` and dynamic `Eigen::MatrixXd`
versions were used, but the comparison to Aither was made using the heap version. Aither uses dynamic allocation for all matrices
`squareMatrix(5)`, `squareMatrix(2)`, but uses static allocation for vectors `genArray`. In Aither, all vectors are of size 7 which is the
maximum number of equations solved.


## Results
The timing results for each of the tests are shown below in seconds. There is no data for the vector multiplication with a scalar and addtion test for Aither
for the size 2 vector because in Aither all vectors are the same length. Therefore the test would be the same as for the larger vector.

| Matrix Type       | Size | Matrix-Matrix Multiplication | Matrix-Vector Multiplication | Matrix Scale & Addition | Vector Scale & Addition | Matrix Inverse |
|---                |---   |---                           |---                           |---                      |---                      |---             |
|`Eigen::MatrixXd`  |5x5   |2.87540                       |1.04917                       |2.30247e-1               |4.68587e-2               |9.86778         |
|`Eigen::Matrix2d`  |2x2   |6.05223e-3                    |3.09800e-3                    |4.20000e-8               |4.70000e-8               |3.15865e-3      |
|`Eigen::MatrixXd`  |2x2   |1.81787                       |6.09033e-1                    |2.65002e-1               |1.67215e-2               |3.70098         |
|`squareMatrix`     |5x5   |1.91885                       |2.39314e-1                    |1.81551                  |2.13851e-2               |6.99836         |
|`squareMatrix`     |2x2   |5.08120e-1                    |1.09769e-1                    |9.79889e-1               |N/A                      |1.45957         |

## Conclusions
Surprisingly Aither outperforms Eigen in all tests with the exception of the matrix multiplication with a scalar and addition test. The statically allocated
Eigen matrix and vector classes far outperform both the dynamically allocated Eigen classes and the Aither classes. This is expected as memory access is faster
to the stack than it is to the heap. Looking at the Eigen [benchmarks](http://eigen.tuxfamily.org/index.php?title=Benchmark), the best performance is expected
for larger matrices than the ones tested here. The results of these tests indicate that it would not be worthwhile to integrate Eigen into Aither. The source
code for these tests can be found [here](https://github.com/mnucci32/eigenVsAither).

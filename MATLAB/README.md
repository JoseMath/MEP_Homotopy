# FiberHomotopy.m
This is the main function.

## Inputs:

- A: k by (k+1) cell array containing all matrix coefficients
- opts: stuct, 'tracker': {'default','fast','conservative'}, 'EndpointN': int.

## Outputs:

- EigenValue: n1*…*nk by k cell, where each column is a lambda variable group
- EigenVector:  n1*…*nk by k cell, where each row is an eigenvector solution (x1,…,xk)
- SEigenValue: eigenvalues of the start system (the k GEPs)
- SEigenVector: eigenvectors of the start system (the k GEPs)
- TimeEachPath: timing of each path
- NumNewtonEachPath: the number of total Newton iterations in each path
- LastNewton: number of Newton iteration in the last step
- JacG: Jacobian matrix of the target system
- C: linear constraints on eigenvectors


# EvaluateHi.m
Sub function; the EvaluateHi function evaluate the Jacobian matrix with the most recent eigenpair


# RBackwardError.m
The function that computes the relative backward error.


## Inputs: 
The approximate eigenvalue, eigenvector and coefficient matrices

- EigenValue  : n_sols by 1 cell
- EigenVector : n_sols by k cell
- A           : k by k+1 cell, matrix coefficients

## Outputs:

- MaxRelError : n_sols by 1 matrix, maximum backward errors among k equations
- RelError    : n_sols by 1 matrix, backward errors
- Position    : n_sols by 1 matrix, where the MaxRelError obtains
- ri          : n\_sols by k matrix, the l2 norm of the residuals
- thetai      : n\_sols by k matrix, ||A\_{i0}||+sum\_k |lambda\_j|*||A\_{ij}||


# extract\_intrinsic\_part.m

The function that extracts intrinsic generalized eigenvalues and vectors from a GEP.
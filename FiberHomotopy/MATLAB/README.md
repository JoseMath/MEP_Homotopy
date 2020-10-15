# FiberHomotopy.m
This is the main function.

## Inputs:

- A: k by (k+1) cell array containing all matrix coefficients

## Outputs:

- EigenValue: n1*…*nk by k cell, where each column is a lambda variable group
- EigenVector:  n1*…*nk by k cell, where each row is an eigenvector solution (x1,…,xk)
- SEigenValue: eigenvalues of the start system (the k GEPs)
- SEigenVector: eigenvectors of the start system (the k GEPs)
- TimeEachPath: timing of each path
- LastNewton: number of Newton iteration in the last step
- JacG: Jacobian matrix of the target system
- C: linear constraints on eigenvectors
- Newtoniteration: total number of Newton iterations


## Example:

```
% Randomly generated system with k = 3, n1=n2=n3=3 (can just copy paste the code to try it)
k = 3;
ni = 3;
n = ni*ones(1,k);
Amatrix = rands(k*ni,(k+1)*ni)+1i*randn(k*ni,(k+1)*ni);
A = mat2cell(Amatrix,n,ni*ones(1,k+1));

[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,NewtonIteration] = FiberHomotopy(A); % naive method with default number of Newton iteration in the last step
[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,NewtonIteration,LastTValue]=restrictedFiberHomotopy(A,m,1000,.00001,.000001);
```

# EvaluateHi.m
Sub function; the EvaluateHi function evaluate the Jacobian matrix with the most recent eigenpair


# RBackwardError.m
The function that compute the relative backward error

## Inputs: 
The approximate eigenvalue, eigenvector and coefficient matrices

## Example:
```
	RBE = BackwardError(EigenValue,EigenVector,A);
```
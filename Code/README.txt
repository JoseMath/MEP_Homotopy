# FiberHomotopy.m
# The main function
	Inputs: 
	# A: k by (k+1) cell array containing all matrix coefficients
	# (optional) Number of Newton iteration in the last step: the default value is k*(max(n)+5)
	# (optional) Method: 1 — structured method, 0 — naive method (default)
	
	Outputs:
	# EigenValue: n1*…*nk by k cell, where each column is a lambda variable group  
	# EigenVector:  n1*…*nk by k cell, where each row is an eigenvector solution (x1,…,xk)
	# SEigenValue: eigenvalues of the start system (the k GEPs)
	# SEigenVector: eigenvectors of the start system (the k GEPs)
	# TimeEachPath: timing of each path
	# LastNewton: number of Newton iteration in the last step
	# JacG: Jacobian matrix of the target system
	# C: linear constraints on eigenvectors
	# Newtoniteration: total number of Newton iterations

	Example:
	% Randomly generated system with k = 3, n1=n2=n3 = 5 (can just copy paste the code to try it)
	k = 3;
	ni = 5;
	n = ni*ones(1,k);
	Amatrix = rands(k*ni,(k+1)*ni)+1i*randn(k*ni,(k+1)*ni);
	A = mat2cell(Amatrix,n,ni*ones(1,k+1));
 
  	[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,NewtonIteration] = FiberHomotopy(A); % naive method with default number of Newton iteration in the last step
	[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,NewtonIteration] = FiberHomotopy(A,50); % naive method with 50 Newton iteration in the last step (recommended; the default number is usually not the best one)
	[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,NewtonIteration] = FiberHomotopy(A,50,1); % structured method with 50 Newton iteration in the last step

#######################

# EvaluateHi.m
# Sub function; the EvaluateHi function evaluate the Jacobian matrix with the most recent eigenpair

#######################

# RBackwardError.m
# The function that compute the relative backward error

	# Inputs: the approximate eigenvalue, eigenvector and coefficient matrices

	# Example:
	RBE = BackwardError(EigenValue,EigenVector,A);




 











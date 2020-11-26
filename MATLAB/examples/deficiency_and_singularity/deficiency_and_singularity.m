addpath(genpath("MutliParEig"))
% Two-parameter eigenvalue problem
% A1*x = (lambda*B1+mu*C1)*x
% A2*x = (lambda*B2+mu*C2)*x

%% Case 1 - singular but not defective
A1 = [1  2;  3  4]; 
B1 = [1  1; -1  1]; 
C1 = [2  1;  5  1];
A2 = [1 -2;  3 -5]; 
B2 = [1 -1; -2  3]; 
C2 = [1 -1;  3  1];

% Delta0 operator determinant
Delta0 = kron(B1,C2) - kron(C1,B2);

% rank of matrix Delta0 is 3 -> Delta0 is singular
ran = rank(Delta0)

% extrinsic (2,2)
% intrinsic (2,2)
A = {A1,B1, C1;A2,B2,C2};
intrinsic(A)


%% Case 2 - defective and singular 
A1 = [0  2;  1  0]; 
B1 = [0  1; 0  0]; 
C1 = [0  1;  0  0];
A2 = [1 0;  0 0]; 
B2 = [1 0; 0  0]; 
C2 = [1 0;  0  0];

% Delta0 operator determinant
Delta0 = kron(B1,C2) - kron(C1,B2);

% rank of matrix Delta0 is 0 -> Delta0 is singular
ran = rank(Delta0)

% extrinsic (2,2)
% intrinsic (1,2)
A = {A1,B1, C1;A2,B2,C2};
intrinsic(A)
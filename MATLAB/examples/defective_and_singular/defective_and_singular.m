addpath(genpath("/Users/dujinhong/Documents/study/Lim/code/MutliParEig_2_5/MutliParEig"))

%% Case 1 - defective but not singular 



%% Case 2 - singular but not defective
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

A = {A1,B1, C1;A2,B2,C2};
intrinsic(A)


%% Case 3 - defective and singular 
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

A = {A1,B1, C1;A2,B2,C2};
intrinsic(A)
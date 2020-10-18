%% Quadratic multiparameter eigenvalue problems (QMEPs)

k = 2; % number of parameters
n = 2; % size of coefficient matrices of the QMEPs

% Q1(lambda,mu) = B1 + lambda B2 + mu B3 + lambda^2 B4 + lambda mu B5 + mu^2 B6
% Q2(lambda,mu) = C1 + lambda C2 + mu C3 + lambda^2 C4 + lambda mu C5 + mu^2 C6
B = randn(n,2*(k+1)*n)+1i*randn(n,2*(k+1)*n);
B = mat2cell(B, n, n*ones(1,2*(k+1)));
C = randn(n,2*(k+1)*n)+1i*randn(n,2*(k+1)*n);
C = mat2cell(C, n, n*ones(1,2*(k+1)));

% linearization to get coefficient matrices A of a two-parameter eigenvalue problem
A = cell(k,k+1);
A{1,1} = [B{1},B{2},B{3};
            zeros(2*n,n),-eye(2*n)];
A{1,2} = [zeros(n),B{4},B{5};
           eye(n) zeros(n,2*n);
           zeros(n,3*n)];
A{1,3} = [zeros(n,2*n),B{6};
            zeros(n,3*n);
            eye(n),zeros(n,2*n)];
A{2,1} = [C{1},C{2},C{3};
            zeros(2*n,n),-eye(2*n)];
A{2,2} = [zeros(n),C{4},C{5};
           eye(n) zeros(n,2*n);
           zeros(n,3*n)];
A{2,3} = [zeros(n,2*n),C{6};
            zeros(n,3*n);
            eye(n),zeros(n,2*n)];


[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,Newtoniteration] = FiberHomotopy(A);
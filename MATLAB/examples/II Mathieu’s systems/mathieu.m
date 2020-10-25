% Mathieu's Systems

addpath(genpath('MutliParEig_2_5/'))
addpath(genpath('FiberHmotopy/'))

%% generate data
% We choose parameters as : n1 = 18, n2 = 38, mode = 1, alpha = 2, beta = 1
% The dimension of the problem would be 18 * 38 = 684

k = 2;
N = 20;
P = (N-2)*(2*N-2);

[A1,B1,C1,A2,B2,C2,~,~,~,~,~,~,~,~] = mathieu_mep(2*N,N,1,2,1);
% K = 1;
A1 = A1+5*B1;
A2 = A2+5*B2;

A{1,1} = A1;
A{1,2} = B1;
A{1,3} = C1;
A{2,1} = A2;
A{2,2} = B2;
A{2,3} = C2;


%% Delta method
tic;
[EigenValueD,EigenVectorD,~] = multipareig(A);
tDelta = toc;
EigenValueD = mat2cell(EigenValueD,ones(1,P),k);
[MaxRelError,RelError,Position,ri,thetai] = RBackwardError(EigenValueD,EigenVectorD,A);
save("res_mathieu_delta.mat")  


%% Fiber Homotopy
tic;
[EigenValueN,EigenVectorN,SEigenValueN,SEigenVectorN,TimeEachPath,NumNewtonEachPath,LastNewton,JacG,C] = FiberHomotopy2(A, struct('tracker', 'conservative'));
tFiber = toc;
[MaxRelError,RelError,Position,ri,thetai] = RBackwardError(EigenValueN,EigenVectorN,A);
save("res_mathieu_fiber.mat")  

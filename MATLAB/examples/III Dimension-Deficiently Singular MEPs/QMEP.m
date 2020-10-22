%% Quadratic multiparameter eigenvalue problems (QMEPs)

addpath(genpath('MutliParEig_2_5/'))
addpath(genpath('FiberHmotopy/'))

%% generate data
rng('default')
k = 2; % number of parameters
n_list = [2 5 10 20 40]; % size of coefficient matrices of the QMEPs
random_matrices = cell(length(n_list), 10);
for i = 1:length(n_list)
    n = n_list(i);
    for n_simulation = 1:10
        % Q1(lambda,mu) = B1 + lambda B2 + mu B3 + lambda^2 B4 + lambda mu B5 + mu^2 B6
        % Q2(lambda,mu) = C1 + lambda C2 + mu C3 + lambda^2 C4 + lambda mu C5 + mu^2 C6
        B = randn(n,2*(k+1)*n)+1i*randn(n,2*(k+1)*n);
        B = mat2cell(B, n, n*ones(1,2*(k+1)));
        C = randn(n,2*(k+1)*n)+1i*randn(n,2*(k+1)*n);
        C = mat2cell(C, n, n*ones(1,2*(k+1)));

        % linearization to get coefficient matrices A of a two-parameter eigenvalue problem
        A = cell(k,k+1);
        A{1,1} = [B{1},B{2},B{3}; zeros(2*n,n),-eye(2*n)];
        A{1,2} = [zeros(n),B{4},B{5}; eye(n) zeros(n,2*n); zeros(n,3*n)];
        A{1,3} = [zeros(n,2*n),B{6}; zeros(n,3*n); eye(n),zeros(n,2*n)];
        A{2,1} = [C{1},C{2},C{3}; zeros(2*n,n),-eye(2*n)];
        A{2,2} = [zeros(n),C{4},C{5}; eye(n) zeros(n,2*n); zeros(n,3*n)];
        A{2,3} = [zeros(n,2*n),C{6}; zeros(n,3*n); eye(n),zeros(n,2*n)];
        
        random_matrices{i, n_simulation} = A;
    end
end
save('random_matrices_QMEP','random_matrices')


%% run the delta method
t = cell(length(n_list),10);
EigenValue = cell(length(n_list),10);
MaxRelError = cell(length(n_list),10);
RelError = cell(length(n_list),10);
Position = cell(length(n_list),10);
ri = cell(length(n_list),10);
thetai = cell(length(n_list),10);

for i = 1:length(n_list)
    res = [];
    n = n_list(i);
    for n_simulation = 1:10
        A = random_matrices{i, n_simulation, 1};
        tic;
        [lambda,mu,eta,X1,X2,X3,~,~,~] = threepareig(A{1,1},A{1,2},A{1,3},A{1,4},A{2,1},A{2,2},A{2,3},A{2,4},A{3,1},A{3,2},A{3,3},A{3,4});
        t{i,n_simulation} = toc;
        
        EigenValue{i,n_simulation} = [lambda,mu,eta]; EigenValue{i,n_simulation} = mat2cell(EigenValue{i,n_simulation}, ones(length(lambda),1), k);        
        EigenVector = [X1.',X2.',X3.']; EigenVector = mat2cell(EigenVector, ones(length(lambda),1), ones(k,1)*n); EigenVector=cellfun(@transpose,EigenVector,'UniformOutput',false);
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue{i,n_simulation},EigenVector,A);                     
    end
end
    
save("res_random_matrices_QMEP_delta.mat")        


%% run the fiber product homotopy method
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

t = cell(length(n_list),10);
MaxRelError = cell(length(n_list),10);
RelError = cell(length(n_list),10);
Position = cell(length(n_list),10);
ri = cell(length(n_list),10);
thetai = cell(length(n_list),10);
NumNewtonEachPath = cell(length(n_list),10);
TimeEachPath = cell(length(n_list),10);

for i = 1:length(n_list)
    res = [];
    n = n_list(i);
    for n_simulation = 1:10
        A = random_matrices{i, n_simulation, 1};
        tic;
        [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath{i,n_simulation},NumNewtonEachPath{i,n_simulation},LastNewton,~,~] = FiberHomotopy(A, struct('tracker', 'conservative'));
        t{i,n_simulation} = toc;
        
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue,EigenVector,A);                     
    end
end
    
save("res_random_matrices_QMEP_fiber.mat")        
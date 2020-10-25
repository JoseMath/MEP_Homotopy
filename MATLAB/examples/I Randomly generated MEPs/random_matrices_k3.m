addpath(genpath('MutliParEig_2_5/'))
addpath(genpath('FiberHmotopy/'))

%% generate data
rng('default')
n_list = [10 15 20 25 30];
k = 3;
random_matrices = cell(length(n_list), 10);
for i = 1:length(n_list)
    n = n_list(i);
    for n_simulation = 1:10
        A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
        A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));
        random_matrices{i, n_simulation} = A;
    end
end
save('random_matrices_k3','random_matrices')


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
        
        EigenValue{i,n_simulation} = [lambda,mu,eta]; EigenValue{i,n_simulation} = mat2cell(EigenValue{i,n_simulation}, ones(size(lambda,1),1), k);        
        EigenVector = [X1.',X2.',X3.']; EigenVector = mat2cell(EigenVector, ones(size(lambda,1),1), ones(k,1)*n); EigenVector=cellfun(@transpose,EigenVector,'UniformOutput',false);
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue{i,n_simulation},EigenVector,A);                     
    end
end
    
save("res_random_matrices_k3_delta.mat")        


%% run the fiber product homotopy method
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

t = cell(length(n_list),10);
EigenValue = cell(length(n_list),10);
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
        [EigenValue{i,n_simulation},EigenVector,SEigenValue,SEigenVector,TimeEachPath{i,n_simulation},NumNewtonEachPath{i,n_simulation},LastNewton,~,~] = FiberHomotopy(A, struct('tracker', 'fast'));
        t{i,n_simulation} = toc;
        
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue{i,n_simulation},EigenVector,A);                     
    end
end
    
save("res_random_matrices_k3_fiber.mat")        
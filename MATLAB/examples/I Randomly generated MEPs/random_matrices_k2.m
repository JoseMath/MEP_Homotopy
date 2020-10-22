addpath(genpath('MutliParEig_2_5/'))
addpath(genpath('FiberHmotopy/'))

%% generate data
rng('default')
n_list = [5 15 30 50];
k = 2;
random_matrices = cell(length(n_list), 10);
for i = 1:length(n_list)
    n = n_list(i);
    for n_simulation = 1:10
        A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
        A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));
        random_matrices{i, n_simulation} = A;
    end
end
save('random_matrices_k2','random_matrices')


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
        [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath{i,n_simulation},NumNewtonEachPath{i,n_simulation},LastNewton,~,~] = FiberHomotopy(A, struct('tracker', 'fast'));
        t{i,n_simulation} = toc;
        
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue,EigenVector,A);                     
    end
end
    
save("res_random_matrices_k2_fiber.mat")        
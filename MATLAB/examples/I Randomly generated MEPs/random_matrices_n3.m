addpath(genpath('MutliParEig_2_5/'))
addpath(genpath('FiberHmotopy/'))

%% generate data
rng('default')
k_list = 4:9;
n = 3;
random_matrices = cell(length(k_list), 10);
for i = 1:length(k_list)
    k = k_list(i);    
    for n_simulation = 1:10
        A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
        A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));
        random_matrices{i, n_simulation} = A;
    end
end
save('random_matrices_n3','random_matrices')


%% run the delta method
t = cell(length(k_list),10);
MaxRelError = cell(length(k_list),10);
RelError = cell(length(k_list),10);
Position = cell(length(k_list),10);
ri = cell(length(k_list),10);
thetai = cell(length(k_list),10);

for i = 1:length(k_list)
    k = k_list(i);
    for n_simulation = 1:10
        A = random_matrices{i, n_simulation};
        tic;
        [lambda,X,~] = multipareig(A);
        t{i,n_simulation} = toc;
                
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(mat2cell(lambda,ones(length(lambda),1),k),X,A);                     
    end
end
    
save("res_random_matrices_n3_delta.mat")  


%% run the fiber product homotopy method
pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

t = cell(length(k_list),10);
MaxRelError = cell(length(k_list),10);
RelError = cell(length(k_list),10);
Position = cell(length(k_list),10);
ri = cell(length(k_list),10);
thetai = cell(length(k_list),10);
NumNewtonEachPath = cell(length(k_list),10);
TimeEachPath = cell(length(k_list),10);

for i = 1:length(k_list)
    k = k_list(i);
    for n_simulation = 1:10
        A = random_matrices{i, n_simulation, 1};
        tic;
        [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath{i,n_simulation},NumNewtonEachPath{i,n_simulation},LastNewton,~,~] = FiberHomotopy(A, struct('tracker', 'fast'));
        t{i,n_simulation} = toc;
        
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue,EigenVector,A);                     
    end
end
    
% compute the distance between copies of lambda
delta = cell(length(k_list),10);
for i = 1:length(k_list)
    k = k_list(i);
    pn = size(EigenValue{i,1},1);
    
    for j = 1:10
        eigs = cell2mat(EigenValue{i,j});
        diff_eigs = eigs - eigs(:,1);
        diff_eigs = mat2cell(diff_eigs, ones(pn,1)*k, k);
        diff_eigs = cellfun(@(x) max(max(abs(x))), diff_eigs,'UniformOutput',false);
        delta{i, j}= mean(cell2mat(diff_eigs));
    end
end

save("res_random_matrices_n3_fiber.mat")  


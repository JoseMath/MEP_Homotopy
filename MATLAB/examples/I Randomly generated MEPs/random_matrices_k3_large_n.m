addpath(genpath('FiberHmotopy/'))

%% generate data
rng('default')
n_list = [30 70 150];
k = 3;
random_matrices = cell(length(n_list), 100);
for i = 1:length(n_list)
    n = n_list(i);
    for n_simulation = 1:100
        A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
        A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));
        random_matrices{i, n_simulation} = A;
    end
end
save('random_matrices_k3_large_n','random_matrices')


%% run the fiber product homotopy method
t = cell(length(n_list),10);
MaxRelError = cell(length(n_list),10);
RelError = cell(length(n_list),10);
Position = cell(length(n_list),10);
ri = cell(length(n_list),100);
thetai = cell(length(n_list),100);
NumNewtonEachPath = cell(length(n_list),100);
TimeEachPath = cell(length(n_list),100);

for i = 1:length(n_list)    
    n = n_list(i);
    for n_simulation = 1:100
        A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
        A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));
        
        tic;        
        [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath{i,n_simulation},NumNewtonEachPath{i,n_simulation},LastNewton,~,~] = FiberHomotopy_one(A, struct('tracker', 'fast'));
        t{i,n_simulation} = toc;
        
        [MaxRelError{i,n_simulation},RelError{i,n_simulation},Position{i,n_simulation},ri{i,n_simulation},thetai{i,n_simulation}] = RBackwardError(EigenValue,EigenVector,A);                     
    end
end
    
save("Lim/res_random_matrices_k3_nlarge_fiber.mat")        
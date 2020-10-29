%% Example 5 in the paper
k = 3;
n = 2;

rng(62)
A = randn(k*n,(k+1)*n) + 1i*randn(k*n,(k+1)*n);
A = mat2cell(A, n*ones(k,1), n*ones(k+1,1));

[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,Newtoniteration, T, EigenValueT] = FiberHomotopy_track(A, struct('tracker', 'default'));

save("example5.mat")
load("example5.mat")

colors = ['b','r', 'm'];
figure();
hold on
box on;
for path = 1:size(EigenValueT,1)
    t = T{path};
    eigs = cell2mat(EigenValueT{path}); % num_t * (k*k)
    

    i = 1;
    for copy=1:k
        plot(real(eigs((copy-1)*k+i, t>0.95)), imag(eigs((copy-1)*k+i, t>0.95)), 'color', colors(copy))
    end    

end
print(gcf, 'example5.eps' , '-depsc')


path = 3;
figure();
hold on
box on;
t = T{path};
eigs = cell2mat(EigenValueT{path}); % num_t * (k*k)
i = 1;
for copy=1:k
    plot(eigs((copy-1)*k+i,:), 'color', colors(copy))
end   
print(gcf, 'example5_path1.eps' , '-depsc')
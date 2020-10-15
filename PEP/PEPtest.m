d = 2;
n = 200;
Amatrix = randn(n,(d+1)*n)+1i*randn(n,(d+1)*n);
A = mat2cell(Amatrix,n,n*ones(1,d+1));

% A{2} = A{2}*100000;

[EigenValue,EigenVector,TimeEachPath,NewtonIteration,c] = HomotopyPEP(A);

residual = 1:(n*d);
for i = 1:n*d
    Plambda = A{d+1};
    for j = 1:d
        Plambda = Plambda + EigenValue(i)^(d+1-j)*A{j};
    end
    residual(i) = norm(Plambda*EigenVector{i})/norm(EigenVector{i});
end

scale = norm(A{d+1})*ones(1,n*d);
for i = 1:n*d
    for j = 1:d
        scale(i) = scale(i) + abs(EigenValue(i))^(d+1-j)*norm(A{j});
    end
end

RBE_H = residual./scale;
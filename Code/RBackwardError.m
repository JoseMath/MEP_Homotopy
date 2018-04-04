function [MaxError,Error,Position,ri,thetai] = RBackwardError(EigenValue,EigenVector,A)
% prodn = n1*n2*...*nk
% EigenValue: prodn by 1 cell
% EigenVector: prodn by k cell
% A: matrix coefficients
k = size(A,1);
n = 1:k;
for i = 1:k
    n(i) = size(A{i,1},1);
end
prodn = prod(n);
Error = zeros(prodn,k);
MaxError = 1:prodn;
Position = 1:prodn;
ri = zeros(prodn,k);
thetai = zeros(prodn,k);
E = 1:k;
for i = 1:prodn
    EValue = EigenValue{i};
    EVector = EigenVector(i,:);
    for j = 1:k
        r = A{j,1};
        N = norm(A{j,1});
        for s = 1:k
            r = r - A{j,s+1}*EValue(s);
            N = N + norm(A{j,s+1})*abs(EValue(s));
        end
        r = r * EVector{j}/norm(EVector{j});
        ri(i,j) = norm(r);
        thetai(i,j) = N;
        E(j) = norm(r)/N;
        Error(i,j) = E(j);
    end
    [MaxError(i),Position(i)] = max(E);
end
    
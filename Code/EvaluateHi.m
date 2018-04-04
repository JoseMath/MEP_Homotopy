function [JacTL,JacTR] = EvaluateHi(k,A,S,n)
% The EvaluateHi function evaluate the Hi functions at the most recent
% eigenvalues Lambda, which should be a k by 1 cell array. 
% A is the matrix coeffient cell array;
% n is the dimension vector
% C is the linear constraint cell on eigenvectors
% the output would be a 1 by k cell array; the ci (linear constraints on 
% eigenvectors are also attached in the second row)
JacTL = cell(1,k);
JacTR = cell(1,k);
for i = 1:k
    JacTL{i} = A{i,1};
    for j = 1:k
        JacTL{i} = JacTL{i} - A{i,j+1}*S{k+i}(j);
    end
    
    R = zeros(n(i),k);
    for j = 1:k
        R(:,j) = -A{i,j+1}*S{i};
    end
    JacTR{i} = R;
end
end
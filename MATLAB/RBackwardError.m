function [MaxRelError,RelError,Position,ri,thetai] = RBackwardError(EigenValue,EigenVector,A)
    % Params:
    %    EigenValue  : n_sols by 1 cell
    %    EigenVector : n_sols by k cell
    %    A           : k by k+1 cell, matrix coefficients
    % Returns:
    %    MaxRelError : n_sols by 1 matrix, maximum backward errors among k equations
    %    RelError    : n_sols by k matrix, backward errors
    %    Position    : n_sols by 1 matrix, where the MaxRelError obtains
    %    ri          : n_sols by k matrix, the l2 norm of the residuals
    %    thetai      : n_sols by k matrix, ||A_{i0}||+sum_k |lambda_j|*||A_{ij}||
    
    k = size(A,1);
    n_sols = length(EigenValue);

    RelError = zeros(n_sols,k);
    MaxRelError = 1:n_sols;
    Position = 1:n_sols;
    ri = zeros(n_sols,k);
    ri_inf = zeros(n_sols,k);
    thetai = zeros(n_sols,k);
    E = 1:k;
    
    for i = 1:n_sols
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
            ri_inf(i,j) = norm(r, inf);
            ri(i,j) = norm(r);
            thetai(i,j) = N;
            E(j) = norm(r)/N;
            RelError(i,j) = E(j);
        end
        % exclude paths that converge to Inf or do not converge, in terms of infinity norm
        if max(abs(EValue))>1e6 || max(max(abs(cell2mat(EVector))))>1e6 || max(ri_inf(i,:))>1.0
            ri(i,:) = nan;
        end
        [MaxRelError(i),Position(i)] = max(E);
    end
    
    MaxRelError(isnan(ri(:,1))) = nan;
    RelError(isnan(ri(:,1)),:) = nan;
    Position(isnan(ri(:,1))) = nan;
    thetai(isnan(ri(:,1)),:) = nan;
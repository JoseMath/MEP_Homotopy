function [V, D] = extract_intrinsic(GA, GB)
    % Extract the intrinsic generalized eigenvalues and generalized eigenvectors
    % of a GEP (GA, GB).
   
    
    
    [V,D] = eig(GA,GB);
    
    % Remove Inf eigenvalues
    D = diag(D);   
    V(:, isinf(D)) = [];
    D(isinf(D)) = [];
    
    % Sort eigenvalues by magnitude first, then imaginary part.
    [D, sortIndex] = sort(D,'ascend');     
    V = V(:,sortIndex);
   
    % Remove duplicates
    real_tol = 1e-8;
    imag_tol = 1e-8;
    dD_sorted = diff(D);     
    remove_idx = (abs(real(dD_sorted))<real_tol) & (abs(imag(dD_sorted))<imag_tol);    
    D(remove_idx) = [];
    V(:,remove_idx) = [];
return
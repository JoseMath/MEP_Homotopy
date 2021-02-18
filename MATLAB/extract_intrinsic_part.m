function [V, D] = extract_intrinsic_part(GA, GB)
    real_tol = 1e-8;
    imag_tol = 1e-8;
    
    % Extract the intrinsic generalized eigenvalues and generalized eigenvectors
    % of a GEP (GA, GB).   
    
    % Count the infinity generalized eigenvalues
    [~,D] = eig(GB,GA,'qz','vector');
    numD = sum(abs(real(D))<real_tol) & (abs(imag(D))<imag_tol);
    [V,D] = eig(GA,GB,'qz','vector');
    
    % Remove Inf eigenvalues
    D = sort(D);   
    D(end-numD+1:end) = Inf;
    V(:, isinf(D)) = [];
    D(isinf(D)) = [];
    
    % Sort eigenvalues by magnitude first, then imaginary part.
    [D, sortIndex] = sort(D,'ascend');     
    V = V(:,sortIndex);
   
    % Remove duplicates    
    dD_sorted = diff(D);     
    remove_idx = (abs(real(dD_sorted))<real_tol) & (abs(imag(dD_sorted))<imag_tol);    
    D(remove_idx) = [];
    V(:,remove_idx) = [];
return
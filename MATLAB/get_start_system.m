function [SEigenVector,SEigenValue, JacG, JacL, Lc, C, nCell, sCell, m] = get_start_system(A, n, sumn, k)
    %% Compute JacG, JacL, Lc from randomly generated UG, RG, UL, RL
    % G_i(lambda_i) = JacG_i * lambda_i
    % L_i(lambda_i) = JacL_i * lambda_i + Lc_i
    ik = (k-1)*k;
    UG = [eye(ik),zeros(ik,k)]-[zeros(ik,k),eye(ik)];
    RG = randn(ik)+1i*randn(ik);
    R = mat2cell(randn(ik,k)+1i*randn(ik,k),(k-1)*ones(1,k),k);
    UL = blkdiag(R{:});
    RL = randn(ik)+1i*randn(ik);
    JacG = RG*UG;
    JacL = RL*UL;

    % the constant term of L1,...,Lk
    Lc = RL*ones(ik,1);

    % linear constraints on eigenvectors, i.e, c1,...ck
    C = mat2cell(randn(1,sumn)+1i*randn(1,sumn),1,n);
    
    % get the solutions for the start system, i.e, solve the GEPs
    SEigenVector = cell(k,1);
    SEigenValue = cell(k,1);
    s = UL\ones(ik,1);   % k^2 by 1 specific solution
    nCell = cell(k,1);
    for i = 1:k
        nCell{i} = null(R{i});
    end
    sCell = mat2cell(s,k*ones(1,k),1);

    m = zeros(size(n));  
    for i = 1:k
        GA = A{i,1};
        GB = 0;
        for j = 1:k
            GA = GA - A{i,j+1}*sCell{i}(j);
            GB = GB + A{i,j+1}*nCell{i}(j);
        end
        
        [V, D] = extract_intrinsic_part(GA, GB);
        m(i) = length(D);

        ScaleEigenVector = C{i}*V;
        V = V./repmat(ScaleEigenVector,n(i),1);
        SEigenVector{i} = V;   % n(i) by m(i) square matrix
        SEigenValue{i} = D;
    end
end

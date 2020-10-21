function [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,NumNewtonEachPath,LastNewton,JacG,C] = FiberHomotopy(A)
    % This is the parallelizable version of the fiber product homotopy.    
    %
    % Params:
    %   A                   - k by (k+1) cell, the matrix coefficients of the MEP.
    % Returns:
    %   EigenValue          - n_path by k cell, eigenvalues of the target system.
    %   EigenVector         - n_path by k cell, eigenvectors of the target system.
    %   SEigenValue         - k by 1 cell, eigenvalues of the start system.
    %   SEigenVector        - k by 1 cell, eigenvalues of the start system.
    %   TimeEachPath        - 1 by n_path array, time of each path.
    %   NumNewtonEachPath   - 1 by n_path array, total number of Newton steps of each path.
    %   LastNewton          - 1 by n_path array, the number of Newton steps of each path near the target system.
    %   JacG                - k(k-1) by k^2 cell, the coefficient matrix of the G.
    %   C                   - 1 by k cell, the linear constraints on eigenvectors.                     

    %% debugDisp is set to 1 when debugging
    debugDisp = 0;

    %% Set the number of eigenparameters and dimensions.
    % k is the number of eigenparameters 
    k = size(A,1);
    % n is the extrinsic dimension
    n = 1:k;
    for i = 1:k
        n(i) = size(A{i,1},1);
    end
    sumn = sum(n);

    if debugDisp==1
        disp(k);
        disp(n);
    end

    %% EndpointN is the number of Newton iterations at the endpoint. 
    %EndpointN = k*(max(n)+5);
    %EndpointN = max(20,EndpointN);
    EndpointN=5;

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
    CBlock = [blkdiag(C{:}),zeros(k,k^2)];
    CBlockDiag = blkdiag(C{:});

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
        %option 1
        for j = 1:k
            GA = GA - A{i,j+1}*sCell{i}(j);
            GB = GB + A{i,j+1}*nCell{i}(j);
        end
        
        [V, D] = extract_intrinsic_part(GA, GB);
        m(i) = length(D);

        ScaleEigenVector = C{i}*V;
        V = V./repmat(ScaleEigenVector,n(i),1);
        SEigenVector{i} = V;   % n(i) by m(i) square matrix
        SEigenValue{i} = diag(D);
    end

    if debugDisp==1
        disp(m)
    end


    % loop over all possible paths
    loop = cell(k,1);
    for i = 1:k
        loop{i} = 1:m(i);
    end
    Loop = allcomb(loop{:})';
    pn = prod(m);

    hmax = 10^-2;
    hmin = 10^-6;

    Eterm = ones(k,1);
    JD = -JacL+JacG;


    EigenValue = cell(pn,k);
    EigenVector = cell(pn,k);

    TimeEachPath = 1:pn;
    LastNewton = 1:pn;
    NumNewtonEachPath = 1:pn;
    
    
    % Naive Method
    parfor index = 1:size(Loop, 2)
        path = Loop(:,index);
        
        warning('off', 'MATLAB:illConditionedMatrix')
        warning('off', 'MATLAB:singularMatrix')
        warning('off', 'MATLAB:nearlySingularMatrix')
        
        S = cell(2*k,1);
        for i = 1:k
            S{i} = SEigenVector{i}(:,path(i));
            S{k+i} = nCell{i}*SEigenValue{i}(path(i))+sCell{i};
        end

        Newtoniteration = 0;
        t = 0;
        h = hmax;

        [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
        Jac = [JacH;CBlock;zeros(ik,sumn),JacL];
    
        tic;
        while t <= 1 - h && sum(sum(isnan(cell2mat(S))))==0 %-hmax
            vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];

            DeltaEulorStep = Jac\(-vEulor); % "-" or not?

            % update S and t
            S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);
            t = t + h;

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            JacM = (1-t)*JacL+t*JacG;
            ML = (1-t)*Lc;
            NumberNewton = 1;

            % Newton step
            % Naive way
            EVector = cell2mat(S(1:k));
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
            
            JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
            Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
            DeltaNewtonStep = Jac\(-vNewton);
            S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            EVector = cell2mat(S(1:k));
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                
            
            while NumberNewton < 150 && norm(DeltaNewtonStep,inf) >= 10^-9 && sum(sum(isnan(cell2mat(S))))==0
                JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
                Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
                DeltaNewtonStep = Jac\(-vNewton);
                S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

                [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                EVector = cell2mat(S(1:k));
                vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                NumberNewton = NumberNewton + 1;
            end

            if NumberNewton <= 2 && h <= hmax/2
                h = h*2;
            else
                if NumberNewton == 8 && h >= hmin*2
                    h = h/2;
                end
            end
        end
        Newtoniteration = Newtoniteration + NumberNewton;
      
        % one more step when t gets larger than 1
        if t < 1 && sum(sum(isnan(cell2mat(S))))==0
            JacM = (1-t)*JacL+t*JacG;
            Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
            vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];

            DeltaEulorStep = Jac\(-vEulor); % "-" or not?

            % update S and t
            h = 1 - t;
            S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            JacM = JacG;
            ML = 0;
            NumberNewton = 1;
        
            % Newton step
            % Naive way
            EVector = cell2mat(S(1:k));
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
                        
            while NumberNewton < EndpointN && norm(DeltaEulorStep,inf) >= 10^-9 && sum(sum(isnan(cell2mat(S))))==0 %EndpointN
                JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
                Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
                DeltaNewtonStep = Jac\(-vNewton);
                S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

                [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                EVector = cell2mat(S(1:k));
                vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                NumberNewton = NumberNewton + 1;
            end
            LastNewton(index) = NumberNewton;
        end
        Newtoniteration = Newtoniteration + NumberNewton;
        NumNewtonEachPath(index) = Newtoniteration;
        
        TimeEachPath(index) = toc;
        for vg = 1:k
            EigenVector{index,vg} = S{vg};
            EigenValue{index,vg} = S{k+vg};
        end
    end  
end
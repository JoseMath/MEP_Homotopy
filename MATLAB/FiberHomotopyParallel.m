function [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,NumNewtonEachPath,LastNewton,JacG,C] = FiberHomotopy(A, opts)
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

    % Validate number of input parameters
    narginchk(1, 2);

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
    
    % Analyse user supplied options, if any.
    if nargin < 2, opts = []; end
    if ~isfield(opts,'tracker'), opts.tracker = 'default'; end
    switch(opts.tracker)
        case 'default' 
            EndpointN = 20; % the number of Newton iterations at the endpoint.         
            NormRes = Inf; % the kinds of norm for the residual
            EpsNormRes = 1e-8; % maximum value of norm of the residual required to obtain
            MinNumberNewtonRetain = 3; % minimum number of Newton iterations to retain the same step size h
            MaxNumberNewtonRetain = 5; % maximum number of Newton iterations to retain the same step size h
        case 'fast' 
            EndpointN = 5;
            NormRes = 2;
            EpsNormRes = 1e-8;
            MinNumberNewtonRetain = 3;
            MaxNumberNewtonRetain = 5;
        case 'conservative' 
            EndpointN = max(20, k*(max(n)+5));
            NormRes = Inf;
            EpsNormRes = 1e-9;
            MinNumberNewtonRetain = 3;
            MaxNumberNewtonRetain = 8;
        otherwise
            error('Invalid tracker\n');        
    end
    if isfield(opts,'EndpointN'), EndpointN = opts.EndpointN; end

    %% Compute start system
    % In theory we make the assumption that the intrinsic dimension is measured correctly. 
    % But the intrinsic dimension can be over-estimated when the eigenvalues of an associated GEP 
    % at infinity incorrectly classified as eigenvalues with very large modulus. When this happens, 
    % the homotopy method method will track this start point to infinity or to a copy of an eigenpair of the MEP. 
    ik = (k-1)*k;  
    SEigenVector=cell(3,1);SEigenValue=cell(3,1);JacG=cell(3,1);JacL=cell(3,1);Lc=cell(3,1);C=cell(3,1);nCell=cell(3,1);sCell=cell(3,1);m=cell(3,1);
    for i=1:3
        [SEigenVector{i},SEigenValue{i}, JacG{i}, JacL{i}, Lc{i}, C{i}, nCell{i}, sCell{i}, m{i}] = get_start_system(A, n, sumn, k);
    end
    max_m = max(reshape(cell2mat(m), [k,3]), [], 2);
    min_m = min(reshape(cell2mat(m), [k,3]), [], 2);
    if sum(abs(max_m-min_m))~=0       
        [~,max_i] = max(prod(reshape(cell2mat(m), [k,3]), 1));
        disp('Discrepancy on estimations of intrinsic dimension.')
        disp(['Minimum: [',num2str(min_m.'),'] Maximum: [',num2str(max_m.'),'].'])
        disp(['Using [',num2str(m{max_i}),']. This may cause extra divergent paths or paths leading to multiple copies.'])
    else
        max_i = 1;     
    end
    SEigenVector=SEigenVector{max_i}; SEigenValue=SEigenValue{max_i}; JacG=JacG{max_i}; JacL=JacL{max_i}; Lc=Lc{max_i}; C=C{max_i}; nCell=nCell{max_i}; sCell=sCell{max_i}; m = m{max_i};
    CBlock = [blkdiag(C{:}),zeros(k,k^2)];
    CBlockDiag = blkdiag(C{:});

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
        while t <= 1 - h%-hmax
            vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];

            DeltaEulorStep = Jac\(-vEulor);
            
            % update S and t
            S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);
            t = t + h;

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            JacM = (1-t)*JacL+t*JacG;
            ML = (1-t)*Lc;
            NumberNewton = 1;

            % Newton step
            EVector = cell2mat(S(1:k));
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
            
            JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
            Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
            DeltaNewtonStep = Jac\(-vNewton);
            Newtoniteration = Newtoniteration + 1;
            S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            EVector = cell2mat(S(1:k));
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                
            
            while NumberNewton < 150 && norm(DeltaNewtonStep,NormRes) >= EpsNormRes
                JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
                Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
                DeltaNewtonStep = Jac\(-vNewton);
                Newtoniteration = Newtoniteration + 1;
                S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

                [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                EVector = cell2mat(S(1:k));
                vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                NumberNewton = NumberNewton + 1;
            end

            if NumberNewton < MinNumberNewtonRetain && h <= hmax/2
                h = h*2;
            else
                if NumberNewton == MaxNumberNewtonRetain && h >= hmin*2
                    h = h/2;
                end
            end
        end
        Newtoniteration = Newtoniteration + NumberNewton;

        % one more step when t gets larger than 1
        if t < 1 
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
            EVector = cell2mat(S(1:k));
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
                        
            while NumberNewton < EndpointN && norm(DeltaEulorStep,NormRes) >= EpsNormRes
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
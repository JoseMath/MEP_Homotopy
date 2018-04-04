function [EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,LastNewton,JacG,C,Newtoniteration] = FiberHomotopy2(A,varargin)
% the input A should be a k by (k+1) cell array which contains all the
% matrix coefficients

% to get the number of parameters k, and the dimensions n=(n1,...,nk)
k = size(A,1);
n = 1:k;
for i = 1:k
    n(i) = size(A{i,1},1);
end
sumn = sum(n);
EndpointN = k*(max(n)+5);
EndpointN = max(20,EndpointN);

% get UG, RG, UL, RL --> JacG, JacL, Lc
ik = (k-1)*k;
UG = [eye(ik),zeros(ik,k)]-[zeros(ik,k),eye(ik)];
RG = randn(ik)+1i*randn(ik);
% [RG,~,~] = svd(RG); 
R = mat2cell(randn(ik,k)+1i*randn(ik,k),(k-1)*ones(1,k),k);
UL = blkdiag(R{:});
RL = randn(ik)+1i*randn(ik);
% [RL,~,~] = svd(RL); 
JacG = RG*UG;
JacL = RL*UL;

Lc = RL*ones(ik,1);
Newtoniteration = 0;

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


for i = 1:k
    GA = A{i,1};
    GB = 0;
    %option 1
    for j = 1:k
        GA = GA - A{i,j+1}*sCell{i}(j);
        GB = GB + A{i,j+1}*nCell{i}(j);
    end
    
    % option 2: without the for loop; should compare the speed later
    % MC = cell2mat(A{i,2:(k+1)});
    % GA = GA - MC*kron(s(((i-1)*k+1):i*k),eye(n(i)));
    % GB = GB + MC*kron(n(((i-1)*k+1):i*k),eye(n(i)));
    [V,D] = eig(GA,GB);
    ScaleEigenVector = C{i}*V;
    V = V./repmat(ScaleEigenVector,n(i),1);
    SEigenVector{i} = V;   % n(i) by n(i) square matrix
    SEigenValue{i} = diag(D);
end

% track each path
% for all n1*n2*...*nk start solutions (to be completed)
% let's just fix a start point for now
%  S = cell(2*k,1);

% test one path
%     path = n(1);
%     for i = 1:k
%         S{i} = SEigenVector{i}(:,path);
%         S{k+i} = nCell{i}*SEigenValue{i}(path)+sCell{i};
%     end

S = cell(2*k,1);
% loop over all possible paths
loop = cell(k,1);
for i = 1:k
    loop{i} = 1:n(i);
end
Loop = allcomb(loop{:})';
pn = prod(n);
rs = ceil(pn/50);
Sample = randsample(pn,rs);
Loop = Loop(:,Sample);

hmax = 10^-2;
hmin = 10^-6;

Eterm = ones(k,1);
JD = -JacL+JacG;


EigenValue = cell(rs,k);
EigenVector = cell(rs,k);

TimeEachPath = 1:rs;
LastNewton = 1:rs;
index = 1;
if nargin == 1 % Structured method
    for path = Loop
        for i = 1:k
            S{i} = SEigenVector{i}(:,path(i));
            S{k+i} = nCell{i}*SEigenValue{i}(path(i))+sCell{i};
        end
        
        t = 0;
        h = hmax;

        [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
        Jac = [JacH;CBlock;zeros(ik,sumn),JacL];


        % Euler step

        % StopByNumberofNewton = [];
        % NumberofTotalOuterloop = 0;
    

    %if nargin == 1
        tic;
        while t <= 1 - h
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
            % Naive way
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            
            
            Decomp = cell(1,k);
                for i = 1:k
                    Decomp{i} = inverse(JacTL{i});
                end
                Constraints = cell(1,k);
                for i = 1:k
                    Cons = Decomp{i}*JacTR{i};
                    Constraints{i} = -C{i}*Cons;
                end
                J = [JacM;blkdiag(Constraints{:})];
                DeltaNewtonLambda = J\[-RJacM;Eterm];
                DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
                maxDelta = 0;
                for i = 1:k
                    Rside = JacTR{i}*DeltaNewtonLambda{i};
                    Delta = -Decomp{i}*Rside - S{i};
                    maxDelta = max(maxDelta,norm(Delta,inf));
                    S{i} = S{i} + Delta;
                end
                for i = 1:k
                    maxDelta = max(maxDelta,norm(DeltaNewtonLambda{i},inf));
                    S{i+k} = S{i+k}+DeltaNewtonLambda{i};
                end

                [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
                
            
            while NumberNewton < 7 && maxDelta >= 10^-6
    %             JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
    %             Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
    %             DeltaNewtonStep = Jac\(-vNewton);
    %             S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);
    % 
    %             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
    %             EVector = cell2mat(S(1:k));
    %             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
    %             NumberNewton = NumberNewton + 1;
    %         end

    % 
    % 
    % 
    %         structured method
                Decomp = cell(1,k);
                for i = 1:k
                    Decomp{i} = inverse(JacTL{i});
                end
                Constraints = cell(1,k);
                for i = 1:k
                    Cons = Decomp{i}*JacTR{i};
                    Constraints{i} = -C{i}*Cons;
                end
                J = [JacM;blkdiag(Constraints{:})];
                flag = decomposition(J);
                %if rcond(J) < eps
                if isIllConditioned(flag)
                  break
                end
                DeltaNewtonLambda = flag\[-RJacM;Eterm];
                DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
                maxDelta = 0;
                for i = 1:k
                    Rside = JacTR{i}*DeltaNewtonLambda{i};
                    Delta = -Decomp{i}*Rside - S{i};
                    maxDelta = max(maxDelta,norm(Delta,inf));
                    S{i} = S{i} + Delta;
                end
                for i = 1:k
                    maxDelta = max(maxDelta,norm(DeltaNewtonLambda{i},inf));
                    S{i+k} = S{i+k}+DeltaNewtonLambda{i};
                end

                [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
                NumberNewton = NumberNewton + 1;
            end

%             if NumberNewton == 7 && norm(vNewton) >= 10^-9
%                 StopByNumberofNewton = [StopByNumberofNewton,1];
%             else
%                 StopByNumberofNewton = [StopByNumberofNewton,0];
%             end
        % Block method of the Newton step
        % for i = 1:k
        %     vn = null(JacTL{i}.');
            if NumberNewton <= 2 && h <= hmax/2
                h = h*2;
            else
                if NumberNewton == 7 && h >= hmin*2
                    h = h/2;
                end
            end
%         NumberofTotalOuterloop = NumberofTotalOuterloop + 1;
            JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
            Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
        end
    
    
    
    
    
    %     % one more step when t gets larger than 1
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
            % Naive way
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            
            
            Decomp = cell(1,k);
                for i = 1:k
                    Decomp{i} = inverse(JacTL{i});
                end
                Constraints = cell(1,k);
                for i = 1:k
                    Cons = Decomp{i}*JacTR{i};
                    Constraints{i} = -C{i}*Cons;
                end
                J = [JacM;blkdiag(Constraints{:})];
                flag = decomposition(J);
                if ~isIllConditioned(flag)
                    DeltaNewtonLambda = flag\[-RJacM;Eterm];
                    DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
                    maxDelta = 0;
                    for i = 1:k
                        Rside = JacTR{i}*DeltaNewtonLambda{i};
                        Delta = -Decomp{i}*Rside - S{i};
                        maxDelta = max(maxDelta,norm(Delta,inf));
                        S{i} = S{i} + Delta;
                    end
                    for i = 1:k
                        maxDelta = max(maxDelta,norm(DeltaNewtonLambda{i},inf));
                        S{i+k} = S{i+k}+DeltaNewtonLambda{i};
                    end

                    [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                    RJacM = JacM*cell2mat(S((k+1):2*k))-ML;


                    while NumberNewton < 80 && maxDelta >= 10^-9
            %             JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
            %             Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
            %             DeltaNewtonStep = Jac\(-vNewton);
            %             S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);
            % 
            %             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            %             EVector = cell2mat(S(1:k));
            %             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
            %             NumberNewton = NumberNewton + 1;
            %         end

            % 
            % 
            % 
            %         structured method
                        Decomp = cell(1,k);
                        for i = 1:k
                            Decomp{i} = inverse(JacTL{i});
                        end
                        Constraints = cell(1,k);
                        for i = 1:k
                            Cons = Decomp{i}*JacTR{i};
                            Constraints{i} = -C{i}*Cons;
                        end
                        J = [JacM;blkdiag(Constraints{:})];
                        flag = decomposition(J);
                        % if rcond(J) < eps
                        if isIllConditioned(flag) 
                          break
                        end
                        DeltaNewtonLambda = flag\[-RJacM;Eterm];
                        DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
                        maxDelta = 0;
                        for i = 1:k
                            Rside = JacTR{i}*DeltaNewtonLambda{i};
                            Delta = -Decomp{i}*Rside - S{i};
                            maxDelta = max(maxDelta,norm(Delta,inf));
                            S{i} = S{i} + Delta;
                        end
                        for i = 1:k
                            maxDelta = max(maxDelta,norm(DeltaNewtonLambda{i},inf));
                            S{i+k} = S{i+k}+DeltaNewtonLambda{i};
                        end

                        [JacTL,JacTR] = EvaluateHi(k,A,S,n);
                        RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
                        NumberNewton = NumberNewton + 1;
                    end
                end

%         if NumberNewton == 7 && norm(vNewton) >= 10^-9
%             StopByNumberofNewton = [StopByNumberofNewton,1];
%         else
%             StopByNumberofNewton = [StopByNumberofNewton,0];
%         end
%     end

    
    
    
            % shrink step size
%     hmax = 10^-7;
%     hmin = 10^-6;
%     if h > hmax
%         h = hmax;
%     end
%     while t < 1
%         vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];
% 
%         DeltaEulorStep = Jac\(-vEulor); % "-" or not?
% 
%         % update S and t
%         S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);
%         t = t + h;
% 
%         [JacTL,JacTR] = EvaluateHi(k,A,S,n);
%         JacM = (1-t)*JacL+t*JacG;
%         ML = (1-t)*Lc;
%         NumberNewton = 0;
% 
%         % Newton step
%         % Naive way
%         EVector = cell2mat(S(1:k));
%         RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
%         vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
%         while NumberNewton < 7 && norm(vNewton) >= 10^-9
% %             JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
% %             Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
% %             DeltaNewtonStep = Jac\(-vNewton);
% %             S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);
% % 
% %             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
% %             EVector = cell2mat(S(1:k));
% %             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
% %             NumberNewton = NumberNewton + 1;
% %         end
% 
% % 
% % 
% % 
% %         structured method
%             Decomp = cell(1,k);
%             for i = 1:k
%                 Decomp{i} = inverse(JacTL{i});
%             end
%             Constraints = cell(1,k);
%             for i = 1:k
%                 Cons = Decomp{i}*JacTR{i};
%                 Constraints{i} = -C{i}*Cons;
%             end
%             J = [JacM;blkdiag(Constraints{:})];
%             DeltaNewtonLambda = J\[-RJacM;Eterm];
%             DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
%             for i = 1:k
%                 Rside = JacTR{i}*DeltaNewtonLambda{i};
%                 S{i} = -Decomp{i}*Rside;
%             end
%             for i = 1:k
%                 S{i+k} = S{i+k}+DeltaNewtonLambda{i};
%             end
%             
%             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
%             EVector = cell2mat(S(1:k));
%             RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
%             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
%             NumberNewton = NumberNewton + 1;
%         end
% 
%         if NumberNewton == 7 && norm(vNewton) >= 10^-9
%             StopByNumberofNewton = [StopByNumberofNewton,1];
%         else
%             StopByNumberofNewton = [StopByNumberofNewton,0];
%         end
%     % Block method of the Newton step
%     % for i = 1:k
%     %     vn = null(JacTL{i}.');
%         if NumberNewton <= 2 && h <= hmax/2
%             h = h*2;
%         else if NumberNewton == 7 && h >= hmin*2
%                 h = h/2;
%             end
%         end
%         NumberofTotalOuterloop = NumberofTotalOuterloop + 1;
%     end
%     
%     
        end
        TimeEachPath(index) = toc;
        
        % save results
        for vg = 1:k
            EigenVector{index,vg} = S{vg};
            EigenValue{index,vg} = S{k+vg};
        end
        index = index + 1;
        index
    end
    
    
else  % Naive Method
    for path = Loop
        for i = 1:k
            S{i} = SEigenVector{i}(:,path(i));
            S{k+i} = nCell{i}*SEigenValue{i}(path(i))+sCell{i};
        end
        
        t = 0;
        h = hmax;

        [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
        Jac = [JacH;CBlock;zeros(ik,sumn),JacL];
    
        tic;
        while t <= 1 - h%-hmax
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
            Newtoniteration = Newtoniteration + 1;
            S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);

            [JacTL,JacTR] = EvaluateHi(k,A,S,n);
            EVector = cell2mat(S(1:k));
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
                
            
            while NumberNewton < 8 && norm(DeltaNewtonStep,inf) >= 10^-9
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




            % structured method
        %         Decomp = cell(1,k);
        %         for i = 1:k
        %             Decomp{i} = inverse(JacTL{i});
        %         end
        %         Constraints = cell(1,k);
        %         for i = 1:k
        %             Cons = Decomp{i}*JacTR{i};
        %             Constraints{i} = -C{i}*Cons;
        %         end
        %         J = [JacM;blkdiag(Constraints{:})];
        %         DeltaNewtonLambda = J\[-RJacM;Eterm];
        %         DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
        %         for i = 1:k
        %             Rside = JacTR{i}*DeltaNewtonLambda{i};
        %             S{i} = -Decomp{i}*Rside;
        %         end
        %         for i = 1:k
        %             S{i+k} = S{i+k}+DeltaNewtonLambda{i};
        %         end
        %         
        %         [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        %         EVector = cell2mat(S(1:k));
        %         RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
        %         vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
        %         NumberNewton = NumberNewton + 1;
        %     end

%             if NumberNewton == 7 && norm(vNewton) >= 10^-9
%                 StopByNumberofNewton = [StopByNumberofNewton,1];
%             else
%                 StopByNumberofNewton = [StopByNumberofNewton,0];
%             end
    % Block method of the Newton step
    % for i = 1:k
    %     vn = null(JacTL{i}.');
            if NumberNewton <= 2 && h <= hmax/2
                h = h*2;
            else
                if NumberNewton == 8 && h >= hmin*2
                    h = h/2;
                end
            end
%         NumberofTotalOuterloop = NumberofTotalOuterloop + 1;
        end
    
    
    
    
    
    
    
        %     % one more step when t gets larger than 1
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
            % Naive way
            EVector = cell2mat(S(1:k));
            RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
            vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
            
            
            while NumberNewton < 430 && norm(DeltaEulorStep,inf) >= 10^-9  %EndpointN
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
        
        %
        %
        %
        %         structured method
        %         Decomp = cell(1,k);
        %         for i = 1:k
        %             Decomp{i} = inverse(JacTL{i});
        %         end
        %         Constraints = cell(1,k);
        %         for i = 1:k
        %             Cons = Decomp{i}*JacTR{i};
        %             Constraints{i} = -C{i}*Cons;
        %         end
        %         J = [JacM;blkdiag(Constraints{:})];
        %         DeltaNewtonLambda = J\[-RJacM;Eterm];
        %         DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
        %         for i = 1:k
        %             Rside = JacTR{i}*DeltaNewtonLambda{i};
        %             S{i} = -Decomp{i}*Rside;
        %         end
        %         for i = 1:k
        %             S{i+k} = S{i+k}+DeltaNewtonLambda{i};
        %         end
        %
        %         [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        %         EVector = cell2mat(S(1:k));
        %         RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
        %         vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
        %         NumberNewton = NumberNewton + 1;
        %     end
        
%         if NumberNewton == 7 && norm(vNewton) >= 10^-9
%             StopByNumberofNewton = [StopByNumberofNewton,1];
%         else
%             StopByNumberofNewton = [StopByNumberofNewton,0];
%         end
        end
    
    
    
    
    
%     hmax = 10^-7;
%     hmin = 10^-6;
%     if h > hmax
%         h = hmax;
%     end
%     while t < 1
%         vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];
% 
%         DeltaEulorStep = Jac\(-vEulor); % "-" or not?
% 
%         % update S and t
%         S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);
%         t = t + h;
% 
%         [JacTL,JacTR] = EvaluateHi(k,A,S,n);
%         JacM = (1-t)*JacL+t*JacG;
%         ML = (1-t)*Lc;
%         NumberNewton = 0;
% 
%         % Newton step
%         % Naive way
%         EVector = cell2mat(S(1:k));
%         RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
%         vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
%         while NumberNewton < 7 && norm(vNewton) >= 10^-9
%             JacH = [blkdiag(JacTL{:}),blkdiag(JacTR{:})];
%             Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
%             DeltaNewtonStep = Jac\(-vNewton);
%             S = mat2cell(cell2mat(S) + DeltaNewtonStep,[n,k*ones(1,k)],1);
% 
%             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
%             EVector = cell2mat(S(1:k));
%             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacM*cell2mat(S((k+1):2*k))-ML];
%             NumberNewton = NumberNewton + 1;
%         end
% 
% % 
% % 
% % 
% %         structured method
% %             Decomp = cell(1,k);
% %             for i = 1:k
% %                 Decomp{i} = inverse(JacTL{i});
% %             end
% %             Constraints = cell(1,k);
% %             for i = 1:k
% %                 Cons = Decomp{i}*JacTR{i};
% %                 Constraints{i} = -C{i}*Cons;
% %             end
% %             J = [JacM;blkdiag(Constraints{:})];
% %             DeltaNewtonLambda = J\[-RJacM;Eterm];
% %             DeltaNewtonLambda = mat2cell(DeltaNewtonLambda,k*ones(1,k),1);
% %             for i = 1:k
% %                 Rside = JacTR{i}*DeltaNewtonLambda{i};
% %                 S{i} = -Decomp{i}*Rside;
% %             end
% %             for i = 1:k
% %                 S{i+k} = S{i+k}+DeltaNewtonLambda{i};
% %             end
% %             
% %             [JacTL,JacTR] = EvaluateHi(k,A,S,n);
% %             EVector = cell2mat(S(1:k));
% %             RJacM = JacM*cell2mat(S((k+1):2*k))-ML;
% %             vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
% %             NumberNewton = NumberNewton + 1;
% %         end
% 
%         if NumberNewton == 7 && norm(vNewton) >= 10^-9
%             StopByNumberofNewton = [StopByNumberofNewton,1];
%         else
%             StopByNumberofNewton = [StopByNumberofNewton,0];
%         end
%     % Block method of the Newton step
%     % for i = 1:k
%     %     vn = null(JacTL{i}.');
%         if NumberNewton <= 2 && h <= hmax/2
%             h = h*2;
%         else
%             if NumberNewton == 7 && h >= hmin*2
%                 h = h/2;
%             end
%         end
%         NumberofTotalOuterloop = NumberofTotalOuterloop + 1;
%     end
        TimeEachPath(index) = toc;
        for vg = 1:k
            EigenVector{index,vg} = S{vg};
            EigenValue{index,vg} = S{k+vg};
        end
        index = index + 1;
        index
    end  
end





    
    
    
% [JacTL,~] = EvaluateHi(k,A,S,n);
% EVector = cell2mat(S(1:k));
% vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;JacG*cell2mat(S((k+1):2*k))];
% 
% residual = cell(k+1,k);
% for i = 1:k
%     for j = 1:k
%         Residual = A{j,1};
%         for q = 1:k
%         	Residual = Residual - A{j,q+1}*S{k+i}(q);
%         end
%         residual{j,i} = Residual*S{j};
%     end
%     residual{k+1,i} = C{i}*S{i}-1;
% end
%         
% % Error = cell(k,1);
% % [JacTL,JacTR] = EvaluateHi(k,A,S,n,C);
% % for i = 1:k
% %     Error{i} = JacTL{i}*S{i};
% % end
%     
% 
% error = norm(vNewton);
end
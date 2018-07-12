function [EigenValue,EigenVector,TimeEachPath,NewtonIteration,c] = HomotopyPEP(A)
% the input A should be a k by (k+1) cell array which contains all the
% matrix coefficients

d = size(A,2)-1;  % to get the degree of the PEP problem
n = size(A{1},1); % to get the dimension of the matrices

% linear constraints on the eigenvector, i.e, c
c = randn(1,n) + 1i*randn(1,n);

% generate the random start system
% D = cell(1,d+1);
% for i = 1:(d+1)
%     D{i} = diag(randn(1,n) + 1i*randn(1,n));
% end


D = randn(n,d+1) + 1i*randn(n,d+1);

% NormA = 1:(d+1);
% for i = 1:(d+1)
%     NormA(i) = norm(A{i},'fro');
% end
% D = D*diag(NormA)*exp(1i*rand*2*pi);

Dcell = cell(1,d+1);
for i = 1:(d+1)
   Dcell{i} = diag(D(:,i));
end

% get the solutions for the start system
SEigenVector = diag(1./c);
SEigenValue = cell(n,1);

for i = 1:n
   SEigenValue{i} = roots(D(i,:));
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

% loop over all possible paths
% loop = cell(k,1);
% for i = 1:k
%     loop{i} = 1:n(i);
% end
% Loop = allcomb(loop{:})';

hmax = 10^-2;
hmin = 10^-6;

% Eterm = ones(k,1);
% JD = -JacL+JacG;
% 
% pn = prod(n);
EigenValue = 1:(n*d);
EigenVector = cell(1,n*d);
% 
TimeEachPath = 1:(n*d);
% LastNewton = 1:pn;
index = 0;
NewtonIteration = 0;

for i = 1:n
    for j = 1:d
        EVT = SEigenVector(:,i);
        EVL = SEigenValue{i}(j);

        t = 0;
        h = hmax;

        [Dlambda,Plambda,TR] = JacPEP(EVT,EVL,t,d,n,Dcell,A);
        Jac = [(1-t)*Dlambda+t*Plambda,TR;c,0];

        tic;
        while t < 1 %-hmax
            vEulor = [(Plambda-Dlambda)*EVT;0];

            DeltaEulorStep = Jac\(-vEulor); 

            % update S and t
            EVT = EVT + h*DeltaEulorStep(1:n);
            EVL = EVL + h*DeltaEulorStep(n+1);
            t = t + h;

            [Dlambda,Plambda,TR] = JacPEP(EVT,EVL,t,d,n,Dcell,A);
            NumberNewton = 1;
            NewtonIteration = NewtonIteration + 1;

            % Newton step
            % Naive way
            vNewton = [((1-t)*Dlambda+t*Plambda)*EVT;c*EVT-1];
            Jac = [(1-t)*Dlambda+t*Plambda,TR;c,0];
            DeltaNewtonStep = Jac\(-vNewton);
            EVT = EVT + DeltaNewtonStep(1:n);
            EVL = EVL + DeltaNewtonStep(n+1);

            [Dlambda,Plambda,TR] = JacPEP(EVT,EVL,t,d,n,Dcell,A);
            vNewton = [((1-t)*Dlambda+t*Plambda)*EVT;c*EVT-1];


            while NumberNewton < 40 && norm(DeltaNewtonStep,inf) >= 10^-8
            	Jac = [(1-t)*Dlambda+t*Plambda,TR;c,0];
                DeltaNewtonStep = Jac\(-vNewton);
                EVT = EVT + DeltaNewtonStep(1:n);
                EVL = EVL + DeltaNewtonStep(n+1);

                [Dlambda,Plambda,TR] = JacPEP(EVT,EVL,t,d,n,Dcell,A);
                vNewton = [((1-t)*Dlambda+t*Plambda)*EVT;c*EVT-1];
                NumberNewton = NumberNewton + 1;
                NewtonIteration = NewtonIteration + 1;
            end


            if NumberNewton <= 2 && h <= hmax/2
                h = h*2;
            else
                if NumberNewton == 40 && h >= hmin*2
                    h = h/2;
                end
            end
%         NumberofTotalOuterloop = NumberofTotalOuterloop + 1;
        end







    %     % one more step when t gets larger than 1
%     if t < 1 
%         JacM = (1-t)*JacL+t*JacG;
%         Jac = [JacH;CBlock;zeros(ik,sumn),JacM];
%         vEulor = [zeros(sumn+k,1);JD*cell2mat(S((k+1):(2*k)))+Lc];
% 
%         DeltaEulorStep = Jac\(-vEulor); % "-" or not?
% 
%         % update S and t
%         h = 1 - t;
%         S = mat2cell(cell2mat(S) + h*DeltaEulorStep,[n,k*ones(1,k)],1);
% 
%         [JacTL,JacTR] = EvaluateHi(k,A,S,n);
%         JacM = JacG;
%         ML = 0;
%         NumberNewton = 1;
% 
%         % Newton step
%         % Naive way
%         EVector = cell2mat(S(1:k));
%         RJacM = JacM*cell2mat(S((k+1):(2*k)))-ML;
%         vNewton = [blkdiag(JacTL{:})*EVector;CBlockDiag*EVector-Eterm;RJacM];
% 
% 
%         while NumberNewton < 90 && norm(DeltaEulorStep,inf) >= 10^-9  %EndpointN
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
%         LastNewton(index) = NumberNewton;

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
    TimeEachPath(index+1) = toc;
    EigenVector{index+1} = EVT;
    EigenValue(index+1) = EVL;
    index = index + 1;
    index
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
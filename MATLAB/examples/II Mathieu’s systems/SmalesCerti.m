function [alpha,f,mu,vNewton] = SmalesCerti(EigenValue,EigenVector,A,C,JacG)
    k = size(A,1);
    n = 1:k;
    for i = 1:k
        n(i) = size(A{i,1},1);
    end
    prodn = prod(n);
    sumn = sum(n);

    alpha = 1:prodn;
    CBlock = [blkdiag(C{:}),zeros(k,k^2)];
    CBlockDiag = blkdiag(C{:});
    ik = (k-1)*k;
    Eterm = ones(k,1);

    f = norm(cell2mat(C))^2+norm(JacG,'fro')^2;
    for j = 1:k
        f = f + (norm(cell2mat(A(j,:)),'fro')^2)/2;
    end
    f = sqrt(f);
    
%     S = cell(2*k,1);
%     for i = 1:k
%         S{i} = SEigenVector{i}(:,path(i));
%         S{k+i} = nCell{i}*SEigenValue{i}(path(i))+sCell{i};
%     end

    for i = 1:prodn
        EVector = EigenVector(i,:);
        EVector = EVector';
        S = [EVector;EigenValue(i,1)];
        for s = 2:k
            S = [S;EigenValue(i,s)];
        end
        [JacTL,JacTR] = EvaluateHi(k,A,S,n);
        Df = [blkdiag(JacTL{:}),blkdiag(JacTR{:});CBlock;zeros(ik,sumn),JacG];
        vNewton = [blkdiag(JacTL{:})*cell2mat(EVector);CBlockDiag*cell2mat(EVector)-Eterm;JacG*cell2mat(S((k+1):2*k))];
        
        norm_z = norm(cell2mat(S));
        Delta_d = diag([sqrt(2*(1+norm_z^2))*ones(1,sumn),ones(1,k^2)]);
        mu = norm(Df\Delta_d)*f;
        gamma = mu*sqrt(2/(1+norm_z^2));
        alpha(i) = gamma*norm(Df\vNewton);
    end
end
    
   
    
    
function writeBertiniInput(Str1,Str2,A,L,G,varargin)
% Str1: directory, for example 'C:/MyFolder/'
% Str2: name of the text file, for example 'file.txt'
% A: coefficient matrices
% L: a list of k (k-1)-by-k matrices (L1,L2,...,Lk), constraints on lambdas
% so that Li(lambdai) = 1
% G: a list of (k-1) k-by-k matrices (G1,G2,...,G(k-1)), constraints on
% lambda(i+1) - lambda(i)

id = fullfile(Str1, Str2);
fid = fopen(id,'wt');
fprintf(fid, ['%% This input file was written with the Bertini.m2 Macaulay2 package.' '\n\n']);
fprintf(fid, ['CONFIG' '\n\n']);
if nargin == 5
    fprintf(fid, ['ParameterHomotopy : 2 ;' '\n']);
    fprintf(fid, ['PrintPathProgress : 100 ;' '\n\n']);
end
fprintf(fid, ['%%%%%%ENDCONFIG;' '\n']);
fprintf(fid, ['END;' '\n\n']);
fprintf(fid, ['INPUT;' '\n\n']);

k = size(A,1);
Dim = [];
for i = 1:k
    dim = size(A{i,1},1);
    Dim = [Dim,dim];
    fprintf(fid, ['hom_variable_group x' num2str(i) 'v1']);
    for j = 2:dim
        fprintf(fid, [', ' 'x' num2str(i) 'v' num2str(j)]);
    end
    fprintf(fid, [ ' ' ';' '\n\n']);
end
fprintf(fid, ['variable_group' ' ' 'lam1v1' ]);
for i = 1:k
    if i == 1
        for j = 2:k
            fprintf(fid, [ ', ' 'lam1v' num2str(j) ]);
        end
    else
        for j = 1:k
            fprintf(fid, [ ', ' 'lam' num2str(i) 'v' num2str(j) ]);
        end
    end
end
fprintf(fid, [ ' ' ';' '\n\n\n\n']);
fprintf(fid, [ 'parameter fpT ;' '\n\n\n\n']);
fprintf(fid, 'constant ii');
for i = 1:k
    fprintf(fid, [ ', ' 'lam' num2str(i) 'v0' ]);
end
fprintf(fid, [ ' ;' '\n']);
fprintf(fid, [ 'ii = I; ' '\n']);
for i = 1:k
    fprintf(fid, [ 'lam' num2str(i) 'v0 = 1 ;' '\n']);
end
fprintf(fid, '\n');
fprintf(fid, 'function jade0');
s = k + sum(Dim);
for i = 1:(s-1)
    fprintf(fid, [', ' 'jade' num2str(i)]);
end
fprintf(fid, [' ;' '\n\n']);
for i = 1:k
    for t = 1:k-1
        fprintf(fid, ['linearL' num2str(t+(i-1)*(k-1)) ' = ((lam' num2str(i) 'v1)*(' num2str(L{i}(t,1),'%10.15e') '+' num2str(L{i}(t,1),'%10.15e') '*ii)']);
        for j = 2:k
            fprintf(fid, ['+(lam' num2str(i) 'v' num2str(j) ')*(' num2str(L{i}(t,j),'%10.15e') '+' num2str(L{i}(t,j),'%10.15e') '*ii)']);
        end
        fprintf(fid, [')-1 ;' '\n\n']);
    end
end
for i = 1:k-1
    for t = 1:k
        fprintf(fid, ['linearG' num2str(t+(i-1)*k) '= ((' num2str(G{i}(t,1),'%10.15e') '+' num2str(G{i}(t,1),'%10.15e') '*ii)*(lam' num2str(i) 'v1-lam' num2str(i+1) 'v1)' ]);
        for j = 2:k
            fprintf(fid, ['+(' num2str(G{i}(t,j),'%10.15e') '+' num2str(G{i}(t,j),'%10.15e') '*ii)*(lam' num2str(i) 'v' num2str(j) '-lam' num2str(i+1) 'v' num2str(j) ')']);
        end
        fprintf(fid, [') ;' '\n\n']);
    end
end
fprintf(fid, '\n');
for i = 1:k*(k-1)
    fprintf(fid, ['jade' num2str(i-1) ' = (1-fpT)*linearL' num2str(i) '+(fpT)*linearG' num2str(i) ' ; ' '\n\n']);
end
Cumsum = cumsum(Dim);
Cumsum(end) = [];
Cumsum = [0,Cumsum];
for i = 1:k
    dimi = Dim(i);
    for t = 1:dimi
        fprintf(fid, ['jade' num2str(k*(k-1)+Cumsum(i)+t-1) '= (x' num2str(i-1) 'v1)*((' num2str(real(A{i,1}(t,1)),'%10.15e') '+' num2str(imag(A{i,1}(t,1)),'%10.15e') '*ii)*lam' num2str(i) 'v0']);
        for j = 2:k+1
            fprintf(fid, ['+(' num2str(real(A{i,j}(t,1)),'%10.15e') '+' num2str(imag(A{i,j}(t,1)),'%10.15e') '*ii)*lam' num2str(i) 'v' num2str(j-1)]);
        end
        fprintf(fid, ')');
        for p = 2:dimi
            fprintf(fid, ['+(x' num2str(i-1) 'v' num2str(p) ')*((' num2str(real(A{i,1}(t,p)),'%10.15e') '+' num2str(imag(A{i,1}(t,p)),'%10.15e') '*ii)*lam' num2str(i) 'v0']);
            for j = 2:k+1
                fprintf(fid, ['+(' num2str(real(A{i,j}(t,p)),'%10.15e') '+' num2str(imag(A{i,j}(t,p)),'%10.15e') '*ii)*lam' num2str(i) 'v' num2str(j-1)]);
            end
            fprintf(fid, ')');
        end
        fprintf(fid, [' ;' '\n\n']);
    end
end
fprintf(fid, '\n');
fprintf(fid, 'END;');
fclose(fid);
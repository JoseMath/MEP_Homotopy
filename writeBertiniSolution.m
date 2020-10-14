function writeBertiniSolution(Str1,Str2,Csol)
% Str1: directory, for example 'C:/MyFolder/'
% Str2: name of the text file, for example 'file.txt'
% Csol: solution matrix, where columns are solutions (x1,x2,...,xk,lambda1,...,lambdak)
% List: extrinsic dimensions of the MEP (size of polynomial matrices [n1,...,nk])
id = fullfile(Str1, Str2);
fid = fopen(id,'wt');
NumberofSol = size(Csol,2);
NumberofRow = size(Csol,1);
strNumberofSol = num2str(NumberofSol);
fprintf(fid, [ strNumberofSol '\n\n']);

for i = 1:NumberofSol
    for j = 1:NumberofRow
        fprintf(fid, [num2str(real(Csol(j,i)),'%10.15e') ' ' num2str(imag(Csol(j,i)),'%10.15e') '\n']);
    end
    fprintf(fid, '\n');
end
fclose(fid);
end

% fid = fopen( '/Users/yilingyou/Downloads/MatlabOutput.txt', 'wt');
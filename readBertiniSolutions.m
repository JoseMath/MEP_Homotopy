function Csol = readBertiniSolutions(Str1,Str2,List)
% Str1: directory, for example 'C:/MyFolder/'
% Str2: name of the text file, for example 'file.txt'
% List: extrinsic dimensions of the MEP (size of polynomial matrices [n1,...,nk])

% return
% Csol: columns are solutions (x1,x2,...,xk,lambda1,...,lambdak)
k = length(List);
fileID = fopen(strcat(Str1,Str2));
C = textscan(fileID,'%f');

C = C{1};
C(1) = [];
CRe = C(1:2:end);
CIm = C(2:2:end);
Csol = CRe + 1j*CIm;

N = sum(List) + k*2;
Csol = reshape(Csol,N,[]);
end
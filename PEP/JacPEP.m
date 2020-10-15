function [Dlambda,Plambda,TR] = JacPEP(EVT,EVL,t,d,n,Dcell,A)
Dlambda = Dcell{d+1};
Plambda = A{d+1};
Dprime = zeros(n);
Pprime = zeros(n);
for i = 1:d
    Dlambda = Dlambda + EVL^(d+1-i)*Dcell{i};
    Plambda = Plambda + EVL^(d+1-i)*A{i};
    Dprime = Dprime + (d+1-i)*EVL^(d-i)*Dcell{i};
    Pprime = Pprime + (d+1-i)*EVL^(d-i)*A{i};
end
TR = ((1-t)*Dprime+t*Pprime)*EVT;
end
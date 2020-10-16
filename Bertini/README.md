# Write Bertini File From MATLAB.

```MATLAB
    Str1='/Users/jo/Desktop/Dump/'
    Str2='input_target'
	k = 3;	ni = 3;	n = ni*ones(1,k);    m=2*ones(1,k);
	Amatrix = rands(k*ni,(k+1)*ni)+1i*randn(k*ni,(k+1)*ni);
	A = mat2cell(Amatrix,n,ni*ones(1,k+1));
    L1=rands(k-1,k);    L2=rands(k-1,k);    L3=rands(k-1,k);
    L={L1,L2,L3}
    G1=rands(k,k);
    G2=rands(k,k);
    G={G1,G2}
    writeBertiniInput(Str1,Str2,A,L,G)
```
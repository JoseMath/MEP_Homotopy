k=2 % Number of parameters
theN=[3,3] %Dimension of MEP
assert(length(theN)==k)  
%%
A=cell(k,k+1)
for i = 1:k
  for j=1:k+1
    A{i,j}=randn(theN(i),theN(i));
  end    
end
A
%%
FiberHomotopy2(A)
%DEMO_MODEL_UPDATE   demo model updating problem with two parameters from Cottin
%
% For a given matrices A, B, and C of the same size and values s1 and s2
% we are looking for lambda and mu such that A - lambda*B - mu*C has 
% eigenvalues s1 and s2 (we do not care about the remaining eigenvalues)
% 
% We can write this as a singular two-parameter eigenvalue problem
% 
% (A - s1 I) x = lambda B x + mu C x 
% (A - s2 I) y = lambda B y + mu C y 
%
% In the generic case there are n(n-1)/2 solutions.
%
% References:
%  - N. Cottin. Dynamic model updating � a multiparameter eigenvalue problem, 
%    Mech. Syst. Signal Pr. 15 (2001) 649�665.
%  - A. Muhic, B. Plestenjak, On the singular two-parameter eigenvalue problem, 
%    Electron. J. Linear Algebra 18 (2009) 420�437.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(1);
end

n = 50; % computation is feasible for approx. n < 50

A = rand(n);
B = rand(n);
C = rand(n);
I = eye(n);

% prescribed eigenvalues are 2 and 3
s1 = 2;
s2 = 3;

%% Solve the singular two-parameter eigenvalue problem by MultiParEig package
tic
opts.singular = 1; 
[lambda,mu] = twopareig(A-s1*I,B,C,A-s2*I,B,C,opts);
disp(toc)

A = {A-s1*I,B,C;A-s2*I,B,C};
k = 2;

% calculate the backward error
EigenValue = [lambda,mu]; EigenValue = mat2cell(EigenValue, ones(length(lambda),1), k);        
EigenVector = [X1',X2']; EigenVector = mat2cell(EigenVector, ones(length(lambda),1), ones(k,1)*n); EigenVector=cellfun(@transpose,EigenVector,'UniformOutput',false);
[MaxRelError_delta,RelError_delta,Position_delta,ri_delta,thetai_delta] = RBackwardError(EigenValue,EigenVector,A);  


%% Solve the singular two-parameter eigenvalue problem by fiber product homotopy method
tic
[EigenValue,EigenVector,SEigenValue,SEigenVector,TimeEachPath,NumNewtonEachPath,LastNewton,JacG,C] = FiberHomotopy(A, struct('tracker', 'conservative'));
disp(toc)
% calculate the backward error
[MaxRelError_fiber,RelError_fiber,Position_fiber,ri_fiber,thetai_fiber] = RBackwardError(EigenValue,EigenVector,A);
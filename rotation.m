function [eval, evec] = rotation(eta, D, D2, W0, k)

A = Amat_incomp(eta, D, D2, W0, k);
B = Bmat_incomp(eta, D, D2, k);

% solve generalized EV problem
[evec,eval] = eig(A,B);
eval = diag(eval);

clear A B


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BUILD MATRIX A 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = Amat_incomp(eta, D, D2, W0, k)

n = length(eta);
A = zeros(n);   % 5 columns 
clear i

dW0=D*W0;
d2W0 =D2*W0; 

%%%%%%%%%% GOVERNING EQUATIONS %%%%%%%%%%

% psi
A(1:n, 1:n) = -k*diag(W0)*(D2 - k^2*eye(n)) - k*(diag(d2W0) + 0.2*eye(n));


%%%%%%%%%% BOUNDARY CONDITIONS %%%%%%%%%%
%Dirichlet
A(1,:)=0;
A(n,:)=0;
A(1,1)=1;
A(n,n)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BUILD MATRIX B 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = Bmat_incomp(eta, D, D2, k)

n = length(eta);
B = zeros(n);
clear i

%%%%%%%%%% GOVERNING EQUATIONS %%%%%%%%%%

% psi
B(1:n, 1:n) = (D2 - k^2*eye(n));

%%%%% BOUNDARY CONDITIONS %%%%%%
% B must be zero wherever we have a Neumann or Dirichlet BC in A
B(1,:)=0;
B(n,:)=0;




function [M_est U_est V_est obj] = l2_rpca_mask_alm_fast(M,W,r,lambda,U,V,maxIterIN, rho)
% Ricardo Silveira Cabral
% Human Sensing Lab, Carnegie Mellon University
% If you use this code in your research, please be so kind as to cite:
% Unifying Nuclear Norm and Bilinear Factorization Approaches for Low-Rank Matrix Decomposition
% Ricardo S. Cabral, Fernando De la Torre, Joao P. Costeira, Alexandre Bernardino
% ICCV 2013


%% Robust low-rank matrix approximation with missing data and outliers
% min |W.*(M-E)|_1 + lambda*|V|_*
% s.t., E = UV, U'*U = I
%
%Input:
%   M: m*n data matrix
%   W: m*n indicator matrix, with '1' means 'observed', and '0' 'missing'.
%   r: the rank of r
%   lambda: the weighting factor of the trace-norm regularization, 1e-3 in default.
%   rho: increasing ratio of the penalty parameter mu, usually 1.05.
%   maxInterIN: maximum iteration number of inner loop, usually 100.
%   signM: if M>= 0, then signM = 1, otherwise, signM = 0;
%Output:
%   M_est: m*n full matrix, such that M_est = U_est*V_est
%   U_est: m*r matrix
%   V_est: r*n matrix
%   L1_error: the L1-norm error of observed data only.

%% In-default parameters
[m n] = size(M); %matrix dimension
if nargin < 7
  rho = 1.05;
  
  maxIterIN = 5000;
end

if nargin < 5
  U = randn(m,r);
  V = randn(n,r);
  lambda = 1e-3;
end
if nargin < 3
  disp('Please input the data matrix M, the indicator W and the rank r, and try again.');
end

% works 5000 1e-6 1.05
% works 5000 1e-5 1.05
maxIterOUT = 2000;
max_mu = 1e20;
mu = 1e-3;%mu = 1e-5; % GOOD FOR SFM
M_norm = norm(M,'fro');
tol = 1e-9*M_norm;
%tol = 1e-10
cW = ones(size(W)) - W; %the complement of W.
display = 1; %display progress
%% Initializing optimization variables as zeros
E = randn(m,n);
Y = zeros(m,n); %lagrange multiplier

%% caching
lr = lambda*eye(r);
%% Start main outer loop
iter_OUT = 0;
while iter_OUT < maxIterOUT
  iter_OUT = iter_OUT + 1;
  
  itr_IN = 0;
  obj_pre = 1e20;
  %start inner loop
  while itr_IN < maxIterIN
    %update U
    tmp = mu*E + Y;
    U = (tmp)*V/(lr + mu*(V'*V));
    
    %update V
    V = (tmp)'*U/(lr + mu*(U'*U));
    
    %update E
    UV = U*V';
    temp1 = UV - Y/mu;
    E = 1/(mu+2)*(2*M+mu*temp1).*W + temp1.*cW;
    
    %evaluate current objective
    leq = E - UV;
    obj_cur = sum(sum(abs(W.*(M-E)).^2)) + lambda/2*(norm(U,'fro')^2 + norm(V,'fro')^2) + sum(sum(Y.*leq)) + mu/2*norm(leq,'fro')^2;
    
    %check convergence of inner loop
    if abs(obj_cur - obj_pre) <= 1e-12*abs(obj_pre)
      break;
    else
      obj_pre = obj_cur;
      itr_IN = itr_IN + 1;
    end
  end
 % itr_IN
  stopC = norm(leq,'fro');
  if display && (iter_OUT==1 || mod(iter_OUT,50)==0 || stopC<tol)
    obj = sum(sum(abs(W.*(M-U*V')).^2)) + lambda/2*(norm(U,'fro')^2 + norm(V,'fro')^2);
    
    disp(['iter ' num2str(iter_OUT) ',mu=' num2str(mu,'%2.1e') ...
      ',obj=' num2str(obj) ',stopALM=' num2str(stopC,'%2.3e')]);
  end
  if stopC<tol
    break;
  else
    %update lagrage multiplier
    Y = Y + mu*leq;
    %update penalty parameter
    mu = min(max_mu,mu*rho);
  end
end

%% Denormalization
U_est = U;
V_est = V;

%U_est = sqrt(scale)*U; V_est = sqrt(scale)*V;
%L1_error = sum(sum(abs(W.*(scale*M-M_est))));
M_est = U_est*V_est';
obj = sum(sum(abs(W.*(M-M_est)).^2)) + lambda/2*(norm(U,'fro') + norm(V,'fro'));

end

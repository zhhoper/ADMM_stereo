% written by debingzhang
% if you have any questions, please fell free to contact
% debingzhangchina@gmail.com

% version 1.0, 2012.11.16




function [ X,iterations] = admmAXB( A,B,X,M,known,rho)

% my_admm
% solve  minimize ||X||_*-trace(A*X*B')

MAX_ITER = 200;

[m,n] = size(M);
%r = size(A,1);

U = X;
V = zeros(m,n);

%eigX = svd(X);

AB = A'*B;
%BA = B'*A;
eeppss = 0.0001;
%localpsnr = [];
for k = 1:MAX_ITER
    % X-update
    tem = U - V/rho;    lastX = X;
    [u,sigma,v] = svd(tem);
    X = u*max(sigma-1/rho,0)*v';
    if(norm(X-lastX,'fro')/norm(M,'fro')<eeppss)
        break;
    end
    
    %U-update
    lastU = U;
    U = (V+AB+rho*X)/rho;
    U = M.*known + U.*(ones(size(U))-known);
    if(norm(U-lastU,'fro')/norm(M,'fro')<eeppss)
        break;
    end
    
    %V-update
    if(norm(X-U,'fro')/norm(M,'fro')<eeppss)
        break;
    end
    V = V+rho*(X-U);
    %localpsnr = [localpsnr PSNR(Xfull,X,inmissing)]; 
    
end

iterations = k;

end

%function obj = objective(X,A,B,BA)
%    sigma = svd(X);
%    obj = sum(sigma) - trace(X*BA);
%end


function [XX,iterations] = my_apgl(A,B,X,M,known,eps,lambda)%,Xfull,inmissing)
    % min  g(x)+h(x)    g(x) = -trace(AXB)+\lambda*||X_\Omega - M_\Omega||_F^2, h(x) = ||X||_*
    if (~exist('eps','var'))
        eps = 0.1;
    end
    if (~exist('lambda','var'))
        lambda = 0.00000000000000000001;
    end

    AB = A'*B;BA = B'*A;

    lastX = X;
    Y = lastX;

    tlast = 1;t = 1;
    XX = X;
    %eeppss = 1e-4;
    %localpsnr = [];
    for k=1:200
        %k
        %whos Y t AB M known
        [u,sigma,v] = svd(Y+t*(AB-lambda*(Y-M).*known));

        X = u*max(sigma-t,0)*v';

        t = (1 + sqrt(1+4*tlast*tlast))/2;

        Y = X + (tlast-1)/(t)*(X - lastX);
        
        XX = X;
        lastX = X;
        tlast = t;
%         if(norm(lastX-X,'fro')/norm(M,'fro'))<eeppss
%             %if(k>=10)
%             break;
%             %end
%         end
       %localpsnr = [localpsnr PSNR(Xfull,X,inmissing)];
       
       history.len = k;
       history.t = t;
       history.objval(k) = trace_norm(X) - trace(X*BA) + lambda/2*norm((X-M).*known,'fro')^2;
      % history.objval(k);
     % history.realval(k) = real_norm(X,size(A,1));
        if(k>=2 && -(history.objval(k)-history.objval(k-1))<eps)
            %if(k>=10)
                break;
            %end
        end
        
        
    end
    iterations = k;
end

function obj = trace_norm(X)
    if(size(X,1)>2000)
        sigma = svds(X,100);
    else
    sigma = svd(X);
    end
    obj = sum(sigma);
end
% function obj = real_norm(X,r)
%     sigma = svd(X);
%     obj = sum(sigma(r+1,end));
% end
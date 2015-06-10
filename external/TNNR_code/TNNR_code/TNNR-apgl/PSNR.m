function [ psnr ] = PSNR( Xfull,Xrecover,missing )
%PSNR Summary of this function goes here
%   Detailed explanation goes here
%     if(size(Xfull) == size(Xrecover))
%     end
    Xrecover = max(0,Xrecover);
    Xrecover = min(255,Xrecover);
    [m,n,dim] = size(Xrecover);
    MSE = 0;
    for i =1 : dim
        MSE = MSE + norm((Xfull(:,:,i)-Xrecover(:,:,i)).*missing(:,:,i),'fro')^2;
    end
    MSE = MSE/nnz(missing);
    
    psnr = 10*log10(255^2/MSE);
        
        
end


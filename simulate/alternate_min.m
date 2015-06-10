function [Lnew,Snew]=alternate_min(M,L0)

maxit=50;
L=L0;
for i=1:maxit
    
    Snew=(L'*L)\(L'*M);
    obj(1)=0.5*norm(M-L*Snew,'fro')^2;
    Lnew=(M*Snew')/(Snew*Snew');
    L=Lnew;
    obj(i+1)=0.5*norm(M-Lnew*Snew,'fro')^2;
end
figure; plot(obj);
end
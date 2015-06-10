function [G, ZNew] = getGBR(ZGT, ZComplete)
% returns the GBR transformation matrix between depth fields that transforms M into MGT

[m, n] = size(ZGT);


A = [createCoordArray(1:m,1:n),ZComplete(:), ones(m*n,1)];

B = ZGT(:);

res = A \ B;
gbr = res(1:3);

if gbr(3,1) == 0
    gbr = [0, 0, 1]';
end

constant = res(4);

G = eye(3,3);
G(3,:) = gbr';

ZNew = constant + reshape(A(:,1:3)*gbr, m,n);


end

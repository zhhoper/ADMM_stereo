function [G, ZNew] = getGBR(ZGT, ZComplete, mask)
% returns the GBR transformation matrix between depth fields that transforms M into MGT

[m, n] = size(ZGT);

tmpCoordinate = createCoordArray(1:m,1:n);
tmpCoordinate = tmpCoordinate(mask(:), :);
num = sum(sum(mask));

A = [tmpCoordinate, ZComplete(:), ones(num,1)];

B = ZGT(:);
B = B(mask(:));

res = A \ B;
gbr = res(1:3);

if gbr(3,1) == 0
    gbr = [0, 0, 1]';
end

constant = res(4);

G = eye(3,3);
G(3,:) = gbr';

tmpA = A(:, 1:3)*gbr;
tmpA = vec2mat_mask(tmpA, mask);
ZNew = constant + reshape(tmpA, m,n);


end

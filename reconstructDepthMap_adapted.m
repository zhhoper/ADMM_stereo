function [Z, C, D] = reconstructDepthMap_adapted(S, imgsize)
addpath('simulate');
rows = imgsize(1);
cols = imgsize(2);

pixels = rows*cols;

[Sx, Sy] = getGradientField(S,imgsize);
Sx = Sx(1:end-1,:);
Sy = Sy(:, 1:end-1);

% x direction
[indx, indy] = meshgrid(1:rows-1, 1:cols);
C1 = zeros((rows-1)*cols, pixels);
indx = indx';
indy = indy';
tind1 = sub2ind(imgsize, indx(:), indy(:));
tind2 = sub2ind(imgsize, indx(:)+1, indy(:));
tind3 = sub2ind(size(C1), (1:(rows-1)*cols)', tind1);
tind4 = sub2ind(size(C1), (1:(rows-1)*cols)', tind2);
C1(tind3) = -1;
C1(tind4) = 1;
tSx = Sx(~isnan(Sx));
D1 = tSx(:);


% y direction
[indx, indy] = meshgrid(1:rows, 1:cols-1);
C2 = zeros(rows*(cols-1), pixels);
indx = indx';
indy = indy';
tind1 = sub2ind(imgsize, indx(:), indy(:));
tind2 = sub2ind(imgsize, indx(:), indy(:)+1);
tind3 = sub2ind(size(C2), (1:rows*(cols-1))', tind1);
tind4 = sub2ind(size(C2), (1:rows*(cols-1))', tind2);
C2(tind3) = -1;
C2(tind4) = 1;
tSy = Sy(~isnan(Sy));
D2 = tSy(:);
C = [C1;C2];
D = [D1;D2];
C(end+1,1) = 1;
D(end+1) = 0;
Z = C\D;

% Solve Least squares with regularization constant lambda
%             lambda = max(eig(C' * C)) * 1e-3;
%             Z = inv(C' * C + lambda * eye(size(C,2))) * C' * D;

end

function [Sx,Sy] = getGradientField(S,imgsize)
S = normc(S);
Sx = -reshape(S(1,:) ./ S(3,:), imgsize);
Sy = -reshape(S(2,:) ./ S(3,:), imgsize);

end
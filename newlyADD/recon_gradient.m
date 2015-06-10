function Z = recon_gradient(gradx, grady, mask, imgsize)

rows = imgsize(1);
cols = imgsize(2);

[xEdge, yEdge] = findEdges(mask);
validX = logical(mask.*(~xEdge));
validY = logical(mask.*(~yEdge));
validX = validX(1:end-1, :);
validY = validY(:, 1:end-1);
indX = validX(:);
indY = validY(:);

pixels = rows*cols;
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
C1 = C1(indX,:);

sx = gradx(1:end-1,:);
sy = grady(:, 1:end-1);

D1 = sx(:);
D1 = D1(indX);

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
C2 = C2(indY,:);

D2 = sy(:);
D2 = D2(indY);

C = [C1;C2];
D = [D1;D2];
C(end+1,1) = 1;
D(end+1) = 0;

Z = C\D;
end
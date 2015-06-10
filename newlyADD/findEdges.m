function [maskEdgeX, maskEdgeY] = findEdges(mask)
% maskEdge = findEdges(mask)
% 
% This function is used to find the edges based on mask
% Definition of edge: x: if current pixel is valid, but x+1 is not valid
% y: if current pixel is valid, but y+1 is not

[row, col] = size(mask);
% for x, first shift the image along x axis
tmpMask = zeros(row, col);
for i = 1 : row-1
    tmpMask(i,:) = mask(i+1,:);
end

maskEdgeX = xor(mask, tmpMask).*mask;

% for y, first shift the image along y axis
tmpMask = zeros(row, col);
for i = 1 : col-1
    tmpMask(:,i) = mask(:, i+1);
end

maskEdgeY = xor(mask, tmpMask).*mask;
end
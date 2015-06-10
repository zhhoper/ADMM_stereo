function [h] = showNormalsOnly(S, size)

rows=size(1);
cols = size(2);

% coords = createCoordArray(linspace(-2,2,rows),linspace(-2, 2,cols));
coords = createCoordArray(1:rows,1:cols);
X = coords(:,1);
Y = coords(:,2);
Z = ones(rows, cols);
h = figure('Name','First Rank 3 approx. using SVD');
S3 = 0.5* normc(S);
quiver3(X,Y,Z(:),S3(1,:)',S3(2,:)', S3(3,:)',0);
xlabel('x');
ylabel('y');
axis equal;
end

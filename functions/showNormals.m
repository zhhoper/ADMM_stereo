function [h] = showNormals(varargin)
S = varargin{1};
Z = varargin{2};

[rows, cols] = size(Z);
%[X,Y] = meshgrid(linspace(2,-2,cols),linspace(2, -2,rows));
%b = reshape(S', rows,cols,3); %%%% FIX this reversing!!!
%coords = createCoordArray(linspace(-rows/2,2,rows),linspace(-2, 2,cols));
coords = createCoordArray(1:rows,1:cols);
X = coords(:,1);
Y = coords(:,2);
if nargin==4
    visible = varargin{4};
else 
    visible = 'on';
end

if nargin==3
    title = varargin{3};
    h = figure('Name',title, 'visible', visible);
else 
    h = figure('visible',visible);
end
S = 0.5* my_normc(S);
quiver3(X,Y,Z(:),S(:,1),S(:,2), S(:,3),2);
xlabel('x');
ylabel('y');
hold on
surf(reshape(X,[rows,cols]),reshape(Y, [rows, cols]),Z, 'EdgeColor','none')%,'FaceColor', 'interp');
axis equal;
hold off;

end

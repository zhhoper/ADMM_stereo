function [ coords ] = createCoordArray( rows,cols )
%Creates a array of points that go down a whole row, column for column

[X, Y] = meshgrid(cols,rows);
coords = [Y(:),X(:)];



end


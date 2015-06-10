% test for findEdges

% construct mask x
x = zeros(5,5);
x(2,2:3) = 1;
x(3, 2:4) = 1;
x(4, 2:3) = 1;

[xEdge, yEdge] = findEdges(x);
ccc = 0;
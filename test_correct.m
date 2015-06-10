% test correct_position
[X, Y] = meshgrid(1:3,1:2);
X = X - 1;
Y = Y - 1;

[tx, ty] = correct_position(X, Y, 3, 2);
cc
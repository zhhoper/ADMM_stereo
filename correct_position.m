function [x,y] = correct_position(x,y,row, col)
% when using circulate boundary condition, this function can help correct
% positions
x(x < 1) = x(x < 1) + row;
y(y < 1) = y(y < 1) + col;
x(x > row) = x(x > row) - row;
y(y > col) = y(y > col) - col;
end
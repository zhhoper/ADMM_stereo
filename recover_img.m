function x = recover_img(dx, dy)
% recover image from gradient x and y

[row, col] = size(dx);

%% get the gradient
[Gx,~ ] = gradient_xy(dx);
[~, Gy] = gradient_xy(dy);

G = Gx + Gy;

%% get the matrix
[X, Y] = meshgrid(1:row, 1:col);
% center
center.x = X+1;
center.y = Y+1;
[center.x, center.y] = correct_position(center.x, center.y, row, col);
% upper
upper.x = X;
upper.y = Y+1;
[upper.x, upper.y] = correct_position(upper.x, upper.y, row, col);
% down
down.x = X+2;
down.y = Y+1;
[down.x, down.y] = correct_position(down.x, down.y, row, col);
% left
left.x = X+1;
left.y = Y;
[left.x, left.y] = correct_position(left.x, left.y, row, col);
% right
right.x = X+1;
right.y = Y+2;
[right.x, right.y] = correct_position(right.x, right.y, row, col);

totalNum = row*col;
A = sparse([], [], [], totalNum, totalNum, 5*totalNum);

ind_center = sub2ind([row, col], center.x, center.y);
A(sub2ind(size(A), 1:totalNum, ind_center(:)')) = -4;

ind_upper = sub2ind([row, col], upper.x, upper.y);
A(sub2ind(size(A), 1:totalNum, ind_upper(:)')) = 1;

ind_down = sub2ind([row, col], down.x, down.y);
A(sub2ind(size(A), 1:totalNum, ind_down(:)')) = 1;

ind_left = sub2ind([row, col], left.x, left.y);
A(sub2ind(size(A), 1:totalNum, ind_left(:)')) = 1;

ind_right = sub2ind([row, col], right.x, right.y);
A(sub2ind(size(A), 1:totalNum, ind_right(:)')) = 1;

%x = linsolve(A, G(:));
x = A\G(:);
x = reshape(x, [row, col]);

function [dx, dy] = gradient_xy(x)
% in this part, we do not use circular boundary condition but reflect the last
% row (column). So the last element in the gradient fields will always be 0
[row, col] = size(x);
tdx = zeros(row, col);
tdx(end,:) = x(end,:);
tdx(1:end-1,:) = x(2:end,:);
dx = tdx - x;

tdy = zeros(row, col);
tdy(:,end) = x(:,end);
tdy(:,1:end-1) = x(:,2:end);
dy = tdy - x;
end
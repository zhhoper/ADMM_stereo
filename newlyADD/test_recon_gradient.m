% test recon_gradient 
clc;
clear;
x = rand(5,5);

% construct mask 
mask = zeros(5,5);
mask(2,2:3) = 1;
mask(3, 2:4) = 1;
mask(4, 2:3) = 1;


[dx, dy] = gradient_xy(x);
dx = dx.*mask;
dy = dy.*mask;

z = recon_gradient(dx, dy, mask, [5,5]);
z = reshape(z, [5,5]);
[tdx, tdy] = gradient_xy(z);
tdx = tdx.*mask;
tdy = tdy.*mask;


clc;
clear all;

num = 40;
x = rand(num,num);
[gx, gy] = gradient_xy(x);
s = [-gx(:)';-gy(:)';ones(1,num*num)];
z = reconstructDepthMap_adapted(s, [num,num]);
tx = reshape(z,[num,num]);

[gxx, gyy] = gradient_xy(tx);
tx = gx - gxx;
ty = gy - gyy;
cccc  =0;
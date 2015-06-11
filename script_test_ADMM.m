tic
clear;  close all; 
addpath('solvers');
addpath('functions');
addpath('external');
addpath('simulate')
addpath('newlyADD');

% generate data
Sim = generate_data(0,1);

% initialzie data
[X, Y, lam, l3, Z, lambda, ini_SGRB] = initialize_ADMM(Sim);
% [~, ~, ~, ~, ~, ~, ini_SGRB] = initialize_ADMM(Sim);
% [X, Y, lam, l3, Z, lambda] = initialize_random_ADMM(Sim);
var = struct;
var.M = Sim.MNoise;
var.W = Sim.INC;
var.X = X;
var.Y = Y;
var.lam = lam;
var.l3 = l3;
var.Z = Z;
var.lambda = lambda;
[var.img_row, var.img_col] = size(Sim.Z);

% what should be the correct parameter?
var.c2 = 10;
var.tau = 1;

result = ADMM_opt(var);

numPixel = var.img_row*var.img_col;

depth = result.var.Z;
light = [result.var.X(3:end, 1:2), result.var.l3];
normal = [-result.var.X(1:2, 3:end)', ones(numPixel,1)];
%normal = [result.var.X(1:2, 3:end)', -ones(numPixel,1)];
albedo = -(sqrt(sum(normal.^2, 2)))./result.var.lam;
normal = my_normr(normal);

%depth=reconstructDepthMap_adapted(-normal',[var.img_row,var.img_col]);

[LGBR,ZGBR,SGBR,RGBR,G] = applyGBR(Sim.Z, Sim.R, depth, light, normal, albedo);

if (LGBR(:,1)'*Sim.L(:,1))<0
    LGBR=-LGBR; SGBR=-SGBR;
end

h = showLights(Sim.L, LGBR, light);
% showNormals(normal, reshape(depth,[var.img_row var.img_col]), 'Surface normals BEFORE GBR');
showNormals(SGBR, ZGBR, 'Surface normals AFTER GBR');
set(h, 'Name', 'Final lightning');
set(h, 'Name', 'Final corresponding mesurements');
 
Mask=ones(var.img_row, var.img_col); %Mask(2:var.img_row-1, 2:var.img_col-1)=1; %Mask that deletes boundary
Mask1=Mask(:);

recZ=reshape(ZGBR,[var.img_row, var.img_col]);
fprintf('Z diff %.4g\n', norm((Sim.Z - recZ).*Mask,'fro'));
fprintf('L diff %.4g\n', norm(Sim.L - LGBR,'fro'));
fprintf('S diff %.4g\n', norm((Sim.S - SGBR).*repmat(Mask1,[1,3]),'fro'));
fprintf('R diff %.4g\n', norm((Sim.R - RGBR).*Mask1));

figure;
subplot(3,3,1),imagesc((reshape(Sim.S(:,1), [var.img_row, var.img_col])));
subplot(3,3,2),imagesc((reshape(ini_SGRB(:,1), [var.img_row, var.img_col])));
subplot(3,3,3),imagesc((reshape(SGBR(:,1), [var.img_row, var.img_col])));

subplot(3,3,4),imagesc((reshape(Sim.S(:,2), [var.img_row, var.img_col])));
subplot(3,3,5),imagesc((reshape(ini_SGRB(:,2), [var.img_row, var.img_col])));
subplot(3,3,6),imagesc((reshape(SGBR(:,2), [var.img_row, var.img_col])));

subplot(3,3,7),imagesc((reshape(Sim.S(:,3), [var.img_row, var.img_col])));
subplot(3,3,8),imagesc((reshape(ini_SGRB(:,3), [var.img_row, var.img_col])));
subplot(3,3,9),imagesc((reshape(SGBR(:,3), [var.img_row, var.img_col])));

figure;
subplot(2,2,1)
semilogy(result.out.obj); title('Convergence plot of overall objective function')
subplot(2,2,2)
semilogy(result.out.obj1); title('Convergence plot of Data consistency term')
subplot(2,2,3)
semilogy(result.out.obj2); title('Convergence plot of Integrability Constraint term')
subplot(2,2,4)
semilogy(result.out.obj2); title('Covergence plot of Constraint violation');

toc
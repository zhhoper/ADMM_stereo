tic
clear;  close all; 
addpath('solvers');
addpath('functions');
addpath('external');
addpath('simulate')
addpath('newlyADD');

% generate data
Sim = generate_data(0,1);

% synthesize a mask
mask = syn_generateMask(Sim.resolutionX, Sim.resolutionY);

% initialzie data
[X, Y, lam, l3, Z, lambda, ini_SGRB, error] = initialize_ADMM(Sim, mask);
% [~, ~, ~, ~, ~, ~, ini_SGRB] = initialize_ADMM(Sim);
% [X, Y, lam, l3, Z, lambda] = initialize_random_ADMM(Sim);
var = struct;
var.M = Sim.MNoise(:, mask(:));
var.W = Sim.INC(:, mask(:));
var.X = X;
var.Y = Y;
var.lam = lam;
var.l3 = l3;
var.Z = Z;
var.lambda = lambda;
var.mask = mask;
[var.img_row, var.img_col] = size(Sim.Z);

% what should be the correct parameter?
var.c2 = 1;
var.tau = 1;

result = ADMM_opt(var);

% numPixel = var.img_row*var.img_col;
numPixel = sum(mask(:));  % number of valid pixels

depth = result.var.Z;
light = [result.var.X(3:end, 1:2), result.var.l3];
normal = [-result.var.X(1:2, 3:end)', ones(numPixel,1)];
%normal = [result.var.X(1:2, 3:end)', -ones(numPixel,1)];
albedo = -(sqrt(sum(normal.^2, 2)))./result.var.lam;
normal = my_normr(normal);

%depth=reconstructDepthMap_adapted(-normal',[var.img_row,var.img_col]);
t_depth = depth(mask(:));
Rmask = Sim.R(mask(:));
% [LGBR,ZGBR,SGBR,RGBR,G] = applyGBR(Sim.Z, Sim.R, depth, light, normal, albedo, mask);
[LGBR,ZGBR,SGBR,RGBR,G] = applyGBR(Sim.Z, Rmask, t_depth, light, normal, albedo, mask);

if (LGBR(:,1)'*Sim.L(:,1))<0
    LGBR=-LGBR; SGBR=-SGBR;
end

h = showLights(Sim.L, LGBR, light);
% showNormals(normal, reshape(depth,[var.img_row var.img_col]), 'Surface normals BEFORE GBR');

[xEdge, yEdge] = findEdges(mask);
tmp_mask = logical(mask.*(~xEdge).*(~yEdge));
xEdgeInd = xEdge(mask(:));
yEdgeInd = yEdge(mask(:));
edgeInd = xEdgeInd | yEdgeInd;
SGBR = SGBR(~edgeInd,:);

S1 = vec2mat_mask(SGBR(:,1), tmp_mask);
S2 = vec2mat_mask(SGBR(:,2), tmp_mask);
S3 = vec2mat_mask(SGBR(:,3), tmp_mask);
tS = [S1(:), S2(:), S3(:)];

showNormals(tS, ZGBR, 'Surface normals AFTER GBR');
set(h, 'Name', 'Final lightning');
set(h, 'Name', 'Final corresponding mesurements');
 

recZ=reshape(ZGBR,[var.img_row, var.img_col]);
fprintf('Z diff %.4g\n', norm((Sim.Z - recZ).*tmp_mask,'fro')/sum(tmp_mask(:)) );
fprintf('L diff %.4g\n', norm(Sim.L - LGBR,'fro'));
fprintf('S diff %.4g\n', norm((Sim.S(tmp_mask(:),:) - SGBR),'fro')/sum(tmp_mask(:)) );
fprintf('R diff %.4g\n', norm((Sim.R(tmp_mask(:)) - RGBR(~edgeInd)))/sum(tmp_mask(:)) );

fprintf('..........................\n');
fprintf('error for baseline\n');
fprintf('Z diff %.4g\n', error.depth );
fprintf('L diff %.4g\n', error.light );
fprintf('S diff %.4g\n', error.normal );
fprintf('R diff %.4g\n', error.albedo );


figure;
subplot(3,3,1),imagesc((reshape(Sim.S(:,1), [var.img_row, var.img_col])));
% subplot(3,3,2),imagesc(vec2mat_mask(ini_SGRB(:,1), mask));
% subplot(3,3,3),imagesc(vec2mat_mask(SGBR(:,1), mask));
subplot(3,3,2),imagesc(vec2mat_mask(ini_SGRB(~edgeInd,1), tmp_mask));
subplot(3,3,3),imagesc(vec2mat_mask(SGBR(:,1), tmp_mask));



subplot(3,3,4),imagesc((reshape(Sim.S(:,2), [var.img_row, var.img_col])));
% subplot(3,3,5),imagesc(vec2mat_mask(ini_SGRB(:,2), mask));
% subplot(3,3,6),imagesc(vec2mat_mask(SGBR(:,2), mask));
subplot(3,3,5),imagesc(vec2mat_mask(ini_SGRB(~edgeInd,2), tmp_mask));
subplot(3,3,6),imagesc(vec2mat_mask(SGBR(:,2), tmp_mask));


subplot(3,3,7),imagesc((reshape(Sim.S(:,3), [var.img_row, var.img_col])));
% subplot(3,3,8),imagesc(vec2mat_mask(ini_SGRB(:,3), mask));
% subplot(3,3,9),imagesc(vec2mat_mask(SGBR(:,3), mask));
subplot(3,3,8),imagesc(vec2mat_mask(ini_SGRB(~edgeInd,3), tmp_mask));
subplot(3,3,9),imagesc(vec2mat_mask(SGBR(:,3), tmp_mask));

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
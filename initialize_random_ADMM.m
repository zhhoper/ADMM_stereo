function [X, Y, lam, l3, Z, lambda] = initialize_random_ADMM(Sim)
% to initialzie the variables for ADMM

[row, col] = size(Sim.Z);
Sol = SimpleSolver(Sim.MNoise, [row, col], 3, struct());
Sol.solve();
[M, L, S, R, Z] = Sol.extract(); %L=-L; S=-S;
[LGBR, ZGBR, SGBR, RGBR, G] =  applyGBR(Sim.Z, Sim.R, Z, L, S, R);

numImg = size(L,1);
numPixel = row*col;


% baseline %

% 
% X = zeros(numImg+2, numPixel+2);
% X(1:2, 1:2) = eye(2);
% X(3:end, 1:2) = -LGBR(:, 1:2);
% %X(1:2, 3:end) = -[SGBR(:,1)'./SGBR(:,3)'; SGBR(:,2)'./SGBR(:,3)'];
% X(1:2, 3:end) = -[SGBR(:,1)'./SGBR(:,3)'; SGBR(:,2)'./SGBR(:,3)'];
% X(3:end, 3:end) = M*diag(lam) + l3*ones(1, numPixel);
% Y = X;
% lambda = zeros(size(X));

% random

lam = (-1 + 100)*rand(row*col,1)-100;
l3 = 50*rand(numImg,1)+190;
Z = 10*rand(row, col)+15;

X = rand(numImg+2, numPixel+2);
X(1:2, 1:2) = eye(2); 

%GT initialization
L1=500*rand(numImg,2)-250; S1=2*rand(2,row*col)-1;

X(1:2, 3:end) = S1; %initialize surface normal only
X(3:end, 1:2) = L1;
X(3:end, 3:end) =L1*S1;%Sim.MNoise*diag(lam) + l3*ones(1, numPixel);

Y = X;

lambda = zeros(size(X));
% Sol = SimpleSolver(Sim.MNoise, [row, col], 3, struct());
% Sol.solve();
%
% [M, L, S, R, Z] = Sol.extract();
% [LGBR, ZGBR, SGBR, RGBR, G] =  applyGBR(Sim.Z, Sim.R, Z, L, S, R);
%
% numImg = size(L,1);
% numPixel = row*col;
%
% lam = RGBR;
% l3 = LGBR(:,3);
% Z = ZGBR;
% X = zeros(numImg+2, numPixel+2);
% X(1:2, 1:2) = eye(2);
% X(3:end, 1:2) = LGBR(:, 1:2);
% X(1:2, 3:end) = -[SGBR(:,1)'./SGBR(:,3)'; SGBR(:,2)'./SGBR(:,3)'];
% X(3:end, 3:end) = M*diag(lam);
% Y = X;
% lambda = zeros(size(X));

end


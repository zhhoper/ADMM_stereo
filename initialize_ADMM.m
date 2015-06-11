function [X, Y, lam, l3, Z, lambda, SGBR, error, final] = initialize_ADMM(Sim, mask)
% to initialzie the variables for ADMM

[row, col] = size(Sim.Z);

Sol = SimpleSolver(Sim.MNoise, [row, col], 3, mask, struct());
Sol.solve();

[M, L, S, R, Z] = Sol.extract();

maskR = Sim.R(mask(:));


[LGBR, ZGBR, SGBR, RGBR, G] =  applyGBR(Sim.Z, maskR, Z, L, S, R, mask);

if (LGBR(:,1)'*Sim.L(:,1))<0
    LGBR=-LGBR; SGBR=-SGBR;
end

numImg = size(L,1);
% numPixel = row*col;

% we will use the pixels within the mask to deal with this
numPixel = sum(mask(:));

lam=-1./(abs(SGBR(:,3)).*RGBR); %lam = RGBR;

l3 = LGBR(:,3);
Z = ZGBR;
X = zeros(numImg+2, numPixel+2);
X(1:2, 1:2) = eye(2);
X(3:end, 1:2) = LGBR(:, 1:2);
%X(1:2, 3:end) = -[SGBR(:,1)'./SGBR(:,3)'; SGBR(:,2)'./SGBR(:,3)'];
X(1:2, 3:end) = -[SGBR(:,1)'./SGBR(:,3)'; SGBR(:,2)'./SGBR(:,3)'];
X(3:end, 3:end) = M*diag(lam) + l3*ones(1, numPixel);
Y = X;
lambda = zeros(size(X));


recZ=reshape(ZGBR,[row, col]);
error.depth = norm((Sim.Z - recZ).*mask,'fro')/numPixel;
error.light = norm(Sim.L - LGBR,'fro');
error.normal = norm((Sim.S(mask(:),:) -SGBR),'fro')/numPixel;
error.albedo = norm((Sim.R(mask(:)) - RGBR))/numPixel;
error.M=norm(Sim.M(:,mask(:)) - M);

final.light=LGBR;
final.depth=ZGBR;
final.normal=SGBR;
final.albedo=RGBR;

end


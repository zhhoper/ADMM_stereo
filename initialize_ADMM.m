function [X, Y, lam, l3, Z, lambda, SGBR,error,final] = initialize_ADMM(Sim)
% to initialzie the variables for ADMM

[row, col] = size(Sim.Z);

% synthesize a mask
mask = syn_generateMask(row, col);

Sol = SimpleSolver(Sim.MNoise, [row, col], 3, struct());
Sol.solve();

[M, L, S, R, Z] = Sol.extract();



[LGBR, ZGBR, SGBR, RGBR, G] =  applyGBR(Sim.Z, Sim.R, Z, L, S, R);

if (LGBR(:,1)'*Sim.L(:,1))<0
    LGBR=-LGBR; SGBR=-SGBR;
end


numImg = size(L,1);
numPixel = row*col;


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

Mask=ones(row, col); %Mask(2:var.img_row-1, 2:var.img_col-1)=1; %Mask that deletes boundary
Mask1=Mask(:);

 recZ=reshape(ZGBR,[row, col]);
error.depth = norm((Sim.Z - recZ),'fro')/numPixel;
error.light = norm(Sim.L - LGBR,'fro');
error.normal = norm((Sim.S - SGBR).*repmat(Mask1,[1,3]),'fro')/numPixel;
error.albedo = norm((Sim.R - RGBR).*Mask1)/numPixel;
error.M=norm(Sim.M - M);


final.light=LGBR;
final.depth=ZGBR;
final.normal=SGBR;
final.albedo=RGBR;

end


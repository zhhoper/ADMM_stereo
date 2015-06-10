close all; clear all; clear
addpath('../solvers'); addpath('../functions'); addpath('../external'); addpath('../simulate');

%seed = round(rand()*1e6);
%% Test Case 1 : no Noise, no missing data
seed = 58597;
rng(seed);

size = [40,40];
noiseStrength = 0.0;
knowledgeRatio = 0.4;
rows = size(1);
cols = size(2);
lights = 20;
pixels = prod(size);

scenario.object = 'peaks';          % random || peaks     ||  ellipse
scenario.lighting = 'random';       % random    ||  fitted
scenario.albedo = 'random';         % number || 'random'


Sim = Simulator(lights,size,scenario,noiseStrength, knowledgeRatio);
Sim.simulate();
showNormals(Sim.S, Sim.Z, 'Ground Truth');
f = showImages(Sim.MNoise, size);
set(f,'Name','Input Images');


p.imgsize = size;

%[M, L, S, R, Z, solresult.breakReason] = PlainSolver(Sim.MNoise, p);

Sol=SimpleSolver(Sim.MNoise,size,3, struct());
Sol.solve();

[M, L, S, R, Z] = Sol.extract(); %L=-L; S=-S;

[LGBR,ZGBR,SGBR,RGBR,G] = applyGBR(Sim.Z, Sim.R, Z, L, S, R);
if (LGBR(:,1)'*Sim.L(:,1))<0
    LGBR=-LGBR; SGBR=-SGBR;
end

h = showLights(Sim.L, LGBR, L);
% showNormals(S, reshape(Z,size), 'Surface normals BEFORE GBR');
showNormals(SGBR, ZGBR, 'Surface normals AFTER GBR');
% set(h, 'Name', 'Final lightning');
% set(h, 'Name', 'Final corresponding mesurements');
fprintf('Z diff %.4g\n', norm(Sim.Z - reshape(ZGBR,size),'fro'));
fprintf('L diff %.4g\n', norm(Sim.L - LGBR,'fro'));
fprintf('S diff %.4g\n', norm(Sim.S - SGBR,'fro'));
fprintf('R diff %.4g\n', norm(Sim.R - RGBR));
fprintf('M diff %.4g\n', norm(Sim.M - M));

knownInds = find(Sim.MNoise>0);
kr = numel(knownInds)/ numel(Sim.MNoise);
fprintf('%.4g knowledge ratio\n',kr);


%% TODO: This should be a GBR, but it isn't...
Sim.S\S

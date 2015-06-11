function Sim = generate_data(nz,know)

seed = 58597;
rng(seed);

imsize = [40,40];
noiseStrength = nz;
knowledgeRatio = know;
rows = imsize(1);
cols = imsize(2);
lights = 20;
pixels = prod(imsize);

scenario.object = 'ellipse';          % random || peaks     ||  ellipse
scenario.lighting = 'random';       % random    ||  fitted
scenario.albedo = 'random';         % number || 'random'


Sim = Simulator(lights,imsize,scenario,noiseStrength, knowledgeRatio);
Sim.simulate();
% showNormals(Sim.S, Sim.Z, 'Ground Truth');
% f = showImages(Sim.MNoise, imsize);
% set(f,'Name','Input Images');
end
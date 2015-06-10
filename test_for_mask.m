
x = rand(10,10);
mask = syn_generateMask(10,10);
x = x.*mask;

y = mat2vec_mask(x, mask);

tx = vec2mat_mask(y, mask);
ccc = 0;

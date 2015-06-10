function [ fig] = plotNoiseDepthError(stats)




%stats = stats(stats.KnowledgeRatio>=0.9,:);
fig = plotStats(stats, 'NoiseStrength', 'errDepth');

xlabel('noise level');
ylabel('recovered depth error per pixel');
title('Noise Level ~ Recovery error');


end







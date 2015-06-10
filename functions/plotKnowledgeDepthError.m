function [ fig] = plotKnowledgeDepthError( stats, noiseStrength)

%stats = stats(stats.NoiseStrength<=noiseStrength,:);

fig = plotStats(stats, 'KnowledgeRatio', 'errDepth','median','linearinterp');
end


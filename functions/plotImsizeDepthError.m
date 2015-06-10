function [ fig] = plotImsizeDepthErr( stats)



stats.Scenario = nominal(categorical(stats.Scenario,'Ordinal',true));
stats.Method = nominal(categorical(stats.Method,'Ordinal',true));
stats.NoiseStrength = nominal(categorical(stats.NoiseStrength,'Ordinal',true));
stats = stats(stats.KnowledgeRatio==1,:);
stats = stats(stats.Method == 'simple',:);

stats{:,'errDepth'} = stats{:,'errDepth'}./stats{:,'Pixels'};


fig = figure;
hold on;


COLORS = ['r', 'g','b','y','m','c','w','k'];
pos = 1;
for c = categories(stats.NoiseStrength)'
    statsC = stats(stats.NoiseStrength == c,:);
    scatter(statsC.Pixels, statsC.errDepth,'x');
    f = fit(statsC{:,3},statsC{:,11},'exp1');
%     plot(f, COLORS(pos));
    pos = pos +1;    
end


legend(categories(stats.NoiseStrength));



xlabel('number of pixels per image');
ylabel('recovered depth error per pixel');
title('Image size ~ Recovery error');


hold off;

end
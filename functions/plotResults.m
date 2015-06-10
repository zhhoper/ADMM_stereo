%% Plot results
close all;

%stats = qextract(res_qtesting);

stats.Method = nominal(categorical(stats.Method,'Ordinal',true));
stats.Scenario = nominal(categorical(stats.Scenario,'Ordinal',true));
stats.breakReason = nominal(categorical(stats.breakReason,'Ordinal',true));

stats = stats(stats.Method ~= 'plain',:);
%stats = stats(stats.breakReason ~= 'stalling',:);

stats{stats.KnowledgeRatio>0.95,6} = 1;
stats.KnowledgeRatio = round(stats.KnowledgeRatio/.01) * .01;

for s = categories(stats.Scenario)'
    statsF = stats(stats.Scenario == s, :);
    kr = 1;
    fig1 = plotNoiseDepthError(statsF(statsF.KnowledgeRatio>=kr,:));
    title(strcat(s, num2str(kr)));
    kr = 0.75;
    fig1 = plotNoiseDepthError(statsF(statsF.KnowledgeRatio==kr,:));
    title(strcat(s, num2str(kr)));
    kr = 0.5;
    fig1 = plotNoiseDepthError(statsF(statsF.KnowledgeRatio==kr,:));
    title(strcat(s, num2str(kr)));   
    nr = 0;
    fig2 = plotKnowledgeDepthError(statsF(statsF.NoiseStrength==nr,:));
    title(strcat(s, num2str(nr)));   
    nr = 0.05;
    fig2 = plotKnowledgeDepthError(statsF(statsF.NoiseStrength==0.05,:));
    title(strcat(s, num2str(nr)));   
    nr = 0.1;
    fig2 = plotKnowledgeDepthError(statsF(statsF.NoiseStrength==0.1,:));
    title(strcat(s, num2str(nr)));   
end
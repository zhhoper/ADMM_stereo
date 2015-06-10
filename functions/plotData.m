function [ fig] = plotStats(stats, nameX, nameY)


stats{:,'errDepth'} = stats{:,'errDepth'}./stats{:,'Pixels'};

names = stats.Properties.VariableNames;
xvar = strcmp(nameX,names);
yvar = find(strcmp(nameY,names));


fig = figure;
hold on;


COLORS = ['r', 'g','b','y','m','c','w','k'];
pos = 1;
plots = [];
for c = categories(stats.Method)'
    statsC = stats(stats.Method == c{1},:);
    scatter(statsC{:,xvar}, statsC{:,yvar},'.', COLORS(pos));
    try 
        means = grpstats([statsC.NoiseStrength, statsC.errDepth], statsC.NoiseStrength, {@median});
        if (numel(means) > 1)
            line = fit(means(:,1), means(:,2),'linearinterp');
            x = min(stats.NoiseStrength):0.001:max(stats.NoiseStrength);
            y = feval(line,x);
            p = plot(x,y, 'DisplayName', c{1}, 'Color', COLORS(pos));
            plots = [plots; p];
            pos = pos +1;
        end
    catch 
        fprintf('There was an error in plotting the means...\n');
    end
        
end


legend(plots);



xlabel('noise level');
ylabel('recovered depth error per pixel');
title('Noise Level ~ Recovery error');


hold off;

end







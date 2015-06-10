%% Plots stats
function [ fig] = plotStats(varargin)
stats = varargin{1};
nameX = varargin{2};
nameY = varargin{3};

usePoints = 'median';
fitFunction = 'linearinterp';
if nargin > 4
    usePoints = varargin{4};
    fitFunction = varargin{5};
end

fig = figure;
hold on;

stats.Scenario = nominal(categorical(stats.Scenario,'Ordinal',true));
stats.Method = nominal(categorical(stats.Method,'Ordinal',true));
stats{:,'errDepth'} = stats{:,'errDepth'}./stats{:,'Pixels'};

names = stats.Properties.VariableNames;
xvar = strcmp(nameX,names);
yvar = strcmp(nameY,names);




COLORS = ['r', 'g','b','m','c','y','w','k'];
pos = 0;
plots = [];
for c = categories(stats.Method)'
    try
    pos = pos +1;
    statsC = stats(stats.Method == c{1},:);
    %scatter3(statsC{:,'KnowledgeRatio'},statsC{:,'NoiseStrength'},statsC{:,yvar});
    %scatter(statsC{:,xvar}, statsC{:,yvar},'.', COLORS(pos));
    switch usePoints
        case 'median'
            points = grpstats([statsC{:,xvar}, statsC{:,yvar}], statsC{:,xvar}, {@median});
        case 'mean'
            points = grpstats([statsC{:,xvar}, statsC{:,yvar}], statsC{:,xvar});
        otherwise
            points = [statsC{:,xvar}, statsC{:,yvar}];
    end
     
        
            line = fit(points(:,1), points(:,2),fitFunction);
            x = min(stats{:,xvar}):0.001:max(stats{:,xvar});
            y = feval(line,x);
            p = plot(x,y, 'DisplayName', c{1}, 'Color', COLORS(pos));
            plots = [plots; p];
        
    catch e
        e
        fprintf('There was an error in interpolation of %s,%s in method %s using %s.\n',nameX, nameY,c{1},usePoints);
    end
       
end


legend(plots);



xlabel(nameX);
ylabel(nameY);



hold off;

end







function [ h1, h2] = showDepthMap(varargin)

Z = varargin{1};
if nargin> 1
    ZOrig = varargin{2};
    lo = min(min(Z(:)'), min(ZOrig(:)'));
    hi = max(max(Z(:)'), max(ZOrig(:)'));
end

visible = 'on';
if nargin >= 6
    visible = varargin{6};
end



m = size(Z,1);
n = size(Z,2);

switch nargin
    case 1
        h1 = figure('Name', 'Recovered Surface', 'visible', visible); surf(Z);
        h2 = 0;
    case 2
        h1 = figure('Name', 'Recovered Surface', 'visible', visible); surf(Z);
        ac = gca;
        set(ac,'zlim',[lo, hi]);
        h2 = figure('Name', 'Ground Truth', 'visible', visible); surf(ZOrig);
        ac = gca;
        set(ac,'zlim',[lo, hi]);

    case {4,5,6}
        R = varargin{3};
        ROrig = varargin{4};
        h1 = figure('Name', 'Recovered Surface', 'visible', visible); surf(Z, reshape(R, m, n));
        ac = gca;
        set(ac,'zlim',[lo, hi]);
        caxis([0, 1]);
        colorbar;

        if nargin >=5 
            text = varargin{5};
            MyBox = uicontrol('style','text');
            set(MyBox,'String',text);
            set(MyBox,'Position',[0, 0, 150, 115]);
        end
        
        h2 = figure('Name', 'Ground Truth', 'visible', visible); surf(ZOrig, reshape(ROrig, m,n));
        ac = gca;
        set(ac,'zlim',[lo, hi]);
        caxis([0, 1]);
        colorbar;
        

end

        
        



end


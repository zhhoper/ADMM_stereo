% Saves the measurements to separate images 
function [fig] = showImages(varargin)
    M = varargin{1};
    s = varargin{2};
    if nargin == 3
        vis = varargin{3};
    else 
        vis = 'on';
    end
    cols = ceil(sqrt(size(M,1)));
    rows = ceil(size(M,1)/ cols);
    fig = figure('Name','Input', 'visible',vis);
    
    M(M<0) = 0;
    range = [max(0,min(min(M))), max(max(M))];
    for i=1:size(M,1)
        IMG = mat2gray(reshape(M(i,1:end), s ), range);
        IMG = rot90(IMG);
        subplot(rows,cols,i);
        imshow(IMG);
    end
end

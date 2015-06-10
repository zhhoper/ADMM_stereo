function [MNew] = completeMatrixToPoints(varargin)
M = varargin{1};
[m, n] = size(M);
if nargin == 3 
    
    rangeX = varargin{2};
    rangeY = varargin{3};
    if numel(rangeX) > 2 
        error('rangeX');
    end
    rows = linspace(rangeX(1), rangeX(2),m);
    cols = linspace(rangeY(1), rangeY(2),n);  
else
    
        
   rows = 1:m;
   cols = 1:n;
end




C = createCoordArray(rows,cols);
MNew = [C, M(:)];
end
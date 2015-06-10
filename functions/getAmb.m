function [A, SNew] = getAmb(SGT, S)
% returns the GBR transformation matrix that transforms M into MGT

A = (S' \ SGT')';
SNew = (A*S);

% A = M \ MGT;
% MNew = M * A;


end

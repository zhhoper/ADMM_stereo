function [ RNew ] = gbrAlbedos( G, R, S)
%gbrAlbedos transform albedo values according to the given GBR
%transformation matrix G and the unresolved surface matrix S

G = inv(G)';

n1 = sum((G*S)'.^2,2);
n2 = sum(S'.^2,2);


RNew = R .* sqrt(n1 ./ n2)';


end


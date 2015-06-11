function [LRes,ZRes,SRes,RRes,G] = applyGBR(ZGT,RGT, ZGBR, LGBR, SGBR, RGBR, mask)

[G, ZRes] = getGBR(ZGT,ZGBR, mask);

LRes = LGBR*G';
SRes = my_normr(SGBR/G);
%SRes = my_normr(G\SGBR')';
RRes = RGBR .* sqrt( sum((SGBR/G)'.^2) ./ sum(SGBR'.^2))';
%RRes = RGBR .* sqrt( sum((G\SGBR').^2) ./ sum(SGBR'.^2))';

a = RRes\RGT;
RRes = RRes * a;
LRes = LRes / a;

end
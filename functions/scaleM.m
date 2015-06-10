function MScaled = scaleM( M, MRef )
%SCALEM Try to move and scale M so that it is as close a possible to MRef

mean1 = mean(M(:)');
meanRef = mean(MRef(:)');
dmean = meanRef-mean1;

var1 = var(M(:)');
varRef = var(MRef(:)');
f = varRef/var1;
MScaled = M * f + dmean;

end


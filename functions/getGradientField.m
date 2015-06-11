function [Sx,Sy] = getGradientField(S, imsize)
    S = my_normc(S);
    Sx = reshape(S(1,:) ./ S(3,:), imsize);
    Sy = reshape(S(2,:) ./ S(3,:), imsize);

end
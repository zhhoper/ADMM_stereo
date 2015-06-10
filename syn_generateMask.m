function mask = syn_generateMask(row, col)
% We try to generate a round mask
mask = zeros(row, col);

max_radius = min(row, col)/2;
cx = row/2;
cy = col/2;

radius = max_radius - 3;  % radius is 3 pixels smaller than the maximum radius

for i = 1 : row
    for j = 1 : col
        dist = sqrt((i - cx)^2 + (j - cy)^2);
        if dist < radius
            mask(i,j) = 1;
        end
    end
end

mask = logical(mask);
end
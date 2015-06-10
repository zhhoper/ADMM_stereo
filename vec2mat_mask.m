function A = vec2mat_mask(vec_data, mask)
% A = vec2mat_mask(vec_data, mask)
%
% This function is used to convert a vector to a matrix

[row, col] = size(mask);
tmpA = zeros(row*col,1);
tmpA(mask(:)) = vec_data;

A = reshape(tmpA, [row, col]);

end
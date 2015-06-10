function vec_data = mat2vec_mask(A, mask)
% vec_data = mat2vec_mask(A, mask)
%
% This function is used to convert a matrix to a vector w.r.t. mask

tmpA = A(:);
tmpMask = mask(:);
vec_data = tmpA(tmpMask);

end
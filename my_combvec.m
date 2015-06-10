function comb = my_combvec(c1, c2)
% this function is used to combine c1 and c2 to generate all possible
% combinations of these two vectors

num1 = size(c1,2);
num2 = size(c2,2);

tmp1 = repmat(c1, 1, num2);
tmp2 = repmat(c2, num1, 1);
comb = [tmp1;tmp2(:)'];
end
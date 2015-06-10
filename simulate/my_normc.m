function n = my_normc(x)
[row, ~] = size(x);

n = x./(repmat(sqrt(sum(x.^2, 1)),[row, 1]));
end
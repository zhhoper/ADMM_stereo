function v = shrinkToSum(v, s)
%% shrink the values in the positive, ordered vector, s.t. sum(v) = s and all 
%% elements in v >= 0
t = v(v>0);
diff = sum(t)-s;
while diff > 0
    t = t-diff/numel(t);
    diff = -sum(t(t<0));
    t = t(t>0);
end

m = size(t,2);
l = size(v,2);

v = [t, zeros(1,l-m)];


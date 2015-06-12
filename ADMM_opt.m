function result = ADMM_opt(var)
% var contains all the variables we need to solve this problem
% var.M contains the observations
% var.W contains the mask
% var.X, var.l3, var.lam, var.Z, var.img_row, var_img_col

max_ite = 600;
threshold = 0.001;
count = 0;
%var_old = var;
diff = inf;

[obj, obj1, obj2, obj3] = objective_fun(var); %obj_t=obj1+obj2+obj3;
fprintf('count: %d, obj %.5e, obj1 %.5e, obj2 %.5e, obj3 %.5e\n', 0, obj, obj1, obj2, obj3);

out = struct;
out.residual = zeros(max_ite,1);
out.obj = zeros(max_ite,1);

while count < max_ite && diff > threshold %
    obj_old = obj;
    var.Z = opt_Z(var.X, var.mask, var.img_row, var.img_col);
     %[~, ~, obj1,~] = objective_fun(var);
     %obj1
    var.lam = opt_lam(var);
    %obj2 = objective_fun(var);
    var.l3 = opt_l3(var);
    %obj3 = objective_fun(var);
    var.X = opt_X(var);
    %obj4 = objective_fun(var);
    var.Y = opt_Y(var);
    %obj5 = objective_fun(var);
    var.lambda = opt_lambda(var);
    %obj6 = objective_fun(var);
    %fprintf('obj value %.5e, %.5e, %.5e, %.5e, %.5e, %.5e \n', obj1, obj2, obj3, obj4, obj5, obj6);
%     residual = norm(var.Z - var_old.Z, 'fro')...
%         + norm(var.lam - var_old.lam, 'fro')...
%         + norm(var.l3 - var_old.l3, 'fro')...
%         + norm(var.X - var_old.X, 'fro')...
%         + norm(var.Y - var_old.Y, 'fro')...
%         + norm(var.lambda - var_old.lambda, 'fro');
    
%     norm_old = norm(var_old.Z,'fro') + norm(var_old.lam, 2)...
%         + norm(var_old.l3, 2) + norm(var_old.X, 'fro')...
%         + norm(var_old.Y, 'fro') + norm(var.lambda, 'fro');
%     
%     residual = residual/norm_old;
%    var_old = var;
    count = count + 1;
%    out.residual(count) = residual;
    [obj, obj1, obj2, obj3] = objective_fun(var);
    %obj_t=obj1+obj2+obj3;
    diff = abs(obj - obj_old);
    %if mod(count,10) == 1
    if count == 20
        tmp = fix(log(obj3/obj1));
        if tmp > 0
            var.tau = var.tau*(tmp+1);
%         elseif tmp < 0
%             var.tau = var.tau/(-tmp + 1);
        end
    end
    out.obj(count) = obj;
    out.obj1(count)=obj1;
    out.obj2(count)=obj2;
    out.obj3(count)=obj3;
    fprintf('count: %d, obj %.5e, obj1 %.5e, obj2 %.5e, obj3 %.5e\n', count, obj, obj1, obj2, obj3);
end

[U, S, V] = svds(var.X, 2);
var.X=U*S*V';

result = struct;
result.var = var;
result.out = out;


end

function [obj, obj1, obj2, obj3] = objective_fun(var)
% M, W, X, l3, lam, Z, Y, lambda, img_row, img_col)
% compute objecitve function
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col, var.c2, var.tau
% NOTE: WE TREATE THE X COORDINATE AS ALONG THE COLUMN, Y COORDINATE AS
% ALONG THE ROW

c2 = var.c2;
tau = var.tau;

%numPixel = var.img_row*var.img_col;   % number of pixels
numPixel = sum(var.mask(:));   % number of valide pixels
[row, col] = size(var.X);

x12_x = var.X(1, 3:col);
x12_y = var.X(2, 3:col);
x22 = var.X(3:row, 3:col);

img_Z = var.Z;
[GZ_x, GZ_y] = gradient_xy(img_Z);
GZ_x = GZ_x(var.mask(:));
GZ_y = GZ_y(var.mask(:));

tmpOne = ones(1, numPixel);
obj1 = 0.5*norm(var.W.*(var.M*diag(var.lam) - x22 + var.l3*tmpOne), 'fro')^2;
obj2 = (norm(GZ_x(:) - x12_x')^2 + norm(GZ_y(:) - x12_y')^2);
obj3 = norm(var.Y - var.X, 'fro')^2;
obj33 = norm(var.Y - var.X +  var.lambda, 'fro')^2;
obj = obj1 + 0.5*c2*obj2 + 0.5*tau*obj33;
end

function Z = opt_Z(X, mask, img_row, img_col)
% optimize w.r.t. Z
% solve possion equation to get the result

% gx = reshape(X(1, 3:end), [img_row, img_col]);
% gy = reshape(X(2, 3:end), [img_row, img_col]);
% Z = recover_img(gx, gy);

len=length(X(1:2,3:end));
S=[-X(1:2,3:end);ones(1,len)];
Z1 = reconstructDepthMap_adapted_mask(S, [img_row,img_col], mask);
%Z1 = reconstructDepthMap(S,[img_row,img_col]);

%Z=reshape(Z1,[img_row,img_col]);
Z = vec2mat_mask(Z1, mask);
end

function lam = opt_lam(var)
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col
% optimize w.r.t. lam
% NOTE: LAM IS A DIAGONAL MATRIX

%numPixel = var.img_row*var.img_col;
numPixel = sum(var.mask(:));  % number of valid pixels
tmpOne = ones(1, numPixel);

X22 = var.X(3:end, 3:end);
a = X22 - var.l3*tmpOne;

tmp1 = var.M'*(var.W.*var.M);
tmp2 = var.M'*(var.W.*a);

% we need to consider about divided by 0
numerator = diag(tmp1'*tmp2);
denominator = diag(tmp1'*tmp1);
tind1 = denominator ~= 0;

% if denominator is 0, we use the old value
lam = var.lam;
lam(tind1) = min(numerator(tind1)./denominator(tind1), -1);

end

function l3 = opt_l3(var)
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col
% optimize w.r.t. l3

a = var.W.*(var.M*diag(var.lam) - var.X(3:end, 3:end));
% numPixel = var.img_row*var.img_col;
numPixel = sum(var.mask(:));  % number of valid pixels
tmpOne = ones(1, numPixel);

numerator = -a*tmpOne';
denominator = sum(var.W, 2);
tind1 = denominator ~= 0;
% for those denominator is 0, we just use the old value
l3 = var.l3;
l3(tind1) = numerator(tind1)./denominator(tind1);

end

function X = opt_X(var)
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col
% optimize w.r.t. X

X = var.X;

% get x2,1
X(3:end, 1:2) = var.Y(3:end, 1:2) + var.lambda(3:end, 1:2);


% get X1,2
[gx, gy] = gradient_xy(var.Z);
[xEdge, yEdge] = findEdges(var.mask);
edge = xEdge | yEdge;

gx(edge(:)) = 0;
gy(edge(:)) = 0;

gx = gx(var.mask(:));
gy = gy(var.mask(:));

X(1:2, 3:end) = (var.tau*var.Y(1:2, 3:end) + var.tau*var.lambda(1:2, 3:end)...
    +var.c2*[gx, gy]')/(var.tau + var.c2);

% get X2,2
% numPixel = var.img_row*var.img_col;
numPixel = sum(var.mask(:));  % number of valid pixels
tmpOne = ones(1, numPixel);
B = var.tau*(var.Y(3:end, 3:end) + var.lambda(3:end, 3:end)) + var.W.*...
    (var.M*diag(var.lam) + var.l3*tmpOne);

c = var.W + var.tau;  % c cannot be 0
X(3:end, 3:end) = B./c;
end

function Y = opt_Y(var)
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col
% optimize w.r.t. Y

tmpY = var.X - var.lambda;
% [U, S, V] = svd(tmpY);
% Y = U(:, 1:3)*S(1:3, 1:3)*V(:, 1:3)';

[U, S, V] = svds(tmpY, 2);
Y = U*S*V';
end

function lambda = opt_lambda(var)
% var.M is the observation
% var.W is the mask
% var.lambda is the lagrangian parameter.
% var.X, var.l3, var.lam, var.Z, var.Z, var.lambda, var.img_row,
% var.img_col
% optimize w.r.t. lambda

lambda = var.lambda + (var.Y - var.X);
end

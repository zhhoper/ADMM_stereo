function [Z, C, D] = reconstructDepthMap(S,imgsize)
addpath('../simulate');
rows = imgsize(1);
cols = imgsize(2);

pixels = rows*cols;
num_equ = (rows-1)*cols + (cols-1)*rows;
C = zeros(num_equ, pixels);
D = zeros(num_equ, 1);
% x direction
pos = 0;
[Sx, Sy] = getGradientField(S,imgsize);
for j=1:cols
    for i=1:rows-1
        
        ind = sub2ind(imgsize,i,j);
        dx = Sx(i,j);
        if not(isnan(dx))
            pos = pos +1;
            D(pos) = dx;
            V = zeros(1, pixels);
            V(1, ind) = -1;
            V(1, sub2ind(imgsize,i+1,j)) = 1;
            C(pos,:) = V;
        end
    end
end
% y direction
for i=1:rows
    for j=1:cols-1
        
        ind = sub2ind(imgsize,i,j);
        dy = Sy(i,j);
        if not(isnan(dy))
            pos = pos +1;
            D(pos) = dy;
            V = zeros(1, pixels);
            V(1, ind) = -1;
            V(1, sub2ind(imgsize,i,j+1)) = 1;
            C(pos,:) = V;
            
        end
    end
end
C(pos+1,1) = 1;
D(pos+1) = 0;


Z = C\D;
% Solve Least squares with regularization constant lambda
%             lambda = max(eig(C' * C)) * 1e-3;
%             Z = inv(C' * C + lambda * eye(size(C,2))) * C' * D;

end

function [Sx,Sy] = getGradientField(S,imgsize)
S = my_normc(S);
Sx = reshape(S(1,:) ./ S(3,:), imgsize);
Sy = reshape(S(2,:) ./ S(3,:), imgsize);

end
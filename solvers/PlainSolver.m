%% This is a very primitive solver for comparison. 
%% It always returns closest plain as a surface.
%%
%% params: imgsize
%% result:  M,L,S,R,Z, breakReason
function [M,L,S,R,Z, breakReason] = PlainSolver(MIn, params) 

    knownInds = MIn > 0;
    if size(knownInds(knownInds==0),1) > 0
        fprintf('Completing matrix using Wiberg\n');
        try 
            [U, V] = damped_wiberg_new(MIn, knownInds, 3);
        catch 
            fprintf('Wiberg faild, using Cabral\n');
            [~,U,V] = l2_rpca_mask_alm_fast(MIn,knownInds, 3);
        end
        M = U * V';
    else 
        M = MIn;
    end


    %% Rank 3 approximation
    [A,B,C] = svd(M);
    B = sqrt(B(1:3,1:3));
    L = A(:,1:3)*B;
    S = B*C(:,1:3)';

    [S,A] = projectToIntegrableSurface(S, params.imgsize);
    L = L/A;
    %% Integrability
    pixels = size(M,2);
    C = zeros(pixels*2,2);    
    C(1:pixels,1) = S(3,:)' ./ S(1,:)';
    C(pixels+1:end,2) = S(3,:)' ./ S(2,:)';
    
    x = C \ -ones(pixels*2,1);
    
    [a,b] = ind2sub(params.imgsize,1:pixels);
    Z = reshape([a;b]'*x,params.imgsize);
    
    
    R = ones(size(M,2),1);
    S = S';
    breakReason = 'done';
end

function [S, A, C] = projectToIntegrableSurface(S, imgsize)
            
npixels = size(S,2);
C = zeros(npixels,6);
for p = 1:size(S,2)
    [x,y] = ind2sub(imgsize, p);
    if (x ~= imgsize(1) && y ~= imgsize(2)) 
        s = S(:,p);
        s_x = S(:,sub2ind(imgsize,x+1,y)) - s;
        s_y = S(:,sub2ind(imgsize,x,y+1)) - s;
        C(p,:) = [cross(s_x,s)', cross(s_y,s)'];
    end
end


[~,~,VC] = svd(C);
x = VC(:,end);


CoFac = zeros(3,3);
CoFac(1,:) = x(1:3);
CoFac(2,:) = x(4:6);
CoFac(3,3) = 1;

A = inv(CoFac)';
S = A*S;

end
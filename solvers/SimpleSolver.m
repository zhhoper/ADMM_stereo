%% SimpleSolver 
% Solves the Photometric stereo problem by first completing the missing
% entries and then using svd to obtain a factorization $M = LS$ 
% Next S is projected onto the set of integrable surface normals and then the
% depthmap Z is reconstructed using the projected, integrable surface normals.

classdef SimpleSolver < handle
    properties 
        M;
        r;
        imgsize;
        params;    
        result;
        MNew;
        LNew;
        SNew
        ZNew;
        LGT = [];
        SGT = [];
        ANew;
    end
    methods
        function Sol = SimpleSolver(M, imgsize, r, params)
            Sol.r = r;
            Sol.M = M;
            Sol.params = params;
            Sol.imgsize = imgsize;
        end
        
        function [] = initialize(varargin)    
            Sol = varargin{1};
            if nargin > 1
                Sol.LGT = varargin{3};
                Sol.SGT = varargin{4};
            end
        end
        
                
        function result = solve(Sol) 
            
            Sol.MNew = Sol.M;
            % Run Wiberg to complete matrix if we have missing data
            knownInds = Sol.M > 0;
            if size(knownInds(knownInds==0),1) > 0
                fprintf('Completing matrix using Wiberg\n');
                %display(Sol.M);
                %display(knownInds);
                try 
                    [U, V] = damped_wiberg_new(Sol.M, knownInds, 3);
                    
                catch 
                    
                    fprintf('Wiberg faild, using Cabral\n');
                    [~,U,V] = l2_rpca_mask_alm_fast(Sol.M,knownInds, 3);
                end
                
                Mtemp=U * V';
                Mtemp(knownInds)=Sol.M(knownInds);
                Sol.MNew= Mtemp;
            end
            %else
            
            
            
            
            M = Sol.MNew;
            [ U, W, V ] = svd(M);
            % remove every sigular value except the biggest three to
            % remove noise and make the Matrix rank 3
            for i = 4:min(size(W))
                W(i,i) = 0;
            end
            % new rank three matrix M3
            M3 = U*W*V';
            Sol.MNew = M3;
            % choose an initial factorization into Light and Shape
            W_root = sqrt(W);
           
            Sol.LNew = U* ( W_root(:, 1:3));
            Sol.SNew = W_root(1:3, :) * V';
            
           
            
            
            %% Solve for a surface
            if numel(Sol.SGT) > 0
                % Use ground truth if given
                [Sol.ANew, Sol.SNew] = getAmb(Sol.SGT,Sol.SNew);            
            else            
                %Resolve Ambiguity up to GBR
                [Sol.SNew, Sol.ANew] = Sol.projectToIntegrableSurface(Sol.SNew);
                Sol.LNew = Sol.LNew /Sol.ANew;
            end
            
            
%             % Uncomment to display the gradient fields that are used to
%             % reconstruct the depth map.
%             x = 1:Sol.imgsize(1);
%             y = 1:Sol.imgsize(2);
%              [Zx, Zy] = Sol.getGradientField(Sol.SNew);
%              figure('Name', 'Zx'); mesh(rot90(Zx,3));
%              xlabel('x');
%              ylabel('y');
%              figure('Name','Zy'); mesh(rot90(Zy,3));
%              xlabel('x');
%              ylabel('y');
%             Sol.ZNew = g2s(Zx', Zy', x',y');
            
            
            Sol.ZNew = Sol.reconstructDepthMap();
            
            result = struct();
            result.breakReason = 'done';
        end
        
        function [M, L, S, R, Z] = extract(Sol)
            M = Sol.MNew;
            L = Sol.LNew;
            SNorms = sqrt(sum(Sol.SNew.^2));
            S = my_normc(Sol.SNew)';
            R = SNorms';
            Z = Sol.ZNew;
        end
 
        function [S, A, C] = projectToIntegrableSurface(Sol,S)
            
            npixels = size(S,2);
            C = zeros(npixels,6);
            for p = 1:size(S,2)
                [x,y] = ind2sub(Sol.imgsize, p);
%                 if (x ~= 1 && y ~= 1 && x ~= Sol.imgsize(1) && y ~= Sol.imgsize(2)) 
                if (x ~= Sol.imgsize(1) && y ~= Sol.imgsize(2)) 
                    s = S(:,p);
%                     s_x = (S(:,sub2ind(Sol.imgsize,x+1,y)) - S(:,sub2ind(Sol.imgsize,x-1,y)))*0.5;
%                     s_y = (S(:,sub2ind(Sol.imgsize,x,y+1)) - S(:,sub2ind(Sol.imgsize,x,y-1)))*0.5;
                    s_x = S(:,sub2ind(Sol.imgsize,x+1,y)) - s;
                    s_y = S(:,sub2ind(Sol.imgsize,x,y+1)) - s;
%                     s_x = S(:,sub2ind(Sol.imgsize,x+1,y));
%                     s_y = S(:,sub2ind(Sol.imgsize,x,y+1));
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
        
        function [Sx,Sy] = getGradientField(Sol, S)
            S = my_normc(S);
            Sx = -reshape(S(1,:) ./ S(3,:), Sol.imgsize);
            Sy = -reshape(S(2,:) ./ S(3,:), Sol.imgsize);
            
        end
        
        
        
        % creates an overdetermined system of linear equations 
        % and solves it using Linear Least Squares
        function Z = reconstructDepthMap(varargin)
            Sol = varargin{1};
            if nargin == 1
                S = Sol.SNew;
            else 
                S = varargin{2};
            end
            
            
            rows = Sol.imgsize(1);
            cols = Sol.imgsize(2);
            
            pixels = rows*cols;
            num_equ = (rows-1)*cols + (cols-1)*rows;
            C = zeros(num_equ, pixels);
            D = zeros(num_equ, 1);
            % x direction
            pos = 0;
            [Sx, Sy] = Sol.getGradientField(S);
            for j=1:cols
                for i=1:rows-1
                    
                    ind = sub2ind(Sol.imgsize,i,j);
                    dx = Sx(i,j);
                    if not(isnan(dx))
                        pos = pos +1;
                        D(pos) = dx;
                        V = zeros(1, pixels);
                        V(1, ind) = -1;
                        V(1, sub2ind(Sol.imgsize,i+1,j)) = 1;
                        C(pos,:) = V;
                    end
                end                
            end
            % y direction
            for i=1:rows
                for j=1:cols-1
                    
                    ind = sub2ind(Sol.imgsize,i,j);                    
                    dy = Sy(i,j);
                    if not(isnan(dy))
                        pos = pos +1;
                        D(pos) = dy;
                        V = zeros(1, pixels);
                        V(1, ind) = -1;
                        V(1, sub2ind(Sol.imgsize,i,j+1)) = 1;
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
        
        
%         function Z_MAT = saveDepthMap(Sol, Z)
%             Z_MAT = reshape(Z,Sol.imgsize );
%             Z_MAT = Z_MAT';
%             %Z_IMG = mat2gray(Z_MAT);
%             %imwrite(Z_IMG, 'solution.jpg'); 
%         end
    end %methods
end %classdef
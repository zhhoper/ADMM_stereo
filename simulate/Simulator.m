%% The simulator class 
% produces surface normals that are pointing outwards and 
% light sources that are shining from behind the object through to get 
% positive values for the Intensities matrix $M = LS$
classdef Simulator < handle
    properties (SetAccess = private)
        numberOfLightSources = 5;
        resolutionX, resolutionY;
        scenario;
        S; % Shape of the object
        M, MNoise; % The measurement matrix (without and with noise added)
        Z;      % Depth Map of the object
        R;      % Albedos of the object
        unknownInds; % The unknown indices
        INC;
       
    end
    properties (SetAccess = public)
        ellipse = [];
        L = [];  % Light Sources
        noiseStrength;
        knowledgeRatio = 1;
        rangeX;%= [-2, 2];
        rangeY;% = [-2, 2];        
    end
    methods
        function Sim = Simulator(lightSources, resolution, scenario, noiseStrength, knowledgeRatio)
            Sim.numberOfLightSources = lightSources;
            Sim.resolutionX = resolution(1);
            Sim.resolutionY = resolution(2);
            Sim.rangeX = [-(Sim.resolutionX-1)/2, (Sim.resolutionX-1)/2];
            Sim.rangeY = [-(Sim.resolutionY-1)/2, (Sim.resolutionY-1)/2];
            Sim.scenario = scenario;
            Sim.noiseStrength = noiseStrength;
            if nargin == 5
                Sim.knowledgeRatio = knowledgeRatio;
            end
        end
        
        %% Main Simulation function
        % Creates all the data for the given scenario, i.e.
        % M,L,S,Z,R,MNoise, INC
        function M = simulate(Sim)
            rows = linspace(Sim.rangeX(1), Sim.rangeX(2),Sim.resolutionX);
            cols = linspace(Sim.rangeY(1), Sim.rangeY(2),Sim.resolutionY);            
            COORDS = createCoordArray(rows,cols);
            m = numel(rows);
            n = numel(cols);
            switch Sim.scenario.object
                
                case 'face'
                    load('/fs/jacobsvideo/Photo-stereo/data1_ps.mat','shape');
                    Zuu=reshape(shape(:,3)/2,[100,100]); Zuu=Zuu-5; 
                    Zuu(Zuu<0)=0;
                    Zuu=Zuu(1:85,1:85);
                    
                    Sim.Z=imresize(Zuu,[Sim.resolutionY,Sim.resolutionY]);                  
                    
                    
                    
                
                case 'new'
                    [X,Y]=meshgrid(1:Sim.resolutionX,1:Sim.resolutionY);%Sim.rangeX(1):Sim.rangeX(2),Sim.rangeY(1):Sim.rangeY(2));
                    
                    Z1=300-((X-20).^2+(Y-20).^2); Z1(Z1<0)=0; 
                    Z1=Z1/5;

                    Z2=50-(X-7).^2-(Y-7).^2; Z2(Z2<0)=0;  
                    Z3=50-(X-33).^2-(Y-7).^2; Z3(Z3<0)=0;  
                    Z4=50-(X-33).^2-(Y-33).^2; Z4(Z4<0)=0;  
                    Z5=50-(X-7).^2-(Y-33).^2; Z5(Z5<0)=0;  
                    Sim.Z=Z1+Z2+Z3+Z4+Z5;
                    
                case 'new1'
                    
                    [X,Y]=meshgrid(1:Sim.resolutionX,1:Sim.resolutionY);
                    
                    Z1=300-((X-15).^2+(Y-10).^2); Z1(Z1<0)=0; 
                    Z1=Z1/5;
                    
                    Z2=300-((X-10).^2+(Y-10).^2); Z2(Z2<0)=0; 
                    Z2=Z2/5;
                    
                    Z3=300-((X-28).^2+(Y-20).^2); Z3(Z3<0)=0; 
                    Z3=Z3/5;
                    
                                       
                    Sim.Z=(Z1+Z2+Z3)/3;
                    
                    
                case 'new2'
                    
                    [X,Y]=meshgrid(1:Sim.resolutionX,1:Sim.resolutionY);
                    
                    Z1=((X-20).^2+(Y-20).^2)-300; Z1(Z1<0)=0; 
                    Z1=Z1/5;
                    
                    Z2=150-((X-10).^2+(Y-10).^2); Z2(Z2<0)=0; 
                    Z2=Z2/2;
                    
                    Z3=150-((X-30).^2+(Y-30).^2); Z3(Z3<0)=0; 
                    Z3=Z3/2;
                    
                    Z4=150-((X-10).^2+(Y-30).^2); Z4(Z4<0)=0; 
                    Z4=Z4/2;
                    
                    Z5=150-((X-30).^2+(Y-10).^2); Z5(Z5<0)=0; 
                    Z5=Z5/2;
                    
                    Sim.Z=(Z1+Z2+Z3+Z4+Z5)/2;
                    
                case 'ellipse'
                    myprintf('Test with everything fixed \n\n');
                    center = diag(myRandom(-10, 10, 3));
                    x = center(1);
                    y = center(2);
                    z = center(3);
                    
                    ca = max((Sim.rangeX(1) - x)^2, (Sim.rangeX(2) - x)^2)*2;
                    cb = max((Sim.rangeY(1) - y)^2, (Sim.rangeY(2) - y)^2)*2;
                    a = myRandom(sqrt(ca),sqrt(ca)*2);
                    b = myRandom(sqrt(cb), sqrt(cb)*2);
                    c = myRandom(30, 75);
                    E = Ellipse(x,y,z, a,b,c);
                    Sim.Z = reshape(E.depth(COORDS),m,n);
                case 'peaks'
                    Zt = peaks(max(m,n)*2);
                    p = round(rand()*m);
                    q = round(rand()*n);
                    Sim.Z = Zt(p+1:p+m,q+1:q+n);
                case 'both'
                    % with random mixture   
                case 'random'
                    Sim.Z = rand(m,n);
                otherwise
                    error('Unknown object param passed in scenario');
            end
            
            if strcmp(class(Sim.scenario.albedo),'double')
                Sim.R = ones(m*n,1) * Sim.scenario.albedo;   % uniform albedo
            elseif strcmp(Sim.scenario.albedo,'random')
                Sim.R = rand(m*n, 1);
            else
                error('Albedo value is not a number and also not random');
            end
            
            Zx = diff(Sim.Z,1);
            Zx = [Zx;Zx(end,:)];
            Zy = diff(Sim.Z,1,2);
            Zy = [Zy,Zy(:,end)];

            Sim.S = [-Zx(:),-Zy(:),ones(m*n,1)];
            Sim.S = my_normr(Sim.S);
%             [nx,ny,nz] = surfnorm(reshape(COORDS(:,1),m,n), reshape(COORDS(:,2),m,n), Sim.Z);
%             Sim.S = -[nx(:),ny(:),nz(:)];
            
%             switch Sim.scenario.albedo
%                 otherwise
%                     Sim.R = 1;
%             end
            
            % if fitted the lightsources are tried to be placed s.t. they don't 
            % cast any shadows in the images.
            switch Sim.scenario.lighting
                case 'fitted'
                    Sim.L = Sim.generateLightSources(Sim.S);
                case 'random'
                    Sim.L = Sim.generateLightSources(0.7);
                otherwise
                    error('Unknown lighting param passed in scenario');
            end
                        
            % Calculate Measurements matrix
            % And set all negative values in the matrix to zero
            Sim.M = Sim.L * Sim.S' * diag(Sim.R);
            M = Sim.M;   
            [Sim.MNoise, Sim.INC] = Sim.addNoise(Sim.M);
        end
        
        % Creates a Lightmatrix with the given number of different light
        % sources. If given a vector of surface normals the generated light
        % sources are tried to be choosen so that every nonzero normal is
        % bigger than zero
        % If a number is given it
        % describes the part of the unit square, the light sources are
        % allowed to lie on. i.e. 1 means only light sources from above
        % and 0 means light sources from anywhere in a halfsphere are
        % allowed.
        function L = generateLightSources(varargin)
            Sim = varargin{1};
            n = Sim.numberOfLightSources;
            L = rand(n,3);
            
            switch numel(varargin{2})
                case 1
                    r = varargin{2};
                    for i = 1:n   
                        x = rand()*2*r - r;
                        ymax = sqrt(r^2-x^2);
                        y = rand()*2*ymax - ymax;
                        dx = [1,0,-x/sqrt(1-x^2 - y^2 )];
                        dy = [0,1,-y/sqrt(1-x^2 - y^2 )];
                        L(i,:) = cross(dx, dy);
                    end
                otherwise              % Try to choose light sources s.t. no shadows exist
                    S = varargin{2};
                    for i=1:n
                        L(i,:) = perceptron(S)';
                        L = my_normr(L);            
                    end
            end
            
            L = my_normr(L);            
            L = L * 255;    
            Sim.L = L;
            
            
        end
        
        %% Add Noise and Missing data
        % adds noise the measurements M by adding to the already existing
        % shadows (missing entries) in M so many, that the desired
        % knowledgeRatio in obtained
        % This function can also be called multiple times if only noise and
        % knowledgeRatio is subject to change, but the generated underlying
        % shape should stay the same.
        function [M, INC] = addNoise(Sim, M)
            if nargin==1
                M = Sim.M;
            end
            N = rand(Sim.numberOfLightSources,Sim.resolutionY*Sim.resolutionX);
            N = 2* (N - 0.5);
            lo = min(min(M));
            hi = max(max(M));
            noiseStrength = Sim.noiseStrength*(hi-lo);
            M = M + N*noiseStrength;
            
            entries = numel(M);
            INC = M>0;
            currInc = sum(sum(INC));
            currMiss = entries - currInc;
            totalMissing = round(entries * (1-Sim.knowledgeRatio));
            
            m = currMiss;
            for i = randperm(entries)
                if m >= totalMissing
                    break;
                end
                if INC(i)
                    m = m+1;
                    INC(i) = 0;
                end
                
            end
            
           
            %% add missing entries
            M = INC.* M;
            Sim.MNoise = M;
        end
        

        
        

    end %methods
end %classdef
        
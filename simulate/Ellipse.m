%% The ellipse class 
% The only thing that is still used is the generation of the depth map
% values. Not really an important class.
classdef Ellipse < handle
   % The following properties can be set only by class methods
   properties (SetAccess = private)
       % store the center point of the ellipse
      x = 0;
      y = 0;
      z = 0;
      % Store the focal values
      a = 0;
      b = 0;
      c = 0;
      myalbedo = {1};
   end
   methods (Static) 
       function E = generateEllipseContainingArea(XInterval, YInterval) 
           center = diag(myRandom(-10, 10, 3));
           x = center(1);
           y = center(2);
           z = center(3);
           ca = max((XInterval(1) - x)^2, (XInterval(2) - x)^2)*2;
           cb = max((YInterval(1) - y)^2, (YInterval(2) - y)^2)*2;
           a = myRandom(sqrt(ca),sqrt(ca)*2);
           b = myRandom(sqrt(cb), sqrt(cb)*2);
           c = myRandom(30, 75);
           E = Ellipse(x,y,z, a,b,c);
       end
       
       function E = generateEllipseInArea(XInterval, YInterval)
           a = (XInterval(2) - XInterval(1)) / 3;
           b = (YInterval(2) - YInterval(1)) / 3;
           x = rand()*a + XInterval(1) + a;
           y = rand()*b + YInterval(1) + b;
           z = rand();
%            r1 = rand() * a + a;
%            r2 = rand() * b + b;
           r1 = max(abs(XInterval(1)-x), abs(x-XInterval(2)));
           r2 = max(abs(YInterval(1)-y), abs(y-YInterval(2)));
           %r3 = 0.5 + rand()*(a+b)/2*3;
           r3 =  rand()*5;
%            x = mean(XInterval);
%            y = mean(YInterval);
%            r1 = (XInterval(2) - XInterval(1)) / 2;
%            r2 = (YInterval(2) - YInterval(1)) / 2;
           E = Ellipse(x,y,z, r1, r2, r3);
%            E = Ellipse(0,0, 0, 1,1,1);
           
       end
   end % Static methods
   methods
       % creates a new ellipse with center point (x,y,z) and focals (a,b,c)
       function E = Ellipse(x,y,z, a,b,c)
           E.x = x;
           E.y = y;
           E.z = z;
           E.a = a;
           E.b = b;
           E.c = c;
       end
       
       % if called with no argument random albedo between 0,1
       % if called with const,  then constant albedo
       % if called with const and vector, then first order function
       % if called with const, vector and matrix, then second order function
       function [] = setAlbedoFunction(varargin) 
           E = varargin{1};
           E.myalbedo = varargin(2:end);
           
           if nargin > 2 && size(E.myalbedo{2},2) ~= 1
               E.myalbedo{2} = E.myalbedo{2}';
           end
       end
       
       
       function A = albedo(varargin)
           E = varargin{1};
           switch nargin
               case 3 % given x, y coordinates
                   x = varargin{2};
                   y = varargin{3};
                   v = [x, y];
                   switch numel(E.myalbedo)
                       case 0
                           A = rand();
                       case 1
                           A = E.myalbedo{1};
                       case 2
                           A = v * E.myalbedo{2} + E.myalbedo{1};
                           A = 1 / (1 + exp(-A));     % Sigmoid function ensures range is between zero and one
                       case 3
                           A = v * E.myalbedo{3} * v' + v * E.myalbedo{2} + E.myalbedo{1};
                           A = 1 / (1 + exp(-A));     % Sigmoid function ensures range is between zero and one
                   end
                   
                   %A = sin(A);
               case 2 % given a matrix containing the coordinates
                   C = varargin{2};
                   s = size(C);
                   if s(2) ~= 2
                       C = C';
                   end
                   s = size(C);
                   A = zeros(s(1),1);
                   for i=1:s(1)
                       x = C(i,1);
                       y = C(i,2);
                       A(i) = E.albedo(x,y);
                   end
           end
       end
       
       % returns the z/depth value at a given point (x,y)
       function z = depth(varargin)
           switch nargin
               case 3
                   E = varargin{1};
                   x = varargin{2};
                   y = varargin{3};
                   
                   z  =E.z + E.c * sqrt( 1 - (x-E.x)^2 / E.a^2 - (y-E.y)^2 / E.b^2);
                   if ~isreal(z)
                       z = 0;
                   end
               case 2 
                   E = varargin{1};
                   C = varargin{2};
                   s = size(C);
                   if s(2) ~= 2
                       C = C';
                   end
                   s = size(C);
                   z = zeros(s(1),1);
                   for i=1:s(1)
                       z(i) = E.depth(C(i,1),C(i,2));
                   end
           end                      
               
       end
            
      
       % returns the normal vector at a given point (x,y) or a Matrix
       % filled with points. The returned normal is zero, if the given
       % point does not correspond to a point on the ellipse.
       function n = normal(varargin)
           switch nargin
               case 3
                   E = varargin{1};
                   x = varargin{2};
                   y = varargin{3};
                   xd = E.depth(x+1,y) - E.depth(x,y);
                   yd = E.depth(x,y+1) - E.depth(x,y);
                   n = -my_normr([xd,yd,-1]);
                   
%                    root = sqrt( 1- (x-E.x)^2 / E.a^2 - (y-E.y)^2 / E.b^2 );
%            
%                    dx = -E.c/E.a^2 * (x-E.x) / root;               
%                    dy = -E.c/E.b^2 * (y-E.y) / root;
%                    if not(isreal(dx) && isreal(dy)) || any(isnan(dx)) || any(isnan(dy))
%                        n = [0,0,0];
%                    else
%                        n = cross([1,0,dx],[0,1,dy]);
%                        n = normr(n);
% 
%                        a = E.albedo(x,y);
%                        if (a >0)
%                            n = n * a;
%                        else
%                            n = 0;
%                        end
% 
%                    end
               case 2
                   E = varargin{1};
                   C = varargin{2};
                   s = size(C);
                   if s(2) ~= 2
                       C = C';
                   end
                   s = size(C);
                   n = zeros(s(1),3);
                   for i=1:s(1)
                       x = C(i,1);
                       y = C(i,2);
                       n(i,:) = E.normal(x,y);
                   end
                   
           end           
       end % function normal
   end % methods
end % classdef
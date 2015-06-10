function [h] = showLights(varargin)
L1 = my_normr(varargin{1});

r = 3;                 %# A radius value

L = L1;
h = figure;
quiver3(-L(:,1)*r, -L(:,2)*r, -L(:,3)*r,L(:,1), L(:,2), L(:,3),0, 'LineWidth', 5, 'Color', 'g');
hold on
if nargin>=2
    L2 = my_normr(varargin{2});
    L = L2;
    quiver3(-L(:,1)*r, -L(:,2)*r, -L(:,3)*r,L(:,1), L(:,2), L(:,3),0, 'LineWidth', 5, 'Color', 'r');
end
if nargin>=3
    L3 = my_normr(varargin{3});
    L = L3;
    quiver3(-L(:,1)*r, -L(:,2)*r, -L(:,3)*r,L(:,1), L(:,2), L(:,3),0, 'LineWidth', 5, 'Color', 'b');
end


xlabel('x');
ylabel('y');


[x,y,z] = sphere;      %# Makes a 21-by-21 point sphere
x = x(11:end,:);       %# Keep top 11 x points
y = y(11:end,:);       %# Keep top 11 y points
z = z(11:end,:);       %# Keep top 11 z points
legend('Ground Truth','Final solution', 'Initial solution');
m = mesh(r.*x,r.*y,r.*z);  %# Plot the surface
set(m,'facecolor','none');
set(m,'edgecolor', [0.5, 0.5, 0.5]);
axis equal;            %# Make the scaling on the x, y, and z axes equal
%axis equal;
hold off;

end

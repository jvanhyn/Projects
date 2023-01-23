function phaseplane(g,y1_interval,y2_interval,density)
% phaseplane(g, [y1min,y1max], [y2min,y2max], density)
%   
%   Plots the phase plane for a two-equation system of first order ODEs
%       Y' = g(t,Y), where Y is the vector [y1; y2].
%   The number of points chosen in the ranges y1min < y1 < y1max and
%       y2min < y2 < y2max is equal to density.
%
%   g is an @ function.
%   density is a positive integer.
%
% Example:  y1' = cos(y2)
%           y2' = -y1
%   Plot a phase plane for y1 in [-5,5], y2 in [-6,4],
%   on a 15x15 grid:
%
%   g = @(t,Y) [cos(Y(2)); -Y(1)];
%   phaseplane(g,[-5,5],[-6,4],15)

y1 = linspace(y1_interval(1),y1_interval(2),density);
y2 = linspace(y2_interval(1),y2_interval(2),density);
[x,y] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));
t=0;
for i = 1:numel(x)
    dY = g(t,[x(i); y(i)]);
    u(i) = dY(1);
    v(i) = dY(2);
end
quiver(x,y,u,v,'Color',[0.13 0.45 0.81]); figure(gcf);
axis tight equal square
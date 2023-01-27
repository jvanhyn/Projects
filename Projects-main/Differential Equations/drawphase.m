function drawphase(g,tmax,y1start,y2start)
% drawphase(g, tmax, y1start, y2start)
%   
%   Plots the phase path for a solution of the two-equation system of
%       first order ODEs Y' = g(t,Y), where Y is the vector [y1; y2].
%   The arguments y1start and y2start give initial values for y1 and y2;
%       that is, they specify y1 and y2 when t=0. The argument tmax is the
%       final time, so that t runs from 0 to tmax.
%
%   g is an @ function.
%
% Example:  y1' = cos(y2)
%           y2' = -y1
%   We want the phase path over the first 50 seconds that starts at the
%   point (y1,y2) = (1,-3).
%
%   g = @(t,Y) [cos(Y(2)); -Y(1)];
%   drawphase(g, 50, 1, -3)

[ts,ys] = ode45(g,[0,tmax],[y1start;y2start]);
axis manual
    plot(ys(:,1),ys(:,2),'LineWidth',1)
    plot(ys(1,1),ys(1,2),'k*','LineWidth',1) % initial point
    plot(ys(end,1),ys(end,2),'ko','LineWidth',1)
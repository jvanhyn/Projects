function drawode(f,t_interval,t0,y0)
% drawode(f, [tmin,tmax], t0, y0)
%   
%   Plot a solution for the first order ODE y' = f(t,y)
%   using t-values from start time tmin to end time tmax
%       with initial value y(t0) = y0
%
%   f is an @ function
%     or the name of an m-file with quotes.
%   t0 need not actually be between tmin and tmax.
%
% Example: y' = (1 - t^2)*exp(y) for -4 < t < 4; y(1)=-1
%
%   f = @(t,y) (1 - t^2)*exp(y)
%   drawode(f, [-4,4], 1, -1)

tmin=t_interval(1);
tmax=t_interval(2);
[ta,ya] = ode45(f,[t0,tmin],y0);
axis manual
hold on
plot(ta,ya,'LineWidth',1)
ax = gca;    % Color matching for R2014b+; remove % to enable
if (ax.ColorOrderIndex - 1) == 0
    ax.ColorOrderIndex = 7;
else
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
end
[ta,ya] = ode45(f,[t0,tmax],y0);
plot(ta,ya,'LineWidth',1)
hold off
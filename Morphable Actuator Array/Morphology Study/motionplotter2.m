clear; clc;

% Meshing
minT = 0; maxT = pi; resT = 100;
minR = 0; maxR = 5 ; resR = 10;
[T,R] = meshgrid(linspace(minT,maxT,resT),linspace(minR,maxR,resR)); 

% Time Vector
stept = 100;
limt = 10;
mint = 0;
t = linspace(mint,limt,stept);

% Wave Properties
w = 11;
zscale = 1;
A = zscale/maxR;


%% Calculate Surface
x = cell(size(t));
y = cell(size(t));
z = cell(size(t));

for k = 1:(length(t)-1)
    phi = pi/4*(sin(w*T+t(k))).*(sin(T).^2);
    x{k} =  [R.*cos(phi).*cos(T) R.*cos(phi).*cos(T)]; 
    y{k} =  [R.*cos(phi).*sin(T) -R.*cos(phi).*sin(T)];
    z{k} =  [A*(sqrt(2))*R.*sin(phi) A*(sqrt(2))*R.*sin(phi)];
end
%% Display Curvature 
figure;
surf(x{1},y{1},z{1});
zlim([-5 5]);
xlim([-5 5]);
ylim([-5 5]);
%% Animate Kinematics 

% Display Propperties 
C = 'red';

% Define Figure
fig = figure;
num = fig.Number;
ax = axes;

% Plot Surface 
for k = 1:(length(t)-1)

    if ishandle(num) == false
      break;
    end 

    surf(ax,x{k},y{k},z{k});
%    axis off
    set(gca,'color',[0 0 0])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])

    zlim(ax,[-5 5])
    xlim(ax,[-5 5])
    ylim(ax,[-5 5])
    drawnow
    pause(0.2)
end
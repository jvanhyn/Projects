%% 2.1 


% Given the complex nature of this function, its much
% easier to see a solution like this, especially when
% computational methods are hard, or impossible.  


%% 2.2
clear 
clc
f = @(x,y) (exp(-x)-y)*(exp(-x)+2+y);
slopefield(f, [-10,10],[-10,10],20,'k',1)
hold on
drawode(f, [-10,10], 2, 3)
drawode(f, [-10,10], 1, -3)
hold off

% If the initial value was (0, -.25), 
% the function would be linear, the 
% solution would otherwise be far more % complex



%% 2.3
clear 
clc
f = @(x,y) x + 2*y ;
slopefield(f, [-5,5],[-5,5], 50)
hold on
drawode(f, [-5,5], 0, -0.25)
drawode(f, [-5,5], 0, -0.33)
drawode(f, [-5,5], 0, -1)
drawode(f, [-5,5], 0, -0.24)
hold off
%% 2.4
clear 
clc
f = @(x,y) y-1 ;
slopefield(f, [-6,6],[-6,6], 50)
hold on
drawode(f, [-6,6], 1,1)
drawode(f, [-6,6], 1,0.9)
drawode(f, [-6,6], 1,1.1)
drawode(f, [-6,6], 1,0.6)
drawode(f, [-6,6], 1,0.8)
drawode(f, [-6,6], 1,1.2)
hold off
%% 2.5
f = @(t,y) 2*(5-y);
slopefield(f, [0,10],[-8,10], 50)
hold on
drawode(f, [0,10], 1,1)
drawode(f, [0,10], 1,6)
drawode(f, [0,10], 1,77)
hold off

% A seems to represent the ambient tempurature. 
%% 2.6
clear
clc
f = @(t,y) 0.4*(293.706-y);

slopefield(f, [-10,25],[0,400], 50)
hold on
drawode(f, [-10 ,25], 0, 252.039)
hold off

% dy/dt = 0.4*(51-y) , y(0) = -6
% A should be 41 degrees, because that is the abient tempurature in the
% fridge. If we consider 39 degreees F to be defrosted, the the chicken should
% be defrosted after t = 3.75. I we had the ambient tempurature at 69
% degrees F, the chicken would defrost after t = 2.4.

%% 2.7
g = @(t,Y) [2; -3];
tf = 1;
phaseplane(g, [-5,5], [-5,5], 25)
hold on
drawphase(g, tf, 1, 3)
drawphase(g, tf, 1, -1)
drawphase(g, tf, -1, 1)
drawphase(g, tf,  1, 1)
hold off

% Changing x' and y' change the slope of each vector.

%% 2.8
g = @(t,Y) [Y(2); -Y(1)];
tf = 50;
phaseplane(g, [-5,5], [-5,5], 25)
hold on
drawphase(g, tf, 1, 3)
drawphase(g, tf, 1, 2)
drawphase(g, tf, 1, 1)
hold off

%% 2.9
f = @(x,y) x^2 + y^3 ;
slopefield(f, [-5,5],[-5,5], 30)
drawode(f, [-5,5], 1, 1)
drawode(f, [-5,5], 1, 2)
drawode(f, [-5,5], 1, -1)
figure 
g = @(t,Y) [1; Y(1)^2 + Y(2)^3];
phaseplane(g, [-5,5], [-5,5], 10)

hold on
drawphase(g, tf, 1, 3)
drawphase(g, tf, -1, -.5)
drawphase(g, tf, -1, 1)
drawphase(g, tf, -1, -1)
hold off

%% 2.9
g = @(t,Y) [Y(1)*(1-Y(2)); Y(2)*(Y(1)-1)];
phaseplane(g, [-1,5], [-1,5], 20)
hold on
tf = 50;
drawphase(g, tf, 2, .5)
drawphase(g, tf, 2, 1)
drawphase(g, tf, 2, 4)
hold off




% Meshing
clc 
clear
minT = 0; maxT = 2*pi; resT = 100;
minR = 0; maxR = 5 ; resR = 10;
[T,R] = meshgrid(linspace(minT,maxT,resT),linspace(minR,maxR,resR)); 
CX = R.*cos(T);  % Circular X
CY = R.*sin(T);  % Circular Y
Z0 = 5;
%%
minX = -5; maxX = 5; resX = 100;
minY = -5; maxY = 5; resY = 100;
[X,Y] = meshgrid(linspace(minX,maxX,resX),linspace(minY,maxY,resY)); 
Z0 = 5;
%% Define Surface
F = @(X,Y) Z0-sqrt(X.^2+Y.^2);
CS = F(CX,CY);
S = F(X,Y);
%% Place Initial Point on Cone
p1 = [5 0];
p1x = p1()
%% Render Cone
figure
surf(X,Y, S)
% hold on
% plot3(p1(1),p1(2),F(p1(1),p1(2)), '.k', 'MarkerSize', 50);
% hold off
%% Normal Vector at point
surfnorm(CX,CY,CS)
[NX,NY,NZ]= surfnorm(CX,CY,CS);
NX
NY
NZx


clear 
close all
clc
% Simulation
%% Input
n = 100;
t = linspace(0,60,n)';

dtu = 30;
toff = dtu/60*n;

Ustp = @(t) 2000*(t >= 0) ;
u = Ustp(t) - Ustp(t-dtu);

%% Simulation
k = -33;                            % arbitrary 'real' k
b = 0.05;                           % arbitrary 'real' b

G31_P = tf(k*[0 0 1],[1 b 0]);      % transfer function from voltage to position
G31_V = tf(k*[0 1 0],[1 b 0]);      % transfer function from voltage to velocity 
G31_A = tf(k*[1 0 0],[1 b 0]);      % transfer function from voltage to acceleration 

y1 = lsim(G31_P, u, t);             % simluated position
y2 = lsim(G31_V, u, t);             % simluated velocity
y3 = lsim(G31_A, u, t);             % simluated acceleration


%% Fit to Data
f1 = fit(t(1:toff),y2(1:toff),'poly1');      % Linear Fit for K
f2 = fit(t(toff:end),y2(toff:end),'exp1');   % Exponential Fit for B

K = f1.p1;                                   % calculated K
B = f2.b;                                   % calculated B

%% Plot
figure(2)
sgtitle("Gyrosope Axis 1")

subplot(4,1,1)
plot(t,u,"LineWidth",2)
grid on
grid minor
title('Input')

subplot(4,1,2)
plot(t,y1,"LineWidth",2)
grid on
grid minor
title("Position")

subplot(4,1,3)
plot(t,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")

subplot(4,1,4)
plot(t,y3,"LineWidth",2)
grid on
grid minor
title("Acceleration")

%% Display

disp("Output Parameters:")
disp("K = " + K)
disp("B = " + B)
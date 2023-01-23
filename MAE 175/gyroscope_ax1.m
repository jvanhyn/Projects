clear 
close all
load("stepdata.mat");
clc

%% DATA
S = (step1 + step2 + step3 + step4 + step5)./5;

u = S(:,3);
y1 = S(:,4);
t1 = S(:,2);
y2 = diff(S(:,4));
t2 = t1(1:end-1);
y3 = diff(y2);
t3 = t1(1:end-2);

[M, I] = min(y2);
ton = t2(1:I);
toff = t2(I:end);
%% Fit to Data
f1 = fit(ton,y2(1:I),'poly1');              % Linear Fit for K
f2 = fit(toff,y2(I:end),'exp1');            % Exponential Fit for B

K0 = f1.p1;                                   % calculated K
beta = -f2.b;                                   % calculated B

%% Plot
figure(1)
sgtitle("Gyrosope Axis 1")

subplot(4,1,1)
plot(t1,u,"LineWidth",2)
grid on
grid minor
title('Input')

subplot(4,1,2)
plot(t1,y1,"LineWidth",2)
grid on
grid minor
title("Position")

subplot(4,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")

subplot(4,1,4)
plot(t3,y3,"LineWidth",2)
grid on
grid minor
title("Acceleration")
hold off

%% Display

disp("Output Parameters:")
disp("K = " + K0)
disp("B = " + beta)

%% Model 
G31 = tf(K0,[1 beta 0]);

%% Controller 
K = -100;
Ki = -10;
Kd = -10;

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(Kd*[1 0],1);

D1 = P+I+D
GD = G31*D1;

figure(2)
rlocus(G31*D1)

figure(3)
margin(GD)

H = GD/(1+GD);

[Gm,Pm,Wcg,Wcp] = margin(GD);
[mag,phase,wout] = bode(GD,{0,5000});
S = stepinfo(H,'RiseTimeLimits',[0 1]);

eta = Pm/100
RT = S.RiseTime
MP = S.Overshoot

step(H)

%% Simulation
t = 0.009*(1:678);

dtu = 3;
Ustp = @(t) 2000*(t >= 0);

u = Ustp(t) - Ustp(t-dtu);

G31_P = tf(K0*[0 0 1],[1 beta 0]);      % transfer function from voltage to position
G31_V = tf(K0*[0 1 0],[1 beta 0]);      % transfer function from voltage to velocity 
G31_A = tf(K0*[1 0 0],[1 beta 0]);      % transfer function from voltage to acceleration 

y1 = lsim(G31_P, u, t);             % simluated position
y2 = lsim(G31_V, u, t);             % simluated velocity
y3 = lsim(G31_A, u, t);         % simluated acceleration

figure(4)
sgtitle("Gyrosope Axis 1")

subplot(4,1,1)
plot(t,u,"LineWidth",2)
grid on
grid minor
title('Input')
xlim([0 6]);

subplot(4,1,2)
plot(t,y1,"LineWidth",2)
grid on
grid minor
title("Position")
xlim([0 6]);

subplot(4,1,3)
plot(t,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6]);

subplot(4,1,4)
plot(t,y3,"LineWidth",2)
grid on
grid minor
title("Acceleration")
xlim([0 6]);

%%
%% Simulation
t = 0.009*(1:678);

dtu = 3;
Ustp = @(t) 2000*(t >= 0);

u = Ustp(t) - Ustp(t-dtu);

G31_P = tf(K0*[0 0 1],[1 beta 0]);      % transfer function from voltage to position
G31_V = tf(K0*[0 1 0],[1 beta 0]);      % transfer function from voltage to velocity 
G31_A = tf(K0*[1 0 0],[1 beta 0]);      % transfer function from voltage to acceleration 

y1 = lsim(H*tf([0 0 1],1), u, t);             % simluated position
y2 = lsim(H*tf([0 1 0],1), u, t);             % simluated velocity
%y3 = lsim(H*tf([1 0 0],1), u, t);         % simluated acceleration

figure(5)
sgtitle("Gyrosope Axis 1")

subplot(3,1,1)
plot(t,u,"LineWidth",2)
grid on
grid minor
title('Input')
xlim([0 6]);

subplot(3,1,2)
plot(t,y1,"LineWidth",2)
grid on
grid minor
title("Position")
xlim([0 6]);

subplot(3,1,3)
plot(t,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6]);




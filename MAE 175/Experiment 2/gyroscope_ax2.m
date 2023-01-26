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
y3 = lowpass(diff(y2),0.1);
t3 = t1(1:end-2);

%% Figure 1
figure(1)
sgtitle("LAB DATA - AX2", "FontName","Cambria")

subplot(3,1,1)
plot(t1,u,"LineWidth",2)
grid on
grid minor
title('Input')
xlim([0 6])
subplot(3,1,2)
plot(t1,y1,"LineWidth",2)
grid on
grid minor
title("Position")
xlim([0 6])

subplot(3,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6])
%% Predictions 
% yon = y2(1:)
% [M I] = max(yw);
% RT = t2(I);

%% Simulation Prediction

Khw = 9.7656e-03;        % Hardware Gain
t = 0.009*(0:677);       % Simulation Time
dtu = 3;                 % Step Length

wn=2*pi*2;
beta=0.1;
K=750;


Ustp = @(x) (x > 0) - 0*(x>dtu);          % Step function
u = (Ustp(t) - Ustp(t-dtu));            % Simulation Input


G22 = tf(K*wn^2,[1 2*beta*wn wn^2]);      % Transfer function from Voltage to Position


y1 = lsim(G22, u, t);           % Simluated position using LinSim
t1 = t;   

y2 = diff(y1);                          % Velocity
t2 = t1(1:end-1);

y3 = diff(y1,2);                        % Acceleration
t3 = t1(1:end-2);


figure(2)
sgtitle("SYM DATA - AX2 PREDICTION", "FontName","Cambria")

subplot(3,1,1)
plot(t1,u,"LineWidth",2)
grid on
grid minor
title('Input')
xlim([0 6])
subplot(3,1,2)
plot(t1,y1,"LineWidth",2)
xlim([0 6])
grid on
grid minor
title("Position")
subplot(3,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6])
%% 
[M, I] = max(y1);
trise= t(I);
m = mean(y1(1:334));
[m,i] = min(abs(((y1-m)/M)-0.01));
tsettle = t(i);


hold on 
subplot(3,1,2)
xline(tsettle)
xline(trise)
hold off

wn = 1.8/trise*2;
beta = 4.6/tsettle/wn;
K = y1(344);


%% 77777
G22 = tf(K*wn^2,[1 2*beta*wn wn^2]);      % Transfer function from Voltage to Position


y1 = lsim(G22, u, t);           % Simluated position using LinSim
t1 = t;   

y2 = diff(y1);                          % Velocity
t2 = t1(1:end-1);

y3 = diff(y1,2);                        % Acceleration
t3 = t1(1:end-2);

%% Figure 2
figure(3)
sgtitle("SYM DATA - AX2 VALIDATION", "FontName","Cambria")

subplot(3,1,1)
plot(t1,u,"LineWidth",2)
grid on
grid minor
title('Input')
xlim([0 6])
subplot(3,1,2)
plot(t1,y1,"LineWidth",2)
grid on
grid minor
title("Position")
xlim([0 6])

subplot(3,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6])

st = step(G22);
S = stepinfo(G22,'RiseTimeLimits',[0 1]);
eta = Pm/100;
RT = S.RiseTime;
MP = S.Overshoot;
ST = S.SettlingTime;

hold on 
subplot(3,1,2)
xline(RT)
xline(ST)
hold off

%% Controller Design 
figure()
rlocus(G22)

a = 1;
b = 11;

% D1 = tf(1,[1 0]);
% D2 = tf([1 a],1);
% D3 = tf([1 b],1);
% 
% D4 = D1*D2*D3;
% 
% GD1 = G22*D4;
figure()
rlocus(GD1)

figure
margin(GD1)
%%

K = 0.0001;
Ki = 1
Kd = 1000

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(Kd*[1 0],1);


D1 = (P+I+D);
GD = G22*D1;
% 
figure(1)
hold on
margin(D1)

figure(1)
margin(G22)
% hold off
figure(2)
margin(GD)


H = GD/(1+GD);

[Gm,Pm,Wcg,Wcp] = margin(GD);
[mag,phase,wout] = bode(GD,{0,5000});
S = stepinfo(H,'RiseTimeLimits',[0 1]);

eta = Pm/100;
RT = S.RiseTime;
MP = S.Overshoot
figure(5)
hold on
step(H)
hold off

%% Simulation
u = sin(t);
y1 = lsim(H*tf([0 0 1],1), u, t);             % simluated position
y2 = lsim(H*tf([0 1 0],1), u, t);             % simluated velocity

figure(4)
sgtitle("SYM DATA - AX2 PID", "FontName","Cambria")

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

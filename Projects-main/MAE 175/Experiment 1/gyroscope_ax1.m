clear 
close all
load("week2data.mat");
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

[M, I] = min(y2);
ton = t2(1:I);
toff = t2(I:end);
%% Fit to Data
f1 = fit(ton,y2(1:I),'poly1');              % Linear Fit for K

K0 = f1.p1;                                 % calculated K
beta = 0;                                   % guess at B

%% Plot
figure(1)
sgtitle("LAB DATA - AX1", "FontName","Cambria")

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
ylim([-4 0]*10^4)
xlim([0 6])


subplot(3,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
ylim([-110 0])
xlim([0 6])

% subplot(4,1,4)
% plot(t3,y3,"LineWidth",2)
% grid on
% grid minor
% title("Acceleration")
% ylim([-0.5 0.5])
% xlim([0 6])


%% Display

disp("Output Parameters:")
disp("K = " + K0)
disp("B = " + beta)

%% Model 
G31 = tf(K0,[1 beta 0]);

%% Validation
Khw = 9.7656e-03;        % Hardware Gain
t = 0.009*(0:677);       % Simulation Time
dtu = 3;                 % Step Length

Ustp = @(x) (x > 0) - (x>dtu);          % Step function
u = (Ustp(t) - Ustp(t-dtu));            % Simulation Input


G31_P = tf(K0*[0 0 1],[1 beta 0]);      % Transfer function from Voltage to Position


y1 = 1/Khw*lsim(G31_P, u, t);           % Simluated position using LinSim
t1 = t;   

y2 = diff(y1);                          % Velocity
t2 = t1(1:end-1);

y3 = diff(y1,2);                        % Acceleration
t3 = t1(1:end-2);

%% Figure 2
figure(2)
sgtitle("SYM DATA - AX1", "FontName","Cambria")

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
ylim([-4 0]*10^4)


subplot(3,1,3)
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title("Velocity")
xlim([0 6])


% subplot(4,1,4)
% plot(t3,y3,"LineWidth",2)
% grid on
% grid minor
% title("Acceleration")
% xlim([0 6])
% ylim([-0.5 0.5])
% hold off 

%% Controller 
K = -100;
Ki = -10;
Kd = -10;

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(Kd*[1 0],1);

D1 = (P+I+D);
GD = G31*D1;


H = GD/(1+GD);

[Gm,Pm,Wcg,Wcp] = margin(GD);
[mag,phase,wout] = bode(GD,{0,5000});
S = stepinfo(H,'RiseTimeLimits',[0 1]);

eta = Pm/100;
RT = S.RiseTime;
MP = S.Overshoot


%% Simulation

y1 = lsim(H*tf([0 0 1],1), u, t);             % simluated position
y2 = lsim(H*tf([0 1 0],1), u, t);             % simluated velocity
%y3 = lsim(H*tf([1 0 0],1), u, t);         % simluated acceleration

figure(4)
sgtitle("SYM DATA - AX1 PID", "FontName","Cambria")

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

%%
figure()
subplot(1,3,[1 2]);
margin(GD)

subplot(1,3,3)
step(H)


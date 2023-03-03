clear 
close all
load("stepdata.mat");
load("PIDdata.mat")
clc

%% DATA
S = (step1 + step2 + step3 + step4 + step5)./5;

u1 = S(:,3);
y1 = S(:,4);
t1 = S(:,2);
dy1 = diff(S(:,4));
dt1 = t1(1:end-1);
y3 = lowpass(diff(dy1),0.1);
t3 = t1(1:end-2);

[M, I] = min(dy1);
ton = dt1(1:I);
toff = dt1(I:end);
%% Fit to Data
f1 = fit(ton,dy1(1:I),'poly1');              % Linear Fit for K
K0 = f1.p1;                                  % calculated K
beta = 0.0001;                               % guess at Beta (very small)
%% Display
disp("Output Parameters:")
disp("K = " + K0)
disp("B = " + beta)
%% Model 
G31 = tf(K0,[1 beta 0]);
%% Validation / Simulation
t2 = linspace(min(t1),max(t1),length(t1));        % Simulation Time
dtu = 3;                                          % Dwell Time 

Ustp = @(x) (x > 0);                              % Step function
us = (Ustp(t2) - Ustp(t2-dtu));                   % Simulation Input
u2 = u1/8000*pi;                                  % Conversion from encoder counts to radians

G31_P = tf(K0*[0 0 1],[1 beta 0]);                % Transfer function from Encoder Count to Position

ys = lsim(G31_P, u2, t2);                         % Simluated position using LinSim

dys = diff(ys);                                   % Simulated Velocity
dt2 = t2(1:end-1);                                % Corresponding Time Vector                                                

%% Hardware Gain
Khw = K0/(min(dys)/dtu);                 % Estimated Hardware Gain
y2 = Khw*ys;                             % Model Position
dy2 = Khw*dys;                           % Model Velocity

%% Figure 1
fig1 = figure(1)
tcl1 = tiledlayout(3,1);

nexttile
plot(t1,u1,"LineWidth",2)
grid on
grid minor
title('Input',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
xpos = get(yl,'Pos');

nexttile
plot(t1,y1,"LineWidth",2)
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos(1) pos(2) pos(3)]);


nexttile
plot(dt1,dy1,"LineWidth",2)
grid on
grid minor
title('Velocity' ,"FontName","Cambria")
ylabel('Counts / s','Fontweight','bold','FontName',"Cambria")
yl = ylabel('Counts / s','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos(1) pos(2) pos(3)]);


title(tcl1,'G31 – ENCODER3 DATA', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig1,"fig1.fig");
%% Figure 2
lbl = "K_0 = " + K0;
figure(2)
close
fig2 = figure(2);
hold on
plot(dt1,dy1,"LineWidth",2)

plot(t1,f1.p1*t1,"r-.","LineWidth",2)
grid on
grid minor
xlim([0,dtu])
legend("AX 3",lbl)
ylabel("Counts / Second","Fontweight","bold","FontName","Cambria")
xlabel("Time (s)","Fontweight","bold","FontName","Cambria")
title("Experiment 1 – Parameter Fit ")
subtitle("Velocity vs Time")
hold off
%% Figure 3
clear tcl1 
close
fig3 = figure(3)
tcl1 = tiledlayout(3,1);

nexttile
plot(t2,u1,"LineWidth",2)
grid on
grid minor
title('Input',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
xpos2 = get(yl,'Pos');

nexttile
plot(t2,y2,"LineWidth",2)
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos2(1) pos(2) pos(3)]);


nexttile
plot(dt2,dy2,"LineWidth",2)
grid on
grid minor
title('Velocity' ,"FontName","Cambria")
ylabel('Counts / s','Fontweight','bold','FontName',"Cambria")
yl = ylabel('Counts / s','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos2(1) pos(2) pos(3)]);


title(tcl1,'G31 – Simulated Responce', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig3,"fig3.fig");
%% Figure 4
figure(4)
close
clear yl
fig4 = figure(4)
tcl1 = tiledlayout(2,1);

nexttile
hold on
plot(t1,y1,"LineWidth",2)
plot(t2,y2,"LineWidth",2)
xlim([0 3])
hold off
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
xpos3 = get(yl,'Pos');
xlim([0 3])


nexttile
hold on
plot(dt1,dy1,"LineWidth",2)
plot(dt2,dy2,"LineWidth",2)
xlim([0 3])
hold off
grid on
grid minor
title('Velocity',"FontName","Cambria")
yl = ylabel('Counts / Second','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos3(1) pos(2) pos(3)]);


title(tcl1,'G31 – Model Validation', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
legend(tcl1.Children(end),'AX 3','Model')
savefig(fig4,"fig4.fig");

%% Controller 
K = -30;
ki = 0.1; % Integral Gain Fraction
kd = 0.1; % Derivative Gain Fraction 


P = tf(1,1);
I = tf(ki,[1 0]);
D = tf(kd*[1 0],1);

PD =  (P+D);
PID = (P+I+D);

GD1 = G31_P*K;
GD2 = G31_P*PD;
GD3 = G31_P*PID;

%% Model Analysis
figure(5)
rlocusplot(G31);
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
%% Controller Design 1
rlocusplot(GD1)
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
%% Controller Design 2
rlocusplot(GD2)
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
%% Controller Design 3
rlocusplot(GD3)
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
%% Performance 1
figure(9)
margin(GD1)
%% Performance 1
margin(K*GD3)
%% Closed Loop
H = GD3/(1+K*GD3);

%% Pole Placement

%% Closed Loop Responce 

figure(6)
rlocus(GD)



%%
% [Gm,Pm,Wcg,Wcp] = margin(GD);
% [mag,phase,wout] = bode(GD,{0,5000});
% S = stepinfo(H,'RiseTimeLimits',[0 1]);
% 
% eta = Pm/100;
% RT = S.RiseTime;
% MP = S.Overshoot
% 
% figure(4)
% subplot(1,3,[1 2]);
% margin(GD)
% 
% subplot(1,3,3)
% step(H)
% 

%%

S = (control1 + control2 + control3)./3;
u1 = S(:,3);
y1 = S(:,4);
t1 = S(:,2);
dy1 = diff(S(:,4));
t2 = t1(1:end-1);
y3 = lowpass(diff(dy1),0.1);
t3 = t1(1:end-2);

[M, I] = min(dy1);
ton = t2(1:I);
toff = t2(I:end);

figure(5)
sgtitle("LAB DATA - PID AX1", "FontName","Cambria")

subplot(3,1,1)
plot(t1,u1,"LineWidth",2)
grid on
grid minor
title('Input')

subplot(3,1,2)
plot(t1,y1,"LineWidth",2)
grid on
grid minor
title("Position")


subplot(3,1,3)
plot(t2,dy1,"LineWidth",2)
grid on
grid minor
title("Velocity")

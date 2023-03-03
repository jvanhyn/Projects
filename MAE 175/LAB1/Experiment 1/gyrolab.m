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
dy1 = diff(S(:,4))./mean(diff(t1));
dt1 = t1(1:end-1);

%% Trimming
[M, i] = min(dy1);
ton = dt1(1:i);
toff = dt1(i:end);
%% Fit to Data
f1 = fit(ton,dy1(1:i),'poly1');              % Linear Fit for K
K0 = f1.p1/2000;                                  % calculated K
beta =  0.01;                               % guess at Beta (very small)
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
u2 = u1;                                   % Conversion from encoder counts to radians

y2 = lsim(G31, u2, t2);                           % Simluated position using LinSim

dy2 = diff(y2)./diff(t2');                                   % Simulated Velocity
dt2 = t2(1:end-1);                                % Corresponding Time Vector                                                


%% Controller 
Khw = 1;
K = -0.1;
ki = 0.01; % Integral Gain Fraction
kd = 25; % Derivative Gain Fraction 

Ki = ki*K; % Integral Gain 
Kd = kd*K; % Derivative Gain 

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(Kd*[1 0],1);

PD =  (P+D);
PID = Khw*(P+I+D);

GD1 = G31;
GD2 = G31*PD;
GD3 = G31*PID;

%% Performance 
figure(9)
margin(GD3)
%% Close Loop Step Responce
H = feedback(GD3,1);
stp = lsim(H,u1,t2);
S = stepinfo(H)
dstp = diff(stp);
step(H)

%% Model Analysis
figure(5)
rlocusplot(G31);
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
%% Controller Design 1
rlocusplot(-GD1)
  axIm = findall(gcf,'String','Imaginary Axis (seconds^{-1})');
    axRe = findall(gcf,'String','Real Axis (seconds^{-1})');
    set(axIm,'String','Imaginary Axis');
    set(axRe,'String','Real Axis');
    %%
    hold on
    k = -0.1;
    plot(pole(k*GD1/(1+k*GD1)),'x')
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

%% Results 
C = (control1 + control2 + control3)./3;

u3 = C(:,3);
y3 = C(:,4);
t3 = C(:,2);
dy3 = diff(C(:,4));
dt3 = t3(1:end-1);





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


title(tcl1,'\bf G31 – Step Responce', "FontName","Cambria")
subtitle(tcl1,'ENCODER3 Data vs. Time', "FontName","Cambria")
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
title("G31 – Parameter Fit ","FontName","Cambria")
subtitle("Velocity vs Time","FontName","Cambria")
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


title(tcl1,'\bfG31 – Simulated Step Responce', "FontName","Cambria")
subtitle(tcl1,'ENCODER3 Data vs. Time', "FontName","Cambria")
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
set(yl,'Pos',[xpos3(1) - 0.1 pos(2) pos(3)]);


title(tcl1,'\bf G31 – Model Validation', "FontName","Cambria")
subtitle(tcl1,'ENCODER3 Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
legend(tcl1.Children(end),'ENCODER3','Model')
savefig(fig4,"fig4.fig");
%%
figure(5)
plot(t2,stp,'LineWidth',2)
grid on
grid minor
xlabel('Time (s)','FontWeight','bold')
ylabel('Counts','FontWeight','bold')
title("G31 - Simulated Controller Step Responce","FontName","Cambria")
subtitle("Position Vs Time","FontName","Cambria")

%%
figure(6)
close
clear yl
fig4 = figure(5)
tcl1 = tiledlayout(3,1);

nexttile
hold on
plot(t3,u3,"LineWidth",2)
hold off
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
xpos3 = get(yl,'Pos');

nexttile
hold on
plot(t3,y3,"LineWidth",2)
hold off
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos3(1) pos(2) pos(3)]);

nexttile
hold on
plot(dt3,dy3,"LineWidth",2)
hold off
grid on
grid minor
title('Velocity',"FontName","Cambria")
yl = ylabel('Counts / Second','Fontweight','bold','FontName',"Cambria");
pos = get(yl,'Pos');
set(yl,'Pos',[xpos3(1) pos(2) pos(3)]);

title(tcl1,'\bf G31 – Closed-Loop Control Results', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig4,"fig4.fig");
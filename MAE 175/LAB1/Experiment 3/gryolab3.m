%% DATA
clear 
close all
load("week3data.mat");
clc

S = week3part1(1:226,:,:,:);
u1 = S(:,3);
y1 = smooth(S(:,4));
t1 = S(:,2);
dy1 = diff(S(:,4))./diff(t1);
dt1 = t1(1:end-1);


%% Fit to Data
f1 = fit(t1,y1,'poly1');              
y1 = y1-f1.p1*t1;

S = stepinfo(dy1,dt1,"RiseTimeThreshhold",[0 1],'SettlingTimeThreshold',0.2);
ST = S.SettlingTime;
yss = 281;

[pks,locs] = findpeaks(dy1,dt1);
I = 4;
pks = pks(1:I);
locs = locs(1:I);
[pks1,locs1] = findpeaks(y1,t1);
pks1 = pks1(1:I);
locs1 = locs1(1:I);

wd = 2*pi*(length(pks)-1)/(locs(end)-locs(1));

bwn = 1/(locs1(end)-locs1(1))*log((pks1(1)-yss)/(pks1(end)-yss));

wn = sqrt(wd^2+bwn^2);
z2 = bwn/wn;
K1 = yss;


%% Model 
G22 = 1/2000*tf(K1*wn^2,[1 2*z2*wn wn^2]); 
%% Validation / Simulation
t2 = linspace(0, 3, 1000);        % Simulation Time
dtu = 3;                                          % Dwell Time 

Ustp = @(x) (x > 0);                              % Step function
us = (Ustp(t2) - Ustp(t2-dtu));                   % Simulation Input
u2 = 2000*us;                                   % Conversion from encoder counts to radians

y2 = lsim(G22, u2, t2);                           % Simluated position using LinSim

dy2 = diff(y2)./diff(t2');                                   % Simulated Velocity
                              % Corresponding Time Vector                                                
dt2 = t2(1:end-1);  


%% Controller 
K = 100;
ki = 1; % Integral Gain Fraction
kd = 10; % Derivative Gain Fraction 

Ki = ki*K; % Integral Gain 
Kd = kd*K; % Derivative Gain 

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(kd*[1 0],1);

PD =  (P+D);
PID = (P+I+D);

GD1 = G22;
GD2 = G22*PD;
GD3 = G22*PID;

%% Close Loop
H = feedback(GD3,1);
step(H,t2)
%% Performance 
figure(9)
margin(GD3)
%%

%% Step Respnce
stp = lsim(H,u2,t2);
dstp = diff(stp);
plot(t2,stp)
xlim([0 3])


%% Model Analysis
figure(5)
rlocusplot(G22);
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
  axIm = findall(gcf,'St2ring','Imaginary Axis (seconds^{-1})');
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
C = week2part2;
u3 = C(:,3);
y3 = C(:,4);
t3 = C(:,2);
dy3 = diff(C(:,4));
dt3 = t3(1:end-1);





%% Figure 1
fig1 = figure(1);
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

figure(2)
close
fig2 = figure(2);
hold on
plot(dt1,dy1,"LineWidth",2)
%plot(locs,pks,"r.","LineWidth",1)
yyaxis right
plot(dt1,sin(wd*dt1),"r-.","LineWidth",1)
ylim([-5 5])
yyaxis left
plot(dt1,281*(wn/sqrt(1-z2^2))*exp(-bwn*(dt1-locs(1))),"g-.","LineWidth",1)
grid on
grid minor
xlim([0,dtu])
l = legend('$$\bf AX 3$$',"$$\bf sin( \omega_d t)$$","$$\bf K({\omega_n} / \sqrt{1-\zeta^2}) e^{-\beta \omega_n t}$$")
set(l,'Interpreter','latex') 
ylabel("Counts / Second","Fontweight","bold","FontName","Cambria")
xlabel("Time (s)","Fontweight","bold","FontName","Cambria")
title("Experiment 2 – Parameter Fit ")
subtitle("Velocity vs Time")
hold off
savefig(fig2,"fig2.fig");

%% Figure 3
figure(3)
clear tcl1 
close
fig3 = figure(3)
tcl1 = tiledlayout(3,1);

nexttile
plot(t2,u2,"LineWidth",2)
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


title(tcl1,'G22 – Simulated Responce', "FontName","Cambria")
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


title(tcl1,'G22 – Model Validation', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
legend(tcl1.Children(end),'AX 3','Model')
savefig(fig4,"fig4.fig");
%%
figure(5)
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

title(tcl1,'G22 – Model Validation', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig5,"fig5.fig");
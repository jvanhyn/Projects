%% Import Data
clear; close all; clc; 


load("STEPRESPONCES")
load("STEPRESPONCES1")
load("STEPRESPONCES2")
load("STEPRESPONCES3")

t_PP = (STEPRESPONCES(:,1));
t_YY = (STEPRESPONCESS2(:,1));

dt_PP = diff(t_PP)
dt_YY = diff(t_YY)
STEP_PP = STEPRESPONCES(:,4)+40;
STEP_PY = STEPRESPONCESS1(:,4);
STEP_YP = STEPRESPONCESS2(:,4)+40;
STEP_YY = STEPRESPONCESS3(:,4);

dSTEP_PP = smooth(diff(smooth(STEPRESPONCES(:,4),30))./dt_PP,30);
dSTEP_PY = smooth(diff(smooth(STEPRESPONCESS1(:,4),30))./dt_PP,30);
dSTEP_YP = smooth(diff(smooth(STEPRESPONCESS2(:,4),30))./dt_YY,30);
dSTEP_YY = smooth(diff(smooth(STEPRESPONCESS3(:,4),30))./dt_YY,30);


%%
m1 = 50;
m2 = 1;
m3 = 10;
m4 = 1; 
n1 = 190;
n2 = 180;
n3 = 70;
u = 10;

g = fittype('a*x-b*exp(-c*x)+d');
f1 = fit(t_PP(m1:n1),STEP_PP(m1:n1),g,'StartPoint',[[ones(size(t_PP(m1:n1))), -exp(-t_PP(m1:n1))]\STEP_PP(m1:n1); 1;m1]);
f2 = fit(t_PP(m2:n1),STEP_PY(m2:n1),g,'StartPoint',[[ones(size(t_PP(m2:n1))), -exp(-t_PP(m2:n1))]\STEP_PY(m2:n1); 1;m2]);
f3 = fit(t_YY(m3:n2),STEP_YP(m3:n2),g,'StartPoint',[[ones(size(t_YY(m3:n2))), -exp(-t_YY(m3:n2))]\STEP_YP(m3:n2); 1;m3]);
f4 = fit(t_YY(m4:n3),STEP_YY(m4:n3),g,'StartPoint',[[ones(size(t_YY(m4:n3))), -exp(-t_YY(m4:n3))]\STEP_YY(m4:n3); 1;m4]);

B_P = f1.c

B_Y = f4.c

A_PP = f1.a*B_P/100
A_PY = f2.a*B_P/100
A_YP = f3.a*B_Y/100
A_YY = f4.a*B_Y/100


syms x
eq1 = f1.a*x-f1.b*exp(-f1.c*x)+f1.d;
eq2 = f2.a*x-f2.b*exp(-f2.c*x)+f2.d;
eq3 = f3.a*x-f3.b*exp(-f3.c*x)+f3.d;
eq4 = f4.a*x-f4.b*exp(-f4.c*x)+f4.d;

deq1 = diff(f1.a*x-f1.b*exp(-f1.c*x)+f1.d);
deq2 = diff(f2.a*x-f2.b*exp(-f2.c*x)+f2.d);
deq3 = diff(f3.a*x-f3.b*exp(-f3.c*x)+f3.d);
deq4 = diff(f4.a*x-f4.b*exp(-f4.c*x)+f4.d);

%% 
fig1 = figure(1);
tcl1 = tiledlayout(2,2);

n = 200

nexttile
hold on
plot(t_PP(1:n),STEP_PP(1:n),"LineWidth",2)
fplot(@(x) f1.a*x-f1.b*exp(-f1.c*x)+f1.d,'o')
grid on
grid minor
title('Step On Pitch Motor',"FontName","Cambria")
subtitle('Pitch')
yl = ylabel('Pitch Angle \rm (deg)','Fontweight','bold','FontName',"Cambria");
xlim([0 1])
ylim([-1 60])

nexttile
ax = gca;
hold on
plot(t_YY(1:n),STEP_YP(1:n),"LineWidth",2)
fplot(@(x) f3.a*x-f3.b*exp(-f3.c*x)+f3.d,'o')
xlim([0 1])
ylim([-1 60])
grid on
grid minor
title('Step On Yaw Motor',"FontName","Cambria")
subtitle("Pitch","FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");

nexttile
hold on
plot(t_PP(1:n),STEP_PY(1:n),"LineWidth",2)
fplot(@(x) f2.a*x-f2.b*exp(-f2.c*x)+f2.d,'o')
xlim([0 1])
ylim([-1 60])
grid on
grid minor
subtitle('Yaw',"FontName","Cambria")
yl = ylabel('Yaw Angle \rm (deg)','Fontweight','bold','FontName',"Cambria");


nexttile
hold on
plot(t_YY(1:n),STEP_YY(1:n),"LineWidth",2)
fplot(@(x) f4.a*x-f4.b*exp(-f4.c*x)+f4.d,'o')
xlim([0 1])
ylim([-1 60])
grid on
grid minor
subtitle("Yaw","FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");


title(tcl1,'\bf Step Responce - Position', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig1,"fig1.fig");
legend(ax,"Measured Position","Data Fit")

%% 
fig2 = figure(2)
tcl1 = tiledlayout(2,2);

n = 200

nexttile
hold on
plot(t_PP(1:n),dSTEP_PP(1:n),"LineWidth",2)
grid on
grid minor
title('Step on Pitch Motor',"FontName","Cambria")
subtitle("Pitch","FontName","Cambria")
yl = ylabel('Yaw Angle \rm (deg)','Fontweight','bold','FontName',"Cambria");

nexttile
hold on
plot(t_YY(1:n),dSTEP_PY(1:n),"LineWidth",2)
grid on
grid minor
subtitle("Yaw","FontName","Cambria")

yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");

nexttile
plot(t_PP(1:n),dSTEP_PY(1:n),"LineWidth",2)
grid on
grid minor
title('Step on Yaw Motor',"FontName","Cambria")
subtitle("Pitch","FontName","Cambria")

yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");

nexttile
plot(t_YY(1:n),dSTEP_YY(1:n),"LineWidth",2)
grid on
grid minor
subtitle("Yaw","FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");


title(tcl1,'\bf Step Responce - Velocity', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig1,"fig1.fig");
%%
m = 50;
g = fittype('a*x-b*exp(-c*x)+d');
f1 = fit(t_PP(m:200),STEP_PP(m:200),g,'StartPoint',[[ones(size(t_PP(m:200))), -exp(-t_PP(m:200))]\STEP_PP(m:200); 1;m]);
f2 = fit(t_PP(m:200),STEP_PY(m:200),g,'StartPoint',[[ones(size(t_PP(m:200))), -exp(-t_PP(m:200))]\STEP_PY(m:200); 1;m]);


% f2 = fit(tYY(1255:end),dYY(1255:end),g,'StartPoint',[[ones(size(tYY)), -exp(-tYY)]\dYY; 1]);
%% Modeling
b_p  = 9.2818;
b_y  = 3.4930;

a_pp = 2.3553;
a_py = 0.0789;
a_yp = 0.2406;
a_yy = 0.7909;

A=[0 0 1 0
   0 0 0 1
   0 0 -b_p 0
   0 0 0 -b_y];
B=[0 0
   0 0
   a_pp a_py
   a_yp a_yy];
C=[1 0 0 0
   0 1 0 0];
D=[0 0
   0 0];

G=ss(A,B,C,D);

Gpp = tf(a_pp,[1 b_p 0]);
Gpy = tf(a_py,[1 b_p 0]);
Gyp = tf(a_yp,[1 b_y 0]);
Gyy = tf(a_yy,[1 b_y 0]);

S = [Gpp, Gpy; Gyp, Gyy];
P = tf('s')*S;

t = (linspace(0,1,1000))';
u_step = (t>0);
u = u_step;
u1 = u_step.*[1 0];
u2 = u_step.*[0 1];
uss1 = u_step.*[1 0 0 0 0 0];
uss2 = u_step.*[0 1 0 0 0 0];


Y1 = lsim(P,u1,t);
Y2 = lsim(P,u2,t);

%% Parameter Estimation 
dYP = smooth(diff(smooth(LAB2part1S1(:,4),20))./diff(LAB2part1S1(:,3)),30)*pi/180/20;
dYY = -smooth(diff(LAB2part1S2(:,4)),50)/.01*pi/180/20;
tYP = squeeze(LAB2part1S1(1:end-1,3) - 4.5);
tYY = (linspace(0,tYP(end),length(tYP))-4.5)';

dPP = smooth(diff(smooth(LAB2part1S1(:,4),20))./diff(LAB2part1S1(:,3)),30)*pi/180/20;
dPY = -smooth(diff(LAB2part1S2(:,4)),50)/.01*pi/180/20;
tPP = squeeze(LAB2part1S1(1:end-1,3,:) - 4.5);
tPY = (linspace(0,tYP(end),length(tYP))-4.5)';
%%

g = fittype('a-b*exp(-c*x)');
f1 = fit(tYP(1255:end),dYP(1255:end),g,'StartPoint',[[ones(size(tYP)), -exp(-tYP)]\dYP; 1]);
f2 = fit(tYY(1255:end),dYY(1255:end),g,'StartPoint',[[ones(size(tYY)), -exp(-tYY)]\dYY; 1]);

%%
g = fittype('a-b*exp(-c*x)');
f1 = fit(t,Y1(:,1),g,'StartPoint',[[ones(size(t)), -exp(-t)]\Y1(:,1); 1]);
f2 = fit(t,Y1(:,2),g,'StartPoint',[[ones(size(t)), -exp(-t)]\Y1(:,1); 1]);
f3 = fit(t,Y2(:,1),g,'StartPoint',[[ones(size(t)), -exp(-t)]\Y1(:,1); 1]);
f4 = fit(t,Y2(:,2),g,'StartPoint',[[ones(size(t)), -exp(-t)]\Y1(:,1); 1]);



%%
f1 = fit(t,Y1(:,1)-max(Y1(:,1)),'exp1');
f2 = fit(t,Y1(:,2)-max(Y1(:,2)),'exp1');
f3 = fit(t,Y2(:,1)-max(Y2(:,1)),'exp1');
f4 = fit(t,Y2(:,2)-max(Y2(:,2)),'exp1');

dY1 = diff(Y1(:,1))./diff(t);
dY2 = diff(Y2(:,2))./diff(t);

A_PP = dY1(2,1)
A_YY = dY2(2,1)

B_P = -f1.b
B_Y = -f4.b
A_PY = -f2.a*B_Y
A_YP = -f3.a*B_P


%% PID Pitch Control

% Initaial Yaw PID
Kyp = 1;  % Proportional Gail
Kyi = 0;  % Integral Gain 
Kyd = 0;  % Derivative Gain 

PID_YAW = (tf(Kyp,1) + tf(Kyi,[1 0]) + tf(Kyd*[50 0],[1 50]));

% Pitch PID
CP = 1.1;
Kpp = 25;  % Proportional Gail
Kpi = 0.1;  % Integral Gain 
Kpd = 3;  % Derivative Gain 

PID_PITCH = CP*(tf(Kpp,1) + tf(Kpi,[1 0]) + tf(Kpd*[50 0],[1 50]));

% Closed Loop Equivilant Transfer Function
Tp = feedback(G,PID_YAW,2,2);
Tp = Tp(1,1);

Ty = feedback(G,PID_PITCH,1,1);
Ty = Ty(2,2);

% Closed Loop Control
GDp = Tp*PID_PITCH;
GDy = Ty*PID_YAW;
Hp = feedback(GDp,1);
Hy = feedback(GDy,1);


% Performance
YCp_info = stepinfo(Hp);
PRT1 = YCp_info.RiseTime;
POS1 = YCp_info.Overshoot;
PST1 = YCp_info.SettlingTime;

% Pitch Tuning Analysis 
[Gm_p1,Pm_p1,Wcg_p1,Wcp_p1] = margin(GDp);
s_p1 = step(Hp);
stepplot(Hp)


CP*[Kpp Kpi Kpd]
%% PID Yaw Control

% Revised Yaw PID
CY = 1.8;
Kyp = 10;  % Proportional Gail
Kyi = 0.01;  % Integral Gain 
Kyd = 3;  % Derivative Gain 
CY*[Kyp Kyi Kyd]
PID_YAW = CY*(tf(Kyp,1) + tf(Kyi,[1 0]) + tf(Kyd*[50 0],[1 50]));

% Pitch PID
% Designed Above

% Closed Loop Equivilant Transfer Function
Tp = feedback(G,PID_YAW,2,2);
Tp = Tp(1,1);

Ty = feedback(G,PID_PITCH,1,1);
Ty = Ty(2,2);

% Closed Loop Control
GDp = Tp*PID_PITCH;
GDy = Ty*PID_YAW;
Hp = feedback(GDp,1);
Hy = feedback(GDy,1);

% Performance
YCp_info = stepinfo(Hp);
PRT = YCp_info.RiseTime;
POS = YCp_info.Overshoot;
PST = YCp_info.SettlingTime;

YCy_info = stepinfo(Hp);
YRT = YCp_info.RiseTime;
YOS = YCp_info.Overshoot;
YST = YCp_info.SettlingTime;

% Yaw Tuning Analysis
[MAGp2,PHASEp2,Wp2] = bode(GDp);
[Gm_p2,Pm_p2,Wcg_p2,Wcp_p2] = margin(GDp);
margin(GDp)
s_p2 = step(Hp);
[MAGy2,PHASEy2,Wy2] = bode(GDy);
[Gm_y2,Pm_y2,Wcg_y2,Wcp_y2] = margin(GDp);
s_y2 = step(Hy);

YAW_TUNING = CY*[Kyp Kyi Kyd];
PITCH_TUNING = CP*[Kpp Kpi Kpd];

[NUMp, DENp] = ss2tf(GDp.A,GDp.B,GDp.C,GDp.D);
[NUMy, DENy] = ss2tf(GDy.A,GDy.B,GDy.C,GDy.D);

%% State Feedback




%% Figure 1
fig1 = figure(1)
tcl1 = tiledlayout(2,2);

nexttile
plot(t,Y1(:,1),"LineWidth",2)
grid on
grid minor
title('\bf Step on Pitch Motor', "FontName","Cambria")

yl = ylabel('Pitch Velocity \rm (rad/s)','Fontweight','bold','FontName',"Cambria");
ylim([0 0.3])

nexttile
plot(t,Y2(:,1),"LineWidth",2)
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
title('\bf Step on Yaw Motor', "FontName","Cambria")

% yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([0 0.3])

nexttile
plot(t,Y1(:,2),"LineWidth",2)
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
xlabel('Time \rm (s)','FontWeight','bold')
yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([0 0.3])


nexttile
plot(t,Y2(:,2),"LineWidth",2)
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
% yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([0 0.3])
xlabel('Time \rm (s)','FontWeight','bold')



title(tcl1,'\bf Simulted Step Responce', "FontName","Cambria")
% subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
% xlabel(tcl1,'Time \rm (s)','FontWeight','bold')
subtitle(tcl1,'Dynamic Helicopter Model', "FontName","Cambria")
savefig(fig1,"fig1.fig");
%% Figure 2
fig2 = figure(2)
tcl1 = tiledlayout(2,2);

nexttile
plot(tPP,dPP,"LineWidth",2)
grid on
grid minor
title('\bf Step on Pitch Motor', "FontName","Cambria")
yl = ylabel('Pitch Velocity \rm (rad/s)','Fontweight','bold','FontName',"Cambria");
ylim([-0.1 0.3])

nexttile
plot(tYP,dYP,"LineWidth",2)
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
title('\bf Step on Yaw Motor', "FontName","Cambria")

% yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([-0.1 0.3])

nexttile
plot(tPY,dPY,"LineWidth",2)
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
xlabel('Time \rm (s)','FontWeight','bold')
yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([0 0.3])


nexttile
plot(tYY(1255:end),dYY(1255:end),"LineWidth",2)
hold on
fplot(@(x) f2.a-f2.b*exp(-f2.c*x))
grid on
grid minor
%title('Position',"FontName","Cambria")
%title('\rm Yaw', "FontName","Cambria")
% yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
ylim([0 0.3])
xlabel('Time \rm (s)','FontWeight','bold')




title(tcl1,'\bf Simulted Step Responce', "FontName","Cambria")
% subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
% xlabel(tcl1,'Time \rm (s)','FontWeight','bold')
subtitle(tcl1,'Dynamic Helicopter Model', "FontName","Cambria")
savefig(fig2,"fig2.fig");
%% Figure 2
fig2 = figure(2)
tcl2 = tiledlayout(2,2);
% 
nexttile
plot(squeeze(Wp2),squeeze(MAGp2),"LineWidth",2)
grid on
grid minor
% 
% nexttile
% (squeeze(Wy2),squeeze(MAGy2),"LineWidth",2)
% grid on
% grid minor
% 
% % yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
% ylim([0 0.3])
% 
% nexttile
% plot(t,Y1(:,2),"LineWidth",2)
% grid on
% grid minor
% %title('Position',"FontName","Cambria")
% %title('\rm Yaw', "FontName","Cambria")
% xlabel('Time \rm (s)','FontWeight','bold')
% yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
% ylim([0 0.3])
% 
% 
% nexttile
% plot(t,Y2(:,2),"LineWidth",2)
% grid on
% grid minor
% %title('Position',"FontName","Cambria")
% %title('\rm Yaw', "FontName","Cambria")
% % yl = ylabel('\bf Yaw Velocity \rm (rad/s)','FontName',"Cambria");
% ylim([0 0.3])
% xlabel('Time \rm (s)','FontWeight','bold')
% 
% 

% title(tcl1,'\bf Simulted Step Responce', "FontName","Cambria")
% subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
% xlabel(tcl1,'Time \rm (s)','FontWeight','bold')
% subtitle(tcl1,'Dynamic Helicopter Model');
%% Figure 2
fig1 = figure(2)
tcl1 = tiledlayout(2,1);

nexttile
plot(t,Y2(:,1),"LineWidth",2)
grid on
grid minor
title('Input',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
ylim([0 0.3])

nexttile
plot(t,Y2(:,2),"LineWidth",2)
grid on
grid minor
title('Position',"FontName","Cambria")
yl = ylabel('Counts','Fontweight','bold','FontName',"Cambria");
ylim([0 0.3])



title(tcl1,'\bf Yaw Motor – Step Responce', "FontName","Cambria")
subtitle(tcl1,'Data vs. Time', "FontName","Cambria")
xlabel(tcl1,'Time (s)','FontWeight','bold')
savefig(fig1,"fig1.fig");
%% Figure 3

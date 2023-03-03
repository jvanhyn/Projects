clear;close all; clc;

%variables
WS.lin = linspace(0,500); %kg/m2
g = 9.81;
% WS = CL_max/2*1/2*rho*Vstall^2;
%% Parameters
%Speeds %m/s
V_stall = 38;                       %%% Revise
V_takeoff = 1.1*V_stall;            %lecture 3 slide 47
V_climb = V_takeoff;                %%%% REVISE
V_maneuver = 38.0556 + 5.14444;
V_cruise = 108.056;
V_approach = 1.3*V_stall;

%Air Densities
rho_sl= 1.2250;                     %kg/m^3
rho_ceil = 1.0556;                  %kg/m^2, 5000ft

%Dynamic Pressures
q_stall = q(rho_sl, V_stall);
q_maneuver = q(rho_sl, V_maneuver);
q_climb = q(rho_ceil, V_climb);
q_cruise = q(rho_ceil, V_cruise);    %Ns/m^2 : 5000 ft @ 220 mph

%Power
P = 2*820e3;                 %Watts - 2 engines

%Gross Weight
W0 = 9000;                   %kg

%load factor
n=3.5;

%lift-drag ratio
LD = 15;
e = 0.8;                      %oswald
AR = 16;                      %aspect ratio
G = 0.0329158;                %200 ft/nmi, climb gradient
% G = .2; %test
% gamma = 15; %climb angle, degrees

c1 = 0.8635; d1 = 0.5632;
Swet = 10^c1*W0^d1;
sweep = 5;                %degrees
TOP = 300;                %takeoff parameter - 3000 ft runway
eta = .7;                 %%double check

%Coefficients
C_D0 = 0.008;

C_lmax = 1.5;
C_Lmax = 0.9*C_lmax*cosd(sweep);
C_Ltakeoff = C_Lmax/1.21;

% Landing
S_landing = 0.507*mps2knots(V_stall)^2;
S_a = 450;
%% Equations

TW.Ceiling = 1/LD;

TW.Maneuver = (C_D0*q_maneuver./WS.lin + n^2/(pi*e*AR*q_maneuver).*WS.lin);

TW.Climb = G + (WS.lin)./(pi*e*AR*q_climb) + C_D0*q_climb./WS.lin;

WS.Stall = 1/2*C_Lmax*rho_ceil*V_stall^2 / g;

TW.Takeoff = WS.lin/(TOP*C_Ltakeoff) * (eta/V_takeoff);

WS.Landing = (S_landing - S_a)*C_Lmax/80;

%% Figure Plot

figure(1)
hold on
set(gca,"FontName","Cambria")

title('T/W vs W/S')
axis([0 500 0 1])
xlabel('W/S')
ylabel('T/W')
yline(TW.Ceiling,'b','LineWidth',2)
xline(WS.Stall,'LineWidth',2)
plot(WS.lin,TW.Maneuver,'r','LineWidth',2)
plot(WS.lin,TW.Climb, 'm','LineWidth',2)
plot(WS.lin, TW.Takeoff, 'g','LineWidth',2)
xline(WS.Landing, 'y','LineWidth',2)

%optimal point
opter = TW.Maneuver;
opty = min(opter);
optx = WS.lin((opter) == opty);
plot(optx,opty,'^k','LineWidth',3)
%
legend('Ceiling', 'Stall', 'Maneuver','Climb', 'Takeoff','Landing','Optimal Point')
grid on
hold off
%% Range / Weight Sensitivity
Wf = linspace(0,W0/2);
LD = linspace(10,18,5);
n_t = 0.35;
e_t = 43.15*10^6;

[LD, Wf] = meshgrid(LD,Wf);


R.WCD = n_t.*e_t./g.*LD.*log(1./(1-Wf./W0))/10000;

figure(2)
plot(Wf/W0,R.WCD,'LineWidth',3)
set(gca,"FontName","Cambria")
title("Range (km) vs Fuel Weight Fraction ")
hleg = legend('10', '12', '13','16', '18','Location','NW')
htitle = get(hleg,'Title');
set(htitle,'String','Lift to Drag Ratio')
grid on
xlabel("Wf/W0")
ylabel("Range (Km)")


%% Mission Profile
close
p0 = 0;
p1 = 0;
p2 = 1520;
p3 = 1520;
p4 = 600;
p5 = 0;
p6 = 0;
pl = 600;
p = [p0 p1 p2 p3 p4 pl p5 p6];

figure(3)
set(gca,"FontName","Cambria")
plot([0 2 3 5 6 7 8 10],p,"LineWidth",3)
hold on
plot([6 7],[600 600],"LineWidth",3)
title("Mission Profile")
ylim([0 2000])
ylabel("Altitude (m)")
xlim([0 10])
xline(0,'--','Takeoff')
xline(2,'--','Climb')
xline(3,'--','Cruise')
xline(5,'--','Decent')
xline(6,'--','Loiter')
xline(7,'--','Approach')
xline(8,'--','Landing')
x = gca; % Get handle
%xticks('')
%set(x,'XTickLabel',{'a','b','c','d'})

xticklabels({'','Short Takeoff','','',"50km - 800km",'','','','','Short Landing'})
xtickangle(0)
%gtext('V_Cr = 98 m/s')

grid on
%% Functions
function dynamicpressure = q(rho,V) 
dynamicpressure = 1/2 * rho .* V .^2;
end

function mps = mph2mps(mph)
mps = mph/2.237;
end

function knots = mps2knots(mps)
knots = mps*1.94384;
end

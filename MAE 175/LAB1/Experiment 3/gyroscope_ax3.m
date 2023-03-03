close
wn=2*pi*2;beta=0.1;K=-950;
G42 = tf(K*wn^2,[1 2*beta*wn wn^2 0]);

dtu = 3;
ts = linspace(0,6,1000);

Ustep = @(x) (x>0);
us =  Ustep(ts) - Ustep(ts-dtu);

y1 = lsim(G42,us,ts);
dy1 = diff(y1);
dt1 = ts(1:end-1);

figure(1)
sgtitle("LAB DATA - AX2", "FontName","Cambria")

subplot(3,1,1)
plot(ts,us,"LineWidth",2)
grid on
grid minor
title('Input')

subplot(3,1,2)
hold on
plot(ts,y1,"LineWidth",2)
grid on
grid minor
title("Position")

subplot(3,1,3)
plot(dt1,dy1,"LineWidth",2)
grid on
grid minor
title("Velocity")
%%
clc
close
Kp = -0.01;
Ki = 0;
Kd = 0;

P = tf(K,1);
I = tf(Ki,[1 0]);
D = tf(Kd*[1 0],1);


P = tf(K,1);
D = tf(Kd*[1 0],1);

D3 = P + D + I;
GD = G42*D3;

H = GD/(1+GD);
% % figure(8)
% % step(H)
% close
% figure(4)
% title('Locus and Resulting Pole Placement')
% hold on
% rlocus(GD)
% plot(pole(H),'x')
% hold off
% figure(5)
step(H)


% figure(5)
% margin(D3)
% xlim([1e-5 1e5])
% figure(6)
% margin(G42)
% xlim([1e-5 1e5])
% figure(7)
% margin(GD)
% xlim([1e-5 1e5])
%%

 

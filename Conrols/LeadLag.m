clear
close all 
clear
syms s z p  K w
G(s) = 1/s^3
D(s) = K*(s+z)/(s+p)
GD(s) = G(s)*D(s)
z = w/sqrt(10)
p = w*sqrt(10)
s = w;
GD = simplify(subs(GD))

K = solve(1==GD,K)

clear

rtime_min = 1;
w = 1.8/rtime_min;
j = 1;
p1 = w*sqrt(10);
z1 = w/sqrt(10);
p2 = w / 100
z2 = w / 10

K1 = sqrt(10)*w^2;
K2 = 1;

G = tf(j,[1 0 0]);
D1 = tf(K1*[1 z1],[1 p1]);
D2 = tf(K2*[1 z2],[1 p2]);
GD = G*D1*D2;
H = GD/(1+GD);

figure(1)
rlocus(GD)

figure(2)
margin(GD)

[Gm,Pm,Wcg,Wcp] = margin(GD);
[mag,phase,wout] = bode(GD);

figure(3)
subplot(3,1,1);
y1 = log10(reshape(mag,[],1));
x = log10(wout);
plot(x,y1)
yline(0)
xline(log10(p1))
xline(log10(z1))
xline(log10(Wcp))

subplot(3,1,2); 
y2 = reshape(phase,[],1);
plot(x,y2);
xline(log10(Wcp))

subplot(3,1,3)
step(H)
S = stepinfo(H,'RiseTimeLimits',[0 1]);
yline(1)
xline(S.RiseTime)
yline(S.Peak)
xline(S.SettlingTime)


eta = Pm/100
RT = S.RiseTime
MP = S.Overshoot

rt = rtime_min
mp = exp(-pi*eta*sqrt(1-eta^2))*100




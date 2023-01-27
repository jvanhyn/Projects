close all 
clear all
p = 20
z = 2
e = 0.1
K = 1000
G = tf(1,[1 0 z])
D1 = K*tf([1 0 z-e],[1 2*p p])
GD1 = G*D1
H = GD1/(1+GD1);

figure
rlocus(GD1)

figure
margin(GD1)

figure 
step(H)
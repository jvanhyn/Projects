%% Problem 1 
M1 = 2;
delta = 17;
beta = 48;
[M2, P2_P1, rho, T2, Mun, Mdn] = ssObliqueShock(M1,delta,beta,1.4)
gamma =1.4;
Mu = 1.37;
Md = 1.96;
M3 = 1.96;
p = ((1+(gamma-1)/2*Mu^2)/(1+(gamma-1)/2*Md^2))^(gamma/(gamma-1));
gamma = 1.4;
T = p^((gamma-1)/gamma);
v = @(M) sqrt((gamma+1)/(gamma-1)) * atan(sqrt((gamma-1)/(gamma+1)*(M^2-1))) - atan(sqrt(M^2-1));
theta =(v(M3) - v(M2))/pi*180+delta
%% Problem 2
clc
M1 = 3;
delta = [10 15 23.34];
betta1 = [27.3827 32 38];
betta2 = [31.795 40.5 54];

[M2, P2_P1] = ssObliqueShock(M1,delta,betta1,1.4);
[M3, P3_P2] = ssObliqueShock(M2,delta,betta2,1.4);

gamma = 1.4;

P01_P1 = (1+(gamma-1)*M1.^2./2).^((gamma)/(gamma-1));
P03_P3 = (1+(gamma-1)*M3.^2./2).^((gamma)/(gamma-1));

P3_P1 = P2_P1.*P3_P2;

P3_P01 = P3_P1./P01_P1;
P03_P01 = P3_P01.*P03_P3;

PP =  1 - P03_P01;


S = gamma*log(P03_P01.^((1-gamma)/gamma))





%%
M1 = 3.2;
delta = 9;
beta1 = 25.21;
beta2 = 29;
[M2, P2_P1, rho, T2_T1] = ssObliqueShock(M1,delta,beta1,gamma)
[M3, P3_P2, rho, T3_T2] = ssObliqueShock(M2,delta,beta2,gamma)
P3_P4 = 1.82;


T2_T1
T3_T1 = T2_T1 * T3_T2
T3_T4 = P3_P4^(gamma/(gamma-1));
T4_T1 = T3_T1/ T3_T4   



%% Problem 2
clc
M1 = 3;
delta = [10 15 20];
betta1 = [27.3827 32 38];
betta2 = [31.795 40.5 54];

[M2, P2_P1, rho, T2, Mun, Mdn] = ssObliqueShock(M1,delta,betta1,1.4);
[M3, P3_P2, rho, T3, Mun, Mdn] = ssObliqueShock(M2,delta,betta2,1.4);

gamma = 1.4;

P01_P1 = (1+(gamma-1)*M1.^2./2).^((gamma)/(gamma-1))
P3_P1 = P2_P1.*P3_P2;
P3_P01 = P3_P1*(1/P01_P1);
P03_P1 = (1+(gamma-1)*M3.^2./2).^((gamma)/(gamma-1))
P03_P01 = P3_P01.*P03_P1;
PP =  1 - P03_P01 

T01_T1 = (1+(gamma-1)*M3.^2./2).^((gamma)/(gamma-1))


s = gamma*log(/(p3./p1).^((gamma-1)/gamma))
M2
M3

%%
[M p rho T Mun Mdn] = ssObliqueShock([3.2 2.7],9,[25 29],1.4)

T(1)/T(2)
T(1)/T(2)*(1/1.82).^((gamma-1)./gamma)



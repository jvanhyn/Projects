%% Midterm 1
% Jonathan Van Hyning
% February 15, 2023

T1 = @(g) g;
%% Problem 1
% What are the properties of Maxtrix A? Is it invertible? 
A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
B = inv(B)

%% Problem 2
% TB2A = T1(g)T2(e)
% TA2B = T2(-e)T1(-g)

%% Problem 3
% Conver ECEF to Geodedic and determine norther or southern hemisphere
X = 4282
Y = -4029 
Z = -2473 
RE = 637813;

R = ECEF2Geo([X, Y, Z], RE);
N = R(1) > 0
S = R(1) < 0

%% Problem 4
% Convert distance in Geo to NED
RE = 6378137;

lat1 = 32.7462;
long1 = -117.1950;
h1 = 1000;

lat2 = 32.7336;
long2 = -117.1897;
h2 = 5;

R1 = [lat1 long1 h1]';
R2 = [lat2 long2 h2]';
R3 = R1-R2;

R3_P =  Geo2ECEP(R3, RE);
NED = [ 0 0 1; 0 1 0; -1 0 0]*R3_P

N = NED(1)
%% Problem 5
b = -74948;
c = 299.792458e6;
tb = -b/c

%% Problem 6
V = 110;
xsi1 = 139;

V_NED = [cosd(xsi1) sind(xsi1) 0]*V;

W = 35;
xsi2 = 40;
W_NED = [cosd(xsi1) sind(xsi1) 0]*W;

V_Ground = V_NED + W_NED;

%% Problem 7
%% Probelm 8
% Convert Quaterion to Eueler Angle
b0 = 0.5036; bx = 0.8421; by = 0.1608; bz = -0.1069;
phi = atan2d(2*(bz*by+bx*b0),bz^2+b0^2-bx^2-by^2);
%% Problem 9
ax = 0; ay = 0; az = -10.8; p = 0; q = 0.6; r = 2.8;
V = 40;
g = -9.81;

theta = atan2d(ax,sqrt((ax-r*V)^2+(az+q*V)^2))
phi = atan2d(-(ay-r*V),-(az+q*V))

psi =(q*sind(phi)+r*cosd(phi))/cosd(theta)
%% Problem 10
u = 419; v = 0; w = 23;
phi = 0;
theta = 10;
psi = 270;

dh= -(-u*sind(theta)+v*sind(phi)*cosd(theta)+w*cosd(phi)*cosd(theta))*60


%% Functions 
function R_GEO = ECEF2Geo(R_ECEP, RE)
X = R_ECEP(1);
Y = R_ECEP(2);
Z = R_ECEP(3);
lat = atan2d(Z,sqrt(X^2+Y^2));
long = atan2d(Y,X);
h = sqrt(X^2+Y^2+Z^2) - RE;
R_GEO = [lat;long; h];
end

function R_ECEP = Geo2ECEP(R_GEO, RE)
lat = R_GEO(1);
long = R_GEO(2);
h = R_GEO(3);
X = (RE+h)*cosd(lat)*cosd(long);
Y = (RE+h)*cosd(lat)*sind(long);
Z = (RE+h)*sind(lat);
R_ECEP = [X;Y; Z];
end
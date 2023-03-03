% utility to simulate and compare the angular velocities of the helicopter model 
% due to a step on yaw motor
%
% Written by R.A. de Callafon, Dept. of MAE, UCSD (2015-2023)
% Report errors in this software to <callafon@ucsd.edu>

% read pitch angle data file due to step on yaw motor
data_file = 'Pitch_Angle_Yaw_Motor_Step1.txt';
disp(['Reading ' data_file '...'])
f=fopen(data_file,'r');
tline = fgetl(f);
k=1;
while tline~=-1,
    tline = fgetl(f);
    if tline~=-1
        [dummy1,dummy2,dummy3,p_data(k)]=strread(tline,'%f\t%f\t%f\t%f\t');
        k=k+1;
    end 
end
fclose(f);

% adjust time vector to be correct on the basis of 200Hz sampling
% as saved time vector in data file only contains integer values!
N=length(p_data);
t_p_data=0:1/200:(N-1)/200;

% compute velocity data on the basis of 200Hz sampling
vp_data = 0*p_data;
vp_data(2:end) = (p_data(2:end)-p_data(1:end-1))/(1/200);
vp_data(1) = vp_data(2);

% read yaw angle data file due to step on yaw motor
data_file = 'Yaw_Angle_Yaw_Motor_Step1.txt';
disp(['Reading ' data_file '...'])
f=fopen(data_file,'r');
tline = fgetl(f);
k=1;
while tline~=-1,
    tline = fgetl(f);
    if tline~=-1
        [dummy1,dummy2,dummy3,y_data(k)]=strread(tline,'%f\t%f\t%f\t%f\t');
        k=k+1;
    end 
end
fclose(f);

% adjust time vector to be correct on the basis of 200Hz sampling
% as saved time vector in data file only contains integer values!
N=length(y_data);
t_y_data=0:1/200:(N-1)/200;

% compute velocity data on the basis of 200Hz sampling
vy_data = 0*y_data;
vy_data(2:end) = (y_data(2:end)-y_data(1:end-1))/(1/200);
vy_data(1) = vy_data(2);

% read values of parameters of helicopte rmodel contained in parameters.m
parameters

% resulting state-space matrices
A=[0 0 1 0
   0 0 0 1
   0 0 -b_p 0
   0 0 0 -b_y];
B=[0 0
   0 0
   a_pp a_py
   a_yp a_yy];
C=[zeros(2,2) eye(2)];
D=zeros(2,2);
G_velocity=ss(A,B,C,D);

% typical rotor dynamics of pitch motor
wn=5;b=0.8;
GM1=tf(wn^2,[1 2*b*wn wn^2]);
% typical rotor dynamics of yaw motor
wn=8;b=0.9;
GM2=tf(wn^2,[1 2*b*wn wn^2]);

% add rotor dynamics
G_velocity_tot = G_velocity*[GM1 0;0 GM2];

% simulate step response of velocity due to step on yaw motor
t_sim = [0:1/200:4-1/200];
N_sim = length(t_sim);
start_of_step = 10; % in number of samples of t_sim
yaw_motor = [zeros(start_of_step,1);ones(N_sim-start_of_step,1)];
pitch_motor = [zeros(N_sim,1)];
v_sim=lsim(180/pi*G_velocity_tot,[pitch_motor yaw_motor],t_sim);  % times 180/pi to get output in deg/s

% compare data (make sure you plot either pitch of yaw angular velocity
figure
subplot(2,1,1)
l=plot(t_sim,v_sim(:,1),'b',t_p_data,vp_data,'r');
set(l,'linewidth',2);
title('pitch angular velocity due to step on yaw motor')
legend('simulated pitch angular velocity','measured pitch angular velocity')
ylabel('angular velocity [deg/s]')
xlabel('time [sec]')
subplot(2,1,2)
l=plot(t_sim,v_sim(:,2),'b',t_y_data,vy_data,'r');
set(l,'linewidth',2);
title('yaw angular velocity due to step on yaw motor')
legend('simulated yaw angular velocity','measured yaw angular velocity')
ylabel('angular velocity [deg/s]')
xlabel('time [sec]')

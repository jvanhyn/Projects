% Parameters for the state space mode of the 2D Helicopter Experiment
% of MAE175a, Dept. of MAE, UCSD
%
% Written by R.A. de Callafon, Dept. of MAE, UCSD (2015-2023)
% Report errors in this software to <callafon@ucsd.edu>

% Estimate of physical parameters of Helicopter
% Kpp = 0.203;     % [V/N]
% Kyy = 0.0720;    % [V/N]
% Kpy = 0.00680;   % [V/N]
% Kyp = 0.0219;    % [V/N]
% Jp  = 0.03836;    % [kg*m^2]
% Jy  = 0.04321;    % [kg*m^2]
% Bp  = 0.800;     % [N*m*s/rad]
% By  = 0.318;     % [N*m*s/rad]
% Mh  = 1.39;      % [kg]
% Lcm = 0.1855;     % [m]
% JTp=Jp+Mh*Lcm^2;  % [kg*m^2]
% JTy=Jy+Mh*Lcm^2;  % [kg*m^2]

% Computation of MIMO transfer function model parameters for Helicopter Model
% b_p  = Bp/JTp;
% b_y  = By/JTy;
% a_pp = Kpp/JTp;
% a_py = Kpy/JTp;
% a_yp = Kyp/JTy;
% a_yy = Kyy/JTy;

% Add your estimate of MIMO transfer function model parameters for Helicopter Model
% (initial guess based on physical parameters above is given here)
b_p  = 9.2818;
b_y  = 3.4930;
a_pp = 2.3553;
a_py = 0.0789;
a_yp = 0.2406;
a_yy = 0.7909;

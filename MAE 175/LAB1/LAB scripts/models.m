% MODELS - Belongs to MAELAB.M
% Template file for specifying SISO Gyroscope transfer function models

% Written by R.A. de Callafon, Dept. of MAE, UCSD (2001-2021)
% Report errors in this software to <callafon@ucsd.edu>

% model from motor 1 to encoder 3 must be specified by G31
K= -32 ; beta=0.5;
G31 = tf(K,[1 beta 0]);

% model from motor 2 to encoder 2 must be specified by G22
wn=2*pi*2;beta=0.1;K=750;
G22 = tf(K*wn^2,[1 2*beta*wn wn^2]);

% model from motor 2 to encoder 4 must be specified by G42
wn=2*pi*2;beta=0.1;K=-950;
G42 = tf(K*wn^2,[1 2*beta*wn wn^2 0]);

% NOTES:
% - any model from motor m to encoder n must be specified by Gnm
% - models must be specified by a tf or ss operation
%%


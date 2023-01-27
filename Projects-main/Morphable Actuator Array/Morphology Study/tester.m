close all
clear

L0=50e-3;%actuator's length
w=15e-3;%actuator's width
n=50;%segment number


 tp0 = 5e-6;     %PVDF's thickness (<=50um)
 tc0 = 0.016*tp0;  %Cu's thickness (<=800nm)
 tx0 = 6e-6;      %PET's thickness (>=6um)
 ts0 = 5e-6;      %silicone's thickness 
 tm0 = tp0;        %Elastic PVDF's thicknes

dtp = 5e-6;
dtc = 50e-9;
dtx = 1e-5;
dL = 10e-3;

tp_i = tp0;
tx_i = tx0;
tc_i = tc0;
L_i = L0;

steps = 10;
length_var = zeros(steps);
M = cell(1,3)

for i = 1:steps
    
    [Lx{i},Ly{i}] = layer_variation(L_i,tp0,tc0,tx0);
    [PETx{i},PETy{i}] = layer_variation(L0,tp0,tc0,tx_i);
    [PVDFx{i},PVDFy{i}] = layer_variation(L0,tp_i,tc_i,tx0);

    tp_i = tp_i+dtp*i;
    tc_i = 0.016*tp_i;
    tx_i = tx_i+dtx;
    L_i = L_i+dL*i;
    
end


length_model = [Lx,Ly];
PET_model = [PETx,PETy];
PVDF_model = [PVDFx,PVDFx];

for i= 1:steps
%     figure(1)
%     hold on
%     plot(Lx{i},-Ly{i}, "LineWidth",2)
%     title("Length Variation")
    
%     figure(2)
%     hold on
%     plot(PETx{i},-PETy{i}, "LineWidth",2)
%     title("PET Variation")
    
%     figure(3)
%     hold on
%     plot(PVDFx{i},-PVDFy{i}, "LineWidth",2)
%     title("PVDF Variation")
%     
end


%%
L0=10e-3;%actuator's length
w=15e-3;%actuator's width
n=50;%segment number

    tp0 = 25e-6;     %PVDF's thickness (<=50um)
    tc0 = 0.016*tp0;  %Cu's thickness (<=800nm)
    tx0 = 6e-6;
    
    [PVDFx,PVDFy] = layer_variation(L0,tp0,tc0,tx0);
    plot(PVDFx,-PVDFy, "LineWidth",2)
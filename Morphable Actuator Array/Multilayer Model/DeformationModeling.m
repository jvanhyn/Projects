clear;
close all
%% Description 
%{
Describes the deformation of a multilayer laminated sheet with embeded
inplane piezo-electrec actuators. Each case represents a
different experiment where different paramenters are varied. The data from
each test is saved and exprted as .xlsx files. 
%}
%% Actuator Properties
L = 50e-3;    %actuator's length
w = 15e-3;    %actuator's width
n = 500;       %segment number
lseg = L/n;   %segement length

%% Material Properties 
Ep = 2.6448e9;  %PVDF's Young's modulus
Em = Ep;        %Elastic PVDF's Young's modulus
Ec = 110e9;     %Cu's Young's modulus
Es = 25e6;      %Silicone's Young's modulus
Ex = 7.25e9;    %PET's Young's modulus

E = [Ep Em Ec Es Ex];

%% Defalt Parameters 
% These values are are used as the defalt paramenters in each experiment. 

% tp = 25e-6;       %PVDF's thickness (<=50um)
% tm = tp;          %Elastic PVDF's thickness
% d31 = 13.58e-12;  %PVDF's d31

% tc = 0.016*tp;    %Cu's thickness   (<=800nm)
% tx = 6e-6;        %PET's thickness  (>=6um)
% ts = 5e-6;        %Silicone's thickness 

% V = 6e7*tp;       %Max voltage applied to the PVDF (<=3000V)

%% Iterative Parameters 

pmin = 1e-7;    %Min PVDF's thickness (<=50um )
xmin = 6e-6;    %Min PET's thickness  (>=6um  )
cmin = 80e-9;   %Min Cu's thickness   (<=800nm)

pstep = 1e-7;   %Delta tp 
xstep = 5e-6;   %Delta tx
cstep = 50e-9;  %Delta tc

steps = 50;     %Total number of iterations

%% Model 

for k = 1:4
    for j = 1:steps
        switch k                        
            case 1  %Varying PVDF thickness, all others at defalt 
   
            tp = pstep*(j-1) + pmin;    %Variable PVDF's thickness
            tm = tp;                  
            d31 = 13.58e-12;          
            
            tc = 0.016*tp;              
            tx = 6e-6;                 
            ts = 5e-6;
            
            V = 6e7*tp;                 
       
            case 2 %Varying Cu thickness, all others at defalt
           
            tp = 25e-6;                 
            tm = tp;                    
            d31 = 13.58e-12;          
            
            tc = cstep*(j-1) + cmin;    %Variable Cu's thickness 
            
            tx = 6e-6;                 
            ts = 5e-6;                     
            V = 6e7*tp;                 
            
            case 3 %Varying PET thickness, all others at defalt 
            
            tp = 25e-6;                 
            tm = tp;                    
            d31 = 13.58e-12;
            
            tc = 0.016*tp;
            tx = xstep*(j-1) + xmin;    %Variable PET's thickness
            ts = 5e-6;
            
            V = 6e7*tp;                  
            
            case 4 %Varying PVDF thickness with modified d31, all others at defalt
                
            tp = pstep*(j-1) + pmin;    %Variable PVDF's thickness
            tm = tp;                    
            d31 = 30e-12;               %PVDF's Modified d31
            
            tc = 0.016*tp;              
            tx = 6e-6;                 
            ts = 5e-6;
            
            V = 6e7*tp;                
            
        end
        
       [XL,YL,X,Y] = LargeDisplacementModel(Ex,tx,Es,ts,Ep,tp,Ec,tc,Em,tm,V,L,n,w,d31);


        switch k
            case 1
            PVDF_thickness(j) = tp;
            PVDFy(:,j) = YL';
            PVDFx(:,j) = XL';
            
            case 2
            Cu_thickness(j) = tc;
            Cuy(:,j) = YL';
            Cux(:,j) = XL';
           
            case 3
            PET_thickness(j) = tx;
            PETy(:,j) = YL';
            PETx(:,j) = XL';
           
            case 4
            PVDF_thickness2(j) = tp;
            PVDFy2(:,j) = YL';
            PVDFx2(:,j) = XL';
           
        end
    end
end

%% Export Data  
spread_start1 = 'A3';

offsety = 2;
offsetx = 0;

a = num2xlcol(steps+offsetx);
b = num2str(n+offsety);

spread_range = strcat(spread_start1,':',a,b);

writematrix(PVDF_thickness,'PVDF.xlsx','Sheet','X','Range',strcat('A1:',a,'1'));
writematrix(PVDF_thickness,'PVDF.xlsx','Sheet','Y','Range',strcat('A1:',a,'1'));
writematrix(PVDFy,'PVDF.xlsx','Sheet','Y','Range',spread_range)
writematrix(PVDFx,'PVDF.xlsx','Sheet','X','Range',spread_range)

writematrix(PVDF_thickness2,'PVDF2.xlsx','Sheet','Y','Range',strcat('A1:',a,'1'));
writematrix(PVDF_thickness2,'PVDF2.xlsx','Sheet','X','Range',strcat('A1:',a,'1'));
writematrix(PVDFy2,'PVDF2.xlsx','Sheet','Y','Range',spread_range)
writematrix(PVDFx2,'PVDF2.xlsx','Sheet','X','Range',spread_range)

writematrix(PET_thickness,'PET.xlsx','Sheet','Y','Range',strcat('A1:',a,'1'));
writematrix(PET_thickness,'PET.xlsx','Sheet','X','Range',strcat('A1:',a,'1'));
writematrix(PETy,'PET.xlsx','Sheet','Y','Range',spread_range)
writematrix(PETx,'PET.xlsx','Sheet','X','Range',spread_range)

writematrix(Cu_thickness,'CU.xlsx','Sheet','Y','Range',strcat('A1:',a,'1'));
writematrix(Cu_thickness,'CU.xlsx','Sheet','X','Range',strcat('A1:',a,'1'));
writematrix(Cuy,'CU.xlsx','Sheet','Y','Range',spread_range)
writematrix(Cux,'CU.xlsx','Sheet','X','Range',spread_range)

%% Plot Tip Displacement
figure(1)
plot(PVDF_thickness,abs(max(PVDFy)))  
figure(2)
plot(Cu_thickness,abs(max(Cuy)))
figure(3)
plot(PET_thickness,abs(max(PETy)))
figure(4)
plot(PVDF_thickness2,abs(max(PVDFy2)))
%% Plot Normalized Deformation
figure(5)
plot(PVDFx,PVDFy)
figure(6)
plot(PETx,PETy)
figure(7)
plot(Cux,Cuy)
figure(8)
plot(PVDFx2,PVDFy2)
%% Save Data
save("PVDF_DATA2","PVDFy2", "PVDFx2", "PVDF_thickness2")
save("PVDF_DATA","PVDFy", "PVDFx", "PVDF_thickness")
save("PET_DATA","PETy", "PETx", "PET_thickness")
save("Cu_DATA","Cuy", "Cux", "Cu_thickness")

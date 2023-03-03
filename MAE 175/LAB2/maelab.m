% MatLab template file for MAE175A 2DOF Helicopter Experiment to help you
% (1) simulate open loop (uncontrolled) behavior of a 2D Helicopter Model
% (2) design controllers (sequential PIDs or via rltool) 
% (3) simulate closed-loop behavior of 2D Helicopter and controllers
% Notes:
% - Edit model parameters in the (default) parameter file PARAMETERS.M
%
% Written by R.A. de Callafon, Dept. of MAE, UCSD (2015-2023)
% Report errors in this software to <callafon@ucsd.edu>

% PLEASE DO NOT MAKE MODIFICATIONS TO THIS FILE
% IF NEEDED, EDIT THE PARAMETERS OF YOUR MODEL IN THE DEFAULT 
% PARAMETERS FILE CALLED PARAMETERS.M
% SEE MAE175A LAB MANUAL OR LECTURE NOTES FOR DETAILS


clc
drawnow;
format compact

disp('- Welcome to MAELAB for 2DOF Quanser Helicopter, v3.2, 2020')

% Ask for the filename with the parameters:
cantfindit=1;
while cantfindit,
    cantfindit=0;
    typo=1;
    while typo==1,
        typo=0;
        if exist('parfile')==1,
            parfilep=input(['- Name of file (between ''s) with model parameters [''' parfile ''']: ']);
            if ~isempty(parfilep),
                parfile=parfilep;
            end
        else
            parfile=input('- Name of file  (between ''s) with model parameters [''parameters'']: ');
            if isempty(parfile),
                parfile='parameters';
            end       
        end
        if abs(parfile(1))<65,
            typo=1;
            clear parfilep parfile
            disp(char(7))
            disp(['==> Sorry, incorrect filename, choosing default name ''parameters''']);            
        end
    end
    if exist('parfile')==1,
        lfname=length(parfile);
        if parfile(lfname-1)~='.', parfile=[parfile '.m']; end
        eval(['testje=exist(''' parfile ''');']);
        if testje~=2,
            disp(['==> Sorry, cannot find ' parfile ]);
            cantfindit=1;
        end
        lfname=length(parfile);
        parfile=parfile(1:lfname-2);
    end
end

% Evoke parameter file
disp(['Evoking ' parfile '.m to update your model parameters']);
eval(parfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The state space model given below is based on the model given in the
% Quanser Laboratory Manual, where we have:
% - two inputs in Voltage respectively for PITCH motor and YAW motor
% - two outputs in rad respectively PITCH and YAW rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of state space model
% JTp=Jp+Mh*Lcm^2;
% JTy=Jy+Mh*Lcm^2;
% b_p  = Bp/JTp;
% b_y  = Bp/JTp;
% a_pp = Kpp/JTp;
% a_py = Kpy/JTp;
% a_yp = Kyp/JTy;
% a_yy = Kyy/JTy;
A=[0 0 1 0
   0 0 0 1
   0 0 -b_p 0
   0 0 0 -b_y];
B=[0 0
   0 0
   a_pp a_py
   a_yp a_yy];
C=[1 0 0 0
   0 1 0 0];
D=[0 0
   0 0];

% Multivariable Model in SS format:
G=ss(A,B,C,D);

% Add simple motor dynamics (motor current settling) and 
% rotor dynamics (rotor speed settling), both ignored in 
% the state space model of the 2D helicopter model
% PITCH motor + rotor
% wn=3.5;b=0.707; % (old values)
% GM1=tf(1,[0.01 1])*tf(wn^2,[1 2*b*wn wn^2]);
wn=5;b=0.8;
GM1=tf(wn^2,[1 2*b*wn wn^2]);
% YAW motor + rotor. Since motor is smaller, current settles faster, but
% rotor dynamics is slower
% wn=5;b=0.9; (old values)
% GM2=tf(1,[0.001 1])*tf(wn^2,[1 2*b*wn wn^2]);
wn=8;b=0.9;
GM2=tf(wn^2,[1 2*b*wn wn^2]);

% include motor dynamics in the model
G=G*[GM1 0;0 GM2];

% default frequency in Bode plots
f=logspace(-2,2,100);

% simulation over 1 minute, samping at 1kHz
Ts=1e-2;
T=[0:Ts:60-Ts];
Np=length(T);

% default pitch and yaw frequency and amplitude for simulation
% frequencies are in Hz
if exist('pitchamplitude')~=1,
    pitchamplitude=5;
    pitchfrequency=0.05;
    yawamplitude=20;
    yawfrequency=0.04;
end

% Some default (starting) controller parameters for PID control
if exist('Kp_P')~=1,
    Kp_P=10;
    Ki_P=0;
    Kd_P=1;

    Kp_Y=1;
    Ki_Y=0;
    Kd_Y=0;
end

% Some default (starting) state weighting for state feedback
if exist('Q11')~=1,
    Q11=200;
    Q22=150;
    Q33=100;
    Q44=200;
    Q55=50;
    Q66=50;
end


menu_iteration=1;
while menu_iteration==1,

    disp(' ');
    disp('-- MAIN MENU --')
    disp('(0) Velocity step response of multivariable model ');   
    disp('(1) Bode plot of multivariable model ');   
    disp('(2) Design PITCH PID control (YAW control fixed)');
    disp('(3) Design YAW PID control (PITCH control fixed)');
    disp('(4) Design multivariable (state feedback) control');
    disp('(5) Design decoupling multivariable (state feedback) control');
    disp('(6) Bode plot of multivariable controller ');   
    disp('(7) Quick (linear) closed-loop simulation');
    disp('(8) Detailed (discrete-time, constrained) closed-loop simulation');
    disp('(9) Exit this menu');


    option=-1;
    while option==-1,
       option=input('Enter your choice: ');
       if isempty(option), option=-1; end
       option=round(max(option));
       if option>9, option=-1; end
       if option<0, option=-1; end
    end
    
if option==0,

    [v_step,t_step] = step(ss(G.A,G.B,180/pi*[zeros(2,2) eye(2) zeros(2,4)],G.D),3);
    figure(1)
    set(1,'Name','Velocity Step Response of Multivariable Plant')
    subplot(2,2,1);
    l=plot(t_step,squeeze(v_step(:,1,1)),'b');
    axis([0 3 0 30]);grid
    set(l,'linewidth',2);
    title('step on pitch motor')
    ylabel('pitch angular velocity [deg/s]')
    subplot(2,2,2);
    l=plot(t_step,squeeze(v_step(:,1,2)),'r');
    axis([0 3 0 30]);grid
    set(l,'linewidth',2);
    title('step on yaw motor')
    subplot(2,2,3);
    l=plot(t_step,squeeze(v_step(:,2,1)),'r');
    axis([0 3 0 30]);grid
    set(l,'linewidth',2);
    ylabel('yaw angular velocity [deg/s]')
    xlabel('time [sec]')
    subplot(2,2,4);
    l=plot(t_step,squeeze(v_step(:,2,2)),'b');
    axis([0 3 0 30]);grid
    set(l,'linewidth',2);
    xlabel('time [sec]')           
    
elseif option==1,
    
    [m11,p11]=bode(G(1,1),2*pi*f);
    [m12,p12]=bode(G(1,2),2*pi*f);
    [m21,p21]=bode(G(2,1),2*pi*f);
    [m22,p22]=bode(G(2,2),2*pi*f);
    
    figure(3)
    set(3,'Name','Phase Bode plot of Multivariable Plant')
    subplot(2,2,1);
    l=semilogx(f,squeeze(p11),'r');
    set(l,'linewidth',1.5);   
    axis([0.01 100 -400 0]);grid
    title('PITCH input');
    ylabel('PITCH output [phase, deg]')
    subplot(2,2,2);
    l=semilogx(f,squeeze(p12),'r');
    set(l,'linewidth',1.5);   
    axis([0.01 100 -400 0]);grid
    title('YAW input');    
    subplot(2,2,3);
    l=semilogx(f,squeeze(p21),'r');
    set(l,'linewidth',1.5);       
    axis([0.01 100 -400 0]);grid
    ylabel('YAW output [phase, deg]')
    xlabel('frequency [Hz]');    
    subplot(2,2,4);
    l=semilogx(f,squeeze(p22),'r');
    set(l,'linewidth',1.5);       
    axis([0.01 100 -400 0]);grid
    xlabel('frequency [Hz]');    
    
    
    figure(2)
    set(2,'Name','Amplitude Bode plot of Multivariable Plant')
    subplot(2,2,1);
    l=loglog(f,squeeze(m11),'r');
    set(l,'linewidth',1.5);   
    axis([0.01 100 1e-6 1e1]);grid
    title('PITCH input');
    ylabel('PITCH output [gain]')
    subplot(2,2,2);
    l=loglog(f,squeeze(m12),'r');
    set(l,'linewidth',1.5);   
    axis([0.01 100 1e-6 1e1]);grid
    title('YAW input');    
    subplot(2,2,3);
    l=loglog(f,squeeze(m21),'r');
    set(l,'linewidth',1.5);       
    axis([0.01 100 1e-6 1e1]);grid
    ylabel('YAW output [gain]')
    xlabel('frequency [Hz]');    
    subplot(2,2,4);
    l=loglog(f,squeeze(m22),'r');
    set(l,'linewidth',1.5);       
    axis([0.01 100 1e-6 1e1]);grid
    xlabel('frequency [Hz]');    
    

elseif option==2,
        
    disp('* First enter fixed YAW PID controller variables');
    
    stability=0;    
    while stability==0,

        if exist('Kp_Y')==1,
            kpp=input(['  Give YAW proportional gain Kp [' num2str(Kp_Y) ']: ']);
            if ~isempty(kpp),
                Kp_Y=kpp;
            end
        else
            Kp_Y=input('  Give YAW proportional gain Kp: ');
            while isempty(Kp_Y),
                clear Kp_Y
                Kp_Y=input('  Give YAW proportional gain Kp: ');
            end
        end
        if exist('Kd_Y')==1,
            kdp=input(['  Give YAW derivative gain Kd [' num2str(Kd_Y) ']: ']);
            if ~isempty(kdp),
                Kd_Y=kdp;
            end
        else
            Kd_Y=input('  Give YAW derivative gain Kd: ');
            while isempty(Kd_Y),
                clear Kd_Y
                Kd_Y=input('  Give YAW derivative gain Kd: ');
            end
        end
        if exist('Ki_Y')==1,
            kip=input(['  Give YAW integral gain Ki [' num2str(Ki_Y) ']: ']);
            if ~isempty(kip),
                Ki_Y=kip;
            end
        else
            Ki_Y=input('  Give YAW integral gain Ki: ');
            while isempty(Ki_Y),
                clear Ki_Y
                Ki_Y=input('  Give YAW integral gain Ki: ');
            end
        end
        Kp_Y=Kp_Y(1);
        Kd_Y=Kd_Y(1);
        Ki_Y=Ki_Y(1);

        % YAW PID control
        disp('  Resulting YAW PID control:');
        PID_YAW=Kp_Y + Ki_Y*tf(1,[1 0]) + Kd_Y*tf([50 0],[1 50])

        % The equivalent open-loop system
        Teq=feedback(G,PID_YAW,2,2);
        Teq=Teq(1,1);

        % check stability
        stability=1;
        clpoles=eig(Teq);
        if max(real(clpoles))>0,
            stability=0;
            disp('==> EQUIVALENT OPEN-LOOP UNSTABLE WITH THIS YAW CONTROLLER');
            disp('* Reenter fixed YAW PID controller variables');
        end
        
    end

    
    % Compute its open-loop Bode and Nuquist response
    [mag,pha]=bode(Teq,2*pi*f);
    mag=squeeze(mag);
    pha=squeeze(pha);    
    [re,im]=nyquist(Teq,2*pi*f);
    re=squeeze(re);
    im=squeeze(im);


    
    PITCH_UPDATE=1;    
    while (PITCH_UPDATE==1),
        
        disp('* Now update PITCH PID controller variables');

        if exist('Kp_P')==1,
            kpp=input(['  Give PITCH proportional gain Kp [' num2str(Kp_P) ']: ']);
            if ~isempty(kpp),
                Kp_P=kpp;
            end
        else
            Kp_P=input('  Give PITCH proportional gain Kp: ');
            while isempty(Kp_P),
                clear Kp_P
                Kp_P=input('  Give PITCH proportional gain Kp: ');
            end
        end
        if exist('Kd_P')==1,
            kdp=input(['  Give PITCH derivative gain Kd [' num2str(Kd_P) ']: ']);
            if ~isempty(kdp),
                Kd_P=kdp;
            end
        else
            Kd_P=input('  Give PITCH derivative gain Kd: ');
            while isempty(Kd_P),
                clear Kd_P
                Kd_P=input('  Give PITCH derivative gain Kd: ');
            end
        end
        if exist('Ki_P')==1,
            kip=input(['  Give PITCH integral gain Ki [' num2str(Ki_P) ']: ']);
            if ~isempty(kip),
                Ki_P=kip;
            end
        else
            Ki_P=input('  Give PITCH integral gain Ki: ');
            while isempty(Ki_P),
                clear Ki_P
                Ki_P=input('  Give PITCH integral gain Ki: ');
            end
        end
        Kp_P=Kp_P(1);
        Kd_P=Kd_P(1);
        Ki_P=Ki_P(1);
        
        % PITCH PID control
        disp('  Resulting PITCH PID control:');
        PID_PITCH=Kp_P + Ki_P*tf(1,[1 0]) + Kd_P*tf([50 0],[1 50])
        
        % compute bode plot of controller
        [magc,phac]=bode(PID_PITCH,2*pi*f);
        magc=squeeze(magc);
        phac=squeeze(phac);
        % compute bode and nyquist plot of loop gain
        [magl,phal]=bode(Teq*PID_PITCH,2*pi*f);
        magl=squeeze(magl);
        phal=squeeze(phal);
        [rel,iml]=nyquist(Teq*PID_PITCH,2*pi*f);
        rel=squeeze(rel);
        iml=squeeze(iml);
        % compute bode plot of closed-loop transfer function
        [magcl,phacl]=bode(feedback(Teq*PID_PITCH,1),2*pi*f);
        magcl=squeeze(magcl);
        

        % check stability
        stability=1;
        clpoles=eig(feedback(Teq,PID_PITCH));
        if max(real(clpoles))>0,
            stability=0;
            disp('==> CLOSED LOOP UNSTABLE');
        end
        
        figure(4)
        plttitle='Bode and Nyquist plots of sequential Loop Transfer and closed-loop PITCH control';
        set(4,'Name',plttitle);
        clf % use to be clg

        subplot(2,2,1);
        l=loglog(f,magl,f,mag,'b--',f,magc,'r--');
        set(l,'linewidth',1.5);   
        hold on
        l=loglog([f(1) f(length(f))],[1 1],'k:');
        set(l,'linewidth',1.5);   
        hold off
        title('Magnitude Loop Transfer')
        xlabel('frequency [Hz]');
        ylabel('magnitude')
        legend('T*C','T','C')

        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        subplot(2,2,2);
        l=plot(rel,iml);
        set(l,'linewidth',1.5);   
        axis([-4 4 -4 4]);
        hold on
        l=plot(-1,0,'*');
        set(l,'linewidth',1.5);   
        hold off
        title('Nyquist contour Loop Transfer')
        xlabel('real')
        ylabel('imag')

        subplot(2,2,3);
        l=semilogx(f,phal,f,pha,'b--',f,phac,'r--');
        set(l,'linewidth',1.5);   
        hold on
        l=semilogx([f(1) f(length(f))],[-180 -180],'k:');
        set(l,'linewidth',1.5);   
        hold off
        xlabel('frequency [Hz]');
        ylabel('phase [deg]');
        title('Phase Loop Transfer')
        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        subplot(2,2,4);
        l=loglog(f,magcl);
        set(l,'linewidth',1.5);   
        hold on
        l=loglog([f(1) f(length(f))],[1 1],'k:');
        set(l,'linewidth',1.5);   
        hold off
        title('Magnitude closed-loop PITCH control')
        xlabel('frequency [Hz]');
        ylabel('magnitude')
        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        if stability==1,
          text(0.1,0.1,'closed loop stable','color','green','Units','normalized');
        else
          text(0.1,0.1,'closed loop unstable!','color','red','Units','normalized');
        end
        
        yn=[];
        y=1;n=0;Y=1;N=0;
        while isempty(yn),
            yn=input('* Continue updating PITCH PID control? [y/n]: ');
            if ~((yn==1)|(yn==0)),
                yn=[];
            end
        end
        if yn==0,
            PITCH_UPDATE=0;
        end
    end
    
    if stability==1,
        % Resulting Multivariable controller
        K=[PID_PITCH 0;0 PID_YAW];
    else
        clear K
    end
            
elseif option==3,
        
    disp('* First enter fixed PITCH PID controller variables');
    
    stability=0;
    while stability==0,

        if exist('Kp_P')==1,
            kpp=input(['  Give PITCH proportional gain Kp [' num2str(Kp_P) ']: ']);
            if ~isempty(kpp),
                Kp_P=kpp;
            end
        else
            Kp_P=input('  Give PITCH proportional gain Kp: ');
            while isempty(Kp_P),
                clear Kp_P
                Kp_P=input('  Give PITCH proportional gain Kp: ');
            end
        end
        if exist('Kd_P')==1,
            kdp=input(['  Give PITCH derivative gain Kd [' num2str(Kd_P) ']: ']);
            if ~isempty(kdp),
                Kd_P=kdp;
            end
        else
            Kd_P=input('  Give PITCH derivative gain Kd: ');
            while isempty(Kd_P),
                clear Kd_P
                Kd_P=input('  Give PITCH derivative gain Kd: ');
            end
        end
        if exist('Ki_P')==1,
            kip=input(['  Give PITCH integral gain Ki [' num2str(Ki_P) ']: ']);
            if ~isempty(kip),
                Ki_P=kip;
            end
        else
            Ki_P=input('  Give PITCH integral gain Ki: ');
            while isempty(Ki_P),
                clear Ki_P
                Ki_P=input('  Give PITCH integral gain Ki: ');
            end
        end
        Kp_P=Kp_P(1);
        Kd_P=Kd_P(1);
        Ki_P=Ki_P(1);

        % PITCH PID control
        disp('  Resulting PITCH PID control:');
        PID_PITCH=Kp_P + Ki_P*tf(1,[1 0]) + Kd_P*tf([50 0],[1 50])

        % The equivalent open-loop system
        Teq=feedback(G,PID_PITCH,1,1);
        Teq=Teq(2,2);
    
        % check stability
        stability=1;
        clpoles=eig(Teq);
        if max(real(clpoles))>0,
            stability=0;
            disp('==> EQUIVALENT OPEN-LOOP UNSTABLE WITH THIS PITCH CONTROLLER');
            disp('* Reenter fixed PITCH PID controller variables');
        end
        
    end

    % Compute its open-loop Bode and Nuquist response
    [mag,pha]=bode(Teq,2*pi*f);
    mag=squeeze(mag);
    pha=squeeze(pha);    
    [re,im]=nyquist(Teq,2*pi*f);
    re=squeeze(re);
    im=squeeze(im);
    
    YAW_UPDATE=1;    
    while (YAW_UPDATE==1),
        
        disp('* Now update YAW PID controller variables');

        if exist('Kp_Y')==1,
            kpp=input(['  Give YAW proportional gain Kp [' num2str(Kp_Y) ']: ']);
            if ~isempty(kpp),
                Kp_Y=kpp;
            end
        else
            Kp_Y=input('  Give YAW proportional gain Kp: ');
            while isempty(Kp_Y),
                clear Kp_Y
                Kp_Y=input('  Give YAW proportional gain Kp: ');
            end
        end
        if exist('Kd_Y')==1,
            kdp=input(['  Give YAW derivative gain Kd [' num2str(Kd_Y) ']: ']);
            if ~isempty(kdp),
                Kd_Y=kdp;
            end
        else
            Kd_Y=input('  Give YAW derivative gain Kd: ');
            while isempty(Kd_Y),
                clear Kd_Y
                Kd_Y=input('  Give YAW derivative gain Kd: ');
            end
        end
        if exist('Ki_Y')==1,
            kip=input(['  Give YAW integral gain Ki [' num2str(Ki_Y) ']: ']);
            if ~isempty(kip),
                Ki_Y=kip;
            end
        else
            Ki_Y=input('  Give YAW integral gain Ki: ');
            while isempty(Ki_Y),
                clear Ki_Y
                Ki_Y=input('  Give YAW integral gain Ki: ');
            end
        end
        Kp_Y=Kp_Y(1);
        Kd_Y=Kd_Y(1);
        Ki_Y=Ki_Y(1);

        % YAW PID control
        disp('  Resulting YAW PID control:');
        PID_YAW=Kp_Y + Ki_Y*tf(1,[1 0]) + Kd_Y*tf([50 0],[1 50])
        
        
        % compute bode plot of controller
        [magc,phac]=bode(PID_YAW,2*pi*f);
        magc=squeeze(magc);
        phac=squeeze(phac);
        % compute bode and nyquist plot of loop gain
        [magl,phal]=bode(Teq*PID_YAW,2*pi*f);
        magl=squeeze(magl);
        phal=squeeze(phal);
        [rel,iml]=nyquist(Teq*PID_YAW,2*pi*f);
        rel=squeeze(rel);
        iml=squeeze(iml);
        % compute bode plot of closed-loop transfer function
        [magcl,phacl]=bode(feedback(Teq*PID_YAW,1),2*pi*f);
        magcl=squeeze(magcl);

        % check stability
        stability=1;
        clpoles=eig(feedback(Teq,PID_YAW));
        if max(real(clpoles))>0,
            stability=0;
            disp('==> CLOSED LOOP UNSTABLE');
        end
        
        figure(5)
        plttitle='Bode and Nyquist plots of sequential Loop Transfer and closed-loop YAW control';
        set(5,'Name',plttitle);
        clf % use to be clg

        subplot(2,2,1);
        l=loglog(f,magl,f,mag,'b--',f,magc,'r--');
        set(l,'linewidth',1.5);   
        hold on
        l=loglog([f(1) f(length(f))],[1 1],'k:');
        set(l,'linewidth',1.5);   
        hold off
        title('Magnitude Loop Transfer')
        xlabel('frequency [Hz]');
        ylabel('magnitude')
        legend('T*C','T','C')

        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        subplot(2,2,2);
        l=plot(rel,iml,re,im,'r--');
        set(l,'linewidth',1.5);   
        axis([-4 4 -4 4]);
        hold on
        l=plot(-1,0,'*');
        set(l,'linewidth',1.5);   
        hold off
        title('Nyquist contour Loop Transfer')
        xlabel('real')
        ylabel('imag')

        subplot(2,2,3);
        l=semilogx(f,phal,f,pha,'b--',f,phac,'r--');
        set(l,'linewidth',1.5);   
        hold on
        l=semilogx([f(1) f(length(f))],[-180 -180],'k:');
        set(l,'linewidth',1.5);   
        hold off
        xlabel('frequency [Hz]');
        ylabel('phase [deg]');
        title('Phase Loop Transfer')
        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        subplot(2,2,4);
        l=loglog(f,magcl);
        set(l,'linewidth',1.5);   
        hold on
        l=loglog([f(1) f(length(f))],[1 1],'k:');
        set(l,'linewidth',1.5);   
        hold off
        title('Magnitude closed-loop YAW control')
        xlabel('frequency [Hz]');
        ylabel('magnitude')
        axiss=axis;
        axis([f(1) f(length(f)) axiss(3) axiss(4)]);

        if stability==1,
          text(0.1,0.1,'closed loop stable','color','green','Units','normalized');
        else
          text(0.1,0.1,'closed loop unstable!','color','red','Units','normalized');
        end
        
        yn=[];
        y=1;n=0;Y=1;N=0;
        while isempty(yn),
            yn=input('* Continue updating YAW PID control? [y/n]: ');
            if ~((yn==1)|(yn==0)),
                yn=[];
            end
        end
        if yn==0,
            YAW_UPDATE=0;
        end
    end
    
    if stability==1,
        % Resulting Multivariable controller
        K=[PID_PITCH 0;0 PID_YAW];
    else
        clear K
    end
    
elseif option==4,                
    
    % no decoupling...
    W=eye(2);    

    if exist('Q11')==1,
        Q11p=input(['  Give pitch angle weighting Q(1,1) [' num2str(Q11) ']: ']);
        if ~isempty(Q11p),
            Q11=Q11p;
        end
    else
        Q11=input('  Give pitch angle weighting Q(1,1): ');
        if Q11<0,
            Q11=[];
        end
        while isempty(Q11),
            clear Q11
            Q11=input('  Give pitch angle weighting Q(1,1): ');
        end
    end

    if exist('Q22')==1,
        Q22p=input(['  Give yaw angle weighting Q(2,2) [' num2str(Q22) ']: ']);
        if ~isempty(Q22p),
            Q22=Q22p;
        end
    else
        Q22=input('  Give yaw angle weighting Q(2,2): ');
        if Q22<0,
            Q22=[];
        end
        while isempty(Q22),
            clear Q22
            Q22=input('  Give yaw angle weighting Q(2,2): ');
        end
    end

    if exist('Q33')==1,
        Q33p=input(['  Give pitch angular velocity weighting Q(3,3) [' num2str(Q33) ']: ']);
        if ~isempty(Q33p),
            Q33=Q33p;
        end
    else
        Q33=input('  Give pitch angular velocity weighting Q(3,3): ');
        if Q33<0,
            Q33=[];
        end
        while isempty(Q33),
            clear Q33
            Q33=input('  Give pitch angular velocity weighting Q(3,3): ');
        end
    end

    if exist('Q44')==1,
        Q44p=input(['  Give yaw angular velocity weighting Q(4,4) [' num2str(Q44) ']: ']);
        if ~isempty(Q44p),
            Q44=Q44p;
        end
    else
        Q44=input('  Give yaw angular velocity weighting Q(4,4): ');
        if Q44<0,
            Q44=[];
        end
        while isempty(Q44),
            clear Q44
            Q44=input('  Give yaw angular velocity weighting Q(4,4): ');
        end
    end

    if exist('Q55')==1,
        Q55p=input(['  Give pitch angle integration weighting Q(5,5) [' num2str(Q55) ']: ']);
        if ~isempty(Q55p),
            Q55=Q55p;
        end
    else
        Q55=input('  Give pitch angle integration weighting Q(5,5): ');
        if Q55<0,
            Q55=[];
        end
        while isempty(Q55),
            clear Q55
            Q55=input('  Give pitch angle integration weighting Q(5,5): ');
        end
    end

    if exist('Q66')==1,
        Q66p=input(['  Give yaw angle integration weighting Q(6,6) [' num2str(Q66) ']: ']);
        if ~isempty(Q66p),
            Q66=Q66p;
        end
    else
        Q66=input('  Give yaw angle integration weighting Q(6,6): ');
        if Q66<0,
            Q66=[];
        end
        while isempty(Q66),
            clear Q66
            Q66=input('  Give yaw angle integration weighting Q(6,6): ');
        end
    end

    % extend state space model with two additional integrators
    Ae=[A zeros(4,2);eye(2) zeros(2,4)];
    Be=[B; zeros(2,2)];
    Ce=[C zeros(2,2)];
    De=D;
    Pe=ss(Ae,Be,Ce,De);
            
    % definition of Q and R
    Q=diag([Q11 Q22 Q33 Q44 Q55 Q66]);
    R=eye(2);
    % P.S. a good choice for Q might be
    % Q = diag([100   100    10    10    10    10])
    % and include the static input decoupling
    
    % compute K via LQR (including the static input weighting)
    % Note that with static input weighting, the only thing that 
    % changes is B matrix (D matrix is zero), 
    % so we can still do state space design! So compute K via LQR
    Pd=Pe*W;
    Klqr=lqr(Pd,Q,R);
    
    % bring (unit) weighting back into Klqr
    Klqr=W*Klqr
        
    % compute resulting multivariable controller
    % (note: Klqr was already computed for negative feedback)
    K=Klqr(1:2,1:2) + Klqr(1:2,3:4)*[eye(2)*tf([50 0],[1 50])] + Klqr(1:2,5:6)*[eye(2)*tf(1,[1 0])];

elseif option==5,                
    
    % If you go back to the state space model without motor dynamics
    % Gss=ss(A,B,C,D);
    % and knowing the shape of the multivariable model we can choose
    % [num,den]=tfdata(Gss);
    % W=[num{2,2}(3) -num{1,2}(3);-num{2,1}(3) num{1,1}(3)];
    % so that G*W will be completely decoupled! Unfortunately, this
    % analysis does not include the mtor dynamics... 
    % Alternatively, compute the gain of the total system (including
    % motor dynamics) at an appropriate frequency used for control 
    % and compute scaling on the basis of that. Here we choose at
    % 0.05Hz:
    M=bode(G,0.05*2*pi);
    W=[M(2,2) -M(1,2);-M(2,1) M(1,1)];

    if exist('Q11')==1,
        Q11p=input(['  Give pitch angle weighting Q(1,1) [' num2str(Q11) ']: ']);
        if ~isempty(Q11p),
            Q11=Q11p;
        end
    else
        Q11=input('  Give pitch angle weighting Q(1,1): ');
        if Q11<0,
            Q11=[];
        end
        while isempty(Q11),
            clear Q11
            Q11=input('  Give pitch angle weighting Q(1,1): ');
        end
    end

    if exist('Q22')==1,
        Q22p=input(['  Give yaw angle weighting Q(2,2) [' num2str(Q22) ']: ']);
        if ~isempty(Q22p),
            Q22=Q22p;
        end
    else
        Q22=input('  Give yaw angle weighting Q(2,2): ');
        if Q22<0,
            Q22=[];
        end
        while isempty(Q22),
            clear Q22
            Q22=input('  Give yaw angle weighting Q(2,2): ');
        end
    end

    if exist('Q33')==1,
        Q33p=input(['  Give pitch angular velocity weighting Q(3,3) [' num2str(Q33) ']: ']);
        if ~isempty(Q33p),
            Q33=Q33p;
        end
    else
        Q33=input('  Give pitch angular velocity weighting Q(3,3): ');
        if Q33<0,
            Q33=[];
        end
        while isempty(Q33),
            clear Q33
            Q33=input('  Give pitch angular velocity weighting Q(3,3): ');
        end
    end

    if exist('Q44')==1,
        Q44p=input(['  Give yaw angular velocity weighting Q(4,4) [' num2str(Q44) ']: ']);
        if ~isempty(Q44p),
            Q44=Q44p;
        end
    else
        Q44=input('  Give yaw angular velocity weighting Q(4,4): ');
        if Q44<0,
            Q44=[];
        end
        while isempty(Q44),
            clear Q44
            Q44=input('  Give yaw angular velocity weighting Q(4,4): ');
        end
    end

    if exist('Q55')==1,
        Q55p=input(['  Give pitch angle integration weighting Q(5,5) [' num2str(Q55) ']: ']);
        if ~isempty(Q55p),
            Q55=Q55p;
        end
    else
        Q55=input('  Give pitch angle integration weighting Q(5,5): ');
        if Q55<0,
            Q55=[];
        end
        while isempty(Q55),
            clear Q55
            Q55=input('  Give pitch angle integration weighting Q(5,5): ');
        end
    end

    if exist('Q66')==1,
        Q66p=input(['  Give yaw angle integration weighting Q(6,6) [' num2str(Q66) ']: ']);
        if ~isempty(Q66p),
            Q66=Q66p;
        end
    else
        Q66=input('  Give yaw angle integration weighting Q(6,6): ');
        if Q66<0,
            Q66=[];
        end
        while isempty(Q66),
            clear Q66
            Q66=input('  Give yaw angle integration weighting Q(6,6): ');
        end
    end

    % extend state space model with two additional integrators
    Ae=[A zeros(4,2);eye(2) zeros(2,4)];
    Be=[B; zeros(2,2)];
    Ce=[C zeros(2,2)];
    De=D;
    Pe=ss(Ae,Be,Ce,De);
        
    % definition of Q and R
    Q=diag([Q11 Q22 Q33 Q44 Q55 Q66]);
    R=eye(2);
    % P.S. a good choice for Q might be
    % Q = diag([100   100    10    10    10    10])
    % and include the static input decoupling
    
    % compute K via LQR (including the static input weighting)
    % Note that with static input weighting, the only thing that 
    % changes is B matrix (D matrix is zero), 
    % so we can still do state space design! So compute K via LQR
    Pd=Pe*W;
    % ensure it is really decoupled
    index1=find(Pd.B(:,1)<1e-4);
    index2=find(Pd.B(:,2)<1e-4);
    Pd.B(index1,1)=0*Pd.B(index1,1);
    Pd.B(index2,2)=0*Pd.B(index2,2);
    Klqr=lqr(Pd,Q,R);
    
    % bring weighting back into Klqr
    Klqr=W*Klqr
            
    % compute resulting multivariable controller
    % (note: Klqr was already computed for negative feedback)
    K=Klqr(1:2,1:2) + Klqr(1:2,3:4)*[eye(2)*tf([50 0],[1 50])] + Klqr(1:2,5:6)*[eye(2)*tf(1,[1 0])];
    
elseif option==6,
    
    if exist('K')~=1,
        disp('* Cannot compute Bode plot of controller!');
        disp('* First design stabilizing (multivariable) controller.');
        disp('  Press any key to return to the main menu...');
        pause
    else    
    
        [m11,p11]=bode(K(1,1),2*pi*f);
        [m12,p12]=bode(K(1,2),2*pi*f);
        [m21,p21]=bode(K(2,1),2*pi*f);
        [m22,p22]=bode(K(2,2),2*pi*f);

        figure(7)
        set(7,'Name','Phase Bode plot of Multivariable Controller')
        subplot(2,2,1);
        l=semilogx(f,squeeze(p11),'r');
        set(l,'linewidth',1.5);   
        axis([0.01 100 -300 300]);grid
        title('PITCH error');
        ylabel('PITCH control [phase, deg]')
        subplot(2,2,2);
        l=semilogx(f,squeeze(p12),'r');
        set(l,'linewidth',1.5);   
        axis([0.01 100 -300 300]);grid
        title('YAW error');    
        subplot(2,2,3);
        l=semilogx(f,squeeze(p21),'r');
        set(l,'linewidth',1.5);       
        axis([0.01 100 -300 300]);grid
        ylabel('YAW control [phase, deg]')
        xlabel('frequency [Hz]');    
        subplot(2,2,4);
        l=semilogx(f,squeeze(p22),'r');
        set(l,'linewidth',1.5);       
        axis([0.01 100 -300 300]);grid
        xlabel('frequency [Hz]');    


        figure(6)
        set(6,'Name','Amplitude Bode plot of Multivariable Controller')
        subplot(2,2,1);
        l=loglog(f,squeeze(m11),'r');
        set(l,'linewidth',1.5);   
        axis([0.01 100 1e-1 1e3]);grid
        title('PITCH error');
        ylabel('PITCH control [gain]')
        subplot(2,2,2);
        l=loglog(f,squeeze(m12),'r');
        set(l,'linewidth',1.5);   
        axis([0.01 100 1e-1 1e3]);grid
        title('YAW error');    
        subplot(2,2,3);
        l=loglog(f,squeeze(m21),'r');
        set(l,'linewidth',1.5);       
        axis([0.01 100 1e-1 1e3]);grid
        ylabel('YAW control [gain]')
        xlabel('frequency [Hz]');    
        subplot(2,2,4);
        l=loglog(f,squeeze(m22),'r');
        set(l,'linewidth',1.5);       
        axis([0.01 100 1e-1 1e3]);grid
        xlabel('frequency [Hz]');    
        
    end
    
elseif option==7,                
    
    if exist('K')~=1,
        disp('* Cannot do a closed-loop simulation!');
        disp('* First design stabilizing PITCH and YAW (multivariable) controllers.');
        disp('  Press any key to return to the main menu...');
        pause
    else

        % specification of feedforward and reference signals
        % the feedforward in this case is simply a force disturbance
        % on both the pitch (due to off balance) and the yaw 
        % (due to continuous rotation of pitch rotor)
        % neccessitating the need for integral control
        v=[-2*ones(Np,1) 1*ones(Np,1)];
        
        if exist('pitchamplitude')==1,
        pitchamplitudep=input(['Enter pitch amplitude [' num2str(pitchamplitude) ' deg]: ']);
            if ~isempty(pitchamplitudep),
                pitchamplitude=pitchamplitudep;
            end
        else
            pitchamplitude=input('Enter pitch amplitude [deg]: ');
            while isempty(pitchamplitude),
                clear pitchamplitude
                pitchamplitude=input('Enter pitch amplitude [deg]: ');
            end
        end

        if exist('pitchfrequency')==1,
        pitchfrequencyp=input(['Enter pitch frequency [' num2str(pitchfrequency) ' Hz]: ']);
            if ~isempty(pitchfrequencyp),
                pitchfrequency=pitchfrequencyp;
            end
        else
            pitchfrequency=input('Enter pitch frequency [Hz]: ');
            while isempty(pitchfrequency),
                clear pitchfrequency
                pitchfrequency=input('Enter pitch frequency [Hz]: ');
            end
        end

        if exist('yawamplitude')==1,
        yawamplitudep=input(['Enter yaw amplitude [' num2str(yawamplitude) ' deg]: ']);
            if ~isempty(yawamplitudep),
                yawamplitude=yawamplitudep;
            end
        else
            yawamplitude=input('Enter yaw amplitude [deg]: ');
            while isempty(yawamplitude),
                clear yawamplitude
                yawamplitude=input('Enter yaw amplitude [deg]: ');
            end
        end

        if exist('yawfrequency')==1,
        yawfrequencyp=input(['Enter yaw frequency [' num2str(yawfrequency) ' Hz]: ']);
            if ~isempty(yawfrequencyp),
                yawfrequency=yawfrequencyp;
            end
        else
            yawfrequency=input('Enter yaw frequency [Hz]: ');
            while isempty(yawfrequency),
                clear yawfrequency
                yawfrequency=input('Enter yaw frequency [Hz]: ');
            end
        end

        % step wise changes in PITCH and YAW in rad!
        r=2*pi/360*[pitchamplitude*sign(sin(2*pi*pitchfrequency*T)') yawamplitude*sign(sin(2*pi*yawfrequency*T)')];


        %% PURE LINEAR SIMULATION

        % Compute feedback system with as inputs:
        % v = feedforward on plant input 
        % r = two references (PITCH and YAW reference)
        % and outputs:
        % z = control outputs = plant inputs
        % y = plant outputs (PITCH and YAW position)
        % Computation is done with feedback function

        %SYS1=[I K;G G*K] = [I 0;0 G]*[0 I;I I]*[I 0;0 K]
        SYS1=[eye(2) zeros(2,2);zeros(2,2) G]*[zeros(2,2) eye(2);eye(2) eye(2)]*[eye(2) zeros(2,2);zeros(2,2) K];
        Tcl=feedback(SYS1,eye(2),[3 4],[3 4]);

        % pure linear simulation
        uy=lsim(Tcl,[v r],T);
        u=uy(:,1:2);
        y=uy(:,3:4);

        figure(8)
        set(8,'Name','Quick (Linear) Closed-Loop Simulation Results');
        % plot output [rad] and reference [rad] in degrees
        subplot(2,1,1);
        l=plot(T,y(:,1)*360/2/pi,'r',T,r(:,1)*360/2/pi,'b',T,y(:,2)*360/2/pi,'r--',T,r(:,2)*360/2/pi,'b--');
        set(l,'linewidth',1.5);
        legend('PITCH output','PITCH reference','YAW output','YAW reference');
        title('Quick (Linear) Closed-Loop Simulation Results')
        ylabel('pitch and yaw [deg]')
        xlabel('time [sec]')
        subplot(2,1,2);
        l=plot(T,u(:,1),'m',T,u(:,2),'m--');
        set(l,'linewidth',1.5);
        legend('PITCH control','YAW control')
        ylabel('pitch and yaw control [V]')
        xlabel('time [sec]')
    end
    
elseif option==8,                

    if exist('K')~=1,
        disp('* Cannot do a closed-loop simulation!');
        disp('* First design stabilizing PITCH and YAW (multivariable) controllers.');
        disp('  Press any key to return to the main menu...');
        pause
    else

        % specification of feedforward and reference signals
        % the feedforward in this case is simply a force disturbance
        % on both the pitch (due to off balance) and the yaw 
        % (due to continuous rotation of pitch rotor)
        % neccessitating the need for integral control
        v=[-2*ones(Np,1) 1*ones(Np,1)];
        
        if exist('pitchamplitude')==1,
        pitchamplitudep=input(['Enter pitch amplitude [' num2str(pitchamplitude) ' deg]: ']);
            if ~isempty(pitchamplitudep),
                pitchamplitude=pitchamplitudep;
            end
        else
            pitchamplitude=input('Enter pitch amplitude [deg]: ');
            while isempty(pitchamplitude),
                clear pitchamplitude
                pitchamplitude=input('Enter pitch amplitude [deg]: ');
            end
        end

        if exist('pitchfrequency')==1,
        pitchfrequencyp=input(['Enter pitch frequency [' num2str(pitchfrequency) ' Hz]: ']);
            if ~isempty(pitchfrequencyp),
                pitchfrequency=pitchfrequencyp;
            end
        else
            pitchfrequency=input('Enter pitch frequency [Hz]: ');
            while isempty(pitchfrequency),
                clear pitchfrequency
                pitchfrequency=input('Enter pitch frequency [Hz]: ');
            end
        end

        if exist('yawamplitude')==1,
        yawamplitudep=input(['Enter yaw amplitude [' num2str(yawamplitude) ' deg]: ']);
            if ~isempty(yawamplitudep),
                yawamplitude=yawamplitudep;
            end
        else
            yawamplitude=input('Enter yaw amplitude [deg]: ');
            while isempty(yawamplitude),
                clear yawamplitude
                yawamplitude=input('Enter yaw amplitude [deg]: ');
            end
        end

        if exist('yawfrequency')==1,
        yawfrequencyp=input(['Enter yaw frequency [' num2str(yawfrequency) ' Hz]: ']);
            if ~isempty(yawfrequencyp),
                yawfrequency=yawfrequencyp;
            end
        else
            yawfrequency=input('Enter yaw frequency [Hz]: ');
            while isempty(yawfrequency),
                clear yawfrequency
                yawfrequency=input('Enter yaw frequency [Hz]: ');
            end
        end

        % step wise changes in PITCH and YAW in rad!
        r=2*pi/360*[pitchamplitude*sign(sin(2*pi*pitchfrequency*T)') yawamplitude*sign(sin(2*pi*yawfrequency*T)')];

        %% NON-LINEAR SIMULATION (due to constaints on size of input signal)

        % first discretize G and K
        Gd=c2d(G,Ts,'zoh');
        Kd=ss(c2d(K,Ts,'zoh'));

        % start computing, knowing G.D=0
        yd=zeros(Np,2);
        ud=zeros(Np,2);
        xk=zeros(max(size(Gd.A)),1);
        xck=zeros(max(size(Kd.A)),1);
        fprintf('%s','Simulating: 0%,');
        for k=1:Np,
            % compute plant output
            yd(k,:)=[Gd.C*xk]';
            % knowing plant output, we can compute controller output = plant input
            ud(k,:)=[Kd.C*xck + Kd.D*(r(k,:)-yd(k,:))']' + v(k,:);
            % non-linearity due to bound on control signal of 20V
            ud(k,:)=min(abs(ud(k,:)),20).*sign(ud(k,:));
            % knowing plant input = controller output, we can update plant state
            xk=Gd.A*xk + Gd.B*[ud(k,:)'];
            % we can also update controller state
            xck=Kd.A*xck + Kd.B*(r(k,:)-yd(k,:))';
            if round(10*k/Np)==10*k/Np,
                fprintf('%d%%,',round(100*k/Np));
            end
        end
        disp('done...');


        figure(9)
        set(9,'Name','Detailed (Discrete-Time, Constained Input) Closed-Loop Simulation Results')
        % plot output [rad] and reference [rad] in degrees
        subplot(2,1,1);
        l=plot(T,yd(:,1)*360/2/pi,'r',T,r(:,1)*360/2/pi,'b',T,yd(:,2)*360/2/pi,'r--',T,r(:,2)*360/2/pi,'b--');
        set(l,'linewidth',1.5);
        legend('PITCH output','PITCH reference','YAW output','YAW reference');
        title('Detailed (Discrete-Time, Constained Input) Closed-Loop Simulation Results')
        ylabel('pitch and yaw [deg]')
        xlabel('time [sec]')
        subplot(2,1,2);
        l=plot(T,ud(:,1),'m',T,ud(:,2),'m--');
        set(l,'linewidth',1.5);
        legend('PITCH control','YAW control')
        ylabel('pitch and yaw control [V]')
        xlabel('time [sec]')
    end

elseif option==9,

        menu_iteration=0;
        disp('-- Exiting MAELAB --');

    end

end

clear N PITCH_UPDATE YAW_UPDATE SYS1 Tcl Teq Y axiss option v_step t_step
clear cantfindit clpoles f 
clear im iml k kpp kdp kip l lfname
clear m11 m12 m21 m22 mag magcl magl menu_iteration
clear n otion p11 p12 p21 p22 parfile pha phal phacl
clear wn b
clear pitchamplitudep pitchfrequencyp yawamplitudep yawfrequencyp
clear Q11p Q22p Q33p Q44p Q55p Q66p 
clear plttitle r re rel stability testje typo uy v xk xck yn

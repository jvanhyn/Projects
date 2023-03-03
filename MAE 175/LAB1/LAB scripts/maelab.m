% MatLab template file for MAE175A to help you
% (1) simulate open loop (uncontrolled) behavior of a SISO Gyrocope Model
% (2) design controller (internal PID or via rltool) and simulate 
% closed-loop behavior of SISO Gyroscope model and PID controller
% Notes:
% - Edit model transfer functions (default) models file MODELS.M
% - When using rltool or sistool for control design, export the controller 
%   to workspace using the variable name C

% Written by R.A. de Callafon, Dept. of MAE, UCSD (2001-2021)
% Report errors in this software to <callafon@ucsd.edu>

% PLEASE DO NOT MAKE MODIFCATIONS TO THIS FILE
% EDIT THE PARAMETERS OF YOUR MODEL IN THE DEFAULT FILE MODELS.M
% SEE MAE175A LAB MANUAL FOR DETAILS

clc
drawnow;
format compact

disp('- Welcome to MAELAB, v3.2, 2001-2020')

% ask for motor input
typo=1;
while typo==1,
    if exist('motorin')==1,
        motorinep=input(['Specify motor input (1,2) [' num2str(motorin) ']: ']);
        if ~isempty(motorinep),
            motorin=motorinep;
        end
    else
        motorin=input(['Specify motor input (1,2): ']);
        while isempty(motorin),
            clear motorin
            motorin=input(['Specify motor input (1,2): ']);
        end
    end
    if ~((motorin==1)|(motorin==2)),
        clear motorin
    else
        typo=0;
    end
end    

% ask for encoder input
typo=1;
while typo==1,
    if exist('encoderout')==1,
        encoderoutep=input(['Specify encoder output (1,2,3,4) [' num2str(encoderout) ']: ']);
        if ~isempty(encoderoutep),
            encoderout=encoderoutep;
        end
    else
        encoderout=input(['Specify encoder output (1,2,3,4): ']);
        while isempty(encoderout),
            clear encoderout
            encoderout=input(['Specify encoder output (1,2,3,4): ']);
        end
    end
    if ~((encoderout==1)|(encoderout==2)|(encoderout==3)|(encoderout==4)),
        clear encoderout
    else
        typo=0;
    end
end    

% Ask for the filename with the parameters:
cantfindit=1;
while cantfindit,
    cantfindit=0;
    typo=1;
    while typo==1,
        typo=0;
        if exist('parfile')==1,
            parfilep=input(['Name of file with model transfer functions [''' parfile ''']: ']);
            if ~isempty(parfilep),
                parfile=parfilep;
            end
        else
            parfile=input('Name of file with model transfer functions [''models'']: ');
            if isempty(parfile),
                parfile='models';
            end       
        end
        if abs(parfile(1))<65,
            typo=1;
            clear parfilep parfile
            disp(char(7))
            disp(['==> Sorry, incorrect filename, choosing default name ''models''']);            
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
disp(['Evoking ' parfile '.m to update your transfer function models']);
eval(parfile);
        
% check if the transfer function from motorin to encoderout has been defined
eval(['check=(exist(''G' num2str(encoderout) num2str(motorin) ''')==1);']);
if check==0,
    error(['==> Transfer function G' num2str(encoderout) num2str(motorin) ' not defined in ' parfile '.m!'])    
end    

% model is given by a numerator `numg' and denominator `deng'
eval(['[numg,deng] = tfdata(G' num2str(encoderout) num2str(motorin) ');'])
numg=numg{1};
deng=deng{1};

% to include dynamics of the DC-motor (actuator), the total denominator becomes
deng0=deng;
deng=conv(deng,[1/209 1]);

% to include the hardware gain, the total numerator becomes
numg0=numg;
Khw=9.7656e-3; % hardware gain, only plays a role in controller implementation
numg=Khw*numg;


% For Matlab 5.3 and higher in case students like to use root locus method
% for control design:
G=tf(numg,deng); % includes dynamics of mechanics + motor + hardware gain


% In order to simulate Bode plots,
% a frequency vector is specified here:
f=logspace(-2,2,500);      % in Herz, logarithmically spaced



% creating figures

for figs=1:5,
   figure(figs);
   set(figs,'Visible','off');
end

%
% Now we enter the menu.
% Please search throughout the template file for the different
% options to see what is going on
%

menu_iteration=1;
while menu_iteration==1,

   disp(' ');   
   disp(['- Choose an option for model G(s) with input = motor ' num2str(motorin) ' and output = encoder ' num2str(encoderout) ' :']);
   disp('(1) Bode plot of open loop transfer function G(s)');   
   disp('(2) Simulate/compare open loop step response of G(s)');
   disp('(3) Simulate/compare open loop impulse response of G(s)');
   disp('(4) Design/evaluate feedback controller C(s) for G(s)');
   disp('(5) Simulate closed loop step response');
   disp('(6) Exit this menu (to re-specify parameters of your model)');


option=0;
while option==0,
   option=input('Enter your choice: ');
   if isempty(option), option=0; end
   option=round(max(option));
   if option>6, option=0; end
   if option<1, option=0; end
end

if option==1,
    
   disp('- Pole locations and undamped frequency/damping of resonance modes:');
   [Wn,Z]=damp(deng0);
   disp(['  Frequency [Hz]  |  Damping Ratio [0-1]']);
   disp(['  --------------------------------------']);
   format short e
   k=1;
   while k<=length(Wn),
        if Z(k)>0,
            if Z(k)==1,
                fprintf('   %9.4f      |    %9.4f        (stable real valued pole)',Wn(k)/2/pi,Z(k));disp(' ');
            else
                fprintf('   %9.4f      |    %9.4f        (stable complex conjugate pole pair)',Wn(k)/2/pi,Z(k));disp(' ');
            end     
        else
            if Z(k)==-1,
                fprintf('   %9.4f      |    %9.4f        (marginal or unstable real valued pole)',Wn(k)/2/pi,Z(k));disp(' ');
            else
                fprintf('   %9.4f      |    %9.4f        (unstable complex conjugate pole pair)',Wn(k)/2/pi,Z(k));disp(' ');
            end     
        end
        k=k+1;
        if k<=length(Wn),
            if Wn(k)==Wn(k-1),
                k=k+1;
            end
        end
   end    
   disp(['  --------------------------------------']);

   [Wn,Z]=damp(numg0);
   if isempty(Wn)
       disp('- No zero locations or `anti-resonance'' modes in the model G(s)');
   else
       disp('- Zero locations and undamped frequency/damping of `anti-resonance'' modes:');
       disp(['  Frequency [Hz]  |  Damping Ratio [0-1]']);
       disp(['  --------------------------------------']);
       format short e
       k=1;
       while k<=length(Wn),
           if Z(k)>0,
                if Z(k)==1,
                    fprintf('   %9.4f      |    %9.4f        (LHP real valued zero)',Wn(k)/2/pi,Z(k));disp(' ');
                else
                    fprintf('   %9.4f      |    %9.4f        (LHP complex conjugate zero pair)',Wn(k)/2/pi,Z(k));disp(' ');
                end     
           else
                if Z(k)==-1,
                    fprintf('   %9.4f      |    %9.4f        (marginal or RHP real valued zero)',Wn(k)/2/pi,Z(k));disp(' ');
                else
                    fprintf('   %9.4f      |    %9.4f        (RHP complex conjugate zero pair)',Wn(k)/2/pi,Z(k));disp(' ');
                end     
            end
            k=k+1;
            if k<=length(Wn),
                if Wn(k)==Wn(k-1),
                    k=k+1;
                end
            end
       end    
       disp(['  --------------------------------------']);
   end       
        
   %
   % Here we compute frequency response
   %
   [mag,pha]=bode(numg,deng,2*pi*f);
   
   %
   % Here we plot the Bode response response
   %
   figure(1)
   plttitle=['Bode plot of model G(s) [motor ' num2str(motorin) ' -> encoder ' num2str(encoderout) '] with motor dynamics and hardware gain'];
   set(1,'Userdata',plttitle');
   set(1,'Visible','on');
   clf % use to be clg
   subplot(2,1,1);
   l=loglog(f,mag);
   set(l,'linewidth',1.5);
   title(plttitle)
   xlabel('frequency [Hz]');
   ylabel('magnitude')
   grid
   subplot(2,1,2);
   l=semilogx(f,pha);
   set(l,'linewidth',1.5);   
   xlabel('frequency [Hz]');
   ylabel('phase [deg]');
   grid
   disp(['- Figure 1: ' plttitle])   
   
end


if option==2,

   if exist('stepsize')==1,
      stepsizep=input(['Enter (open loop) step size [' num2str(stepsize) ' Volts]: ']);
      if ~isempty(stepsizep),
          stepsize=stepsizep;
      end
   else
      stepsize=input('Enter (open loop) step size [Volts]: ');
      while isempty(stepsize),
          clear stepsize
          stepsize=input('Enter (open loop) step size [Volts]: ');
      end
   end
   if exist('dwelltime')==1,
      dwelltimep=input(['Enter dwell time [' num2str(dwelltime) ' msec]: ']);
      if ~isempty(dwelltimep),
          dwelltime=dwelltimep;
      end
   else
      dwelltime=input('Enter dwell time [msec]: ');
      while isempty(dwelltime),
          clear dwelltime
          dwelltime=input('Enter dwell time [msec]: ');
      end
   end
   stepsize=abs(stepsize(1));
   dwelltime=abs(dwelltime(1));

   typo=1;
   while typo==1,
      typo=0;
      if exist('sfilename')~=1,
        sfilename=input('File that contains open-loop experiment [ENTER for none]: ');
      else
        osfilename=sfilename; 
        sfilename=input(['File that contains open-loop experiment [ENTER for ''' sfilename ''' or 0 for none]: ']);
        if length(sfilename)==0,
          sfilename=osfilename;
        end
      end  
      if length(sfilename)~=0,
          if abs(sfilename(1))==0,
              clear sfilename osfilename              
          elseif abs(sfilename(1))<65
              typo=1;
              if exist('osfilename')==1,
                  sfilename=osfilename;
              else
                  clear sfilename
              end
          end
      else
          clear sfilename
      end
      if exist('sfilename')==1,
        dot_index=find('.'==sfilename);        
        % default, look for the m file with data if no extension was specified
        if isempty(dot_index),
            sfilename=[sfilename '.m'];
            eval(['testje=exist(''' sfilename ''');']);
            if testje~=2,
                disp(['==> Sorry, cannot find Matlab m-file ' sfilename ]);
                typo=1;
                clear sfilename osfilename              
            end
        else
            eval(['testje=exist(''' sfilename ''');']);
            if testje~=2,
                disp(['==> Sorry, cannot find ECP raw data file ' sfilename ]);
                typo=1;
                clear sfilename osfilename              
            end
        end
      end            
   end        

   if exist('sfilename')==1,
       
      lfname=length(sfilename);
      if (sfilename(lfname)=='m')&(sfilename(lfname-1)=='.'),
          % we are working with an m-file
          disp(['- Evoking ' sfilename ' to define t and y from open loop experiment']);
          % remove the '.m' from the filename again
          lfname=length(sfilename);
          sfilename=sfilename(1:lfname-2);
          clear t y
          eval(sfilename);
          if (exist('y')~=1)|(exist('t')~=1),
            error('incorrect file format: no output y or time vector t defined!');
            t=[];y=[];
          end      
      else
          % we must be working with a ECP raw data file
          encoderpos=['Encoder ' num2str(encoderout) ' Pos'];
          data=readecp(sfilename,'Time',encoderpos);
          t=data(:,1);
          y=data(:,2);
          clear data
      end

   else
       
       t=[];y=[];
   
   end
   
   disp('- Computing step response of G(s)');
   %
   % Here we compute step response
   %
   tstep=linspace(0,2*1e-3*dwelltime,900)';
   ustep=[stepsize*ones(450,1);0*ones(450,1)];   
   % the step response without the hardware gain Khw that only plays a role in
   % feedback controller implementation
   ystep=lsim(numg0,deng,ustep,tstep);

   %
   % Here we plot the time response
   %
   figure(2)
   plttitle=['Step response of model G(s) [motor ' num2str(motorin) ' -> encoder ' num2str(encoderout) '] with motor dynamics'];
   set(2,'Userdata',plttitle');
   set(2,'Visible','on');
   clf % use to be clg

   l=plot(tstep,ystep,t,y,'r:');
   set(l,'linewidth',1.5);
   title(plttitle);
   if length(y)>0,
       legend('simulated step response','measured step response')
   else
       legend('simulated step response')
   end
   xlabel('time [sec]');
   ylabel('output signals [counts]')
   disp(['- Figure 2: ' plttitle]);

end


if option==3,
    
   if exist('imptime')==1,
      imptimep=input(['Enter length of impulse response [' num2str(imptime) ' msec]: ']);
      if ~isempty(imptimep),
          imptime=imptimep;
      end
   else
      imptime=input('Enter length of impulse response [msec]: ');
      while isempty(imptime),
          clear imptime
          imptime=input('Enter length of impulse response [msec]: ');
      end
   end
   if exist('deltime')==1,
      deltimep=input(['Enter delay of impulse excitation [' num2str(deltime) ' msec]: ']);
      if ~isempty(deltimep),
          deltime=deltimep;
      end
   else
      deltime=input('Enter delay of impulse excitation [msec]: ');
      while isempty(deltime),
          clear deltime
          deltime=input('Enter delay of impulse excitation [msec]: ');
      end
   end    
   if exist('impscale')==1,
      impscalep=input(['Enter scaling of impulse response [' num2str(impscale) ']: ']);
      if ~isempty(impscalep),
          impscale=impscalep;
      end
   else
      impscale=input('Enter scaling of impulse response [1]: ');
      if isempty(impscale),
          impscale=1;
      end
   end

   
   
   imptime=abs(imptime(1));
   deltime=abs(deltime(1));
   if deltime>imptime,
       deltime=imptime;
   end
   
   typo=1;
   while typo==1,
      typo=0;
      if exist('ifilename')~=1,
        ifilename=input('Name of file (between ''s) that contains open-loop impulse experiment [ENTER for none]: ');
      else
        oifilename=ifilename; 
        ifilename=input(['Name of file (between ''s) that contains open-loop impulse experiment [ENTER for ''' ifilename ''' or 0 for none]: ']);
        if length(ifilename)==0,
          ifilename=oifilename;
        end
      end  
      if length(ifilename)~=0,
          if abs(ifilename(1))==0,
              clear ifilename oifilename              
          elseif abs(ifilename(1))<65
              typo=1;
              if exist('oifilename')==1,
                  ifilename=oifilename;
              else
                  clear ifilename
              end
          end
      else
          clear ifilename
      end
      if exist('ifilename')==1,
        dot_index=find('.'==ifilename);        
        % default, look for the m file with data if no extension was specified
        if isempty(dot_index),
            ifilename=[ifilename '.m'];
            eval(['testje=exist(''' ifilename ''');']);
            if testje~=2,
                disp(['==> Sorry, cannot find Matlab m-file ' ifilename ]);
                typo=1;
                clear ifilename oifilename              
            end
        else
            eval(['testje=exist(''' ifilename ''');']);
            if testje~=2,
                disp(['==> Sorry, cannot find ECP raw data file ' ifilename ]);
                typo=1;
                clear ifilename oifilename              
            end
        end
      end            
   end        
   
   if exist('ifilename')==1,
       
      lfname=length(ifilename);
      if (ifilename(lfname)=='m')&(ifilename(lfname-1)=='.'),
          % we are working with an m-file
          disp(['- Evoking ' ifilename ' to define t and y from open loop experiment']);
          % remove the '.m' from the filename again
          lfname=length(ifilename);
          ifilename=ifilename(1:lfname-2);
          clear t y
          eval(ifilename);
          if (exist('y')~=1)|(exist('t')~=1),
            error('incorrect file format: no output y or time vector t defined!');
            t=[];y=[];
          end      
      else
          % we must be working with a ECP raw data file
          encoderpos=['Encoder ' num2str(encoderout) ' Pos'];
          data=readecp(ifilename,'Time',encoderpos);
          t=data(:,1);
          y=data(:,2);
          clear data
      end
      
   else
       
       t=[];y=[];
   
   end
   
   disp('- Computing (scaled) impulse response of G(s)');
   %
   % Here we compute impulse response
   %
   timp=linspace(0,1e-3*imptime,900)';  
   % the impulse response without the hardware gain Khw that only plays a role in
   % feedback controller implementation
   yimp=impulse(numg0,deng,timp)*impscale;
   index=find(timp<=deltime*1e-3);
   yimp=[zeros(length(index),1);yimp(1:length(yimp)-length(index))];
   %
   % Here we plot the time response
   %
   figure(2)
   plttitle=['Impulse response of model G(s) [motor ' num2str(motorin) ' -> encoder ' num2str(encoderout) '] with motor dynamics'];
   set(2,'Userdata',plttitle');
   set(2,'Visible','on');
   clf % use to be clg

   l=plot(timp,yimp,t,y,'r:');
   set(l,'linewidth',1.5);   
   if length(y)>0,
       legend('simulated impulse response','measured impulse response')
   else
       legend('simulated impulse response')
   end   
   title(plttitle);
   xlabel('time [sec]');
   ylabel('output signals [counts]')
   disp(['- Figure 2: ' plttitle]);

end




if option==4,
    
   disp('- Model G(s) with motor dynamics also available in the Matlab workspace via the');
   disp('  variable G and can be used with any of your own favorite control design tools.');
   disp('  Here we design a P, PD or PID controller only specifically for ECP.');    
    
    usethiscontrol=0;
    if exist('C')==1,
        usethiscontrol=2;
        while ((usethiscontrol~=0)&(usethiscontrol~=1)),
            usethiscontrol=input('Controller C found in memory: use this controller? [1/0]: ');
            if isempty(usethiscontrol), 
                usethiscontrol=2; 
            end
        end
    end
    
    % Note: tau is needed and introduced to make controller realizable
    %       for the differentiating part
    tau=2*pi*0.0001;
    
    % Ts is the specifed sampling time of the discrete-time controller implementation
    Ts = 0.00442;
    
    
if usethiscontrol==0,        

   disp('- Enter the control design variables for C(s) = kp + kd*s + ki/s');
   if exist('kp')==1,
      kpp=input(['Give proportional gain kp [' num2str(kp) ']: ']);
      if ~isempty(kpp),
          kp=kpp;
      end
   else
      kp=input('Give proportional gain kp: ');
      while isempty(kp),
          clear kp
          kp=input('Give proportional gain kp: ');
      end
   end
   if exist('kd')==1,
      kdp=input(['Give derivative gain kd [' num2str(kd) ']: ']);
      if ~isempty(kdp),
          kd=kdp;
      end
   else
      kd=input('Give derivative gain kd: ');
      while isempty(kd),
          clear kd
          kd=input('Give derivative gain kd: ');
      end
   end
   if exist('ki')==1,
      kip=input(['Give integral gain ki [' num2str(ki) ']: ']);
      if ~isempty(kip),
          ki=kip;
      end
   else
      ki=input('Give integral gain ki: ');
      while isempty(ki),
          clear ki
          ki=input('Give integral gain ki: ');
      end
   end
   kp=kp(1);
   kd=kd(1);
   ki=ki(1);
   

   % General form of controller.
   numc=[kd kp ki];
   denc=[tau 1 0];

   % Clearly, if ki = 0, we can remove the pole and zero at 0 from the
   % the general form of the controller, since we have a PD controller

   if ki==0,
      numc=[kd kp];
      denc=[tau 1];
   end
   
   % For matlab 5.3 and higher
   Cpid=tf(numc,denc);
   
else
    
    % get the information from the controller in workspace
    [numc,denc]=tfdata(C);
    numc=numc{1};
    denc=denc{1};
        
    % For matlab 5.3 and higher (update the controller)
    Cpid=tf(numc,denc);    
   
end   
   
   disp('- Controller transfer function:');
   Cpid
   
if denc(1)==0,
   disp('- Inproper controller! Add additional pole(s) to make controller transfer function (stricly) proper.');
   usethiscontrol=0;
   clear stability
else    
      
   disp('- Computing controller and closed loop Bode plots');   

   % computing Bode plot of controller
   [magc,phac]=bode(numc,denc,2*pi*f);

   % computing loop gain
   [numl,denl]=series(numg,deng,numc,denc);
   % computing Bode plot of loop gain
   [magl,phal]=bode(numl,denl,2*pi*f);
   [rel,iml]=nyquist(numl,denl,2*pi*f);

   % computing closed loop output connection (using negative feedback connection)
   [numcl,dencl]=cloop(numl,denl,-1);
   % computing closed loop input connection (using negative feedback connection)
   dens=dencl;
   nums=conv(numc,deng);
   % computing Bode plot of closed loop transfer function
   [magcl,phacl]=bode(numcl,dencl,2*pi*f);

   % check stability
   stability=1;
   clpoles=roots(dencl);
   if max(real(clpoles))>0,
       stability=0;
       disp('==> CLOSED LOOP UNSTABLE');
   end


   
   %
   % Here we plot the Bode plot of the controller
   %
   figure(3)
   plttitle='Bode plot of controller C(s)';
   set(3,'Userdata',plttitle');
   set(3,'Visible','on');
   clf % use to be clg

   subplot(2,1,1);
   l=loglog(f,magc);
   set(l,'linewidth',1.5);
   title(plttitle)
   xlabel('frequency [Hz]');
   ylabel('magnitude')

   subplot(2,1,2);
   l=semilogx(f,phac);
   set(l,'linewidth',1.5);   
   xlabel('frequency [Hz]');
   ylabel('phase [deg]');
   disp(['- Figure 3: ' plttitle]) 

   %
   % Here we plot the Nyquist plot and closed-loop transfer function
   %
   figure(4)
   plttitle='Frequency plots of loopgain L(s)';
   set(4,'Userdata',plttitle');
   set(4,'Visible','on');
   clf % use to be clg

   subplot(2,2,1);
   l=loglog(f,magl);
   set(l,'linewidth',1.5);   
   hold on
   l=loglog([f(1) f(length(f))],[1 1],'k:');
   set(l,'linewidth',1.5);   
   hold off
   title('magnitude L(s)')
   xlabel('frequency [Hz]');
   ylabel('magnitude')
   
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
   title('Nyquist contour L(s)')
   xlabel('real')
   ylabel('imag')
   text(-5.5,5,'L(s)=G(s)*C(s)','HorizontalAlignment','center')

   subplot(2,2,3);
   l=semilogx(f,phal);
   set(l,'linewidth',1.5);   
   hold on
   l=semilogx([f(1) f(length(f))],[-180 -180],'k:');
   set(l,'linewidth',1.5);   
   hold off
   xlabel('frequency [Hz]');
   ylabel('phase [deg]');
   title('phase L(s)')
   axiss=axis;
   axis([f(1) f(length(f)) axiss(3) axiss(4)]);

   subplot(2,2,4);
   l=loglog(f,magcl);
   set(l,'linewidth',1.5);   
   hold on
   l=loglog([f(1) f(length(f))],[1 1],'k:');
   set(l,'linewidth',1.5);   
   hold off
   title('magnitude L(s)/(1+L(s))  [r(s) -> y(s)]')
   xlabel('frequency [Hz]');
   ylabel('magnitude')
   axiss=axis;
   axis([f(1) f(length(f)) axiss(3) axiss(4)]);

   if stability==1,
      text(0.1,0.1,'closed loop stable','color','green','Units','normalized');
   else
      text(0.1,0.1,'closed loop unstable!','color','red','Units','normalized');
   end

   disp(['- Figure 4: ' plttitle])

end   
   
end

if option==5,    

   if exist('stability')==0,
      disp('==> First run option 4 to design/evaluate the controller ');
   else
       
             % Compute the ZOH discrete-time equivalent of G WITHOUT hardware gain
          % and WITH motor dynamics
          Gd=c2d(tf(numg0,deng),Ts,'zoh');
          [numGd,denGd]=tfdata(Gd);
          numGd=numGd{1};denGd=denGd{1};
          
          % Compute the discrete time equivalent of the continuous-time controller       
          
          if usethiscontrol==0,
              
              % Interesting fact: 
              % the (ideal) PID controller given by
              %  Cpid=tf([kd kp ki],[tau 1 0]);   % 0 < tau << 1
              % can be approximated reasonably well in discrete-time with
              % the recursive formula
              %  error=error+e(k)
              %  u(k) = kp*e(k) + kd/Ts*(e(k)-e(k-1)) + ki*Ts*error
              % which again can be shown to be equivalent to a digital PID
              % controller
              %  Cpidd=tf([kp+kd/Ts+ki*Ts -kp-2*kd/Ts kd/Ts],[1 -1 0],Ts);
              % plotting bode(Cpid,Cpidd) shows the resemblance and is
              % often times better than a FOH approximation c2d(Cpid,Ts,'tustin') 
              % that has better phase properties but creates a pole at -1 causing
              % high frequency gain and excessive oscillation in the control signal
              % or a ZOH approximation c2d(Cpid,Ts,'zoh') that creates a similar 
              % phase but make notches move, change the amplitude or will
              % even be uncomputable for an ideal PID controller
                            
              numcd=[kp+kd/Ts+ki*Ts -kp-2*kd/Ts kd/Ts];
              dencd=[1 -1 0];
              
              % It should be noted that the actual applied discrete-time controller in ECP
              % is done based on a Kp and Ki on r-y and Kd on ydiff only and is done as follows:
              % u(k)=kp*[r(k-1)-y(k-1)]-kd/Ts*[y(k-1)-y(k-2)]-ki*Ts*sum_{m=1}^{k-1}[y(m)-r(m)]
              % Note the extra one step time delay in the control algorithm!
              
              % add the extra time delay to make controller realizable in
              % discrete time in the ECP software
              
              numcd=[0 numcd];
              dencd=[dencd 0];
                            
          else
              
              % Note we used matched poles and zeros, as ZOH creates a complete
              % different controller if you have notches with more phase delay
              % and TUSTIN leads to a controller with a high frequency pole
              % (gain) causing high frequency oscillations in the control
              % signal. It seems that matches poles/zeros is the best for this
              % particular application
                        
              [numcd,dencd]=tfdata(c2d(Cpid,Ts,'matched'));
              numcd=numcd{1};dencd=dencd{1};
              
              % add the extra time delay to make controller realizable in
              % discrete time in the ECP software
              
              if numcd~=0,
              
                  numcd=[0 numcd];
                  dencd=[dencd 0];
                  
              end
                            
          end
              
          % display discrete-time controller implementation being used
          disp(['- Real-time controller implementation with e(k) = r(k) - y(k) and Ts = ' num2str(Ts) 'sec']);
          rtstr=[];
          if usethiscontrol==0,
              if (kp~=0),
                if kp>0,
                    rtstr=[rtstr '+'];
                end
                rtstr=[rtstr sprintf('%0.5g',kp) '*e(k-1)'];
              end
              if (kd~=0),                            
                if kd>0,
                    rtstr=[rtstr '+'];
                end
                rtstr=[rtstr sprintf('%0.5g',kd/Ts) '*[y(k-1)-y(k-2)]'];
              end
              if (ki~=0),                            
                if ki>0,
                    rtstr=[rtstr '+'];
                end
                rtstr=[rtstr sprintf('%0.5g',ki*Ts) '*sum_{t=0}^{k-1} e(t)'];
              end
          else
              for k=1:length(numcd),
                  if numcd(k)~=0,
                      if numcd(k)>0,
                          rtstr=[rtstr '+'];
                      end
                      if abs(numcd(k))==1,
                          rtstr=[rtstr 'e(k-' num2str(k-1) ')'];
                      else
                          rtstr=[rtstr sprintf('%0.5g',numcd(k)) '*e(k-' num2str(k-1) ')'];
                      end
                  end
              end
              for k=2:length(dencd),
                  if dencd(k)~=0,
                        if -dencd(k)>0,
                            rtstr=[rtstr '+'];
                        end    
                        if abs(dencd(k))==1,
                            rtstr=[rtstr 'u(k-' num2str(k-1) ')'];
                        else
                            rtstr=[rtstr sprintf('%0.5g',-dencd(k)) '*u(k-' num2str(k-1) ')'];
                        end
                  end
              end
          end
          if isempty(rtstr),
                rtstr='  u(k) = 0';                  
          else    
              rtstr=['  u(k) = ' rtstr];                  
          end
          
          disp(rtstr)        
       
        disp('- Checking discrete-time closed-loop stability of controller');             
          
        % recheck stability in case you changed the encoder output...       
        % computing loop gain
        [numl,denl]=series(numGd,denGd,Khw*numcd,dencd);
        
        % computing closed loop output connection (using negative feedback connection)
        [numcl,dencl]=cloop(numl,denl,-1);

        % check stability
        stability=1;
        clpoles=roots(dencl);
        if max(abs(clpoles))>1.01,
            abs(clpoles)
            stability=0;            
        end
    
      if stability==1,
          if exist('clstepsize')==1,
              clstepsizep=input(['Enter (closed loop) step size [' num2str(clstepsize) ' counts]: ']);
              if ~isempty(clstepsizep),
                  clstepsize=clstepsizep;
              end
          else
              clstepsize=input('Enter (closed loop) step size [counts]: ');
              while isempty(clstepsize),
                  clear clstepsize
                  clstepsize=input('Enter (closed loop) step size [counts]: ');
              end
          end
          if exist('cldwelltime')==1,
              cldwelltimep=input(['Enter dwell time [' num2str(cldwelltime) ' msec]: ']);
              if ~isempty(cldwelltimep),
                  cldwelltime=cldwelltimep;
              end
          else
              cldwelltime=input('Enter dwell time [msec]: ');
              while isempty(cldwelltime),
                  clear cldwelltime
                  cldwelltime=input('Enter dwell time [msec]: ');
              end
          end
          clstepsize=abs(clstepsize(1));
          cldwelltime=abs(cldwelltime(1));
   
          cltstep=[0:Ts:1e-3*cldwelltime]';
          clustep=[clstepsize*ones(floor(length(cltstep)/2),1);0*ones(length(cltstep)-floor(length(cltstep)/2),1)]; 
          
          
          
          typo=1;
          while typo==1,
            typo=0;
            if exist('s2filename')~=1,
                s2filename=input('File that contains closed-loop experiment [ENTER for none]: ');
            else
                os2filename=s2filename; 
                s2filename=input(['File that contains closed-loop experiment [ENTER for ''' s2filename ''' or 0 for none]: ']);
                if length(s2filename)==0,
                    s2filename=os2filename;
                end
            end  
            if length(s2filename)~=0,
                if abs(s2filename(1))==0,
                    clear s2filename os2filename              
                elseif abs(s2filename(1))<65
                    typo=1;
                    if exist('os2filename')==1,
                        s2filename=os2filename;
                    else
                        clear s2filename
                    end
                end
            else
                clear s2filename
            end
            if exist('s2filename')==1,
                dot_index=find('.'==s2filename);        
                % default, look for the m file with data if no extension was specified
                if isempty(dot_index),
                    s2filename=[s2filename '.m'];
                    eval(['testje=exist(''' s2filename ''');']);
                    if testje~=2,
                        disp(['==> Sorry, cannot find Matlab m-file ' s2filename ]);
                        typo=1;
                        clear s2filename os2filename              
                    end
                else
                    eval(['testje=exist(''' s2filename ''');']);
                    if testje~=2,
                        disp(['==> Sorry, cannot find ECP raw data file ' s2filename ]);
                        typo=1;
                        clear s2filename os2filename              
                    end
                end
            end            
          end

          if exist('s2filename')==1,
            lfname=length(s2filename);
            if (s2filename(lfname)=='m')&(s2filename(lfname-1)=='.'),
                % we are working with an m-file
                disp(['- Evoking ' s2filename ' to define t, u and y from step experiment']);
                % remove the '.m' from the filename again
                lfname=length(s2filename);
                s2filename=s2filename(1:lfname-2);
                clear t u y
                eval(s2filename);
                if (exist('t')~=1)|(exist('u')~=1)|(exist('y')~=1),
                    error('incorrect file format: no time vector t, input u and output y defined!');
                    t=[];y=[];
                end      
            else
                % we must be working with a ECP raw data file
                encoderpos=['Encoder ' num2str(encoderout) ' Pos'];
                motorvoltage=['Control Effort ' num2str(motorin)];
                data=readecp(s2filename,'Time',motorvoltage,encoderpos);
                t=data(:,1);
                u=data(:,2);
                y=data(:,3);
                clear data
            end
          else
            t=[];u=[];y=[];   
          end
          
          % Now simulate in discrete time with input constraints using ZOH equivalent of model G

          % add some zeros on input and output for initial conditions
          N=length(cltstep);
          nz=max([length(numGd) length(denGd) length(numcd) length(dencd)]);
          ystepcl=zeros(nz-1+N,1);
          ustepcl=zeros(nz-1+N,1);
          % used for the INTERNAL state (past outputs) of the controller:
          ustepcl_internal=zeros(nz-1+N,1);  
          % it is important to make a distinction in saturation of state
          % variables (past outputs) and an actual (static) output
          % saturation, as they lead to different simulation results!
          % one is with static nonlinearity at the output, the other one is
          % with a static nonlinearity in the internal feedback
          clustep=[zeros(nz-1,1);clustep];
          % reset integrator state
          integral_error=0;
          % specify control saturation limit
          maxu=10;
          for k=nz:N+nz-1,
              % compute discrete output, assuming we have AT LEAST ONE STEP DELAY in Gd (numGd(1)=0) and denGd(1)=1
              ystepcl(k)=numGd(2:length(numGd))*ustepcl(k-1:-1:k-length(numGd)+1)-denGd(2:length(denGd))*ystepcl(k-1:-1:k-length(denGd)+1);
              % compute control input 
              if usethiscontrol==0,
                  % use the straightforward application of
                  %  error=error+e(k-1)
                  %  u(k) = kp*e(k-1) + kd/Ts*(e(k-1)-e(k-2)) + ki*Ts*error
                  % where an extra step time delay has been added for
                  % real-time implementation
                  integral_error = integral_error + clustep(k-1)-ystepcl(k-1);
                  %ustepcl_internal(k) = Khw*kp*[clustep(k-1)-ystepcl(k-1)] + Khw*kd/Ts*[clustep(k-1)-ystepcl(k-1)-(clustep(k-2)-ystepcl(k-2))] + Khw*ki*Ts*integral_error;                                        
                  % use the formula below if you are using velocity feedback or kd only on output, not on error signal
                  ustepcl_internal(k) = Khw*kp*[clustep(k-1)-ystepcl(k-1)] + Khw*kd/Ts*[-ystepcl(k-1)+ystepcl(k-2)] + Khw*ki*Ts*integral_error;                                        
                  
              else
                  % more advanced controller specified by a discrete numerator numcd and denominator dencd
                  ustepcl_internal(k) = Khw*numcd*[clustep(k:-1:k-length(numcd)+1)-ystepcl(k:-1:k-length(numcd)+1)] - dencd(2:length(dencd))*ustepcl_internal(k-1:-1:k-length(dencd)+1);
              end
              % saturation
              if abs(ustepcl_internal(k))>maxu, 
                  ustepcl(k)=sign(ustepcl_internal(k))*maxu;
              else
                  ustepcl(k)=ustepcl_internal(k);
              end
          end
          % remove zeros again
          clustep=clustep(nz:N+nz-1);
          ystepcl=ystepcl(nz:N+nz-1);
          ustepcl=ustepcl(nz:N+nz-1);

          % this would be the continuous time closed loop step response of output and input
          %ystepcl=lsim(numcl,dencl,clustep,cltstep);
          %ustepcl=Khw*lsim(nums,dens,clustep,cltstep);
          
          figure(5)

          plttitle=['Closed-loop step of model G(s) [motor ' num2str(motorin) ' -> encoder ' num2str(encoderout) '] with motor dynamics'];

          set(5,'Userdata',plttitle');
          set(5,'Visible','on');
          clf % use to be clg
          subplot(2,1,1);
          l=plot(cltstep,clustep,'g-',cltstep,ystepcl,'b-',t,y,'r:');
          set(l,'linewidth',1.5);          
          if length(y)>0,
                legend('commanded reference','simulated closed-loop step','measured closed-loop step');
          else
                legend('commanded reference','simulated closed-loop step');
          end
          title(plttitle);
          ylabel('output & reference [counts]')
          subplot(2,1,2);
          l=plot(cltstep,ustepcl,t,u,'r:');
          set(l,'linewidth',1.5);          
          xlabel('time [sec]');
          ylabel('control signal [V]')
          if usethiscontrol==0,
            title(['Control signal for kp = ' num2str(kp) ', kd = ' num2str(kd) ', ki = ' num2str(ki)]);
          else
            title(['Control signal for user defined controller C'])
          end
          
          %axis([cltstep(1) cltstep(length(cltstep)) -5 5]);

          disp(['- Figure 5: ' plttitle])          
      else
          disp('==> CLOSED LOOP SYSTEM UNSTABLE');
          disp(['    Redesign controller around encoder ' num2str(encoderout) ' with option 4']);
          set(5,'Visible','off')
      end
   end
end

if option==6,

   disp('- Exiting MAELAB, rerun to read model parameters...');
   menu_iteration=0;

end

end

clear typo cantfindit deltatimep deltimep imptimep index lfname figs numg deng f
clear option parfilep testje mag pha plttitle dwelltimep stepsizep
clear axiss numc denc numl denl
clear rel iml magc magcl magl phac phacl phal tau
clear kpp kdp kip



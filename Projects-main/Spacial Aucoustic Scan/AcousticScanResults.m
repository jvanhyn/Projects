
close all
clear
clc
format long G
load('acousticscan738647612.mat')
%% System Details
fs = 1/dt;
%% Analysing Data at Control Point 
% Examining the scan at a single point will tell us how to proceed 
% with analysing points across the spacial scan. If we record the signal
% responce at this point five times and average the result, we are able
% to get a clearer picture. When comparing this average to one raw 
% recorded signal it becomes obvious can tell that there is an extremely
% high noise to signal ratio. Some signal process will be required to 
% extract precice information about the' "signal's" 'arrival time at the
% microphone.

figure(1)
hold on
plot(t*10^3,wave_sig)                               % Raw first signal 
plot(t*10^3,recMatrix_sig,'-','LineWidth',2)        % Average across frist 5 samples
legend('Final Measured Signal','Averaged Signal')
title('Magnitude vs Time')
ylabel('Signal Magnitude')
xlabel('Time (ms)')
hold off

txt1 = ['Examining the scan at a single point will tell us how to proceed ' ...
    'with analysing points across the spacial scan. If we record the signal' ...
    'responce at this point five times and average the result, we are able' ...
    'to get a clearer picture. When comparing this average to one raw ' ...
    'recorded signal it becomes obvious can tell that there is an extremely' ...
    'high noise to signal ratio. Some signal process will be required to ' ...
    'extract precice information about the' "signal's" 'arrival time at the' ...
    'microphone.'];

%% Distance to Control Point 

x0 = 41.5*10^-2;        % X-axis offset at (0,0)
y0 = 9.5*10^-2;         % Y-axis offset at (0,0)
t0 = 1*10^-3;           % Time delay signal emitted

D = sqrt(x0^2+y0^2);    % Distance from speaker to Microphone at (0,0)

dx = 23.5*10^-2;        % Total X Travel
dy = 9.5*10^-2;         % Total Y Travel
x = linspace(0,dx,30);  % X coordinates
y = linspace(0,dy,15);  % Y coordinates

Dx = sqrt((x0-x).^2+y0.^2);
Dy = sqrt(x0.^2+(y0-y).^2);


%% Filtering the Signal
% A bandpass filter is used to extract frequencies near 5kHz.

mean_sig_pb = bandpass(recMatrix_sig,[4500,6000], fs); 

figure(2)
bandpass(recMatrix_sig,[4500,6000], fs); 

%% Detecting Control Point Arrival Time
% The signal is normalized and a threshold trigger must be tuned to find
% the arival time

mean_sig_pbn = abs(mean_sig_pb/max(mean_sig_pb));
thresh1 = 0.25;
atime = 0.0;
ai = 0;

for i = 1:length(mean_sig_pbn)
    if mean_sig_pbn(i) - thresh1 >= 0 
        atime = t(i);
        ai = i;
        break;
    end
end 

figure(3)
hold on
plot(t,mean_sig_pbn,'-','LineWidth',1)
xline(atime,'--','LineWidth',1.5)
yline(mean_sig_pbn(ai),'--','LineWidth',1.5)
yline(thresh1,'g')
grid on
hold off

%% Calculating Control Speed of Sound

sound_speed_control = D/(atime-t0)
error = D/(atime-t0-dt)-D/(atime-t0)

%% Analysing the 2D Scan
load('acousticscan738647646.mat')

time0 = 101;        % index of time delay t0
timef = 500;        % 
window = time0:timef;
twindow = t(window);
sig = recMatrix_sig(window,:,:);



%%
dthresh = 0;
thresh = 0.7 + dthresh*ones(1,15);

onset = 50;
offset = 100;

clear atime ai s spb 
for k = 1:30
    for j = 1:15
        [atime,ai,s] = FrequencyArrivalTime(t,sig(:,k,j),thresh(k),onset,offset);
        trig(k,j) = atime;
        trig_index(k,j) = ai;
        sbp(:,k,j) = s;
        
    end 
end

arrivalTime = trig;

%%
figure(5)
heatmap(trig)
%%
row = 1;
col = 11;

figure(6)
hold on
plot(twindow,sbp(:,row,col),'-','LineWidth',1)
xline(arrivalTime(row,col)+t0,'--','LineWidth',1.5)
yline(sbp(trig_index(row,col),row,col),'--','LineWidth',1.5)
yline(thresh(row),'g')
grid on
hold off

%% 


f1 = fit(Dx(2:end)',arrivalTime(2:end,1),'poly1');
l1 = f1.p1*Dx(2:end) + f1.p2;
l2 = (Dx(2:end))/sound_speed_control + f1.p1*Dx(2) + f1.p2 - (Dx(2))/sound_speed_control 

buff = 0.01;
xneg = 0.001*ones(1,29);
xpos = 0.001*ones(1,29);
yneg = 0.0001*ones(1,29);
ypos = 0.0001*ones(1,29);


figure(4)
hold on
errorbar(Dx(2:end),arrivalTime(2:end,1),yneg,ypos,xneg,xpos,'.')
plot(Dx(2:end),l1,'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
plot(Dx(2:end),l2)
%%

for k = 1:30
    for j = 1:15
        [s,w,time] = spectrogram(sig(:,k,j),35,29,128,fs,'yaxis');

        sp(:,k,j) = smooth(abs(s(7,:)));
        spd = diff(sp);

        trigger1 = 0;
        trigger2 = 0;

        amps1(:,k,j) = (recMatrix_sig(:,k,j))/abs(max(recMatrix_sig(:,k,j)));
        amps2(:,k,j) = smooth(amps1(:,k,j));
        amps3(:,k,j) = sp(:,k,j)/max(sp(:,k,j));

        for i = 1:(numel(time))
            if (sp(i,k,j)/max(sp(:,k,j)) - thresh) >= 0
                trigger1 = time(i);
                i;
                break
            elseif i == (numel(time)-1)
                disp('threshold not met')
                fail(k) = i;
            end
        end

        for i = 1:(numel(time)-1)
            if (spd(i,k,j)/max(spd(:,k,j))- thresh) >= 0
                trigger2 = time(i);
                i;
                break
            elseif i == (numel(time)-1)
                disp('threshold not met')
                faild(k) = i;
            end
        end

        trig1(k,j) = trigger1;
        trig2(k,j) = trigger2;
    end
end
%%
heatmap(trig2)

%%
xstep = (0.7833333333)*10^-2;
ystep = (0.7333333333)*10^-2;

x0 = 41.5*10^-2;
y0 = 9.5*10^-2;

dx = 23.5*10^-2;
dy = 9.5*10^-2;
x = linspace(0,dx,30);
y = linspace(0,dy,15);
xax = x;
yax = y;

D = mean(sqrt(x0^2+(y0-y).^2));
dt1 = mean(trig1,2);
dt2 = mean(trig2,2);

f1 = fit(x',dt1,'poly1');
f2 = fit(x',dt2,'poly1','Exclude',y<(1.6*10^-3));
c1 = 1/f1.p1;
c2 = 1/f2.p1;
ct = 1/340 ;
l1 = f1.p1*x + f1.p2;
l2 = f2.p1*x + f2.p2;
lt = - ct*(x) + D*ct;


buff = 0.01;
xneg = 0.001*ones(1,30);
xpos = 0.001*ones(1,30);
yneg = 0.0001*ones(1,30);
ypos = 0.0001*ones(1,30);


figure(4)
hold on
errorbar(xax(1:end),dt2(1:end),yneg,ypos,xneg,xpos,'.')
plot(xax,l2,'--','LineWidth',2,'Color',[0.9290, 0.6940, 0.1250])
plot(xax,lt,'LineWidth',2,'Color',[0.6350, 0.0780, 0.1840])
xlabel('Distance From Origin (m)')
ylabel('Time (s)')
grid on
xlim([xax(1)-buff xax(end)+buff])
legend('Recorded Delay Times','Fit to Data','Theoretical Delay in Air')
%% section 5 

%% Section 6
figure(6)
pcolor(xax,t-0.001,amps1(:,:,15));
shading flat
c = colorbar;
set(c,'FontSize', 20);
hold on 
a = plot(xax,lt,'LineWidth',4,'Color',[0.9, 0.3, 0.1]);
legend(a,'Theoretical Time Delay in Air')
xlabel('Distance From Origin (m)')
ylabel('Time (s)')
title('Normalized Signal Amplitude vs. X-distance and Time')
caxis([-1 1])
ylim([0 5*10^-3])
hold off
%% Section 7
figure(7)
pcolor(xax,yax,(reshape(amps1(200,:,:),[30 15]))');
title(['Normalized Spacial Pressure Distribution (' num2str(20/10) 'ms)']);
xlabel('Distance From Origin (m)')
ylabel('Distance From Origin (m)')
c = colorbar;
set(c,'FontSize', 20);
c.Label.String = 'Amplitude'
c.Label.FontSize = 12;
caxis([-1 1])
shading flat;
%% Section 8
filename = 'scanning.gif';
for k=1:50
    h = figure(8);
    pcolor(xax,yax,(reshape(amps1(k*10,:,:),[30 15]))');
    title(['Normalized Spacial Pressure Distribution (' num2str(k/10) 'ms)']);
    xlabel('Distance From Origin (m)')
    ylabel('Distance From Origin (m)')
    c = colorbar;
    set(c,'FontSize', 20);
    c.Label.String = 'Amplitude'
    c.Label.FontSize = 12;
    shading flat;
    caxis([-1 1])
    pause(0.5);
    drawnow;
    frame = getframe(h); % You can either remove the argument h
    % or define a figure above, e.g., h=figure(01);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if k == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
%% 10
xstep*10^-3/(343/5000)

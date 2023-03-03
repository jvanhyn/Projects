%% 
clc, clear, clear all

%%%%%%%%%% STEP INPUTS
% () corresponds to pitch data/sheet
% (2) corresponds to yaw data/sheet

% File/Sheet name
filename = 'LAB2_part1';
sheetname = 'Step_P1';
sheetname2 = 'Step_Y1';

% Range of Excel values
range = 'A2:D3640';
range2 = 'A2:D3640';

data = xlsread(filename, sheetname, range);
data2 = xlsread(filename, sheetname2, range2);

time = data(:,1);
time2 = data2(:,1);

amp_recorded = data(:,4);
amp_desired = data(:,2);
amp_recorded2 = data2(:,4);
amp_desired2 = data2(:,2);

figure(1)
subplot(2,1,1)

plot(time,amp_recorded)
xlabel('time(s)')
ylabel('amplitude(deg)')
title('STEP INPUT (NO CONTROLLER)')
hold on
plot(time,amp_desired)
legend('PITCH measured amplitude', 'PITCH desired amplitude')
xlim([0 10])

subplot(2,1,2)
xlabel('time(s)')
ylabel('amplitude(deg)')
plot(time2,amp_recorded2)
hold on
plot(time2,amp_desired2)
xlabel('time(s)')
ylabel('amplitude(deg)')
legend('YAW measured amplitude', 'YAW desired amplitude')
xlim([0 10])



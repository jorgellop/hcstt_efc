% test_DriftFIU
% 
% Save powermeter every few seconds to asses the drift of the optics 
%
% Jorge Llop - Apr 7, 2018

clear all;
close all;

addpath(genpath('utils'));

label = '_0503';
outDir = ['output',filesep,'test_DriftFIU_Fast',label,filesep];
mkdir(outDir);

hcstt_Initialize(true);

count = 1;
tic
% while 1
for II=1:1000
    disp(['Read num ', num2str(count)])
    int_arr(count) = hcstt_GetIntensityFIU(zeros(12,12),1);
    disp([num2str(int_arr(count)),'W'])
%     pause(13)
    elapsedTime_arr(count) = toc;
    count = count + 1;
end
time_arr = 1:count-1;
time_arr = time_arr*15/60;
figure(2)
plot(elapsedTime_arr/60,int_arr)
xlabel('time (min)')
ylabel('Intensity')
title('FIU Drift Test - Tight DM')
export_fig([outDir,'DriftTestPlot_tightDM',label,'.png'],'-r300');
hcstt_DisconnectDevices();

% test_applySpeckleOnFiber
%
% 
%
% Jorge Llop, May 10, 2018

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(true);
hcstt_NewFlatForDM('ImageSharpening_fminconIt2_Apr1');
totalPower = 1.08e-6;
peakInt = totalPower;
numtry = 10;

amp_arr = linspace(20,30,numtry);

for II=1:numtry
    DM_Command = hcstt_DMMapSin(amp_arr(II),0,2.4,0);
    int_arr(II) = hcstt_GetIntensityFIU(DM_Command,20);
end

figure(100)
plot(1:numtry,int_arr/peakInt)
hcstt_DisconnectDevices();

%%
DM_Command = hcstt_DMMapSin(26.125,0,2.4,0);
load('flat.mat');
flat_SN = DM_Command + flat;
save('utils\NewFlat_SpeckleOnFiber_May10.mat','flat_SN')
hcstt_NewFlatForDM('NewFlat_SpeckleOnFiber_May10')
int = hcstt_GetIntensityFIU(zeros(12,12),20)
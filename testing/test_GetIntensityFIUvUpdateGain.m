% test_GetIntensityFIUvUpdateGain
%
% Test new code to get intensity from NI powermeter by Ben Calvin
%
% Jorge Llop - Aug 22, 2018

hcstt_Initialize(true)
while 1
    disp(hcstt_GetIntensityFIUvUpdateGain(zeros(12,12),1))
end
hcstt_DisconnectDevices
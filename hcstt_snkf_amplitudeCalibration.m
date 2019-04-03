% hcstt_snkf_amplitudeCalibration.m
% 
% Apply sinusoid to the DM to send speckles at a range of amplitut]des to
% the SMF to calibrate amplitude
%
% Rebecca Zhang, and Jorge Llop - Feb 25, 2019

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

label = '_Feb25';

% Initialize all devices
hcstt_Initialize(true);

numamp = 50;
amp_arr = linspace(0,25,numamp);

% [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(x_fib_est,0,info)
actxc_fib =2.5500;
ang_fib = -0.0700;

n_pairs = 5;
probe_phase = linspace(0,2*pi * 5/6, n_pairs);
amp_cal = zeros(numamp,4);
for II=1:numamp
    intensity = zeros(n_pairs,1);

    for n=1:1:n_pairs
        dm_actuators_mat0 = hcstt_DMMapSin(amp_arr(II), ang_fib, actxc_fib, probe_phase(n));
        intensity(n) = hcstt_GetIntensityFIU(dm_actuators_mat0(:),5,0);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
    end

    %this part calculates phase from measurements
    t  = 1:n_pairs;                                                 % Time Vector
    sig  = intensity;                                                 % Signal Vector
    Ts = mean(diff(t));                                         % Sampling Time
    Fs = 1/Ts;                                                  % Sampling Frequency
    Fn = Fs/2;                                                  % Nyquist Frequency
    L  = length(sig);
    fts = fft(sig)/L;                                             % Normalised Fourier Transform
    Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
    Iv = 1:length(Fv);                                          % Index Vector
    amp_fts = abs(fts(Iv))*2;                                   % Spectrum Amplitude
    %phs_fts = angle(fts(Iv));                                   % Spectrum Phase

    int =  hcstt_GetIntensityFIU(zeros(12,12),5,0);
    amp_cal(II,2) = (amp_fts(1)-2*int)/2;
    amp_cal(II,3) = amp_fts(2);

end
amp_cal(:,4) = amp_cal(:,2)+amp_cal(:,3);
amp_cal(:,1) = amp_arr';

plot(1:numel(amp_cal(:,2)),amp_cal(:,2))
save(['output',filesep,'amp_cal_',label],'amp_cal');

hcstt_DisconnectDevices()

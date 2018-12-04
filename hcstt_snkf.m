% hcstt_snkf
%
% Speckle Nulling with a Kalman Filter. Main code extracted from
% kf_null_femto.m by Yinzi Xin
%
% Yinzi Xin - 2018
% Jorge Llop - November 29, 2018

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

label = '_test_Nov29';
outDir = ['output',filesep,'EFC_wFiber_LabDemonstration',label,filesep];
mkdir(outDir);

% Load calibration data
load(['Calibrations',filesep,'amp_cal_new_2p4.mat']);
load(['output',filesep,'calibrateDM_Aug01']); % actxcDM, angDm vs pix on camera info
load(['Calibrations',filesep,'angle_cal']);
load(['Calibrations',filesep,'position_cal']);

% Initialize all devices
hcstt_Initialize(true);

num_steps = 10;%30 in Yinzi's code

n_pairs = 4; %initialize with 4 probes (for phase)
%generate equally spaced probes (4 different phis)
probe_phase = [0, 3.14/2, 3.14, 3.14*1.5];
probe_amp = 3;

steps = 1:num_steps;


%Find position of fiber
x_fib_est = 2.5;
% [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(x_fib_est,0,info);
actxc_fib = 2.5;
ang_fib = 0;

% Take background image
take_background = true;
Ncam = 400;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    backgroundCam = im_cam;
    backgroundSMF = hcstt_GetIntensityFIU(zeros(Nact,Nact),15,0);
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    backgroundCam = zeros(Ncam,Ncam);
    backgroundSMF = 0;
end

%load interaction matrix
gamma = [1,0;0,1];

%covariance initialization
P_0 = eye(2);

%process noise (on x)
w_k = [.1,0;0,.1]*10^-5;
Q_k = w_k*w_k';

%measurement noise
phase_noise = 1*10^-5;
amplitude_noise = 0.5*10^-5;
R_k = [phase_noise^2,0;0,amplitude_noise^2];

%start looping here

%initialize values
x_k_m = [0;5]; %prior estimates %JLlop: from meas
P_k_m = P_0; %prior covariance
u_k = [0;0]; %control input

%initialize logs
measurements = zeros(num_steps,2);
estimations = zeros(num_steps,2);
inputs = zeros(2,num_steps);

for i=steps
    
    intensity = zeros(n_pairs,1);
    
    H_k = eye(2); %observation matrix
    
    for n=1:1:n_pairs
        par_probe = [probe_amp, ang_fib, actxc_fib, probe_phase(n)];
        dm_actuators_mat0 = hcstt_DMMapSin(1, ang_fib, actxc_fib, ph_arr(KK));
%         intensity_raw = DE_supMinimizer_3(par_probe);
        intensity_raw = hcstt_GetIntensityFIU((+dm_actuators_mat0(:) + 0)/1e-9,5,backgroundSMF);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm

        %simulate gaussian detector noise
        noise_sigma = 0.5*10^-15;
        %draw from gaussian
        noise_amp = normrnd(0,noise_sigma);
        intensity(n) = intensity_raw+noise_amp;
    end
    
    pause(.5)

    %calculate phase and amplitude from measurements
    t  = 1:n_pairs;                                                 % Time Vector
    sig  = intensity;                                                 % Signal Vector
    Ts = mean(diff(t));                                         % Sampling Time
    Fs = 1/Ts;                                                  % Sampling Frequency
    Fn = Fs/2;                                                  % Nyquist Frequency
    L  = length(sig);
    fts = fft(sig)/L;                                             % Normalised Fourier Transform
    Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
    Iv = 1:length(Fv);                                          % Index Vector
    amp_fts = abs(fts(1:3))*2;                                   % Spectrum Amplitude
    phs_fts = angle(fts(1:3));                                   % Spectrum Phase
    
    if phs_fts(2)<0
       phs_fts(2) = phs_fts(2)+2*pi;
    end
    
    amp_meas = amp_fts(1)/2;
    
    %try subtracting a little below lowest known reading
    %offset = 1*10^-10;
    %amp_meas = amp_meas-offset;
    
    %look up closest DM amplitude value
    [c, index] = min(abs(amp_meas-amp_cal(:,2)));
    
    %if less than lowest calibrated value, interpolate
    if index <= 2
        amp_solved = amp_meas/amp_cal(2,2)*amp_cal(2,1);
    else
        amp_solved = amp_cal(index,1);
    end
    
    phase_solved = phs_fts(2);
    z_k = [phase_solved; amp_solved]; %measurement

    %Kalman filter
    x_k_p = x_k_m+P_k_m*H_k'*inv(H_k*P_k_m*H_k'+R_k)*(z_k-H_k*x_k_m); %update estimate with measurement
    P_k_p = inv(inv(P_k_m)+H_k'*inv(R_k)*H_k); %update covariance
    
    %control input gain (tune)
    K_p = 1;
    
    %control input calculation
    u_k = [-(x_k_p(1)-estimations(i-1,1)); -K_p*(x_k_p(2)-estimations(i-1,2))]; %New control vector is the change in the null solution for the DM shape
    
    %convert control input and actuate DM
    par_in = [x_k_p(2), ang, act, x_k_p(1)];
    pow = DE_supMinimizer_3(par_in);
    
    x_k_m = x_k_p+gamma*u_k; %new prior based on control input
    
    %heuristic nonlinear "bottoming out" adjustment
    %amp_offset = 2;
    %x_k_m(2) = x_k_m(2)+amp_offset;
    
    P_k_m = P_k_p+Q_k; %new prior covariance
    
    %set probe amp to expected amplitude, if too low, set to minimum
    %discernible amplitude
    probe_amp = x_k_p(2)/2;
    
    if probe_amp > 3
       probe_amp = 3; 
    end
    
    if probe_amp < 0.05       
        probe_amp = 0.05;  
    end
    
%     %update flat map
%     DM_Map = DE_DMMapSin(par_in(1),par_in(2),par_in(3),par_in(4));
%     flat = flat + DM_Map;
%     DM_history(:,:,i) = flat;

    %log for analysis
    %real_amp_meas(i) = amp_solved;
    measurements(i,:) = z_k;
    estimations(i,:) = x_k_p;
    inputs(:,i) = u_k;
    
end

figure(1)
plot(steps,measurements(:,1),steps,estimations(:,1))

figure(2)
plot(steps,measurements(:,2),steps,estimations(:,2))

%     [Par, Pmin] = patternsearch(@DE_supMinimizer_3, par0, [], [], [], [], ... 
%          LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));

%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
%Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

% Save all images as single cube
fitswrite(Idat, Idatnm);

%save DM history
fitswrite(DM_history, DMdatnm);

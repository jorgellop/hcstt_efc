% hcstt_snkf_simulation
%
% Speckle Nulling with a Kalman Filter. Main code extracted from
% kf_null_femto.m by Yinzi Xin
%
% Grady Morrissey - 
% Yinzi Xin - 2018
% Jorge Llop - November 29, 2018

clear all;
% close all;
addpath(genpath('utils'),genpath('export_scripts'));

label = '_Apr01';
outDir = ['output',filesep,'EFC_wFiber_LabDemonstration',label,filesep];
mkdir(outDir);

% Load calibration data
% load(['output',filesep,'sn_amp_calDec05.mat']);
% load(['output',filesep,'amp_cal_new_2p4.mat']);
load(['output',filesep,'amp_cal__Feb25.mat']);
load(['output',filesep,'calibrateDM_Aug01']); % actxcDM, angDm vs pix on camera info
% load(['Calibrations',filesep,'angle_cal']);
% load(['Calibrations',filesep,'position_cal']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normI = 1.96e-6; %Intesity of non-coronaPSF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


num_steps = 10;%30 in Yinzi's code

%generate equally spaced probes (4 different phis)
% probe_phase = [0, 3.14/2, 3.14, 3.14*1.5];
n_pairs = 6;
probe_phase = linspace(0,2*pi,n_pairs);
n_pairs = numel(probe_phase);
% probe_amp = 20;

% Take background image
take_background = false;
Ncam = 400;Nact=12;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
%     im_cam = zeros(400,400);
%     for II=1:15
%         im_camII = hcstt_TakeCamImage(true,false,tint);
%         im_cam = im_cam + im_camII/15;
%         pause(0.1)
%     end
%     backgroundCam = im_cam;
    backgroundSMF = hcstt_GetIntensityFIU(zeros(Nact,Nact),15,0);
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
%     backgroundCam = zeros(Ncam,Ncam);
    backgroundSMF = 0;
end


%% Find position of fiber
x_fib_est =2.625;
ang_fib_est = -0.105;
info.backgroundSMF = backgroundSMF;
% [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(x_fib_est,ang_fib_est,info)
actxc_fib =2.65;
ang_fib = -0.1050;

%% SNKF

num_steps = 60;%30 in Yinzi's code
steps = 1:num_steps;

probe_amp = 20;
numamptry = 6;
amptry_arr = linspace(1,15,numamptry);
%load interaction matrix
gamma = [1,0;0,1];

%covariance initialization
P_0 = eye(2); 

numw = 5;
numphn=numw;
numamn=numw;
w_arr = linspace(-5,-3,numw);
phn_arr = linspace(-5,-3,numw);
amn_arr = linspace(-5,-3,numw);
% w_arr = [-5];
% phn_arr = [-5];
% amn_arr = [-5];
for IIw=1:numw
    if IIw~=1
        x_fib_est = actxc_fib;
        ang_fib_est = ang_fib;
        [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(x_fib_est,ang_fib_est,info)
    end
    for IIphn=1:numphn
        for IIamn=1:numamn
            w=w_arr(IIw);
            phn=phn_arr(IIphn);
            amn=amn_arr(IIamn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process noise (on x)
w_k = [.1,0;0,.1]*10^w         ;
Q_k = w_k*w_k';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%measurement noise
phase_noise = 1*10^phn       ;
amplitude_noise = 0.5*10^amn       ;
R_k = [phase_noise^2,0;0,amplitude_noise^2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start looping here

%initialize values
x_k_m = [0;5]; %prior estimates %JLlop: from meas
x_k_p = [0;0];

P_k_m = P_0; %prior covariance
u_k = [0;0]; %control input
estimation = [0;0];
%initialize logs
measurements = zeros(num_steps,2);
estimations = zeros(num_steps,2);
inputs = zeros(2,num_steps);
intensity_arr = zeros(num_steps+1,1);
intensityKF_arr = zeros(num_steps+1,1);
flat_new = 0;

intensity_arr(1) = hcstt_GetIntensityFIU(zeros(12,12),5,backgroundSMF);
intensityKF_arr(1) = hcstt_GetIntensityFIU(zeros(12,12),5,backgroundSMF);

for II=steps %KF loop (II==index)
    disp(['Step ',num2str(II),'/',num2str(num_steps)])
    int_arr = zeros(n_pairs,1);
    
    
    H_k = eye(2); %observation matrix
    %% Speckle Nulling
    for n=1:1:n_pairs %SN loop
        % par_probe = [probe_amp, ang_fib, actxc_fib, probe_phase(n)];
        dm_actuators_mat0 = hcstt_DMMapSin(probe_amp, ang_fib, actxc_fib, probe_phase(n));
        intensity_raw = hcstt_GetIntensityFIU(dm_actuators_mat0(:),5,backgroundSMF);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm

        %simulate gaussian detector noise
%         noise_sigma = 0.5*10^-15;
        %draw from gaussian
%         noise_amp = normrnd(0,noise_sigma);
        int_arr(n) = intensity_raw;%+noise_amp;
    end
    figure(201)
    plot(probe_phase,int_arr)
    title('Probing')
    xlabel('phase')
    ylabel('Intensity')
    pause(.5)
    
    % Finde phase
    yu = max(int_arr);
    yl = min(int_arr);
    yr = (yu-yl);                               % Range of ‘y’
    yz = int_arr-yu+(yr/2);
    ym = mean(int_arr);                               % Estimate offset
    fit = @(b,x)  b(1).*(sin(x+b(2))) + b(3);    % Function to fit
    fcn = @(b) sum((fit(b,probe_phase') - int_arr).^2);                              % Least-Squares cost function
    s = fminsearch(fcn, [yr;  -1;  ym]);                       % Minimise Least-Squares
    probe_phase_fine = linspace(0,2*pi,100);
    sinefit = fit(s,probe_phase_fine);
    [mi,ind_mi] = min(sinefit);
    
    phase_solved = probe_phase_fine(ind_mi); %SN solution to phase
    
    % Find Amplitude
    int_arr = zeros(numamptry,1);
    for JJ=1:numamptry
        dm_actuators_mat0 = hcstt_DMMapSin(amptry_arr(JJ), ang_fib, actxc_fib, phase_solved);
        intensity_raw = hcstt_GetIntensityFIU(dm_actuators_mat0(:),5,backgroundSMF);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
        int_arr(JJ) = intensity_raw;
    end
    figure(202)
    plot(amptry_arr,int_arr)
    title('Amplitude try')
    xlabel('Amp')
    ylabel('Int')
    yu = max(int_arr);
    [yl,ind_mi] = min(int_arr);
    yr = (yu-yl);                               % Range of ‘y’
    yz = int_arr-yu+(yr/2);
%     ym = mean(int_arr);                               % Estimate offset
    fit = @(b,x)  b(1).*((x-b(3)).^2) + b(2);    % Function to fit
    fcn = @(b) sum((fit(b,amptry_arr') - int_arr).^2);                              % Least-Squares cost function
    s = fminsearch(fcn, [yr;yl ; amptry_arr(ind_mi)]);                       % Minimise Least-Squares
    amptry_fine = linspace(1,15,100);
    quadfit = fit(s,amptry_fine);
    [mi,ind_mi] = min(quadfit);

    amp_solved = amptry_fine(ind_mi); %SN solution to amplitude
    
    z_k = [phase_solved; amp_solved]; %measurement

    %% Kalman filter
    x_k_p = x_k_m+P_k_m*H_k'*inv(H_k*P_k_m*H_k'+R_k)*(z_k-H_k*x_k_m); %update estimate with measurement
    P_k_p = inv(inv(P_k_m)+H_k'*inv(R_k)*H_k); %update covariance         eye(2);%
    
    %control input gain (tune)
    K_p = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %control input calculation
    u_k =[0;0];% [(x_k_p(1)-estimation(1)); K_p*(x_k_p(2)-estimation(2))]; %New control vector is the change in the null solution for the DM shape
%      u_k = [-x_k_p(1); -K_p*x_k_p(2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %convert control input and actuate DM
%     dm_actuators_mat0 = hcstt_DMMapSin(u_k(2), ang_fib, actxc_fib, u_k(1));
%     pow = hcstt_GetIntensityFIU(dm_actuators_mat0(:),5,backgroundSMF);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm

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
%     flat_new = hcstt_DMMapSin(x_k_p(2), ang_fib, actxc_fib, x_k_p(1));
%     hcstt_AddToFlatForDM(flat_new);
    
    %log for analysis
    %real_amp_meas(i) = amp_solved;
    measurements(II,:) = z_k;
    estimations(II,:) = x_k_p;
    estimation = x_k_p;
    inputs(:,II) = u_k;
    
    dm_actuators_mat0 = hcstt_DMMapSin(z_k(2), ang_fib, actxc_fib, z_k(1));
    intensity_arr(II+1) = hcstt_GetIntensityFIU(dm_actuators_mat0,5,backgroundSMF);
    dm_actuators_mat0 = hcstt_DMMapSin(x_k_p(2), ang_fib, actxc_fib, x_k_p(1));
    intensityKF_arr(II+1) = hcstt_GetIntensityFIU(dm_actuators_mat0,5,backgroundSMF);
    
    figure(300)
    plot(0:II,log10(intensity_arr(1:II+1)/normI),0:II,log10(intensityKF_arr(1:II+1)/normI))
    title('SNKF')
    ylabel('Raw Contrast')
    xlabel('Step')
    legend('From measurement','With KF')
    drawnow;
end

figure(300)
plot(0:num_steps,log10(intensity_arr/normI),0:num_steps,log10(intensityKF_arr/normI))
title('SNKF')
ylabel('Raw Contrast')
xlabel('Step')
legend('From measurement','With KF')
export_fig(['output\snkf_Mar01\SNKF_RCvIt_w',num2str(w),'_phn',num2str(phn),'_amn',num2str(amn)]);
figure(304)
plot(0:num_steps,log10(intensity_arr/normI))
title('SN from measurement')
ylabel('Raw Contrast (Log scale)')
xlabel('Step')
figure(301)
plot(steps,measurements(:,1),steps,estimations(:,1))
legend('measurement','estimation')
title('Phase')
figure(302)
plot(steps,measurements(:,2),steps,estimations(:,2))
title('Amplitude')
legend('measurement','estimation')

save([outDir,'snkf_w',num2str(w),'_phn',num2str(phn),'_amn',num2str(amn)],'measurements','intensity_arr','intensityKF_arr','estimations','normI')
        end
    end
end


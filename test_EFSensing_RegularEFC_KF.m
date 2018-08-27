% test_EFSensing_RegularEFC_KF
%
% Test EF sensing for Regular EFC pixel-based on the bench
% 
%
% v2 Jorge Llop - Feb 15, 2018
% modified Milan Roberson - Aug 15, 2018
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false);

tint = 10; % exposure time
info.tint = tint;
sidepix = 20; %side pixels to plot when cropping the images

N = 2^10; % Size on array for model propagation
apRad = 142; % Radius of aperture for model propagation
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD = N/apRad/2; % lambda/D (samples) 
Ncam = 400;
lambda0 = 650e-9;
numOfWavelengths = 1; % monochromatic, use 1
percentBW = 10; % percent bandwidth=(Delta lambda)/lambda*100
BW = percentBW/100; % bandwidth 
if(numOfWavelengths > 1)
    lam_fracs = linspace(1-BW/2,1+BW/2,numOfWavelengths);
else
    lam_fracs = 1;
end
lam_arr = lambda0*lam_fracs;
fiberDiam = 2*0.71; % Fiber diam. (lambda_0/D)
% These are parameters used in the propagation routines
% use_fiber = true;
% normal_EFC = false;

use_fiber = false;
normal_EFC = true;

num_DM_shapes = 5;

label = '_testEFSensing_RegularEFC_Aug03';
outDir = ['output',filesep,'EFC_wFiber_LabDemonstration',label,filesep];
mkdir(outDir);

% Variables that go into the model propagation
info.fiberDiam = fiberDiam;
info.useGPU = false;
info.apRad = apRad;
info.lambdaOverD = lambdaOverD;
info.RHO = RHO;
info.THETA = THETA;
info.N = N;
info.lambda0 = lambda0;
info.lam_arr = lam_arr;
info.numOfWavelengths = numOfWavelengths;
info.useApodizer = false;
info.useGPU = false; 
info.FPM = exp(1i*4*THETA);
info.outDir = outDir;
info.xvals = xvals;
info.yvals = yvals;
info.use_fiber = use_fiber;
info.normal_EFC = normal_EFC;
info.LPM = exp(-(RHO/(0.85*apRad)).^1000);
info.num_DM_shapes = num_DM_shapes;

% Define the aperture for the model propagation
wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;

% Define position of DH
x_fib_pix = +15;
% x_fib=3.5;
% y_fib=0;
% info.x_fib = x_fib;
% info.y_fib = y_fib;
info.x_fib_pix = x_fib_pix;
info.y_fib_pix = 0;

% Take background image
take_background = false;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    background = zeros(Ncam,Ncam);
end
info.background = background;

% Load new flat into the DM from Image Sharpening
hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');

% Find Center of camera image
im_cam = zeros(400,400);
tint_findCenter = 0.3;
for II=1:10
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
    im_cam = im_cam + im_camII/55;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
info.x_cent_cam = x_cent_cam;
info.y_cent_cam = y_cent_cam;
if max(im_camII(:))>240
    disp('Find Center image saturated')
    return
end

% Model of the fiber mode shape
% [THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD ,Y);
% fiberDiam_pix = (fiberDiam*lambdaOverD);
% fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
%         exp(-(RHO_fib/(fiberDiam_pix/2)).^2);
% info.fibermode0 = fibermode0;


Nact = 12;  % Number of DM actuators, Nact^2

% Load the DM influence functions
% [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuators(N,apRad);
[posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,1);
info.posDM_x = posDM_x;
info.posDM_y = posDM_y;
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac ); % Influence function. Need of a model for actual DM
info.ac_spac = ac_spac;
info.infl = infl;

% Initialize DMs, etc. 
surf_DM10 = zeros(N); % Intialize the DM surface to flat
DM1_strokes = zeros(N); % Intialize the DM strokes to flat
us = zeros(1,Nact^2); % Initialize the fractional stroke changes to zero
us_total = zeros(Nact^2,1); % Initialize the fractional stroke changes to zero
poke_amp = 1e-9; % Initialize the poke amplitude
    
% maxits = 33; % maximum number of EFC iterations allowed
Gcount = 0; % counter for number of times G matrix was used
Gcountmin = 10; % Minimum number of times to use a G matrix
curr_coupl_SMF = 1;
recalc_G = true; % Initialize the flag to re-calculate the G matrix
% regvals = logspace(-6,-1,6); % Range of regularization parameters to test
regval = nan; % Initial regularization value
coupl_SMF_in_DH = []; % Array to keep track of dark hole irradiance 

% Compute Q, the indeces of the pixels in the DH
[Xcam,Ycam] = meshgrid(-y_cent_cam+1:Ncam-y_cent_cam,-x_cent_cam+1:Ncam-x_cent_cam); 
q_pix = 1;
info.q_pix = q_pix;
Q = zeros(Ncam,Ncam);
Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
Q = and(Q, Ycam >= -(q_pix) );
Q = and(Q, Ycam <= (q_pix) );
num_Q = numel(find(Q));
info.num_Q = num_Q;
info.Q = Q;

% Plot Q
figure(103)
imagesc(Q(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
title('Q')
axis image

% Compute Q4G. Q4G is Q, ie the DH indeces, in the model array
[Xcam4G,Ycam4G] = meshgrid(-Ncam/2-1:Ncam/2-2); 
[THETAcam4G, RHOcam4G] = cart2pol(Xcam4G,Ycam4G);
Q4G = zeros(Ncam,Ncam);
Q4G = and(Xcam4G >= (x_fib_pix-q_pix ), Xcam4G <=  (x_fib_pix+q_pix ));
Q4G = and(Q4G, Ycam4G >= -(q_pix) );
Q4G = and(Q4G, Ycam4G <= (q_pix) ); % 60deg keystone about x-axis
info.Q4G = Q4G;

% Propagate to image plane
wf2_current0 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, zeros(N,N), true, info);

% Plot model
% figure(105)
% immod = abs(wf2_current).^2;
% immod_crop = immod(N/2-Ncam/2:N/2+Ncam/2-1,N/2-Ncam/2:N/2+Ncam/2-1);
% immod_crop(Q4G) = max(immod(:));
% imagesc(immod_crop(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
% axis image

% Measure the intensity in Q with the flat DM
hcstt_UpdateMultiDM(zeros(Nact,Nact))
im_cam0 = zeros(Ncam,Ncam);
for II=1:10
    im_camII = hcstt_TakeCamImage(true,false,tint)-background;
    im_cam0 = im_cam0 + im_camII/10;
end
int_cam = im_cam0(Q)';

% Image the coronPSF with the flat DM
figure(106)
imagesc(im_cam0(x_cent_cam-sidepix:x_cent_cam+sidepix,y_cent_cam-sidepix:y_cent_cam+sidepix))
axis image
colorbar
title('CAM - flat DM')

one_it = false; % are we trying just once, or doing a parameter search?
if(one_it)
    load('BenchModelNormalization_0803') % load  the normalization parameters
    normPower = normPower_normalization*tint/tint_normalization;%0.00055;%
    info.normPower = normPower;
    wf2_current = wf2_current0 * sqrt(normPower);

    info.p2v_dm_sensing = 5; % Size of poke for the sensing

    Eab =  EFSensing_RegularEFC_labTest(wf2_current,zeros(Nact^2,1),info);
    x_hatRe = Eab(1,:);
    x_hatIm = Eab(2,:);
    int_est = abs(x_hatRe).^2+abs(x_hatIm).^2;
    med_cam = median(int_cam);
    mean_cam = mean(int_cam);
    med_est = median(int_est);
    mean_est = mean(int_est);
    int_cam'
    int_est'
%     int_cam = int_cam'
%     int = int'
%     int_cam_norm = int_cam/max(int_cam);
%     int_norm = int/max(int);
    figure(301)
    plot(1:nnz(Q),int_cam)
    xlabel('Pixel label from the DH')
    ylabel('Intensity')
    title('Measured vs Estimated')
    hold on
    plot(1:nnz(Q),int_est)
    legend('measured','estimated')
    hold off
    figure(302)
    plot(1:nnz(Q),int_cam/max(int_cam))
    xlabel('Pixel label from the DH')
    ylabel('Intensity')
    title('Measured vs Estimated - Normalized')
    hold on
    plot(1:nnz(Q),int_est/max(int_est))
    legend('measured','estimated')
    hold off
else
    numtry = 5;
    x_hatRe = zeros(numtry,num_Q);
    x_hatIm = zeros(numtry,num_Q);
    
    P_mat = eye(2*num_Q)*1e-5; % State estimate error covariance
    Q_mat = eye(2*num_Q)*1e-5; % Process noise covariance (control noise)
    R_mat = eye(num_DM_shapes*num_Q)*1e-5; % Observation noise covariance
    
    load('KF_cov_estimate');
    Q_mat = Qnew;
    R_mat = Rnew;
    P_mat = Qnew;
    
    
    Eab = EFSensing_RegularEFC_labTest(wf2_current, zeros(Nact^2,1), info); % Initial guess
    
    for k = 1:numtry

        fprintf('Iteration: %d ',k);

        info.p2v_dm_sensing = 5;
        normPower = normPower_normalization*tint/tint_normalization;
        info.normPower = normPower;

        wf2_current = wf2_current0 * sqrt(normPower);
        
        [Eab, P_mat] =  EFSensing_RegularEFC_labTest_KF(wf2_current,zeros(Nact^2,1),info, Eab, P_mat, Q_mat, R_mat);
        x_hatRe(k,:) = Eab(1,:);
        x_hatIm(k,:) = Eab(2,:);
        
        figure(1000);
        imagesc(P_mat);
        title('Error Covariance Estimate');
        
        int = abs(x_hatRe(k,:)).^2+abs(x_hatIm(k,:)).^2;
        med_cam = median(int_cam);
        mean_cam = mean(int_cam);
        med_est(k) = median(int);
        mean_est(k) = mean(int);
        im_cam_q = int_cam'
        im_q = int'
    %     int_cam = int_cam'
    %     int = int'
    %     int_cam_norm = int_cam/max(int_cam);
    %     int_norm = int/max(int);
        figure(200)
        plot(1:nnz(Q),int_cam)
        xlabel('Pixel label from the DH')
        ylabel('Intensity')
        hold on
        plot(1:nnz(Q),int)
        hold off
%         export_fig([info.outDir,'Pix_int',num2str(normPower)],'-r300');
        err_arr =  (int-int_cam);
    %     err_re_arr =  (int-int_cam)./int_cam*100;
        std_err_re_arr(k) = std(abs(err_arr));
        mean_err_re_arr(k) = mean(int)-mean(int_cam);
        %
    end
%     figure(200)
%     % plot(normPower_arr(1:k),me_err_re_arr(1:k))
%     plot(normPower_arr(1:k),std_err_re_arr(1:k))
%     figure(300)
%     % plot(normPower_arr(1:k),me_err_re_arr(1:k))
%     plot(normPower_arr(1:k),mean_err_re_arr(1:k))
end
% close all;
% figure(100)
% plot(normPower_arr, me_err_re_arr)
% xlabel('Intensity ')
% ylabel('pixel label from the DH')
% export_fig([info.outDir,'BenchErrorMeanPercentage_vs_normPower'],'-r300');

% close all;
% figure(100)
% plot(p2v_arr, me_err_re_arr)
% xlabel('Amplitude sinusoid [nm]')
% ylabel('Error mean percentge intensity')
% export_fig([info.outDir,'BenchErrorMeanPercentage_vs_amplitude'],'-r300');

% figure(101)
% int_hat = x_hatRe.^2+x_hatIm.^2;
% int_real = x_trueRe.^2+x_trueIm.^2;
% plot(p2v_arr,int_hat)
% xlabel('Amplitude sinusoid [nm]')
% ylabel('Int')
% hold on
% plot(p2v_arr,int_real)
% hold off
% export_fig([info.outDir,'ModelInt_vs_amplitude'],'-r300');

hcstt_DisconnectDevices()

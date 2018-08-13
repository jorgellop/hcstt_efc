% Test EF sensing through fiber vs amplitude of the sinusoid probe on the DM.
% Bench based
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize();

N = 2^10;
apRad = 2^6;
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD = N/apRad/2; % lambda/D (samples) 

lambda0 = 800e-9;
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
use_fiber = true;
normal_EFC = false;
    
label = '_testEFSensing_Jan10';
outDir = ['output',filesep,'EFC_wFiber_LabDemonstration',label,filesep];
mkdir(outDir);

% Position of the fiber on the image plane. (lambda_0/D)
% x_fib = 3;
% y_fib = 0;

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
info.FPM = exp(1i*8*THETA);
info.outDir = outDir;
info.xvals = xvals;
info.yvals = yvals;
info.use_fiber = use_fiber;
info.normal_EFC = normal_EFC;

info.LPM = ones(N,N);
wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;
% wf2 = wf2*sqrt(normPower);

% Total Power to normilize model Gu and intensity from the fiber. Need of
% normalization factor
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'');
[normPower,totalPower,x_fib,y_fib,x_cent_cam,y_cent_cam,x_fib_pix] = hcstt_NormalizationFactor(wfin_noerrors,info,'',8.8e-7,218,204,206,204);

info.x_fib = x_fib;
info.y_fib = y_fib;

info.normPower = normPower;

% Model of the fiber mode shape
[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD ,Y);
fiberDiam_pix = (fiberDiam*lambdaOverD);
fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
        exp(-(RHO_fib/(fiberDiam_pix/2)).^2);
info.fibermode0 = fibermode0;

info.LPM = exp(-(RHO/(0.9*apRad)).^1000);

Nact = 12;  % Number of DM actuators, Nact^2

% Load the DM influence functions
ac_spac = round(2*apRad/Nact);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac ); % Influence function. Need of a model for actual DM
info.ac_spac = ac_spac;
info.infl = infl;

% Initialize DMs, etc. 
surf_DM10 = zeros(N); % Intialize the DM surface to flat
DM1_strokes = zeros(N); % Intialize the DM strokes to flat
DM1_strokesKK = zeros(N,N);
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

[posDM_x,posDM_y] = hcstt_PositionDMActuators(N,apRad);
info.posDM_x = posDM_x;
info.posDM_y = posDM_y;

%Take camera image with flat DM
hcstt_UpdateMultiDM(zeros(Nact^2,1))
im_cam = hcstt_TakeCamImage(false,false);
sz = size(im_cam);
hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_flatDM'], sz );

maxits = 30;
x_hatRe = zeros(1,maxits);
x_hatIm = zeros(1,maxits);
int_bench = zeros(1,maxits);
int_mod = zeros(1,maxits);
p2v_arr = linspace(50,300,maxits);



% Update actuator height with the LMS solution, us
us_total = us_total + us';

% Build DM surface from stroke amplitudes 
count = 0; 
for ix = 1:Nact
    xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
    for iy = 1:Nact
        count = count + 1;
        ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        DM1_strokes(xpos,ypos) = us(count)*poke_amp + DM1_strokes(xpos,ypos);
    end
end
surf_DM10 = conv2(DM1_strokes,infl,'same'); 
wf2_current = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM10, true, info);
wf2_current = wf2_current * normPower;

[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD, Y);

fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
    exp(-(RHO_fib/(fiberDiam_pix/2)).^2);

% Run EFC iterations 
for k = 1:maxits
            
    fprintf('Iteration: %d ',k);
    
    % Perform WF sensing
    poke_amp = p2v_arr(k)*1e-9;
    
    cosfct = cos(2*pi*[1:Nact]/(Nact) * x_fib + 0) ;
    a = ones(Nact,Nact);
    di = diag(cosfct);
    us = a * di; 
    
    dm_actuators_mat = us' * poke_amp;
    
    %Add the actuators heights that were already set on the DM:
    dm_actuators_mat = dm_actuators_mat(:) + us_total;
    
    figure(200)
    imagesc(dm_actuators_mat)
    colorbar
    axis image
    export_fig([info.outDir,'DM_shape_',int2str(k)],'-r300');
    close;
    
    count = 0; 
    DM1_strokes = zeros(N,N);
    for ix = 1:Nact
        for iy = 1:Nact
            count = count + 1;
            xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
            DM1_strokesKK(xpos,ypos) = dm_actuators_mat(count) + DM1_strokes(xpos,ypos);
        end
    end
    surf_DM1 = conv2(DM1_strokesKK,infl,'same');
    wf2_prob_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM1, true, info);
    wf2_prob_noerrors = wf2_prob_noerrors * sqrt(normPower);

    hcstt_UpdateMultiDM(+dm_actuators_mat/1e-9)
    int_bench(k) = hcstt_GetIntensityFIU(+dm_actuators_mat/1e-9);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm

%     if(KK ~= 1) openhimg = true; else openhimg = false;  end
    im_cam = hcstt_TakeCamImage(false,true);
    sz = size(im_cam);
    hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_k',int2str(k)], sz );
    hcstt_test_plotModelImage(wf2_prob_noerrors, [outDir,'ModelImage_k',int2str(k)], info )

    Eab_mod = sum(sum(wf2_prob_noerrors.*fibermode0));
%     Eab =  EFSensing_EFCwFiber_labTest(wf2_current,us_total*poke_amp,info);
    int_mod(k) = imag(Eab_mod)^2+real(Eab_mod)^2;

    % Perform WF sensing
    info.p2v_dm_sensing = p2v_arr(k);
    Eab =  EFSensing_EFCwFiber_labTest(wf2_current,us_total*poke_amp,info);
    x_hatRe(k) = Eab(1);
    x_hatIm(k) = Eab(2);

    
    %
end
figure(100)
plot(p2v_arr,x_hatRe)
xlabel('Amplitude sinusoid [nm]')
ylabel('Re/Im')
hold on
plot(p2v_arr,x_hatIm)
hold off
export_fig([info.outDir,'SensedReIm_vs_amplitude'],'-r300');

figure(101)
int_arr = x_hatRe.^2+x_hatIm.^2;
xlabel('Amplitude sinusoid [nm]')
ylabel('Int')
plot(p2v_arr,int_arr)
export_fig([info.outDir,'SensedInt_vs_amplitude'],'-r300');

figure(102)
int_arr = x_hatRe.^2+x_hatIm.^2;
xlabel('Amplitude sinusoid [nm]')
ylabel('Int')
plot(p2v_arr,int_mod)
export_fig([info.outDir,'IntTotalModel_vs_amplitude'],'-r300');
figure(102)
int_arr = x_hatRe.^2+x_hatIm.^2;
xlabel('Amplitude sinusoid [nm]')
ylabel('Int')
plot(p2v_arr,int_bench)
export_fig([info.outDir,'IntTotalBench_vs_amplitude'],'-r300');

hcstt_DisconnectDevices()
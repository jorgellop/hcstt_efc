% Test EF sensing for Regular EFC pixel-based
% Model based
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize();

N = 2^10;
% apRad = 2^6;
apRad = 114;
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
    
label = '_testEFSensing_RegularEFC_Jan19';
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
% rng(3)
% wfin_werrors = complex(ones(N, N), randn(N,N)*1) ;
% wfin_werrors(RHO > apRad) = 0;
error_map = fitsread('surfErrorMap_OX5.fits')*0.1;
wfin_werrors = fftshift(wfin_noerrors).*fftshift(exp(1i*4*pi*error_map'/lambda0));
wfin_werrors = fftshift(wfin_werrors);

% wf2 = wf2*sqrt(normPower);

% Total Power to normilize model Gu and intensity from the fiber. Need of
% normalization factor
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'');
% [normPower,totalPower,x_fib,y_fib,x_cent_cam,y_cent_cam,x_fib_pix] = ...
%     hcstt_NormalizationFactor_RegularEFC(wfin_noerrors,info,'',7.3e-7,217,202,207,202);
% wfnorm_current_werrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_werrors, zeros(N,N), false, info);
x_cent_cam = 208;
y_cent_cam = 205;
x_fib_pix = +13;
normPower = 1;
totalPower = 1;
x_fib=3.5;
y_fib=0;

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

Ncam = 400;
[Xcam,Ycam] = meshgrid(-x_cent_cam:Ncam-x_cent_cam-1,-y_cent_cam:Ncam-y_cent_cam-1); 
q_pix = 3;
info.q_pix = q_pix;
Q = zeros(Ncam,Ncam);
Q = and(Xcam >= (x_fib_pix-q_pix+1 ), Xcam <=  (x_fib_pix+q_pix ));
Q = and(Q, Ycam >= -(q_pix)+1 );
Q = and(Q, Ycam <= (q_pix) );
num_Q = numel(find(Q));
info.num_Q = num_Q;
info.Q = Q;
%Take camera image with flat DM
% hcstt_UpdateMultiDM(zeros(Nact^2,1))
% im_cam = hcstt_TakeCamImage(false,false);
% sz = size(im_cam);
% hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_flatDM'], sz );




maxits = 1;
x_hatRe = zeros(maxits,num_Q);
x_hatIm = zeros(maxits,num_Q);
p2v_arr = linspace(1,50,maxits);
% p2v_arr = 5;

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
wf2_current = wf2_current * sqrt(normPower);
% wf2_current_werrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_werrors, surf_DM10, true, info);
% wf2_current_werrors = wf2_current_werrors * sqrt(normPower);
% hcstt_test_plotModelImage(wf2_current_werrors, [outDir,'ModelImage'], info )

[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD, Y);

fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
    exp(-(RHO_fib/(fiberDiam_pix/2)).^2);

hcstt_UpdateMultiDM(zeros(Nact,Nact))
int_cam = hcstt_TakeCamImage(false,true);
int_cam = int_cam(Q)';

normPower_arr = 600;%linspace(100,2000,40);

% Run EFC iterations 
for k = 1:numel(normPower_arr)
    
    fprintf('Iteration: %d ',k);
    
    info.p2v_dm_sensing = 20;
    info.normPower = normPower_arr(k);
    Eab =  EFSensing_RegularEFC_labTest(wf2_current,us_total*poke_amp,info);
    x_hatRe(k,:) = Eab(1,:);
    x_hatIm(k,:) = Eab(2,:);
    int = abs(x_hatRe(k,:)).^2+abs(x_hatIm(k,:)).^2;
    err_re_arr =  (int-int_cam)./int_cam*100
    me_err_re_arr(k) = median(abs(err_re_arr))
    %
end
close all;
figure(100)
plot(normPower_arr, me_err_re_arr)
xlabel('Normalization factor ')
ylabel('Error mean percentge intensity')
export_fig([info.outDir,'BenchErrorMeanPercentage_vs_normPower'],'-r300');

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

% test_checkModulation.m
%
% Check if the speckles modulate if we put the EFC solution on the DM
%
% Jorge Llop - Aug 21, 2018

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(true);

N = 1024;

apRad = 142;
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD =  N/apRad/2; 

p2v_dm = 4;

rng(3);
% error_map = fitsread('surfErrorMap_OX5.fits'); %nm

wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;

lambda0 = 635e-9;
numOfWavelengths = 1; % monochromatic, use 1
percentBW = 10; % percent bandwidth=(Delta lambda)/lambda*100
BW = percentBW/100; % bandwidth 
if(numOfWavelengths > 1)
    lam_fracs = linspace(1-BW/2,1+BW/2,numOfWavelengths);
else
    lam_fracs = 1;
end
lam_arr = lambda0*lam_fracs;

% These are parameters used in the propagation routines  
info.useGPU = false;
info.RHO = RHO;
info.THETA = THETA;
info.N = N;
info.lambda0 = lambda0;
info.lam_arr = lam_arr;
info.numOfWavelengths = numOfWavelengths;
info.useApodizer = false;
info.useGPU = false; 
info.FPM = exp(1i*4*THETA);
info.xvals = xvals;
info.yvals = yvals;
% info.use_fiber = use_fiber;
% info.normal_EFC = normal_EFC;
% info.tint = tint;
info.LPM = ones(N,N);

% Model WF with flat DM, WF0, we need this to compute Gu, since Gu:
%Gu = WF_DM - WF0

Nact = 12;    

x_fib = 2.5;
load('output\calibrateDM_Aug01'); % actxcDM, angDm vs pix on camera info
% [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(x_fib,0);
actxc_fib = 2.525;
ang_fib = 0.125;
info.actxc_fib = actxc_fib;
info.ang_fib = ang_fib;
%Calculate position in pixels
if actxc_fib>max(actxcDM_arr) 
    r_fib_pix = interp1(actxcDM_arr,distPix_meas,actxc_fib,'linear','extrap');
elseif actxc_fib<min(actxcDM_arr)
    disp('Fiber too far away from star!')
    return;
else
    r_fib_pix = interp1(actxcDM_arr,distPix_meas,actxc_fib);
end
x_fib_pix = r_fib_pix*cos(ang_fib);
y_fib_pix = r_fib_pix*sin(ang_fib);

fiberDiam = 2*0.71; % Fiber diam. (lambda_0/D)
fiberDiam_pix = (fiberDiam*lambdaOverD);

[THETA_fib,RHO_fib] = cart2pol(X - x_fib_pix ,Y - y_fib_pix);
fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
        exp(-(RHO_fib/(fiberDiam_pix/2)).^2);
info.fibermode0 = fibermode0;

info.fiberDiam = fiberDiam;
info.apRad = apRad;
info.lambdaOverD = lambdaOverD;

info.LPM = exp(-(RHO/(0.85*apRad)).^1000);

[posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,2);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );

totalPower = 1.6e-6;
normPower = totalPower/3.8185e10;
peakInt = totalPower;


num_DM_shapes = 17;
ph_arr = linspace(-pi/2, 3/2*pi, num_DM_shapes);
H_mat = zeros(num_DM_shapes,2);
DeltaI_arr = zeros(num_DM_shapes,1);
% ww = x_fib;
apRad2 = 12;
% poke_amp = p2v_dm*1e-9;
DM1_strokesKK = zeros(N,N);

load('C:\Users\COO_user\Documents\GitHub\hcstt_efc\output\EFC_wFiber_LabDemonstration_Fiber_Aug22v2\data_intvsit_dmshapes__Fiber_Aug22v2_DMconfig5_apRad142')
% % load('test_us_total_Aug22');
hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');
int_flatDM = hcstt_GetIntensityFIU(zeros(Nact,Nact),20);
int_EFCsol = hcstt_GetIntensityFIU(us_total,20);
disp(['Suppression: ',num2str(int_flatDM/int_EFCsol)])

num_poke = 5;
poke_arr = linspace(20,30,num_poke);
for II = 1:num_poke
    poke_amp = poke_arr(II)*1e-9;
    for KK = 1:num_DM_shapes
    %     ph = ph_arr(KK);
    %     cosfct = cos(2*pi*[1:apRad2]/(apRad2) * ww + ph) ;
    %     a = ones(apRad2,apRad2);
    %     di = diag(cosfct);
    %     us = a * di; 
        dm_probcosfct = hcstt_DMMapSin(1, ang_fib, actxc_fib, ph_arr(KK));
        dm_actuators_mat0 = dm_probcosfct * poke_amp;

        %Add the actuators heights that were already set on the DM:
        dm_actuators_mat = dm_actuators_mat0(:) ;%+ us_total;

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

        %Gu is the effect of the DM on the image plane
    %     Gu = wf2_prob_noerrors-wf0;
    %     Gu_re = real(Gu);
    %     Gu_im = imag(Gu);

        % Measure the intensity out of the fiber for the positive probe

        int_plus = hcstt_GetIntensityFIU((+dm_actuators_mat0(:) )/1e-9,5);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
        int_plus2 = hcstt_GetIntensityFIU((+dm_actuators_mat0(:) )/1e-9+us_total,5);  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
        %     
        int_arr(KK) = int_plus;
        int_arr2(KK) = int_plus2;
    end
    figure(301+II)
    plot(1:num_DM_shapes,int_arr)
    title(['Flat DM - poke amp: ',num2str(poke_amp)])
    figure(401+II)
    plot(1:num_DM_shapes,int_arr2)
    title(['EFC DM Solution applied - poke amp: ',num2str(poke_amp)])
end
hcstt_DisconnectDevices();

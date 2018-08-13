% hcstt_SuperKCalibration
% 
% Calibrate SuperK source with camera by taking images for different power
% settings of light source
%
% Jorge Llop - Mar 5, 2018

clear all;
close all;

addpath(genpath('utils'));

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

apRad = 116;
tint_superKNorm = 0.05;

label = '_0803';
outDir = ['output',filesep,'test_SuperKNormalization',label,filesep];
mkdir(outDir);

Nact = 12;
sidepix = 20;
sidepix_pow = 2;

hcstt_Initialize(false)

Ncam = 400;
% Take background image
take_background = true;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint_superKNorm);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    info.background = background;
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    background = zeros(Ncam,Ncam);
end

powerSource_arr = 16:1:18;
for KK=1:numel(powerSource_arr) 
    disp(['Set power to: ', num2str(powerSource_arr(KK))])
    prompt = 'Continue? ';
    x = input( prompt );

    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint_superKNorm)-background;
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    [ma,ind_ma] = max(im_cam(:));
    [ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);
    % ind_ma_I = 200;
    % ind_ma_J =200;
    im_cam_crop = im_cam(ind_ma_I-sidepix_pow:ind_ma_I+sidepix_pow,ind_ma_J-sidepix_pow:ind_ma_J+sidepix_pow);
    powerCam_arr(KK) = sum(im_cam_crop(:));
end
% im_cam_crop = im_cam_crop/max(im_cam_crop(:));
figure(100)
plot(powerSource_arr,powerCam_arr)
title('Calibration SuperK with Camera')
xlabel('SuperK power [%]') % x-axis label
ylabel('Peak power at camera [counts]') % y-axis label



hcstt_DisconnectDevices();
% normPower_normalization = normPower;
% peakInt_normalization = ma;
save('utils\SuperKCamCalibration_0803.mat','powerSource_arr','powerCam_arr','tint_superKNorm')

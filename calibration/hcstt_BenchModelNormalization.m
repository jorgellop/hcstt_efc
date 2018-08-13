% hcstt_BenchModelNormalization
% 
% Find Normalization factors for bench and model
%
%
% Jorge Llop - Mar 4, 2018

clear all;
close all;

addpath(genpath('utils'));

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

apRad = 68;
tint_normalization = 0.05;
powerSetting = 16;

label = '_0803';
outDir = ['output',filesep,'test_Normalization',label,filesep];
mkdir(outDir);

Nact = 12;
sidepix = 30;
sidepix_pow = 15;

hcstt_Initialize(false)

hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');
hcstt_UpdateMultiDM(zeros(12));

Ncam = 400;
% Take background image
take_background = true;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint_normalization);
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

% Find center
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_normalization)-background;
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
p=FastPeakFind(im_cam, 3 , 4 , 2, 2);
ind_ma_I = p(1);
ind_ma_J = p(2);
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I2,ind_ma_J2] = ind2sub(size(im_cam),ind_ma);
% ind_ma_I = 200;
% ind_ma_J =200;
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
% im_cam_crop = im_cam_crop/max(im_cam_crop(:));

figure(100)
imagesc(im_cam_crop)
axis image

%% Normalization factor
DM_Command = zeros(12,12);        
info.Nact = Nact;
FPM = false;

info.normalize = false;

info.posDM_x = 0;
info.posDM_y = 0;
numtry = 1; 
info.apRad = apRad;

im_mod = hcstt_TakeModelImage(DM_Command(:),FPM,info);
im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);

total_power_mod = sum(sum(im_mod(N/2-sidepix_pow+1:N/2+sidepix_pow+1,N/2-sidepix_pow+1:N/2+sidepix_pow+1)));
total_power_cam = sum(sum(im_cam(ind_ma_I-sidepix_pow:ind_ma_I+sidepix_pow,ind_ma_J-sidepix_pow:ind_ma_J+sidepix_pow)));

normPower0 = total_power_cam/total_power_mod;
im_mod = im_mod*normPower0;
im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);

im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
figure(200)
imagesc(im_mod_crop);
axis image

%% Calibrate poke intensity


%% Disconnect Devices
hcstt_DisconnectDevices();

%% Account for the power source setting
load('SuperKCamCalibration_0803')

powerSetting_experiment = 50;
p1 = interp1(powerSource_arr,powerCam_arr,powerSetting);
if powerSetting_experiment>max(powerSource_arr)
    p2 = interp1(powerSource_arr,powerCam_arr,powerSetting_experiment,'linear','extrap');
else
    p2 = interp1(powerSource_arr,powerCam_arr,powerSetting_experiment);
end
normSource = p2/p1;

normPower = normPower0*normSource;
peakInt = ma*normSource;

normPower_normalization = normPower;
peakInt_normalization = peakInt;
save('utils\BenchModelNormalization_0803.mat','normPower_normalization','peakInt_normalization','tint_normalization')

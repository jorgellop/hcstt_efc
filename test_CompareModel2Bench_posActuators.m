% test_CompareModel2Bench_posActuators
% 
% 
%
% Jorge Llop - Mar 9, 2018

clear all;
close all;

addpath(genpath('utils'));

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

apRad = 116;
tint_normalization = 0.1;
powerSetting = 16;

label = '_0309';
outDir = ['output',filesep,'test_Normalization',label,filesep];
mkdir(outDir);

Nact = 12;
sidepix = 10;
sidepix_pow = 5;

hcstt_Initialize(false)

Ncam = 400;
% Take background image
take_background = false;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt )
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint_normalization);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    info.background = background;
    prompt = 'Put back light on. Continue? ';
    x = input( prompt )
else
    background = zeros(Ncam,Ncam);
end

% Find center
im_cam = zeros(400,400);
tint_findCen = 0.05;
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_findCen)-background;
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
% p=FastPeakFind(im_cam, 3 , 4 , 2, 2);
% ind_ma_I = p(1);
% ind_ma_J = p(2);
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);
% ind_ma_I = 200;
% ind_ma_J =200;
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
% im_cam_crop = im_cam_crop/max(im_cam_crop(:));

figure(100)
imagesc(im_cam_crop)
axis image

%% Apply a sinusoid to the Dm and compare model to camera
ho = 200;
angDM = pi/4;
actxc = 3;
DM_Command = hcstt_DMMapSin(ho, angDM, actxc, 0);    

% Take camera image
hcstt_UpdateMultiDM(DM_Command)
im_cam = hcstt_TakeCamImage(true,false,tint_normalization);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);

p=FastPeakFind(im_cam_crop, 5 , 40 , 2, 2);
 
figure(101)
imagesc(im_cam_crop)
axis image

info.Nact = Nact;
FPM = false;

info.normalize = false;

for II=1:1
    [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,1);
    info.posDM_x = posDM_x;
    info.posDM_y = posDM_y;    
    info.apRad = apRad;
    info.ac_spac = ac_spac;
    
    im_mod = hcstt_TakeModelImage(DM_Command(:)*1e-9,FPM,info);
    im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);

    total_power_mod = sum(sum(im_mod(N/2-sidepix_pow+1:N/2+sidepix_pow+1,N/2-sidepix_pow+1:N/2+sidepix_pow+1)));
    total_power_cam = sum(sum(im_cam(ind_ma_I-sidepix_pow:ind_ma_I+sidepix_pow,ind_ma_J-sidepix_pow:ind_ma_J+sidepix_pow)));

    normPower = total_power_cam/total_power_mod;
    im_mod = im_mod*normPower;
    im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);

    im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
    figure(2)
    imagesc(im_mod_crop);
    axis image
    drawnow
end
%% Disconnect Devices
hcstt_DisconnectDevices();

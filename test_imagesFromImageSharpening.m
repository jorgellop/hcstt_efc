% test_imagesFromImageSharpening.m
%
% Take pictures with camera of before and after Image Sharpening
%
% Jorge Llop - Aug 1, 2018
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

label = '_0801';
outDir = ['output',filesep,'test_imagesFromImageSharpening',label,filesep];
mkdir(outDir);

hcstt_Initialize(false);

Ncam = 400;
tint = 0.05;
sidepix = 30;
openhimg = false;

im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,openhimg,tint);
    im_cam = im_cam + im_camII;
    pause(0.1)
end
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);
im_cam_crop0 = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
im_cam_crop0 = im_cam_crop0/15;

% Flat DM
im_cam = hcstt_TakeCamImage(true,openhimg,tint);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(111)
imagesc(im_cam_crop)
colorbar
axis image
title('Flat DM - 0.05sec')
im_cam = hcstt_TakeCamImage(true,openhimg,tint*10);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(112)
imagesc(im_cam_crop)
colorbar
axis image
title('Flat DM - 0.5sec')

% New Flat DM
hcstt_NewFlatForDM('ImageSharpening_fmincon_0801')
im_cam = hcstt_TakeCamImage(true,openhimg,tint);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(113)
imagesc(im_cam_crop)
colorbar
axis image
title('New Flat w/ fmincon (cost fct = diff model) - 0.05sec')
im_cam = hcstt_TakeCamImage(true,openhimg,tint*10);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(114)
imagesc(im_cam_crop)
colorbar
axis image
title('New Flat w/ fmincon (cost fct = diff model) - 0.5sec')

% New Flat DM
hcstt_NewFlatForDM('ImageSharpening_fmincon_0801_v3FromZeroX0')
im_cam = hcstt_TakeCamImage(true,openhimg,tint);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(115)
imagesc(im_cam_crop)
colorbar
axis image
title('New Flat w/ fmincon (cost fct = Strehl) - 0.05sec')
im_cam = hcstt_TakeCamImage(true,openhimg,tint*10);
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(116)
imagesc(im_cam_crop)
colorbar
axis image
title('New Flat w/ fmincon (cost fct = Strehl) - 0.5sec')

hcstt_DisconnectDevices()

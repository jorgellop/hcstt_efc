% test_FindEstimateOfNormPower
%
%
%
% Jorge Llop, Feb 15, 2018

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false);

totalPowerModel = 3.7057e+10;
max_im_mod = 1.2529e+09;
pix_crop = 30;

texp1 = 0.03;
im_cam = hcstt_TakeCamImage(true,false,texp1);
imagesc(im_cam(180:220,180:220))
axis image
colorbar
[ma,ind] = max(im_cam(:));
sz = size(im_cam);
[indx,indy] = ind2sub(sz,ind);
PeakPower = sum(sum(im_cam(indx-5:indx+5,indy-5:indy+5)));
max_im_cam = max(im_cam(:));

texp2 = 4;
im_cam = hcstt_TakeCamImage(true,true,texp2);
sz = size(im_cam);
PeakPower_aux = sum(sum(im_cam(indx-5:indx+5,indy-5:indy+5)));

imagesc(im_cam(180:220,180:220))
axis image
colorbar

totalPower = sum(sum(im_cam(sz(1)/2-pix_crop:sz(1)/2+pix_crop,sz(2)/2-pix_crop:sz(2)/2+pix_crop)))-PeakPower_aux...
    +PeakPower* sqrt(texp2/texp1);
normPower = totalPower/totalPowerModel;

hcstt_DisconnectDevices();
eragergf
% test_findFiberPositionWithArtSpeckle
%
%
%
% Jorge Llop - Aug 31, 2018

addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false)
us = hcstt_DMMapSin(100, 0, 2.5, 0);
hcstt_UpdateMultiDM(+us)
tint = 10;
im_cam = hcstt_TakeCamImage(true,false,tint);
x_cent_cam = 200;
y_cent_cam = 200;
figure(1)
imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
axis image
colorbar
hcstt_DisconnectDevices()

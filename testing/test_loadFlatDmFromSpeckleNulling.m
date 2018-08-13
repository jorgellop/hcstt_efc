% test_loadFlatDmFromSpeckleNulling
%
%
%
% Jorge Llop - Feb 22 

addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false)

tint = 3000;

sidepix = 20;
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint);
    im_cam = im_cam + im_camII;
    pause(0.1)
end
ind_ma_I = 200;
ind_ma_J = 200;

im_cam_crop0 = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
im_cam_crop0 = im_cam_crop0/15;
figure(100)
imagesc(im_cam_crop0);
axis image
drawnow

hcstt_NewFlatForDM('SpeckleNulling_flat_Feb22')
hcstt_UpdateMultiDM(zeros(12,12));
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint);
    im_cam = im_cam + im_camII;
    pause(0.1)
end

im_cam_crop0 = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
im_cam_crop0 = im_cam_crop0/15;
figure(200)
imagesc(im_cam_crop0);
axis image
drawnow

hcstt_DisconnectDevices();
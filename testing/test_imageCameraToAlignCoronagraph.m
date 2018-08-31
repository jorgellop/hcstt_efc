% test_imageCameraToAlignCoronagraph
%
%
%
% Jorge Llop - Aug 31, 2018

addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false)
% us = hcstt_DMMapSin(100, 0, 2.5, 0);
hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');
totalPowerEFCSMF = 2.2331e-06;
hcstt_UpdateMultiDM(+zeros(12,12))
tint = 15;
while 1
    im_cam = hcstt_TakeCamImage(true,false,tint);
    x_cent_cam = 200;
    y_cent_cam = 200;
    figure(1)
    imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
    axis image
    colorbar
    drawnow
    int = hcstt_GetIntensityFIU(zeros(12,12),2,0);
    disp(log10(int/totalPowerEFCSMF))
end
hcstt_DisconnectDevices()

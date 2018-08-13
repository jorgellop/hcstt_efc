% test_CompareModel2Bench_forloop
% 
% Put sine functions on the DM in the model and in the bench and save
% images to compare
%
% Jorge Llop - Dec 14, 2017

clear all;
close all;

addpath(genpath('utils'));

load('angle_cal');
load('position_cal');

label = '_1214';
outDir = ['output',filesep,'test_CompareModel2Bench_sinefct',label,filesep];
mkdir(outDir);

Nact = 12;

x_c = 183;
y_c = 178;

x_s = 197;
y_s = 178;

h0  = 50;       %Amplitude
alp = 0;       %Phase delay

hcstt_Initialize()
n_pha = 8;
pha_arr = linspace(0,2*pi,n_pha);

openhimg = false;
im_cam = hcstt_TakeCamImage(true,openhimg,0.05);

hcstt_DisconnectDevices();
count = 1;

for II = 0:2
    distance_speckle = sqrt((x_s-x_c+II*10)^2+(y_s-y_c+II*10)^2);
    [c index] = min(abs(distance_speckle-position_cal(:,2)));
    act = position_cal(index,1);
    for JJ = 1:n_pha
        pha = pha_arr(JJ);
        angle_speckle = atan((y_s-y_c)/(x_s-x_c))+pha;
        [c index] = min(abs(angle_speckle-angle_cal(:,2)));
        ang = angle_cal(index,1);

        DM_Command = DE_DMMapSinv2(h0,ang,act,alp,1);
        
        info.Nact = Nact;
        info.N = 1024;
        im = hcstt_TakeModelImage(DM_Command(:),info);
        hcstt_test_plotModelImage(im, [outDir,'ModelImage_i',int2str(II),'_j',int2str(JJ)], info );
%         hcstt_UpdateMultiDM(DM_Command)
% 
%         if(count ~= 1) openhimg = true; end
%         im_cam = hcstt_TakeCamImage(false,openhimg);
%         sz = size(im_cam);
%         hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_i',int2str(II),'_j',int2str(JJ)], sz );
%         count = count + 1;
    end
end
hcstt_DisconnectDevices();

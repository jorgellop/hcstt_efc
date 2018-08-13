% hcstt_BenchModelPixelPitchMatching
% 
% Align model and bench to have the same pixel pitch
%
% Jorge Llop - Dec 31, 2018

clear all;
close all;

addpath(genpath('utils'));

load('angle_cal');
load('position_cal');

label = '_0131';
outDir = ['output',filesep,'test_PixPitchMatching',label,filesep];
mkdir(outDir);

Nact = 12;

h0  = 50;       %Amplitude
alp = 0;       %Phase delay

hcstt_Initialize()
n_pha = 8;
pha_arr = linspace(0,2*pi,n_pha);

openhimg = false;
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,openhimg,0.05);
    im_cam = im_cam + im_camII;
    pause(0.1)
end
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);

im_cam_crop = im_cam(ind_ma_I-10:ind_ma_I+10,ind_ma_J-10:ind_ma_J+10);
im_cam_crop = im_cam_crop/max(im_cam_crop(:));

DM_Command = zeros(12,12);        
info.Nact = Nact;
N = 1024;
info.N = N;

FPM = false;
apRad_arr = 64:2:256;
numtry = numel(apRad_arr);
diff_arr = zeros(numtry,1);
for II=1:numtry
    
    apRad = apRad_arr(II);
    info.apRad = apRad;
    
    im_mod = hcstt_TakeModelImage(DM_Command(:),FPM,info);
    im_mod_crop = im_mod(N/2-10:N/2+10,N/2-10:N/2+10);
    
    diff_mat = im_mod_crop - im_cam_crop;
%     imagesc(im_mod_crop);
%     axis image
%     drawnow;
    
    diff = sum(diff_mat(:));
    diff_arr(II) = diff;
end
[mi,ind_mi] = min(abs(diff_arr));
apRad = apRad_arr(ind_mi)
info.apRad = apRad; 
im_mod = hcstt_TakeModelImage(DM_Command(:),FPM,info);
im_mod_crop = im_mod(N/2-10:N/2+10,N/2-10:N/2+10);
figure(100)
imagesc(im_mod_crop);
axis image
figure(200)
imagesc(im_cam_crop);
axis image

% plot(apRad_arr,diff_arr)
hcstt_DisconnectDevices();

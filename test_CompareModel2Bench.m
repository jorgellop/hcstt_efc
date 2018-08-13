clear all;
close all;

addpath(genpath('utils'));

label = '_1213';
outDir = ['output',filesep,'test_CompareModel2Bench',label,filesep];
mkdir(outDir);

Nact = 12;
%Create the u matrix of actuator heights in nm
u_mat = zeros(Nact,Nact);
u_mat(1,5) = 0;

%% HCSTT Model Image
% info.Nact = Nact;
% info.N = 1024;
% im = hcstt_TakeModelImage(u_mat(:),info);
% hcstt_test_plotModelImage(im, [outDir,'ModelImage'], info );
%% HCSTT Bench Image
hcstt_Initialize()

%Update the DM shape
hcstt_UpdateMultiDM(u_mat)

%Take image with the camera
im_cam = hcstt_TakeCamImage(false);
sz = size(im_cam);

%Plot and save image
hcstt_test_plotCamImage(im_cam, [outDir,'CamImage'], sz );

%Disconnect all devices
hcstt_DisconnectDevices();

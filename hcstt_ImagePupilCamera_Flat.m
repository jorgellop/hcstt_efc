close all;
clear all;

%{
hcstt_ImagePupilCamera_Flat

With the camera placed at the pupil plane, take an image of the pupil.

Jorge Llop - Feb 12, 2017
%}

addpath(genpath('utils'),genpath('export_scripts'));
label = '_Feb13';
outDir = ['utils',filesep,'AlignmentDMModel',filesep,'ImagePupilCamera',label,filesep];
mkdir(outDir);

hcstt_Initialize(false);

im_cam = hcstt_TakeCamImage_pupil(true,false,10);
figure(100)
imagesc(im_cam)
axis image
colorbar
drawnow;

full_path = [outDir,'PupilImage'];
hcstt_test_plotCamImage(im_cam, full_path , size(im_cam));
im_pupil = im_cam;
save([outDir,'PupilImage',label,'.mat'],'im_pupil')
hcstt_DisconnectDevices();

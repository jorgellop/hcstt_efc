% test_CompareModel2Bench_forloop
% 
% Poke Each actuator both in the model and in the bench and save images in
% the image plane
%
% Jorge Llop - Dec 14, 2017

clear all;
close all;

addpath(genpath('utils'));

label = '_1214v2';
outDir = ['output',filesep,'test_CompareModel2Bench',label,filesep];
mkdir(outDir);

Nact = 12;

hcstt_Initialize()
count = 1;
openhimg = false;
for II = 1:Nact
    for JJ = 1:Nact
        disp(['iteration_i',int2str(II),'_j',int2str(JJ)])
        
        %Create the u matrix of actuator heights in nm
        u_mat = zeros(Nact,Nact);
        u_mat(II,JJ) = 50;

        %%HCSTT Model Image
        info.Nact = Nact;
        info.N = 1024;
        im = hcstt_TakeModelImage(u_mat(:),info);
        hcstt_test_plotModelImage(im, [outDir,'ModelImage_i',int2str(II),'_j',int2str(JJ)], info );
        
        %%HCSTT Bench Image
        %Update the DM shape
        hcstt_UpdateMultiDM(u_mat)

        %Take image with the camera
        if(count ~= 1) openhimg = true; end
        im_cam = hcstt_TakeCamImage(false,openhimg);
        sz = size(im_cam);

        %Plot and save image
        hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_i',int2str(II),'_j',int2str(JJ)], sz );
        count = count + 1;
    end
end
%Disconnect all devices
hcstt_DisconnectDevices();

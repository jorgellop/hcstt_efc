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

distance_speckle = sqrt((x_s-x_c)^2+(y_s-y_c)^2);
angle_speckle = atan((y_s-y_c)/(x_s-x_c));

%convert, define location of speckle
[c index] = min(abs(distance_speckle-position_cal(:,2)));
act = position_cal(index,1);

[c index] = min(abs(angle_speckle-angle_cal(:,2)));
ang = angle_cal(index,1);
par_probe = [50, ang, act, 0];
h0  = par_probe(1);       %Amplitude
q   = par_probe(2);       %Angle of Sin
x   = par_probe(3);       %Actuators/cycle
alp = par_probe(4);       %Phase delay

hcstt_Initialize()

DM_Command = DE_DMMapSin(h0,q,x,alp);
hcstt_UpdateMultiDM(DM_Command)

im_cam = hcstt_TakeCamImage(false,false);

intensity = DE_supMinimizer_3(par_probe);


% count = 1;
% openhimg = false;
% for II = 1:Nact
%     for JJ = 1:Nact
%         disp(['iteration_i',int2str(II),'_j',int2str(JJ)])
%         
%         %Create the u matrix of actuator heights in nm
%         u_mat = zeros(Nact,Nact);
%         u_mat(II,JJ) = 50;
% 
%         %%HCSTT Model Image
%         info.Nact = Nact;
%         info.N = 1024;
%         im = hcstt_TakeModelImage(u_mat(:),info);
%         hcstt_test_plotModelImage(im, [outDir,'ModelImage_i',int2str(II),'_j',int2str(JJ)], info );
%         
%         %%HCSTT Bench Image
%         %Update the DM shape
%         hcstt_UpdateMultiDM(u_mat)
% 
%         %Take image with the camera
%         if(count ~= 1) openhimg = true; end
%         im_cam = hcstt_TakeCamImage(false,openhimg);
%         sz = size(im_cam);
% 
%         %Plot and save image
%         hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_i',int2str(II),'_j',int2str(JJ)], sz );
%         count = count + 1;
%     end
% end
%Disconnect all devices
hcstt_DisconnectDevices();

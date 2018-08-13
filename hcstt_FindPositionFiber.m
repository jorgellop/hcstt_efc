% hcstt_FindPosiotionFiber
% 
% Find the position of the fiber wrt to the center of the star. It requires
% a first estimate of the postition, from which this code will scan around the
% neighbourhood and find the maximum FIU output power
%
% Jorge Llop - Dec 14, 2017

clear all;
close all;

addpath(genpath('utils'));
addpath(genpath('export_scripts'));

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
n_pha = 10;
n_pos = 10;
pha_arr = linspace(0,2*pi,n_pha); %rad
pos_arr = linspace(0,10,n_pos); %pix

int_mat = zeros(n_pos,n_pha);

count = 1;
openhimg = false;
for II = 1:n_pha
    distance_speckle = sqrt((x_s-x_c+pos_arr(II))^2+(y_s-y_c+pos_arr(II))^2);
    [c index] = min(abs(distance_speckle-position_cal(:,2)));
    act = position_cal(index,1);
    for JJ = 1:n_pha
        angle_speckle = atan((y_s-y_c)/(x_s-x_c))+pha_arr(JJ);
        [c index] = min(abs(angle_speckle-angle_cal(:,2)));
        ang = angle_cal(index,1);

        DM_Command = DE_DMMapSin(h0,ang,act,alp);
        hcstt_UpdateMultiDM(DM_Command)

        if(count ~= 1) openhimg = true; end
        im_cam = hcstt_TakeCamImage(false,openhimg);
        sz = size(im_cam);
        hcstt_test_plotCamImage(im_cam, [outDir,'CamImage_i',int2str(II),'_j',int2str(JJ)], sz );
        
%         int_mat(II,JJ) = 
        count = count + 1;
    end
end
% intensity = DE_supMinimizer_3(par_probe);


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

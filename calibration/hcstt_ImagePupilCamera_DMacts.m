% close all;
clear all;

%{
hcstt_ImagePupilCamera_DMacts

With the camera placed at the pupil plane, take images of each actuator.

Jorge Llop - Feb 12, 2017
%}

addpath(genpath('utils'),genpath('export_scripts'));
label = '_Jul26v2';
outDir = ['utils',filesep,'AlignmentDMModel',filesep,'ImagePupilCamera',label,filesep];
mkdir(outDir);
outDirPupIms = ['utils',filesep,'AlignmentDMModel',filesep,'ImagePupilCamera',label,filesep,'AllActuators'];
mkdir(outDirPupIms);

Nact = 12;
poke_amp = 300;

hcstt_Initialize(false);

%Take darks
% u = zeros(Nact,Nact);
% hcstt_UpdateMultiDM(u);
dark = 0;
% numdark = 5;
% for II = 1:numdark
%     im_cam = hcstt_TakeCamImage_pupil(true,false,1);
%     dark = dark + im_cam;
% end
% dark = dark/numdark;

% u = zeros(Nact,Nact);
% u(2,1:Nact) = poke_amp;
% u(1:Nact,2) = poke_amp;
% u(Nact-1,1:Nact) = poke_amp;
% u(1:Nact,Nact-1) = poke_amp;
% % u(1,1:Nact) = poke_amp;
% % u(1:Nact,1) = poke_amp;
% % u(Nact,1:Nact) = poke_amp;
% % u(1:Nact,Nact) = poke_amp;
% % u(2,1:Nact) = poke_amp;
% % u(1:Nact,2) = poke_amp;
% 
% % u(1,1:Nact) = poke_amp
% hcstt_UpdateMultiDM(u);
% im_cam = hcstt_TakeCamImage_pupil(true,false,10);
% % im_cam = im_cam-dark;
% 
% figure(100)
% imagesc(im_cam)
% axis image
% colorbar
% hcstt_DisconnectDevices();

%%
im_mat = zeros(1280, 1024,Nact^2);
for II = 1:Nact^2
    u = zeros(Nact,Nact);
    u(II) = poke_amp;
%      u=u';
    hcstt_UpdateMultiDM(u);
    
    im_cam = hcstt_TakeCamImage_pupil(true,false,0.06);
    im_mat(:,:,II) = im_cam-dark;
    [ind_x,ind_y] = ind2sub([Nact,Nact],II);
    figure(100)
%     title(['Actuator ',num2str(ind_x),', ',num2str(ind_y)])
    imagesc(im_cam-dark)
    axis image
    colorbar
    drawnow;
    full_path = [outDirPupIms,'PupilIm_Actx',num2str(ind_x),'y',num2str(ind_y)];
%     hcstt_test_plotCamImage(im_cam, full_path , size(im_cam));
disp(['Actuator ',num2str(ind_x),', ',num2str(ind_y)])
end

hcstt_DisconnectDevices();
save([outDir,'im_mat',label,'.mat'],'im_mat')

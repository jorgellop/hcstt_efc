% test_takeImages4CentroidAndGif
% 
% Take a bunch of images with the camera to show jitter
%
% Jorge Llop - Mar 7, 2018

clear all;
close all;

addpath(genpath('utils'),genpath('export_scripts'));

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

apRad = 116;
tint_centroid = 2;
% powerSetting = 16;

label = 'Jun29';
outDir = ['output',filesep,'test_CentroidJitter_DriftFIU_Downstream',label,filesep];
mkdir(outDir);

Nact = 12;
% sidepix = 20;
% sidepix_pow = 5;
%x=730, y=295


hcstt_Initialize(false)

% Find Center of camera image
im_cam = zeros(400,400);
tint_findCenter = 4.0;
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
imagesc(im_cam)
pause
Ncam = 400;
% Take background image
take_background = false;
tic
for JJ=1:1000
    disp(['Iteration',num2str(JJ)])
%     im_cam = hcstt_TakeCamImage(true,false,tint);
    
    im_cam = hcstt_TakeCamImage(true,false,tint_centroid);
    hcstt_test_plotCamImage(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20), [outDir,'CamImage_JitterTest_',num2str(JJ)], [41,41] );

    im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
    p=FastPeakFind(im_cam_crop/max(im_cam_crop(:)),1.2 , 2 , 2, 2);
    p1(JJ) = p(1)-21;
    p2(JJ) = p(2)-21;
    disp([p1(JJ),p2(JJ)])
    elapsedTime = toc;
    save([outDir,'CentroidData_p1p2.mat'],'p1','p2','elapsedTime');
end
fig0 = figure(100);
plot(1:numel(p1),p1)
hold on
plot(1:numel(p2),p2)
hold off
xlabel('Image number')
ylabel('Centroid position')
title(['Centroid position - x and y. Elapsed time: ', num2str(elapsedTime),'sec'])
% legend('Coupling SMF','Coupling MMF');
export_fig([outDir,'Centroid',label,'.png'],'-r300');

min_centroid_x = min(p1)
min_centroid_y = min(p2)
max_centroid_x = max(p1)
max_centroid_y = max(p2)
p1(find(abs(p1) >2)) = []; 
p2(find(abs(p2) >2)) = []; 
stddev_centroid_x = std(p1)
stddev_centroid_y = std(p2)

save([outDir,'CentroidData.mat'],'min_centroid_x','min_centroid_y','max_centroid_x','max_centroid_y','stddev_centroid_x','stddev_centroid_y');

hcstt_DisconnectDevices();
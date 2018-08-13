% hcstt_SN_calibrateAnglePos
% 
% Calibrate angle and position. v2: find 2 peaks by masking central spot.
% Corona in.
%
% Jorge Llop - Feb 1, 2018

clear all;
close all;

addpath(genpath('utils'));

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

label = '_0801';
outDir = ['output',filesep,'test_CalibrateAnglePos',label,filesep];
mkdir(outDir);

sidepix = 20;
tint = 0.2;

hcstt_Initialize(false)

openhimg = false;
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,openhimg,tint);
    im_cam = im_cam + im_camII;
    pause(0.1)
end
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);
im_cam_crop0 = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
im_cam_crop0 = im_cam_crop0/15;
figure(101)
imagesc(im_cam_crop0)
colorbar
axis image
title('Find center')

% Mask
[X,Y] = meshgrid(-sidepix:sidepix); 
[THETA, RHO] = cart2pol(X,Y);
mask = zeros(sidepix*2+1,sidepix*2+1);
mask(RHO>5) = 1;

ho = 70;
ph_delay = pi;
numtry = 11;
numang = 11;
numph = 4;
im_cam_crop = zeros(sidepix*2+1,sidepix*2+1);
actxcDM_arr = linspace(2,3,numtry);
ph_arr = linspace(0,pi,numph);
angDM_arr = linspace(-pi/5,pi/5,numang);
xmeas_arr = zeros(numtry,numang);
ymeas_arr = zeros(numtry,numang);
angmeas_arr = zeros(numtry,numang);
% count=1;
for KK = 1:numang
    for II = 1:numtry
%         ang = 0;
        actxc = actxcDM_arr(II);
        im_cam_crop = zeros(sidepix*2+1,sidepix*2+1);
        for JJ = 1:numph
%             DM_Command = hcstt_DMMapSin(ho, angDM_arr(KK), actxc, ph_arr(JJ));
            DM_Command = hcstt_DMMapSin(ho, 0, actxc, ph_arr(JJ));
            hcstt_UpdateMultiDM(DM_Command)

            im_cam = hcstt_TakeCamImage(true,false,tint);

            im_cam_cropJJ = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
            im_cam_crop = im_cam_crop+im_cam_cropJJ;
            pause(0.1)
        end
%         diff_im = abs( im_cam_crop/numph - im_cam_crop0 );
%         im_cam_crop = im_cam_crop.*mask;
        im_cam_crop = im_cam_crop;
        
        p=FastPeakFind(im_cam_crop/numph, 12 +II, 1 , 2, 2)

        figure(201)
        imagesc(im_cam_crop/numph);
        axis image
        colorbar
        title(['Angle: ',num2str(angDM_arr(KK))])
        drawnow

%         sz = size(p)
        if(numel(p)==4)
            x = (p(3)-p(1))/2
            y = (p(4)-p(2))/2
            xmeas_arr(II,KK) = x;
            ymeas_arr(II,KK) = y;
            an = atan(y/x)
            angmeas_arr(II,KK) = an;
        end
%         count = count + 1;
    end
end
dist_mat = sqrt(xmeas_arr.^2+ymeas_arr.^2);
% dist_meas = zeros(numtry,1);
% ang_meas = zeros(numang,1);
actxcDM_arr0 = actxcDM_arr;
clear actxcDM_arr;
count = 1;
for II = 1:numtry
    aux = dist_mat(II,:);
    if nnz(aux)>0
        aux2 = aux(find(aux));
        distPix_meas(count) = mean(aux2);
        actxcDM_arr(count) = actxcDM_arr0(II);
        count = count+1;
    end
end
angDM_arr0 = angDM_arr;
clear angDM_arr;
count = 1;
for II = 1:numang
    aux = angmeas_arr(:,II);
    if nnz(aux)>0
        aux2 = aux(find(aux));
        angPix_meas(count) = mean(aux2);
        angDM_arr(count) = angDM_arr0(II);
        count = count + 1;
    end
end

save('output\calibrateDM_Aug08','actxcDM_arr','angDM_arr','xmeas_arr','ymeas_arr','angmeas_arr','distPix_meas','angPix_meas');
hcstt_DisconnectDevices();

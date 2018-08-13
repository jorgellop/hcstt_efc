% test_SpeckleNulling
% 
% Try speckle nulling
%
% Jorge Llop - Feb 1, 2018

clear all;
close all;

global cam img Xcr Ycr CExp s drv_inf flat

addpath(genpath('utils'),genpath('export_scripts'));

load('output\calibrateDM_Feb15');
load('BenchModelNormalization_Mar07');


% load('position_cal');

N = 1024;
info.N = N;

label = '_0307';
outDir = ['output',filesep,'test_SpeckleNulling',label,filesep];
mkdir(outDir);

sidepix = 20;
tint = 30;

peakInt = peakInt_normalization*tint/tint_normalization;

hcstt_Initialize(false)

% Find Center
im_cam = zeros(400,400);
tint_findCenter = 1.1;
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[ind_ma_I,ind_ma_J] = ind2sub(size(im_camaux),ind_ma);

% Take initial image
openhimg = false;
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,openhimg,tint);
    im_cam = im_cam + im_camII;
    pause(0.1)
end
im_cam_crop0 = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
im_cam_crop0 = im_cam_crop0/15;

imagesc(im_cam_crop0/peakInt);
axis image
colorbar
drawnow

numit = 10;
intDH_arr = zeros(numit,1);
maxDH_arr = zeros(numit,1);
minDH_arr = zeros(numit,1);

pix_x_max = 3;
pix_x_min = -3;
pix_y_max = 16;%floor(max(distPix_meas));
pix_y_min = 11;%floor(min(distPix_meas));
% DH = im_cam_crop0(pix_x_min+21:pix_x_max+21,pix_y_min+21:pix_y_max+21);
im_cam_cropKK = im_cam_crop0;

aux = zeros(size(im_cam_cropKK));
aux(pix_x_min+sidepix+1:pix_x_max+sidepix+1,pix_y_min+sidepix+1:pix_y_max+sidepix+1) = 1;
ind_DH = find(aux);
extensionDH = 1;
aux(pix_x_min-extensionDH+sidepix+1:pix_x_max+extensionDH+sidepix+1,pix_y_min-extensionDH+sidepix+1:pix_y_max+extensionDH+sidepix+1) = 1;
ind_DHext = find(aux);
for KK = 1:numit
    disp(['Iteration ',num2str(KK)])
    
    %Dark hole
    DH = im_cam_cropKK(ind_DH);
    DHext = im_cam_cropKK(ind_DHext);
    [maDH,ind_maDHext] = max(DHext(:)); %Find max in dark hole to where we will apply the antispeckle
    [ind_ma_IDH,ind_ma_JDH] = ind2sub(size(zeros(pix_x_max-pix_x_min+1+extensionDH*2,pix_y_max-pix_y_min+1+extensionDH*2)),ind_maDHext);

    intKK0 = mean(mean(DH));
    disp('Intensity in DH:')
    disp(intKK0)
    
    %Save mean intensity in the DH
    intDH_arr(KK) = mean(DH(:));
    maxDH_arr(KK) = max(DH(:));
    minDH_arr(KK) = min(DH(:));
    
    pix_x = ind_ma_IDH+20+pix_x_min-extensionDH;
    pix_y = ind_ma_JDH+20+pix_y_min-extensionDH;

    pix_x_c = pix_x-21;
    pix_y_c = pix_y-21;
    dist = sqrt(pix_x_c^2+pix_y_c^2);
    ang = atan(pix_x_c/pix_y_c);

    if dist<min(distPix_meas)
        actxc = interp1(distPix_meas,actxcDM_arr,dist,'cubic','extrap');
        angDM = interp1(angPix_meas,angDM_arr,ang);
    else
        actxc = interp1(distPix_meas,actxcDM_arr,dist);
        angDM = interp1(angPix_meas,angDM_arr,ang);
    end
    numamp = 10;
    ho_arr = linspace(0.1,5,numamp);
    numph = 12;
    ph_arr = linspace(0,2*pi,numph);
    int_arr = zeros(numamp,numph);
    for II = 1:numamp
        for JJ = 1:numph
            DM_Command = hcstt_DMMapSin(ho_arr(II), angDM, actxc, ph_arr(JJ));
            hcstt_UpdateMultiDM(DM_Command)

            im_cam = hcstt_TakeCamImage(true,false,tint);

            im_cam_cropJJ = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);

        %     im_cam_crop = im_cam_crop+im_cam_cropJJ;        imagesc(im_cam_crop);
%             int_arr(II,JJ) = im_cam_cropJJ(pix_x,pix_y);
            int_arr(II,JJ) = mean(im_cam_cropJJ(ind_DH));

            imagesc(im_cam_cropJJ/peakInt)
            axis image
            colorbar
            drawnow
            pause(0.1)
        end
    end
    [mi,ind_mi] = min(int_arr(:));
    [ind_mi_amp,ind_mi_ph] = ind2sub(size(int_arr),ind_mi);


    %Update the new image:
    DM_Command = hcstt_DMMapSin(ho_arr(ind_mi_amp), angDM, actxc, ph_arr(ind_mi_ph));
    hcstt_UpdateMultiDM(DM_Command)
    im_cam = hcstt_TakeCamImage(true,false,tint);
    im_cam_cropKK = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
    
    ho_used_arr(KK) = ho_arr(ind_mi_amp);
    ang_used_arr(KK) = angDM;
    
    maxKK = max(max(im_cam_cropKK(ind_DH)));
    intKK = mean(mean(im_cam_cropKK(ind_DH)));
    
    if(intKK<intKK0)
        flat = flat + DM_Command; %New shape of the Dm to which apply next sinusoid
    else
    	if(maxKK<maDH)
            flat = flat + DM_Command;
        else
            maxKK
            maDH
            disp('Max is higher now than before')
        end
    end
    if(maxKK>240)
        return
    end
    hcstt_UpdateMultiDM(zeros(12,12));
end
numit = KK;
% im_cam_crop0(pix_y,pix_x)
figure(100)
imagesc(im_cam_crop0/peakInt);
axis image
colorbar;
export_fig([outDir,'ImageCam_NoCorr'],'-r300');

figure(101)
imagesc(im_cam_cropKK/peakInt);
axis image
colorbar;
export_fig([outDir,'ImageCam_Corr'],'-r300');

figure(200)
plot(1:numit,intDH_arr(1:numit)/peakInt)
hold on
plot(1:numit,maxDH_arr(1:numit)/peakInt)
plot(1:numit,minDH_arr(1:numit)/peakInt)
hold off
legend('average','max','min')
export_fig([outDir,'IntVsIt'],'-r300');

flat_SN = flat;
save([outDir,'SpeckleNulling_flat.m'],'flat_SN')
hcstt_DisconnectDevices();


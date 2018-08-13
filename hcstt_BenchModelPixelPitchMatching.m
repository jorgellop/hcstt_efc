% hcstt_BenchModelPixelPitchMatching
% 
% Align model and bench to have the same pixel pitch
%
% Jorge Llop - Dec 31, 2018

clear all;
close all;

addpath(genpath('utils'));

method =1;

% load('angle_cal');
% load('position_cal');

N = 1024;
info.N = N;

label = '_0806';
outDir = ['output',filesep,'test_PixPitchMatching',label,filesep];
mkdir(outDir);

Nact = 12;
sidepix = 30;

if method ~= 2
    hcstt_Initialize(false)

    openhimg = false;
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,openhimg,0.05);
        im_cam = im_cam + im_camII;
        pause(0.1)
    end
    [ma,ind_ma] = max(im_cam(:));
    [ind_ma_I,ind_ma_J] = ind2sub(size(im_cam),ind_ma);

    im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
    im_cam_crop = im_cam_crop/max(im_cam_crop(:));

    figure(100)
    imagesc(im_cam_crop)
    axis image
end
DM_Command = zeros(12,12);        
info.Nact = Nact;
FPM = false;

findApRad = false;
%% Find apRad
if(findApRad)
    info.normalize = true;
    info.posDM_x = 0;
    info.posDM_y = 0;
    apRad_arr = 100:2:160;
    numtry = numel(apRad_arr);
    diff_arr = zeros(numtry,1);
    for II=1:numtry

        apRad = apRad_arr(II);
        info.apRad = apRad;

        im_mod = hcstt_TakeModelImage(DM_Command(:),FPM,info);
        im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);

        diff_mat = im_mod_crop - im_cam_crop;
    %     imagesc(im_mod_crop);
    %     axis image
    %     drawnow;

        diff = sum(abs(diff_mat(:)));
        diff_arr(II) = diff;
    end
    [mi,ind_mi] = min(abs(diff_arr));
    apRad = apRad_arr(ind_mi)
    info.apRad = apRad; 
    im_mod = hcstt_TakeModelImage(DM_Command(:  ),FPM,info);
    im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
    figure(100)
    imagesc(im_mod_crop);
    axis image
    figure(200)
    imagesc(im_cam_crop);
    axis image
    drawnow;
    figure(101)
    imagesc(im_mod_crop);
    axis image
    figure(102)
    imagesc(log10(im_mod_crop));
    axis image
    title('Log Scale')
    drawnow;
else
    apRad = 142;    
    info.apRad = apRad;
    info.normalize = true;
    info.posDM_x = 0;
    info.posDM_y = 0;
    im_mod = hcstt_TakeModelImage(DM_Command(:),FPM,info);
    im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
    figure(100)
    imagesc(im_mod_crop);
    axis image
    figure(200)
    imagesc(im_cam_crop);
    axis image
    figure(102)
    imagesc(log10(im_mod_crop));
    axis image
    title('Log Scale')

end
%  hcstt_DisconnectDevices();
%  save([outDir,'BenchModelPixPitchMatch_Feb15.mat'],'apRad')
%% Find positions of DM actuators



ho = 200;
ang = 0;
ph_delay = 0;
ph_delay_arr = linspace(0,pi/2*3,4);
numtry = 5;
FPM = true;
info.apRad = apRad;
if method==2
    load('output\calibrateDM_Feb15');
    numtry2 = 15;
    scaleR_arr = linspace(0.835,0.89,numtry2);
    info.normalize = true;
    diff = zeros(numtry2,numtry);
    for JJ = 1:numtry2
        for II = 1:numtry
            actxc = actxcDM_arr(II);
            im_mod_crop = zeros(sidepix*2+1,sidepix*2+1);
            for KK=1:numel(ph_delay_arr)
                DM_Command = hcstt_DMMapSin(ho, ang, actxc,  ph_delay_arr(KK));
                scaleR = scaleR_arr(JJ);
                [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestR(N,apRad,scaleR);
                info.posDM_x = posDM_x;
                info.posDM_y = posDM_y;
                info.ac_spac = ac_spac;        
                im_mod = hcstt_TakeModelImage(DM_Command(:)*1e-9,FPM,info);
                im_mod_crop = im_mod_crop+im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
            end
           
%             im_mod_crop(sidepix-3:sidepix+5,sidepix-3:sidepix+5)=0;
            im_mod_crop = im_mod_crop/max(im_mod_crop(:));
            p=FastPeakFind(im_mod_crop, 0.2 , 1, 1, 2)
            p1_mod(II)=p(1);
            p2_mod(II)=p(2);
            figure(101)
            imagesc(im_mod_crop);
            axis image
            colorbar
            drawnow
            ddd
        end
    %     diff1 = p1_mod-p1_cam;
    %     diff2 = p2_mod-p2_cam;
    end

elseif method==1
%Camera image
actxc_arr = linspace(2,3,numtry);
ph_arr = linspace(0,3*pi/2,numtry);
im_cam_crop_mat = zeros(sidepix*2+1,sidepix*2+1,numtry);
% im_cam_crop_tot = zeros(sidepix*2+1,sidepix*2+1);
for II = 1:numtry
    actxc = actxc_arr(II);
    im_cam_crop = zeros(sidepix*2+1,sidepix*2+1);
    for KK = 1:numel(ph_delay_arr)
        DM_Command = hcstt_DMMapSin(ho, ang, actxc, ph_delay_arr(KK));
        hcstt_UpdateMultiDM(DM_Command)
        im_cam = hcstt_TakeCamImage(true,openhimg,0.05);
        im_cam_crop = im_cam_crop+ im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
    end
    im_cam_crop = im_cam_crop/max(im_cam_crop(:));
    figure(200)
    imagesc(im_cam_crop);
    axis image
    colorbar
    drawnow
%     p=FastPeakFind(im_cam_crop, 0.3 , 0.8, 2, 2)
%     prompt = 'Input p1';
%     p1 = input(prompt);
%     prompt = 'Input p2';
%     p2 = input(prompt);
%     p1_cam(II)=p1;
%     p2_cam(II)=p2;
    im_cam_crop_mat(:,:,II) = im_cam_crop;
%     im_cam_crop_tot = im_cam_crop_tot+im_cam_crop;
end



numtry2 = 15;
scaleR_arr = linspace(0.835,0.89,numtry2);
info.normalize = true;
diff = zeros(numtry2,numtry);
for JJ = 1:numtry2
    for II = 1:numtry
        actxc = actxc_arr(II);
        im_mod_crop = zeros(sidepix*2+1,sidepix*2+1);
        for KK=1:numel(ph_delay_arr)
            DM_Command = hcstt_DMMapSin(ho, ang, actxc,  ph_delay_arr(KK));
            scaleR = scaleR_arr(JJ);
            [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestR(N,apRad,scaleR);
            info.posDM_x = posDM_x;
            info.posDM_y = posDM_y;
            info.ac_spac = ac_spac;        
            im_mod = hcstt_TakeModelImage(DM_Command(:)*1e-9,FPM,info);
            im_mod_crop = im_mod_crop+im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
        end
        im_mod_crop = im_mod_crop/max(im_mod_crop(:));
        p=FastPeakFind(im_mod_crop, 0.2 , 0.8, 1, 2);
        p1_mod(II)=p(1);
        p2_mod(II)=p(2);
        figure(101)
        imagesc(im_mod_crop);
        axis image
        colorbar
        drawnow
        diff_mat = im_cam_crop_mat(:,:,II)-im_mod_crop;
        diff(JJ,II) = std(diff_mat(:));
        
    end
%     diff1 = p1_mod-p1_cam;
%     diff2 = p2_mod-p2_cam;
end
for II=1:numtry
    [mi,ind_min] = min(diff(:,II));
    scaleR_solution_arr(II) = scaleR_arr(ind_min);
end
%         dist1_mod = sqrt((p(1)-p(3))^2+(p(2)-p(4))^2)
%         dist2_mod = sqrt((p(5)-p(3))^2+(p(6)-p(4))^2)
%         dist_mod = mean([dist1_mod,dist2_mod])
%         dist_mod_arr(KK) = dist_mod;
hcstt_DisconnectDevices();
ac_spac_arr = (2*apRad./act_in_apRad_arr);
end
% save([outDir,'BenchModelPixPitchMatch_and_apRadPupilModel.mat'],'im_cam_crop_mat','dist_mod_arr','dist_cam','apRad','ac_spac_arr')
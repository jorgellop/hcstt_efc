% hcstt_ImageShaperning
%
% Fit a set of zernike polynomials to sharpen the image
%
% Jorge Llop - Mar 28, 2018

clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

global cam img Xcr Ycr CExp s drv_inf flat himg pm_scale

hcstt_Initialize(false);

Ncam = 400;
tint = 0.05;
% tint = 6;
info.tint = tint;
avgframes = 10;

take_background = true;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    info.background = background;
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    background = zeros(Ncam,Ncam);
end

im_cam = zeros(Ncam,Ncam);
tint_findCenter = 0.05;
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);

%% Flat DM
im_cam = zeros(Ncam,Ncam);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint)-background;
    im_cam = im_cam + im_camII/15;
    pause(0.05)
end
im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
cmin = min(im_cam_crop(:));
cmax = max(im_cam_crop(:));
figure(100)
imagesc(im_cam_crop,[cmin cmax])
title('Flat DM')
axis image
colorbar
im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);
max0 = sum(im_cam_crop_peak(:))

Nact = 12;
% x_pupil = (-fix(Nact/2):fix((Nact-1)/2)) ;
x_pupil = -(Nact-1)/2:(Nact-1)/2 ;
[X_pupil,Y_pupil] = meshgrid(x_pupil);
[THETAcam, R_pupil] = cart2pol(X_pupil,Y_pupil);
R_pupil = R_pupil/max(R_pupil(:));

% export_fig('flatDM_6sec.png','-r300');

%% New flat - manual 
% hcstt_NewFlatForDM('ImageSharpening_Spheri_Mar28');
% hcstt_UpdateMultiDM(zeros(Nact,Nact));
% im_cam = zeros(Ncam,Ncam);
% for II=1:15
%     im_camII = hcstt_TakeCamImage(true,false,tint)-background;
%     im_cam = im_cam + im_camII/15;
%     pause(0.05)
% end
% im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
% figure(300)
% imagesc(im_cam_crop,[cmin cmax])
% title('New flat - spherical')
% axis image
% colorbar
% im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);
% max0newflat = sum(im_cam_crop_peak(:))


%% New flat - fmincon
% hcstt_NewFlatForDM('ImageSharpening_fmincon_Mar30v3');
% hcstt_UpdateMultiDM(zeros(Nact,Nact));
% im_cam = zeros(Ncam,Ncam);
% for II=1:15
%     im_camII = hcstt_TakeCamImage(true,false,tint)-background;
%     im_cam = im_cam + im_camII/15;
%     pause(0.05)
% end
% im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
% figure(400)
% imagesc(im_cam_crop,[cmin cmax])
% title('New flat - fmincon v3')
% axis image
% colorbar
% im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);
% max0newflat = sum(im_cam_crop_peak(:))

% export_fig('correctefdfmincon_005sec.png','-r300');

% %% New flat - fmincon
% hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');
% hcstt_UpdateMultiDM(zeros(Nact,Nact));
% im_cam = zeros(Ncam,Ncam);
% for II=1:15
%     im_camII = hcstt_TakeCamImage(true,false,tint)-background;
%     im_cam = im_cam + im_camII/15;
%     pause(0.05)
% end
% im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
% figure(500)
% imagesc(im_cam_crop,[cmin cmax])
% title('New flat - fmincon it 2')
% axis image
% colorbar
% im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);
% max0newflat = sum(im_cam_crop_peak(:))

% export_fig('correctefdfmincon_6sec.png','-r300');

%% fmincon
% numzern = 12;
% % A = eye(numzern,numzern);
% lb = ones(numzern,1)*(-75);
% ub = ones(numzern,1)*75;
% % b = ones(numzern,1)*75;
% x0 = zeros(numzern,1);
% % x0 = [2, -1, 2.75, 12, 33.5, 0.5, 0.4, 0.5, 0.6, 1, 10, 4]; % From model-based image sharpening
% x = fmincon(@fct2min_ZernikeImageSharpening,x0,[],[],[],[],lb,ub);
% 
% Z2 = 2*R_pupil.*sin(THETAcam);
% Z3 = 2*R_pupil.*cos(THETAcam);
% Z4 = sqrt(3)*(2*R_pupil.^2-1);
% Z5 = sqrt(6)*(R_pupil.^2).*sin(2*THETAcam);
% Z6 = sqrt(6)*(R_pupil.^2).*cos(2*THETAcam);
% Z7 = sqrt(8)*(3*R_pupil.^3 - 2*R_pupil).*sin(THETAcam);
% Z8 = sqrt(8)*(3*R_pupil.^3 - 2*R_pupil).*cos(THETAcam);
% Z9 = sqrt(8)*(R_pupil.^3).*sin(3*THETAcam);
% Z10 = sqrt(8)*(R_pupil.^3).*cos(3*THETAcam);
% Z11 =  sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
% Z12 = sqrt(10)*(4*R_pupil.^4 - 3*R_pupil.^2).*cos(2*THETAcam);
% Z13 = sqrt(10)*(4*R_pupil.^4 - 3*R_pupil.^2).*sin(2*THETAcam);
% 
% W = x(1)*Z2 + x(2)*Z3 + x(3)*Z4 + x(4)*Z5 + x(5)*Z6 + x(6)*Z7 + x(7)*Z8 + x(8)*Z9 + ...
%     x(9)*Z10 + x(10)*Z11 + x(11)*Z12 + x(12)*Z13; 
% 
% hcstt_UpdateMultiDM(W)
% im_cam = zeros(Ncam,Ncam);
% for JJ=1:avgframes
%     int_cam0 = hcstt_TakeCamImage(true,false,tint)-background;
%     im_cam = im_cam + int_cam0/avgframes;
% end
% im_camaux = im_cam;
% im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
% [ma,ind_ma] = max(im_camaux(:));
% [x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
% im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
% im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);
% 
% maxfmincon = (sum(im_cam_crop_peak(:)))
% figure(101)
% imagesc(im_cam_crop)
% title(['Correcting fmincon '])
% axis image
% colorbar
% flat_SN = flat + W;
% save('utils\ImageSharpening_fmincon_0801_v3FromZeroX0.mat','flat_SN','x')

%% Spherical
% numtry = 20;
% RMS_SA_arr = linspace(-15,15,numtry);
% for II=1:numtry
%     RMS_SA = RMS_SA_arr(II); %nm
%     W = RMS_SA * sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
% %     W = RMS_SA * sqrt(5) * (R_pupil.^2).*cos(2*THETAcam);
%     hcstt_UpdateMultiDM(W)
%     im_cam = zeros(Ncam,Ncam);
%     for JJ=1:avgframes
%         int_cam0 = hcstt_TakeCamImage(true,false,tint)-background;
%         im_cam = im_cam + int_cam0/avgframes;
%     end
%     figure(101)
%     im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
%     imagesc(im_cam_crop)
%     title(['Correcting spherical - ',num2str(II)])
%     axis image
%     colorbar
%     drawnow
%     [ma,ind_ma] = max(im_cam_crop(:));
%     [x_cent_camII,y_cent_camII] = ind2sub(size(im_cam_crop),ind_ma);
%     im_cam_crop_peak = im_cam_crop(x_cent_camII-1:x_cent_camII+1,y_cent_camII-1:y_cent_camII+1);
%     max_arr(II) = sum(im_cam_crop_peak(:));
% end
% figure(301)
% plot(1:numtry,max_arr)
% [ma, ind_ma] = max(max_arr);
% W_spheri = RMS_SA_arr(ind_ma) * sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
% flat_SN = flat + W_spheri;
% save('utils\ImageSharpening_Spheri_Mar28.mat','flat_SN')

%% Astig x
numtry = 10;
RMS_SA_arr = linspace(-20,20,numtry);
for II=1:numtry
    RMS_SA = RMS_SA_arr(II); %nm
%     W = RMS_SA * sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
    W = RMS_SA * sqrt(5) * (R_pupil.^2).*cos(2*THETAcam);
    hcstt_UpdateMultiDM(W)
    im_cam = zeros(Ncam,Ncam);
    for JJ=1:avgframes
        int_cam0 = hcstt_TakeCamImage(true,false,tint)-background;
        im_cam = im_cam + int_cam0/avgframes;
    end
    figure(101)
    im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
    imagesc(im_cam_crop)
    title(['Correcting astig x - ',num2str(II)])
    axis image
    colorbar
    drawnow
    [ma,ind_ma] = max(im_cam_crop(:));
    [x_cent_camII,y_cent_camII] = ind2sub(size(im_cam_crop),ind_ma);
    im_cam_crop_peak = im_cam_crop(x_cent_camII-1:x_cent_camII+1,y_cent_camII-1:y_cent_camII+1);
    max_arr(II) = sum(im_cam_crop_peak(:));
end
figure(301)
plot(1:numtry,max_arr)
[ma, ind_ma] = max(max_arr);
W_astig_x = RMS_SA_arr(ind_ma) * sqrt(5) * (R_pupil.^2).*cos(2*THETAcam);
flat_SN = flat + W_astig_x;
save('utils\ImageSharpening_Astigx_0802.mat','flat_SN')

%% Astig y
numtry = 10;
RMS_SA_arr = linspace(-20,20,numtry);
for II=1:numtry
    RMS_SA = RMS_SA_arr(II); %nm
%     W = RMS_SA * sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
    W = RMS_SA * sqrt(5) * (R_pupil.^2).*sin(2*THETAcam);
    hcstt_UpdateMultiDM(W)
    im_cam = zeros(Ncam,Ncam);
    for JJ=1:avgframes
        int_cam0 = hcstt_TakeCamImage(true,false,tint)-background;
        im_cam = im_cam + int_cam0/avgframes;
    end
    figure(101)
    im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
    imagesc(im_cam_crop)
    title(['Correcting astig y - ',num2str(II)])
    axis image
    colorbar
    drawnow
    [ma,ind_ma] = max(im_cam_crop(:));
    [x_cent_camII,y_cent_camII] = ind2sub(size(im_cam_crop),ind_ma);
    im_cam_crop_peak = im_cam_crop(x_cent_camII-1:x_cent_camII+1,y_cent_camII-1:y_cent_camII+1);
    max_arr(II) = sum(im_cam_crop_peak(:));
end
figure(301)
plot(1:numtry,max_arr)
[ma, ind_ma] = max(max_arr);
W_astig_y = RMS_SA_arr(ind_ma) * sqrt(5) * (R_pupil.^2).*sin(2*THETAcam);
flat_SN = flat + W_astig_y;
save('utils\ImageSharpening_Astigy_0802.mat','flat_SN')

%% Disconnect
hcstt_DisconnectDevices();







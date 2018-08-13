function apRad = test_fmincon_ApRad_Model()
% Attempt at finding apRad.
%
% Milan Roberson
% 7/31/18


%% Setup
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

global cam img Xcr Ycr CExp s drv_inf flat himg pm_scale

hcstt_Initialize(false);

Ncam = 400;
tint = 0.05;
% tint = 6;
info.tint = tint;
avgframes = 100;
info.apRad = 68;
info.Nact = 12;
info.normalize = true;
info.posDM_x = 0;
info.posDM_y = 0;

cleanup = onCleanup(@() hcstt_DisconnectDevices);


%% Find center of PSF
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
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint);
    im_cam = im_cam + im_camII/avgframes;
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


norm_im = im_cam_crop / max(im_cam_crop(:));
x0 = 68;
ub = 100;
lb = 40;
x = fmincon(@(x)fn2min(x, norm_im, info), x0, [], [], [], [], lb, ub);


apRad = x;

mod_im = model_image(apRad, info);
[ma,ind_ma] = max(mod_im(:));
[x_cent,y_cent] = ind2sub(size(mod_im),ind_ma);
im_mod_crop = mod_im(x_cent-20:x_cent+20,y_cent-20:y_cent+20);
diff = im_mod_crop - norm_im;

imagesc(diff);
colorbar;

diff = fn2min(apRad, norm_im, info)

    function im_mod = model_image(apRad, info)
        Nact = 12;
        Ncam = 400;
        x_pupil = -(Nact-1)/2:(Nact-1)/2 ;
        [X_pupil,Y_pupil] = meshgrid(x_pupil);
        [THETA, R_pupil] = cart2pol(X_pupil,Y_pupil);
        R_pupil = R_pupil/max(R_pupil(:));
        
        info.apRad = apRad;
        im_mod = hcstt_TakeModelImage(zeros(Nact), false, info);
    end

    function diff = fn2min(x, im, info)
        mod_im = model_image(x, info);
        [ma,ind_ma] = max(mod_im(:));
        [x_cent,y_cent] = ind2sub(size(mod_im),ind_ma);
        if (x_cent + 20) > 1024 || (x_cent - 20) < 0 || (y_cent + 20) > 1024 || (y_cent - 20) < 0
            diff = Inf;
            return
        end
        im_mod_crop = mod_im(x_cent-20:x_cent+20,y_cent-20:y_cent+20);
        diff = im_mod_crop - im;
        
        [Xcam, Ycam] = meshgrid(-20:20);
        [~, Rcam] = cart2pol(Xcam, Ycam);
        weights = Rcam + sum(Rcam(:))/numel(Rcam);
        weights = weights.^-2;
        diff = sum(diff(:).^2 .* weights(:));
    end
end
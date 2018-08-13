%function [soln, norm_im] = test_ImageSharpening_Model()
% Attempt at modelling the Zernikes present on the DM from an image of the
% current flat.
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
info.apRad = 67;
info.Nact = 12;
info.normalize = true;

% load('utils/AlignmentDMModel/positionsDMActuators_raw');
% info.posDM_x = pos_x;
% info.posDM_y = pos_y;

info.posDM_x = 0;
info.posDM_y = 0;

cleanup = onCleanup(@() hcstt_DisconnectDevices);


%% Take background image
take_background = true;
if(take_background)
    prompt = 'Take out light. Continue? ';
    input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    info.background = background;
    prompt = 'Put back light on. Continue? ';
    input( prompt );
else
    background = zeros(Ncam,Ncam);
end

%% Find center of PSF
im_cam = zeros(Ncam,Ncam);
tint_findCenter = 0.05;
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter) - background;
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);


%% Flat DM

% Using previous result
hcstt_NewFlatForDM('utils/ImageSharpeningModel_0801_flatv2');
hcstt_UpdateMultiDM(zeros(12));
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
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
norm_im = im_cam_crop / max(im_cam_crop(:));


x0 = zeros(13,1);
% x0 = [-6, 1, 4, 11, 27, -.21, .4, .83, 1, 2.2, 6.8, 3.5, 87];
% x0 = randn(13, 1)*20;
x0(13) = 87;
ub = 100*ones(12, 1);
lb = -ub;
lb(13) = 50;
[x, fval, exitflag, output, l, grad, hess] = fmincon(@(x)fn2min(x, norm_im, info), x0, [], [], [], [], lb, ub);


soln = x;

[mod_im, W] = model_image(soln, info);
[ma,ind_ma] = max(mod_im(:));
[x_cent,y_cent] = ind2sub(size(mod_im),ind_ma);
im_mod_crop = mod_im(x_cent-20:x_cent+20,y_cent-20:y_cent+20);

diff = im_mod_crop - norm_im;

% flat_SN = flat - W*1e9;
% save('utils/ImageSharpeningModel_0801_flat.mat', 'flat_SN')
% save('utils/ImageSharpeningModel_0801.mat')

imagesc(diff);
colorbar;
delete cleanup;

    function [im_mod, W] = model_image(x, info)
        info.apRad = x(13);
        x(1:12) = x(1:12) * 1e-9;
        Nact = 12;
        Ncam = 400;
        x_pupil = -(Nact-1)/2:(Nact-1)/2 ;
        [X_pupil,Y_pupil] = meshgrid(x_pupil);
        [THETA, R_pupil] = cart2pol(X_pupil,Y_pupil);
        R_pupil = R_pupil/max(R_pupil(:));
        
        Z2 = 2*R_pupil.*sin(THETA);
        Z3 = 2*R_pupil.*cos(THETA);
        Z4 = sqrt(3)*(2*R_pupil.^2-1);
        Z5 = sqrt(6)*(R_pupil.^2).*sin(2*THETA);
        Z6 = sqrt(6)*(R_pupil.^2).*cos(2*THETA);
        Z7 = sqrt(8)*(3*R_pupil.^3 - 2*R_pupil).*sin(THETA);
        Z8 = sqrt(8)*(3*R_pupil.^3 - 2*R_pupil).*cos(THETA);
        Z9 = sqrt(8)*(R_pupil.^3).*sin(3*THETA);
        Z10 = sqrt(8)*(R_pupil.^3).*cos(3*THETA);
        Z11 =  sqrt(5) * (6*R_pupil.^4 - 6*R_pupil.^2 + 1);
        Z12 = sqrt(10)*(4*R_pupil.^4 - 3*R_pupil.^2).*cos(2*THETA);
        Z13 = sqrt(10)*(4*R_pupil.^4 - 3*R_pupil.^2).*sin(2*THETA);
        
        W = x(1)*Z2 + x(2)*Z3 + x(3)*Z4 + x(4)*Z5 + x(5)*Z6 + x(6)*Z7 + x(7)*Z8 + x(8)*Z9 + ...
            x(9)*Z10 + x(10)*Z11 + x(11)*Z12 + x(12)*Z13;
        
        im_mod = hcstt_TakeModelImage(W(:), false, info);
    end

    function diff = fn2min(x, im, info)
        [mod_im, W] = model_image(x, info);
        [ma,ind_ma] = max(mod_im(:));
        [x_cent,y_cent] = ind2sub(size(mod_im),ind_ma);
        if (x_cent + 20) > 1024 || (x_cent - 20) < 1 || (y_cent + 20) > 1024 || (y_cent - 20) < 1
            diff = Inf;
            return
        end
        im_mod_crop = mod_im(x_cent-20:x_cent+20,y_cent-20:y_cent+20);
        
        diff = im_mod_crop - im;
        
        diff = sum(diff(:).^2);
    end
%end
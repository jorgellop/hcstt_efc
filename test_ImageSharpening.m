function [soln, norm_im] = test_ImageSharpening()
% Attempt at matching modelled image by placing Zernikes on DM
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
tint = 10;
% tint = 2;
info.tint = tint;
avgframes = 100;
info.apRad = 90;
info.Nact = 12;
info.normalize = false;

% load('utils/AlignmentDMModel/positionsDMActuators_raw');
% info.posDM_x = pos_x;
% info.posDM_y = pos_y;

info.posDM_x = 0;
info.posDM_y = 0;

cleanup = onCleanup(@() hcstt_DisconnectDevices);


mask_R = 5;

% %% Take model ideal image
% 
% im_mod = hcstt_TakeModelImage(zeros(144,1), false, info);
% im_mod = im_mod * 70e-6;
% [m, n] = size(im_mod);
% [I, J] = ndgrid(1:m, 1:n);
% 
% [~,ind_ma] = max(im_mod(:));
% [x_cent_mod,y_cent_mod] = ind2sub(size(im_mod),ind_ma);
% 
% im_mod( (I - x_cent_mod).^2 + (J - y_cent_mod).^2 <= mask_R^2) = 0;
% 
% im_mod_crop = im_mod(x_cent_mod-40:x_cent_mod+40,y_cent_mod-40:y_cent_mod+40);
% 
% figure(100)
% imagesc(im_mod_crop)
% axis image
% colorbar
% % input('continue')


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
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
    im_cam = im_cam + im_camII/avgframes;
    pause(0.05)
end
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
cmin = min(im_cam_crop(:));
cmax = max(im_cam_crop(:));
figure(101)
imagesc(im_cam_crop,[cmin cmax])
title('Flat DM')
axis image
colorbar

[m, n] = size(im_cam);
[I, J] = ndgrid(1:m, 1:n);

% Mask out saturated pixels
im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
figure(102)
imagesc(im_cam_crop)
title('Flat DM, masked')
axis image
colorbar

%% Model-Based Image Sharpening

% Using previous result
hcstt_NewFlatForDM('utils/ImageSharpeningModel_0801_flatv2');
hcstt_UpdateMultiDM(zeros(12));
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
    im_cam = im_cam + im_camII/avgframes;
    pause(0.05)
end
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
cmin = min(im_cam_crop(:));
cmax = max(im_cam_crop(:));
figure(103)
imagesc(im_cam_crop,[cmin cmax])
title('Model-Based Image Sharpening')
axis image
colorbar

[m, n] = size(im_cam);
[I, J] = ndgrid(1:m, 1:n);

% Mask out saturated pixels
im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
figure(104)
imagesc(im_cam_crop)
title('Model-Based Image Sharpening, masked')
axis image
colorbar


%% fmincon Image Sharpening v2

% Using previous result
hcstt_NewFlatForDM('utils/ImageSharpening_fmincon_0801_v2');
hcstt_UpdateMultiDM(zeros(12));
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
    im_cam = im_cam + im_camII/avgframes;
    pause(0.05)
end
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
cmin = min(im_cam_crop(:));
cmax = max(im_cam_crop(:));
figure(105)
imagesc(im_cam_crop,[cmin cmax])
title('fmincon Image Sharpening, starting with MBIS')
axis image
colorbar

[m, n] = size(im_cam);
[I, J] = ndgrid(1:m, 1:n);

% Mask out saturated pixels
im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
figure(106)
imagesc(im_cam_crop)
title('fmincon Image Sharpening, starting with MBIS, masked')
axis image
colorbar


%% fmincon Image Sharpening v3

% Using previous result
hcstt_NewFlatForDM('utils/ImageSharpening_fmincon_0801_v3FromZeroX0');
hcstt_UpdateMultiDM(zeros(12));
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
    im_cam = im_cam + im_camII/avgframes;
    pause(0.05)
end
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
cmin = min(im_cam_crop(:));
cmax = max(im_cam_crop(:));
figure(107)
imagesc(im_cam_crop,[cmin cmax])
title('fmincon Image Sharpening, starting with flat')
axis image
colorbar

[m, n] = size(im_cam);
[I, J] = ndgrid(1:m, 1:n);

% Mask out saturated pixels
im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;
im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
figure(108)
imagesc(im_cam_crop)
title('fmincon Image Sharpening, starting with flat masked')
axis image
colorbar


% input('continue?')




%% Minimization
x0 = zeros(13,1);
x0(13) = 1;
ub = 100*ones(13, 1);
ub(13) = 2;
lb = -ub;
lb(13) = 0.5;
[x, fval, exitflag, output, l, grad, hess] = fmincon(@(x)fn2min(x, im_mod_crop), x0, [], [], [], [], lb, ub);


soln = x;


%% Build Zernikes to apply to DM
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

% Apply Zernikes to DM
hcstt_UpdateMultiDM(W);
pause(0.05);


% Acquire image
im_cam = zeros(Ncam,Ncam);
for II=1:avgframes
    im_camII = hcstt_TakeCamImage(true,false,tint) - background;
    im_cam = im_cam + im_camII/avgframes;
    pause(0.05)
end

% Find centroid of image (not sure if necessary)
im_camaux = im_cam/sum(im_cam(:));
[m, n] = size(im_camaux);
[I, J] = ndgrid(1:m, 1:n);

% Mask out saturated pixels
im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;


diff = im_mod - norm_im;

% flat_SN = flat - W*1e9;
% save('utils/ImageSharpeningModel_0801_flat.mat', 'flat_SN')
% save('utils/ImageSharpeningModel_0801.mat')

imagesc(diff);
colorbar;

    function diff = fn2min(x, im_mod)
        
        Nact = 12;
        
        tint = 10;
        % Build Zernikes to apply to DM
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
        
        % Apply Zernikes to DM
        hcstt_UpdateMultiDM(W);
        pause(0.05);
        
        % Acquire image
        im_cam = zeros(Ncam,Ncam);
        for II=1:avgframes
            im_camII = hcstt_TakeCamImage(true,false,tint) - background;
            im_cam = im_cam + im_camII/avgframes;
            pause(0.05)
        end
        
        % Find centroid of image (not sure if necessary)
        [m, n] = size(im_cam);
        [I, J] = ndgrid(1:m, 1:n);
        
        % Mask out saturated pixels
        im_cam( (I - x_cent_cam).^2 + (J - y_cent_cam).^2 <= mask_R^2) = 0;
        im_cam_crop = im_cam(x_cent_cam-40:x_cent_cam+40,y_cent_cam-40:y_cent_cam+40);
        
        diff = std2(im_mod*x(13) - im_cam_crop);
    end
end
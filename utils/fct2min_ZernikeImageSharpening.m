% fct2min_ZernikeImageSharpening
%
% Load new flat array for DM
%
% Jorge Llop - Mar 29, 2018

function f = fct2min_ZernikeImageSharpening(x)

avgframes = 100;
tint = 0.05;
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

hcstt_UpdateMultiDM(W)
im_cam = zeros(Ncam,Ncam);
for JJ=1:avgframes
    int_cam0 = hcstt_TakeCamImage(true,false,tint);
    im_cam = im_cam + int_cam0/avgframes;
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
im_cam_crop_peak = im_cam_crop(21-1:21+1,21-1:21+1);

f = 1/(sum(im_cam_crop_peak(:)));


end
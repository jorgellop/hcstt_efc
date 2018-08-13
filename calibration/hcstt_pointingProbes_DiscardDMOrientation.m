% test_pointingProbesDMOrientation.m
%
% Send a tilted probe on testbed & model to discard 4 out of 8 DM
% orientations
%
% Jorge Llop - May 3, 2018


N = 1024;
apRad = 142;
Nact = 12;

hcstt_Initialize(false);

x_fib_pix = 12;
% actxc = interp1(distPix_meas,actxcDM_arr,x_fib_pix-0.5);
actxc = 3;
angDM = 0;%pi/4;%interp1(angPix_meas,angDM_arr,ang);
dm_probcosfct = hcstt_DMMapSin(70, angDM, actxc, 0);

sidepix = 20;
numph = 4;
ph_arr = linspace(0,3*pi/2,numph);
%%
tint = 0.05;
background = 0;
% Find center
im_cam = zeros(400,400);
for II=1:15
    im_camII = hcstt_TakeCamImage(true,false,0.05)-background;
    im_cam = im_cam + im_camII/15;
    pause(0.1)
end
% p=FastPeakFind(im_cam, 3 , 4 , 2, 2);
% ind_ma_I = p(1);
% ind_ma_J = p(2);
[ma,ind_ma] = max(im_cam(:));
[ind_ma_I2,ind_ma_J2] = ind2sub(size(im_cam),ind_ma);
ind_ma_I = ind_ma_I2;
ind_ma_J = ind_ma_J2;
im_cam_crop = im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
figure(1)
imagesc(im_cam_crop)
axis image
colorbar
title('Flat DM')
hcstt_UpdateMultiDM(+dm_probcosfct(:))
im_cam = hcstt_TakeCamImage(true,false,tint)-background;
figure(2)
imagesc(im_cam(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix))
axis image

hcstt_DisconnectDevices();
%%
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);

info.useGPU = false;
info.apRad = apRad;
info.lambdaOverD =  N/apRad/2;
info.RHO = RHO;
info.THETA = THETA;
info.N = N;
info.lambda0 = 8e-7;
info.lam_arr = 8e-7;
info.numOfWavelengths = 1;
info.useApodizer = false;
info.useGPU = false; 
info.FPM = exp(1i*8*THETA);
info.xvals = xvals;
info.yvals = yvals;
info.use_fiber = false;
info.normal_EFC = true;

% info.LPM = ones(N,N);
info.LPM = exp(-(RHO/(0.8*apRad)).^1000);

poke = 90;
for II=1:1
    [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,II);
    info.posDM_x = posDM_x;
    info.posDM_y = posDM_y;
%         ac_spac = round(2*apRad/Nact);
    infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );
    wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
    wfin_noerrors(RHO > apRad) = 0;

    wf2_noerrors = zeros(N,N);
    for KK=1:numph
        dm_probcosfct = hcstt_DMMapSin(1, angDM, actxc, ph_arr(KK));
        count = 0; 
        DM1_strokes = zeros(N,N);
        for ix = 1:Nact
            for iy = 1:Nact
                count = count + 1;
                xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
                ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
    %             xpos = round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
    %             ypos = round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
                DM1_strokes(xpos,ypos) = dm_probcosfct(count)*poke*1e-9 + DM1_strokes(xpos,ypos);
    %             DM1_strokes(xpos,ypos) = 5*poke*1e-9 + DM1_strokes(xpos,ypos);
            end
        end
        surf_DM1 = conv2(DM1_strokes,infl,'same');
        wf2_prob_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM1, false, info);
        wf2_noerrors0 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, zeros(N,N), false, info);
        wf2_noerrors = wf2_noerrors + wf2_noerrors0;
    end

    int_mod = abs(wf2_prob_noerrors).^2;
    int_mod0 = abs(wf2_noerrors).^2;

    figure(400+II)
    % imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20) - int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
    imagesc(int_mod(N/2-sidepix:N/2+sidepix,N/2-sidepix:N/2+sidepix))
%     imagesc(log10(int_mod(N/2-sidepix:N/2+sidepix,N/2-sidepix:N/2+sidepix)))
    axis image
    title('Model - im')
    colorbar
    drawnow
    
end
figure(500)
imagesc(surf_DM1)
axis image
colorbar
title('DM shape')

% test_powerFactor4Poke.m
%
%
%
% Jorge Llop - May 3, 2018


N = 1024;
apRad = 120;

hcstt_Initialize(true);

x_fib_pix = 12;
% actxc = interp1(distPix_meas,actxcDM_arr,x_fib_pix-0.5);
actxc = 2.4;
angDM = pi/4;%interp1(angPix_meas,angDM_arr,ang);
dm_probcosfct = hcstt_DMMapSin(30, angDM, actxc, 0);

%%
hcstt_UpdateMultiDM(+dm_probcosfct(:))
im_cam = hcstt_TakeCamImage(true,false,0.5)-background;
figure(2)
imagesc(im_cam(200-20:200+20,200-20:200+20))
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

info.LPM = ones(N,N);
poke = 100;
Nact = 12;
for II=1:4
    [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,II);
    info.posDM_x = posDM_x;
    info.posDM_y = posDM_y;
    infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );
    wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
    wfin_noerrors(RHO > apRad) = 0;

    count = 0; 
    DM1_strokes = zeros(N,N);
    for ix = 1:Nact
        for iy = 1:Nact
            count = count + 1;
            xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
            DM1_strokes(xpos,ypos) = dm_probcosfct(count)*poke*1e-9 + DM1_strokes(xpos,ypos);
        end
    end
    surf_DM1 = conv2(DM1_strokes,infl,'same');
    wf2_prob_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM1, true, info);
    wf2_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, zeros(N,N), true, info);
    % wf2_prob_noerrors = wf2_prob_noerrors ;

    int_mod = abs(wf2_prob_noerrors).^2;
    int_mod0 = abs(wf2_noerrors).^2;

    figure(400+II)
    % imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20) - int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
    imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20))
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

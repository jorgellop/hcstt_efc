% test_pointingProbes
%
% Test the pointing of the probes to match the DH position
%
% Jorge Llop - Feb 21, 2018

clear all;
% close all;

addpath(genpath('utils'),genpath('export_scripts'));
load('output\calibrateDM_Feb15');

hcstt_NewFlatForDM('ImageSharpening_fminconIt2_Apr1');

label = '_May03v4';
outDir = ['output',filesep,'test_pointingProbes',label,filesep];
mkdir(outDir);

hcstt_Initialize(false);

N = 2^10;
% apRad = 2^6;
apRad = 120;

x_fib_pix = 12;
% actxc = interp1(distPix_meas,actxcDM_arr,x_fib_pix-0.5);
actxc = 2.4;
angDM = 0;%pi/4;%interp1(angPix_meas,angDM_arr,ang);

Ncam = 400;
tint_findCenter = 1.1;
im_cam = zeros(Ncam,Ncam);
for II=1:25
    im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
    im_cam = im_cam + im_camII/25;
    pause(0.1)
end
im_camaux = im_cam;
im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
[ma,ind_ma] = max(im_camaux(:));
[x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
info.x_cent_cam = x_cent_cam;
info.y_cent_cam = y_cent_cam;

[Xcam,Ycam] = meshgrid(-x_cent_cam:Ncam-x_cent_cam-1,-(Ncam-y_cent_cam):y_cent_cam-1); 
q_pix = 1;
info.q_pix = q_pix;
Q = zeros(Ncam,Ncam);
Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
Q = and(Q, Ycam >= -(q_pix) );
Q = and(Q, Ycam <= (q_pix) );


%%
% DM shape
config = 1;
% folder = 'C:\Users\COO_user\Documents\MATLAB\EFCwFiber_labTestv6\output\EFC_wFiber_LabDemonstration_May03v2_FlipDM\';
folder = 'C:\Users\COO_user\Documents\MATLAB\EFCwFiber_labTestv6\output\EFC_wFiber_LabDemonstration_May03\';
load([folder,'us0_DMconfig',num2str(config),'.mat'])
% dm_probcosfct = hcstt_DMMapSin(1, angDM, actxc, 0);
dm_probcosfct = us0(:)/2.0;%hcstt_DMMapSin(1, angDM, actxc, 0);
max_us0 = max(abs(us0))

%%
% Camera
%
tint = 0.5;
im_cam0 = hcstt_TakeCamImage(true,false,tint);

% im_cam0(Q) = 0;

figure(100)
imagesc(im_cam0(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
title('Flat DM')
axis image 

poke = 1; %nm
apRad2 = 12;
Nact = 12;
a = ones(apRad2,apRad2);
sincfct1 = sinc([2:Nact-1]/(Nact-1) * q_pix*2 * 1.2);
sincfct1 = [0,sincfct1,0];
di = diag(sincfct1);
dm_probsincfct1 = a * di; 

sincfct2 = sinc([2:Nact-1]/(Nact-1) * q_pix*2 * 1.2);
sincfct2 = [0,sincfct2,0];
di = diag(sincfct2);
dm_probsincfct2 = (a * di)'; 

% dm_actuators_mat0 = dm_probcosfct.*(dm_probsincfct1.*dm_probsincfct2) ;

hcstt_UpdateMultiDM(+dm_probcosfct(:)*poke)
% hcstt_UpdateMultiDM(+dm_actuators_mat0(:)*poke)

im_cam = hcstt_TakeCamImage(true,false,tint);

figure(200)
imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
title('Image Camera EFC solution applied')
axis image 
export_fig([outDir,'CAMImage_EFCSol_config',num2str(config),'.png'],'-r300');

diff_cam = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20)-im_cam0(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);
figure(300)
imagesc(diff_cam)
title('Difference: im-im0')
axis image
colorbar
export_fig([outDir,'CAMImDiff',num2str(config),'.png'],'-r300');
diff_neg = diff_cam*0;
diff_neg(find(diff_cam<0)) = 1;
figure(301)
imagesc(diff_neg)
axis image
title('Ones are negative')
export_fig([outDir,'OnesAreNEgative_config',num2str(config),'.png'],'-r300');


%%
% Model
%
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

[posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,1);
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
phase = angle(wf2_prob_noerrors);
phase0 = angle(wf2_noerrors);

figure(400)
% imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20)-int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20) - int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
axis image
title('Model - im-im0')
colorbar
export_fig([outDir,'ModelImDiff',num2str(config),'.png'],'-r300');
figure(401)
% imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20)-int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20))
axis image
title('Model - im')
colorbar
figure(4011)
imagesc(real(wf2_prob_noerrors(N/2-20:N/2+20,N/2-20:N/2+20)))
axis image
title('Model - real')
colorbar
figure(4012)
imagesc(imag(wf2_prob_noerrors(N/2-20:N/2+20,N/2-20:N/2+20)))
axis image
title('Model - imag')
colorbar

figure(402)
% imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20)-int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
imagesc(phase(N/2-20:N/2+20,N/2-20:N/2+20))
axis image
title('Phase')
colorbar
figure(402)
% imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20)-int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
imagesc(phase(N/2-20:N/2+20,N/2-20:N/2+20)-phase0(N/2-20:N/2+20,N/2-20:N/2+20))
% imagesc(phase0(N/2-20:N/2+20,N/2-20:N/2+20))
axis image
title('Phase Diff')
colorbar


figure(500)
        imagesc(xvals/apRad,yvals/apRad,surf_DM1*1e9);
        colormap(gray(256));hcb=colorbar;
        %caxis([0 10e-9]);
        axis image;%axis off;% 
        axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
        %text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
        hx = xlabel('{\itx} / {\itR}');
        hy = ylabel('{\ity} / {\itR}');
        title('DM1 surface height (nm)');
        set(gca,'FontSize', 10,...
                        'TickDir','out',...
                        'TickLength',[.02 .02]);
                     export_fig([outDir,'DMSurf_config',num2str(config),'.png'],'-r300');

% % imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20)-int_mod0(N/2-20:N/2+20,N/2-20:N/2+20))
% imagesc(int_mod(N/2-20:N/2+20,N/2-20:N/2+20))
% axis image
% title('Model - Sinusoid on DM')
hcstt_DisconnectDevices();

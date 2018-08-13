function [ x_hat ] =  EFSensing_RegularEFC_labTest(wf0,us_total,info)
%Performs the Electric Field extraction for pixel based EFC.


load('output\calibrateDM_Aug08');

apRad = info.apRad;
lambdaOverD = info.lambdaOverD; 

rng(3);


x_fib_pix = info.x_fib_pix;
% x_fib_pix_4interp = x_fib_pix-0.5; %matching the position of the actual sinusoid to the position of the DH
x_fib_pix_4interp = abs(x_fib_pix); %
actxc = interp1(distPix_meas,actxcDM_arr,x_fib_pix_4interp);
angDM = interp1(angPix_meas,angDM_arr,atan(info.y_fib_pix/x_fib_pix));

useGPU = info.useGPU;

RHO = info.RHO ;
THETA = info.THETA;
N = info.N;
lambda0 = info.lambda0 ;
useApodizer = info.useApodizer;
FPM = info.FPM;
LPM = info.LPM ;
outDir = info.outDir;
xvals = info.xvals;
yvals = info.yvals;
numOfWavelengths = info.numOfWavelengths;
Q = info.Q;
Q4G = info.Q4G;
num_Q = info.num_Q;

lam_arr = info.lam_arr ;

% Model WF with flat DM, WF0, we need this to compute Gu, since Gu:
%Gu = WF_DM - WF0
normPower = info.normPower;

Nact = 12;    

tint = info.tint;

ac_spac = info.ac_spac;%round(2*apRad/Nact);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );

posDM_x = info.posDM_x;
posDM_y = info.posDM_y;

x_cent_cam = info.x_cent_cam;
y_cent_cam = info.y_cent_cam;

wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;
background = info.background;

num_DM_shapes = 5;
ph_arr = linspace(0, 3/2*pi, num_DM_shapes);

H_regular_mat = zeros(num_DM_shapes,num_Q,2);
DeltaI_regular_arr = zeros(num_DM_shapes,num_Q);

p2v_dm = info.p2v_dm_sensing;
poke_amp = p2v_dm*1e-9;

for KK = 1:num_DM_shapes
%     ph = ph_arr(KK);
%     cosfct = cos(2*pi*[1:apRad2]/(apRad2) * ww + ph) ;
    a = ones(Nact,Nact);
%     di = diag(cosfct);
%     us = a * di; 
%     dm_probcosfct = us';
    dm_probcosfct = hcstt_DMMapSin(1, angDM, actxc, ph_arr(KK));
    sincfct1 = sinc([-Nact/2:Nact/2-1]/(Nact) *  3.5);
%     a = ones(apRad*2,apRad*2);
    di = diag(sincfct1);
    dm_probsincfct1 = a * di; 
%     [rows,cols] = size(dm_prob); 
%     dm_probsincfct1 = padarray_centered(dm_prob,rows,cols,N);
    
    sincfct2 = sinc([-Nact/2:Nact/2-1]/(Nact) * 3.5);
%     a = ones(apRad*2,apRad*2);
    di = diag(sincfct2);
    dm_probsincfct2 = (a * di)'; 
%     [rows,cols] = size(dm_prob); 
%     dm_probsincfct2 = padarray_centered(dm_prob,rows,cols,N);
    
%     dm_actuators_mat0 = dm_probcosfct.*(dm_probsincfct1.*dm_probsincfct2) * poke_amp;
    dm_actuators_mat0 = dm_probcosfct * poke_amp;
    
    %Add the actuators heights that were already set on the DM:
    dm_actuators_mat = dm_actuators_mat0(:) + us_total;
    
    count = 0; 
    DM1_strokes = zeros(N,N);
    for ix = 1:Nact
        for iy = 1:Nact
            count = count + 1;
            xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
            DM1_strokes(xpos,ypos) = dm_actuators_mat(count) + DM1_strokes(xpos,ypos);
        end
    end
    surf_DM1 = conv2(DM1_strokes,infl,'same');
    wf2_prob_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM1, true, info);
    wf2_prob_noerrors = wf2_prob_noerrors * sqrt(normPower);

    %Gu is the effect of the DM on the image plane
    Gu = wf2_prob_noerrors-wf0;
    Gu_re = real(Gu);
    Gu_im = imag(Gu);
    
    % Measure the intensity  of the fiber for the positive probe
    hcstt_UpdateMultiDM((+dm_actuators_mat0(:) + us_total)/1e-9)
    im_cam = zeros(400,400);
    for II=1:25
        im_camII = hcstt_TakeCamImage(true,false,tint)-background;
        im_cam = im_cam + im_camII/25;
    end
    int_regular_plus = im_cam(Q);
    sz = size(im_cam);

    figure(1111)
    im_camaux=im_cam;
    imagesc(im_camaux(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
    title(['Probing phase ',num2str(KK)])
    axis image
    colorbar
    drawnow;
    figure(1112)
    im_camaux=im_cam;
    im_camaux(Q) = 0;
    imagesc(im_camaux(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
    title(['Probing phase ',num2str(KK)])
    axis image
    colorbar
    drawnow;
    %
    
    % Measure the intensity for the positive probe
    hcstt_UpdateMultiDM((-dm_actuators_mat0(:) + us_total)/1e-9)
    im_cam = zeros(400,400);
    for II=1:25
        im_camII = hcstt_TakeCamImage(true,false,tint)-background;
        im_cam = im_cam + im_camII/25;
    end
    int_regular_minus = im_cam(Q);
    %
    
    %
    DeltaI_regular_arr(KK,:) = int_regular_plus - int_regular_minus;

    % Compute the ith element of the observation matrix H
%     H_mat(KK, :) = 4*[sum(sum(Gu_re.*fibermode0)),sum(sum(Gu_im.*fibermode0))];
    Gu_re2 = Gu_re(N/2-sz(1)/2:N/2+sz(1)/2-1,N/2-sz(2)/2:N/2+sz(2)/2-1);
    Gu_im2 = Gu_im(N/2-sz(1)/2:N/2+sz(1)/2-1,N/2-sz(2)/2:N/2+sz(2)/2-1);
    H_regular_mat(KK, :, :) = 4*[Gu_re2(Q4G),Gu_im2(Q4G)];
    %     
end
x_regular_hat = zeros(2,num_Q);
H = zeros(num_DM_shapes,2);
for II = 1 : num_Q
    H(:,:) = H_regular_mat(:,II,:);
    x_regular_hat(:,II) =  pinv(H)*DeltaI_regular_arr(:,II);
end
x_hat = x_regular_hat;
end

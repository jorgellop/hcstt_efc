function [ x_regular_hat, x_regular_true ] =  EFSensing_RegularEFC_ModelOnly(wf0,wf0_werrors,wfin_werrors,us_total,info)
%Performs the Electric Field extraction for a regular pixel based EFC
%
%All model based

N = info.N;

apRad = info.apRad;
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD = info.lambdaOverD; 

fiberDiam = info.fiberDiam; % Fiber diam. (lambda_0/D)

rng(3);
% error_map = fitsread('surfErrorMap_OX5.fits'); %nm

wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;


x_fib = info.x_fib; % Position of the fiber on the image plane. (lambda_0/D)
y_fib = info.y_fib;

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

lam_arr = info.lam_arr ;

% Model WF with flat DM, WF0, we need this to compute Gu, since Gu:
%Gu = WF_DM - WF0
normPower = info.normPower;

Nact = 12;    

fiberDiam_pix = fiberDiam*lambdaOverD;

[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD, Y);

fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
    exp(-(RHO_fib/(fiberDiam_pix/2)).^2);

ac_spac = round(2*apRad/Nact);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac );

posDM_x = info.posDM_x;
posDM_y = info.posDM_y;

Q = info.Q;
num_Q = info.num_Q;

num_DM_shapes = 4;
ph_arr = linspace(pi/2, 2*pi, num_DM_shapes);
H_mat = zeros(num_DM_shapes,2);
H_regular_mat = zeros(num_DM_shapes,num_Q,2);
DeltaI_arr = zeros(num_DM_shapes,1);
DeltaI_regular_arr = zeros(num_DM_shapes,num_Q);
ww = x_fib;
p2v_dm = info.p2v_dm_sensing;
Nact = 12;
poke_amp = p2v_dm*1e-9;
DM1_strokesKK = zeros(N,N);
for KK = 1:num_DM_shapes
    ph = ph_arr(KK);
    cosfct = cos(2*pi*[1:Nact]/(Nact) * (ww-0.0) + ph) ;
    a = ones(Nact,Nact);
    di = diag(cosfct);
    us = a * di; 
    
    dm_probcosfct = us';
        
    sincfct1 = sinc([1:Nact]/(Nact) * fiberDiam * 1.2);
%     a = ones(apRad*2,apRad*2);
    di = diag(sincfct1);
    dm_probsincfct1 = a * di; 
%     [rows,cols] = size(dm_prob); 
%     dm_probsincfct1 = padarray_centered(dm_prob,rows,cols,N);
    
    sincfct2 = sinc([1:Nact]/(Nact) * fiberDiam * 1.2);
%     a = ones(apRad*2,apRad*2);
    di = diag(sincfct2);
    dm_probsincfct2 = (a * di)'; 
%     [rows,cols] = size(dm_prob); 
%     dm_probsincfct2 = padarray_centered(dm_prob,rows,cols,N);
    
    dm_actuators_mat = dm_probcosfct.*(dm_probsincfct1.*dm_probsincfct2) * poke_amp;

    count = 0; 
    DM1_strokes = zeros(N,N);
    for ix = 1:Nact
        for iy = 1:Nact
            count = count + 1;
            xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
            ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
            DM1_strokesKK(xpos,ypos) = dm_actuators_mat(count) + DM1_strokes(xpos,ypos);
        end
    end
    surf_DM1 = conv2(DM1_strokesKK,infl,'same');
    wf2_prob_noerrors = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM1, true, info);
    wf2_prob_noerrors = wf2_prob_noerrors * sqrt(normPower);
    wf2_prob_werrors_plus = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_werrors, surf_DM1, true, info);
    wf2_prob_werrors_plus = wf2_prob_werrors_plus * sqrt(normPower);
    
    % Measure the intensity out of the fiber for the positive probe
    int_plus = abs(sum(sum(wf2_prob_werrors_plus.*fibermode0)))^2;  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
    int_regular_plus = abs(wf2_prob_werrors_plus(Q)).^2;  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
    
    wf2_prob_werrors_minus = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_werrors, -surf_DM1, true, info);
    wf2_prob_werrors_minus = wf2_prob_werrors_minus * sqrt(normPower);

    int_minus = abs(sum(sum(wf2_prob_werrors_minus.*fibermode0)))^2; % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
    int_regular_minus = abs(wf2_prob_werrors_minus(Q)).^2;  % dm_actuators_mat is a 12^2x1 array with the actuators heights in nm
    %
    
    %Gu is the effect of the DM on the image plane
    Gu = wf2_prob_noerrors-wf0;
    Gu_re = real(Gu);
    Gu_im = imag(Gu);

    %
    DeltaI_arr(KK,1) = int_plus - int_minus;
    DeltaI_regular_arr(KK,:) = int_regular_plus - int_regular_minus;

    % Compute the ith element of the observation matrix H
    H_mat(KK, :) = 4*[sum(sum(Gu_re.*fibermode0)),sum(sum(Gu_im.*fibermode0))];
    H_regular_mat(KK, :, :) = 4*[Gu_re(Q),Gu_im(Q)];
    %     
end

x_true = sum(sum(wf0_werrors.*fibermode0))';
x_regular_true = wf0_werrors(Q);
x_hat = pinv(H_mat)*DeltaI_arr;
% x_regular_hat = pinv(H_regular_mat)*DeltaI_regular_arr;

x_regular_hat = zeros(2,num_Q);
H = zeros(num_DM_shapes,2);
for II = 1 : num_Q
    H(:,:) = H_regular_mat(:,II,:);
    x_regular_hat(:,II) =  pinv(H)*DeltaI_regular_arr(:,II);
end
% 
% wf2_prob_true = prescription_DM1toImage_compact_vFiberCoupling( wf1_werrors, zeros(N,N), true, info);
% wf2_prob_true = wf2_prob_true / sqrt(normI);
% re_true = sum(sum(real(wf2_prob_true).*fibermode0));
% im_true = sum(sum(imag(wf2_prob_true).*fibermode0));
% x_true = [re_true, im_true]
% norm = max(max(abs(wf2_prob).^2))
end

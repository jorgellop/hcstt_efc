% EFC with fiber, code for lab demonstration
% Notes:
%   - Regularization parameter: which one to chose? This choice affects how
%   fast the dark hole is digged, and the damping of the actuators pulse
%   - Normalization of the model WF: the propagated model used to build the
%   G matrix and to do the sensing of the actual WF, has to be normalized
%   in the same way as the actual WF.
clear all;
close all;
addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(false);
hcstt_NewFlatForDM('SpeckleNulling_flat_Feb26');

N = 2^10;
apRad = 118;
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
lambdaOverD = N/apRad/2; % lambda/D (samples) 
tint = 2;

lambda0 = 650e-9;
numOfWavelengths = 1; % monochromatic, use 1
percentBW = 10; % percent bandwidth=(Delta lambda)/lambda*100
BW = percentBW/100; % bandwidth 
if(numOfWavelengths > 1)
    lam_fracs = linspace(1-BW/2,1+BW/2,numOfWavelengths);
else
    lam_fracs = 1;
end
lam_arr = lambda0*lam_fracs;

fiberDiam = 2*0.71; % Fiber diam. (lambda_0/D)
% These are parameters used in the propagation routines
use_fiber = false;
normal_EFC = true;
    
label = '_Feb26v2';
outDir = ['output',filesep,'EFC_RegularEFC_LabDemonstration',label,filesep];
mkdir(outDir);

% Position of the fiber on the image plane. (lambda_0/D)
% x_fib = 3;
% y_fib = 0;

info.fiberDiam = fiberDiam;
info.useGPU = false;
info.apRad = apRad;
info.lambdaOverD = lambdaOverD;
info.RHO = RHO;
info.THETA = THETA;
info.N = N;
info.lambda0 = lambda0;
info.lam_arr = lam_arr;
info.numOfWavelengths = numOfWavelengths;
info.useApodizer = false;
info.useGPU = false; 
info.FPM = exp(1i*8*THETA);
info.outDir = outDir;
info.xvals = xvals;
info.yvals = yvals;
info.use_fiber = use_fiber;
info.normal_EFC = normal_EFC;
info.tint = tint;

info.LPM = ones(N,N);
wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
wfin_noerrors(RHO > apRad) = 0;
% wf2 = wf2*sqrt(normPower);

% Total Power to normilize model Gu and intensity from the fiber. Need of
% normalization factor
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'');
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'set',8.8e-7,218,204,206,204);
% x_fib_pix = +12;
normPower = 0.00015;
totalPower = 1;
x_fib=3.5;
y_fib=0;
x_fib_pix = 14;
info.x_fib_pix = x_fib_pix;
info.x_fib = x_fib;
info.y_fib = y_fib;
info.p2v_dm_sensing = 1;

info.normPower = normPower;

if normal_EFC
    x_cent_cam = 200;
    y_cent_cam = 200;
    Ncam = 400;
    [Xcam,Ycam] = meshgrid(-x_cent_cam:Ncam-x_cent_cam-1,-(Ncam-y_cent_cam):y_cent_cam-1); 
    q_pix = 1;
    info.q_pix = q_pix;
    Q = zeros(Ncam,Ncam);
    Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
    Q = and(Q, Ycam >= -(q_pix) );
    Q = and(Q, Ycam <= (q_pix) );
    num_Q = numel(find(Q));
    info.num_Q = num_Q;
    info.Q = Q;
end

% Model of the fiber mode shape
[THETA_fib,RHO_fib] = cart2pol(X - x_fib * lambdaOverD ,Y);
fiberDiam_pixII = (fiberDiam*lambdaOverD);
fibermode0 = sqrt(2/(pi*(fiberDiam_pixII/2)^2))* ...
        exp(-(RHO_fib/(fiberDiam_pixII/2)).^2);
info.fibermode0 = fibermode0;

info.LPM = exp(-(RHO/(0.9*apRad)).^1000);

Nact = 12;  % Number of DM actuators, Nact^2
[posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuators(N,apRad);
info.posDM_x = posDM_x;
info.posDM_y = posDM_y;

% Load the DM influence functions
% ac_spac = round(2*apRad/Nact);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac ); % Influence function. Need of a model for actual DM
info.ac_spac = ac_spac;
info.infl = infl;

% Initialize DMs, etc. 
surf_DM10 = zeros(N); % Intialize the DM surface to flat
DM1_strokes = zeros(N); % Intialize the DM strokes to flat
us = zeros(Nact^2,1); % Initialize the fractional stroke changes to zero
us_total = zeros(Nact^2,1); % Initialize the fractional stroke changes to zero
poke_amp = 1e-9; % Initialize the poke amplitude
    
maxits = 33; % maximum number of EFC iterations allowed
Gcount = 0; % counter for number of times G matrix was used
Gcountmin = 10; % Minimum number of times to use a G matrix
curr_coupl_SMF = 1;
curr_int = 1e9;
recalc_G = true; % Initialize the flag to re-calculate the G matrix
% regvals = logspace(-6,-1,6); % Range of regularization parameters to test
regval = nan; % Initial regularization value
int_in_DH = []; % Array to keep track of dark hole irradiance 
int_est_in_DH = []; % Array to keep track of dark hole irradiance 

%Take camera image with flat DM
hcstt_UpdateMultiDM(zeros(Nact^2,1))
im_cam = hcstt_TakeCamImage(true,false,tint);
sz_imcam = size(im_cam);
info.sz_imcam = sz_imcam;
hcstt_test_plotCamImage(im_cam(sz_imcam(1)/2-20:sz_imcam(1)/2+20,sz_imcam(1)/2-20:sz_imcam(1)/2+20), [outDir,'CamImage_flatDM'], [41,41] );

% Run EFC iterations 
for k = 1:maxits
    fprintf('Iteration: %d ',k);
        
    % Update actuator height with the LMS solution, us
%     us_total = us_total + us;
    us_total = us_total + us;
    
    % Build DM surface from stroke amplitudes 
    count = 0; 
%     DM1_strokes = zeros(N);
    for ix = 1:Nact
        xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
        for iy = 1:Nact
            count = count + 1;
            ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
            DM1_strokes(xpos,ypos) = us(count)*poke_amp + DM1_strokes(xpos,ypos);
        end
    end
    surf_DM10 = conv2(DM1_strokes,infl,'same');

    %Make of plot of the current DM surface 
    fig0 = figure('visible','off','color','none');
    imagesc(xvals/apRad,yvals/apRad,surf_DM10*1e9);
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
    set(gca,'YDir','normal');
    set(fig0,'units', 'inches', 'Position', [0 0 5 5])
    export_fig([info.outDir,'DM1surf_',num2str(k),'.png']);
    close(fig0);    
        
    % Simualte WF in image plane with current DM shape, this is needed for the WF sensing 
    wf2_current = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM10, true, info);
    wf2_current = wf2_current * sqrt(normPower);
    
    % Perform WF sensing
    if normal_EFC
        Eab =  EFSensing_RegularEFC_labTest(wf2_current,us_total*poke_amp,info);
    else
        Eab =  EFSensing_EFCwFiber_labTest(wf2_current,us_total*poke_amp,info);
    end
    %
    
    % Build WF vector
    if normal_EFC
        Eabreg = [Eab(1,:)';Eab(2,:)'; zeros(Nact^2,1)];
    else
        Eabreg = [Eab(1);Eab(2); zeros(Nact^2,1)];
    end

    % Check progress 
    
    % Update intensity over fiber
%     curr_int_measured = hcstt_GetIntensityFIU(us_total*poke_amp/1e-9,100);  % WATCH OUT: us is a 12^2x1 array with the actuators heights in nm
%     prev_irr = curr_irr;
    curr_int_est = sum(sum(abs(Eab).^2))/numel(Eab(1,:))/numel(lam_fracs);

    if ~normal_EFC
        prev_coupl_SMF = curr_coupl_SMF;
        curr_coupl_SMF = sum(abs(Eab).^2/totalPower)/numel(lam_fracs);
        fprintf(fprintf('Coupl SMF: %g ',curr_coupl_SMF));fprintf('Reg. Val.: %g ',regval);
        coupl_SMF_in_DH = [coupl_SMF_in_DH, curr_coupl_SMF];
        prev_int = prev_coupl_SMF;
        curr_int = curr_coupl_SMF;
    else
        figure(100)
        hcstt_UpdateMultiDM(+us_total)
        im_cam = hcstt_TakeCamImage(true,false,tint);
        imagesc(im_cam(sz_imcam(1)/2-20:sz_imcam(1)/2+20,sz_imcam(1)/2-20:sz_imcam(1)/2+20))
        axis image
        colorbar
        drawnow;
        
        prev_int = curr_int;
        curr_int = mean(im_cam(Q));
    end    
    int_in_DH = [int_in_DH, curr_int];
    int_est_in_DH = [int_est_in_DH, curr_int_est];
    
    fprintf(['Mean MeasIntensity in DH: ', num2str(curr_int)])
    fprintf(['Mean EstIntensity in DH: ', num2str(curr_int_est)])
    
    % Determine whether to continue and/or calculate a new G matrix
    if(prev_int<curr_int)
%         return;
    elseif(Gcount > Gcountmin)
        recalc_G = true;
    end

    % Sets the new DM poke amplitude based on current dark hole irr. 
    poke_amp = 1e-9;%sqrt(curr_irr)*lambda0;
        
    % calculate the G matrix, if needed
    if(recalc_G)
        Gcount = 0;
        disp('Calculating the G matrix for EFC.')
        if k==1
                load('G_it1_Feb22')
        else
            G = calculateGmatrix_vFiberCouplingOneFibxel_broadband( wfin_noerrors, surf_DM10, Q, Nact, ac_spac, poke_amp, infl, lambda0, N , info);
            save(['output\G.mat'],'G');
        end
        recalc_G = false;
    end
    Gcount = Gcount + 1; % Count the times the current G matrix has bee used
            
    Gsplit = [real(G);imag(G)]; % Splits the G-matrix into real and imaginary parts 
    
    % Regalarization value, to be updated each iteration?
    numreg = 11;
    regval_arr = logspace(1e-4,1e1,numreg);
    curr_reg_arr = zeros(1,numreg) + nan;
    for II=1:numreg
    %     regval = 0.1; % To be determined
        regval = regval_arr(II);
        %

        Greg = [Gsplit;regval*eye(Nact^2)]; 

        % Compute the  new actuator heights
        usII = -1*pinv(Greg)*Eabreg; %  
        %
        if max(usII)<20
            hcstt_UpdateMultiDM(+(usII+us_total)*poke_amp/1e-9)
            im_cam = hcstt_TakeCamImage(true,false,tint);
            curr_reg_arr(II) = mean(im_cam(Q));
        end
    end
    [mi,ind_min] = min(curr_reg_arr);
    regval = regval_arr(ind_min);
    
    save([info.outDir,'DMshapes.mat'],'surf_DM10');
    save([info.outDir,'DM1_strokes.mat'],'DM1_strokes');
    
    Greg = [Gsplit;regval*eye(Nact^2)]; % Final G matrix to be used
    % Compute the  new actuator heights
    us = -1*pinv(Greg)*Eabreg; %  
    %
    
    %Check if everything is OK
    us_max = max(us*poke_amp*1e9)
%     prompt = 'Continue (0 to end)? ';
%     cont = input(prompt);
%     if cont==0 
%         return 
%     end
end
    
fig0 = figure(2);
plot(1:k,int_in_DH)
xlabel('Iteration')
ylabel('Mean Intensity in DH')
title(['EFC 3x3 pixels - Mean Intensity vs it'])
% legend('Coupling SMF','Coupling MMF');
export_fig([outDir,'MeanInt_vs_it',label,'.png'],'-r300');

fig0 = figure(3);
plot(1:k,int_est_in_DH)
xlabel('Iteration')
ylabel('Mean EST Intensity in DH')
title(['EFC 3x3 pixels - Mean EST Intensity vs it'])
% legend('Coupling SMF','Coupling MMF');
export_fig([outDir,'MeanESTInt_vs_it',label,'.png'],'-r300');

hcstt_UpdateMultiDM(+us_total)
fig0 = figure(3);
im_cam = hcstt_TakeCamImage(true,false,tint);
hcstt_test_plotCamImage(im_cam(sz_imcam(1)/2-20:sz_imcam(1)/2+20,sz_imcam(1)/2-20:sz_imcam(1)/2+20), [outDir,'CamImage_final'], [41,41] );
hcstt_DisconnectDevices();
% Save actual data
% save([info.outDir,'coupl_SMF_in_DH.mat'],'coupl_SMF_in_DH');

% fig0 = figure(3);
% plot(1:k,log10(Spl_pix_arr))
% hold on
% plot(1:k,log10(Spl_fib_arr))
% hold on
% plot(1:k,log10(S_pix_arr))
% hold on
% plot(1:k,log10(S_fib_arr))
% hold off
% xlabel('Iteration')
% ylabel('log(Throughput)')
% title(['Throughput vs it'])
% legend('Throughput planet MMF','Throughput planet SMF','Throughput speckles MMF','Throughput speckles SMF');
% export_fig([outDir,'FigThroughput_vs_it',label,'.png'],'-r300');
% 
% fig0 = figure(4);
% plot(1:k,log10(Spl_pix_arr))
% hold on
% plot(1:k,log10(Spl_fib_arr))
% hold off
% xlabel('Iteration')
% ylabel('log(Throughput)')
% title(['Throughput vs it'])
% legend('Throughput planet MMF','Throughput planet SMF');
% export_fig([outDir,'FigThroughputPlanet_vs_it',label,'.png'],'-r300');


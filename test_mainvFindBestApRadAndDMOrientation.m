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

use_fiber = false;
normal_EFC = true;
debug = true;

label = '_Aug12';
outDir = ['output',filesep,'EFC_wFiber_LabDemonstration',label,filesep];
mkdir(outDir);

load('output\calibrateDM_Aug01'); % actxcDM, angDm vs pix on camera info

% Initialize all devices
if use_fiber
    hcstt_Initialize(true);
else
    hcstt_Initialize(false);
end

% load('benchModelPixelPitchMatching_scaleR_apRad_May7')
N = 2^10;

[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);
[THETA,RHO] = cart2pol(X,Y);
tint = 25;


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

% These are parameters used in the propagation routines  
info.useGPU = false;
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


% Total Power to normilize model Gu and intensity from the fiber. Need of
% normalization factor
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'');
% [normPower,totalPower,x_fib,y_fib] = hcstt_NormalizationFactor(wfin_noerrors,info,'set',8.8e-7,218,204,206,204);
load('BenchModelNormalization_0803')
normPower = normPower_normalization*tint/tint_normalization;%0.00055;%
peakInt = peakInt_normalization*tint/tint_normalization;
totalPower = 1.16e-6;
if use_fiber
    normPower = totalPower/3.8185e10;
    peakInt = totalPower;
end
x_fib=3.5;
y_fib=0;
x_fib_pix = -11;
info.x_fib_pix = x_fib_pix;
info.y_fib_pix = 0;
info.x_fib = x_fib;
info.y_fib = y_fib;
if normal_EFC
    info.p2v_dm_sensing = 2;
else
    info.p2v_dm_sensing = 10;
end

info.normPower = normPower;
% Take background image
take_background = false;
Ncam = 400;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
    im_cam = zeros(400,400);
    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    background = im_cam;
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    background = zeros(Ncam,Ncam);
end
info.background = background;

% Update shape of DM with SN result
hcstt_NewFlatForDM('ImageSharpeningModel_0801_flatv2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% hcstt_NewFlatForDM('NewFlat_SpeckleOnFiber_May10');
% hcstt_NewFlatForDM('output\NewFlat_testPReviousEFCRun');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

if normal_EFC    
    % Find Center of camera image
    im_cam = zeros(400,400);
    tint_findCenter = 0.3;
    for II=1:25
        im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
        im_cam = im_cam + im_camII/25;
        pause(0.1)
    end
    im_camaux = im_cam;
    im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
    [ma,ind_ma] = max(im_camaux(:));
    [x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
    if max(im_camII(:))>240
        disp('Find Center image saturated')
        return
    end

    % Save image of flat DM without SN
    im_cam = hcstt_TakeCamImage(true,false,tint)-background;
    sz_imcam = size(im_cam);
    info.sz_imcam = sz_imcam;
    hcstt_test_plotCamImage(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20), [outDir,'CamImage_flatDM'], [41,41] );
    
    
    [Xcam,Ycam] = meshgrid(-y_cent_cam+1:Ncam-y_cent_cam,-x_cent_cam+1:Ncam-x_cent_cam); 
    [THETAcam, RHOcam] = cart2pol(Xcam,Ycam);
    q_pix = 1;
    info.q_pix = q_pix;
    Q = zeros(Ncam,Ncam);
    Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
    Q = and(Q, Ycam >= -(q_pix) );
    Q = and(Q, Ycam <= (q_pix) );
    IWA_pix = (x_fib_pix-q_pix );
    OWA_pix = (x_fib_pix+q_pix );
    angleDH = pi/12;
%     Q = and(RHOcam >= IWA_pix, RHOcam <= OWA_pix);
%     Q = and(Q, Xcam > 0);
%     Q = and(Q, abs(THETAcam) < angleDH); % 60deg keystone about x-axis

    
    Q4G = zeros(Ncam,Ncam);
%     [Xcam4G,Ycam4G] = meshgrid(-Ncam/2:Ncam/2-1); 
    [Xcam4G,Ycam4G] = meshgrid(-Ncam/2-1:Ncam/2-2); 
%     [THETAcam4G, RHOcam4G] = cart2pol(Xcam4G,Ycam4G);
    Q4G = zeros(Ncam,Ncam);
    Q4G = and(Xcam4G >= (x_fib_pix-q_pix ), Xcam4G <=  (x_fib_pix+q_pix ));
    Q4G = and(Q4G, Ycam4G >= -(q_pix) );
    Q4G = and(Q4G, Ycam4G <= (q_pix) ); % 60deg keystone about x-axis
%     Q4G = and(RHOcam4G >= IWA_pix, RHOcam4G <= OWA_pix);
%     Q4G = and(Q4G, Xcam4G > 0);
%     Q4G = and(Q4G, abs(THETAcam4G) < angleDH); % 60deg keystone about x-axis

    num_Q = numel(find(Q));
    info.num_Q = num_Q;
    
    info.Q4G = Q4G;
    
else
    x_cent_cam = 200;
    y_cent_cam = 200;
end



counttot = 0;
Nact = 12;  % Number of DM actuators, Nact^2
tic
% num_apRad = numel(apRad_arr);

% info.LPM = exp(-(RHO/(0.85*apRad)).^1000);

for apRII=1:1
    apRad = 142;%apRad_arr(apRII);
%     scaleR = scaleR_arr(apRII);
    lambdaOverD = N/apRad/2; % lambda/D (samples) 

    fiberDiam = 2*0.71; % Fiber diam. (lambda_0/D)
    fiberDiam_pix = (fiberDiam*lambdaOverD);

    info.fiberDiam = fiberDiam;
    info.apRad = apRad;
    info.lambdaOverD = lambdaOverD;

    info.LPM = exp(-(RHO/(0.85*apRad)).^1000);
    
    % Wavefront at entrance for model propagation
    wfin_noerrors = complex(ones(N, N), zeros(N, N)) ;
    wfin_noerrors(RHO > apRad) = 0;

    for posII=1:8


        [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,posII);
        info.posDM_x = posDM_x;
        info.posDM_y = posDM_y;
        info.ac_spac = ac_spac;

        % Load the DM influence functions
        % ac_spac = round(2*apRad/Nact);
        infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac ); % Influence function. Need of a model for actual DM
        info.infl = infl;

        % Initialize DMs, etc. 
        surf_DM10 = zeros(N); % Intialize the DM surface to flat
        DM1_strokes = zeros(N); % Intialize the DM strokes to flat
        us = zeros(Nact^2,1); % Initialize the fractional stroke changes to zero
        us_total = zeros(Nact^2,1); % Initialize the fractional stroke changes to zero
        poke_amp = 1e-9; % Initialize the poke amplitude

        maxits = 30; % maximum number of EFC iterations allowed
        Gcount = 0; % counter for number of times G matrix was used
        Gcountmin = 10; % Minimum number of times to use a G matrix
        curr_coupl_SMF = 1;
        curr_int = 1e9;
        recalc_G = true; % Initialize the flag to re-calculate the G matrix
        % regvals = logspace(-6,-1,6); % Range of regularization parameters to test
        regval = nan; % Initial regularization value
        int_in_DH = []; % Array to keep track of dark hole irradiance 
        int_est_in_DH = []; % Array to keep track of dark hole irradiance 
        coupl_SMF_in_DH = [];

        %Take camera image with flat DM
        hcstt_UpdateMultiDM(zeros(Nact^2,1))
        im_cam = hcstt_TakeCamImage(true,false,tint)-background;
        sz_imcam = size(im_cam);
        info.sz_imcam = sz_imcam;
    %     hcstt_test_plotCamImage(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20), [outDir,'CamImage_initial_DMconfig',num2str(posII)], [41,41] );
        [Xcam,Ycam] = meshgrid(-y_cent_cam+1:Ncam-y_cent_cam,-x_cent_cam+1:Ncam-x_cent_cam); 
        [THETAcam, RHOcam] = cart2pol(Xcam,Ycam);
        Q = zeros(Ncam,Ncam);
        Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
        Q = and(Q, Ycam >= -(q_pix) );
        Q = and(Q, Ycam <= (q_pix) );
        info.Q = Q;
        info.x_cent_cam = x_cent_cam;
        info.y_cent_cam = y_cent_cam;

        % Run EFC iterations 
        for k = 1:maxits
            counttot = counttot+1;
            %%%%%%%%%%%%%%
            %Find Center again just in case it moved...
            if normal_EFC
%                 tint_findCenter = 1.1;
%                 im_cam = zeros(Ncam,Ncam);
%                 for II=1:25
%                     im_camII = hcstt_TakeCamImage(true,false,tint_findCenter);
%                     im_cam = im_cam + im_camII/25;
%                     pause(0.1)
%                 end
%                 im_camaux = im_cam;
%                 im_camaux(190:210,190:210) = im_camaux(190:210,190:210)*1000;
%                 [ma,ind_ma] = max(im_camaux(:));
%                 [x_cent_cam,y_cent_cam] = ind2sub(size(im_camaux),ind_ma);
% 
%                 % Recalculate Q with the new center:
%                 [Xcam,Ycam] = meshgrid(-y_cent_cam+1:Ncam-y_cent_cam,-x_cent_cam+1:Ncam-x_cent_cam); 
%                 [THETAcam, RHOcam] = cart2pol(Xcam,Ycam);
%                 Q = zeros(Ncam,Ncam);
%                 Q = and(Xcam >= (x_fib_pix-q_pix ), Xcam <=  (x_fib_pix+q_pix ));
%                 Q = and(Q, Ycam >= -(q_pix) );
%                 Q = and(Q, Ycam <= (q_pix) );
%                 info.Q = Q;
            else
                %Recalculate position of fiber
%                 [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(2.4,0);
                actxc_fib = 2.4;
                ang_fib = 0;
                info.actxc_fib = actxc_fib;
                info.ang_fib = ang_fib;
                %Save to check drift of optical system
                actxc_fib_check_arr(counttot) = actxc_fib;
                ang_fib_check_arr(counttot) = ang_fib;
                %Recalculate position in pixels
                if actxc_fib>max(actxcDM_arr) 
                    r_fib_pix = interp1(actxcDM_arr,distPix_meas,actxc_fib,'linear','extrap');
                elseif actxc_fib<min(actxcDM_arr)
                    disp('Fiber too far away from star!')
                    return;
                else
                    r_fib_pix = interp1(actxcDM_arr,distPix_meas,actxc_fib);
                end
                x_fib_pix = r_fib_pix*cos(ang_fib);
                y_fib_pix = r_fib_pix*sin(ang_fib);
                % Model of the fiber mode shape
                [THETA_fib,RHO_fib] = cart2pol(X - x_fib_pix ,Y - y_fib_pix);
                fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
                        exp(-(RHO_fib/(fiberDiam_pix/2)).^2);
                info.fibermode0 = fibermode0;

            end
    %         IWA_pix = (x_fib_pix-q_pix );
    %         OWA_pix = (x_fib_pix+q_pix );
    %         angleDH = pi/12;
    %         Q = and(RHOcam >= IWA_pix, RHOcam <= OWA_pix);
    %         Q = and(Q, Xcam > 0);
    %         Q = and(Q, abs(THETAcam) < angleDH); % 60deg keystone about x-axis
            %%%%%%%%%%%%%%%%%%

            fprintf('Iteration: %d ',k);
            fprintf('\n')
            
            % Update actuator height with the LMS solution, us
            us_total = us_total + us;

            % Build DM surface from stroke amplitudes 
            count = 0; 
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
            figure(100);
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
    %         set(gca,'YDir','normal');
    %         set(fig0,'units', 'inches', 'Position', [0 0 5 5])
            if debug
                export_fig([info.outDir,'DM1surf_',num2str(k),'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png']);
            else
                if k == maxits
                    export_fig([info.outDir,'DM1surf_final_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png']);
                end
            end
    %         close(fig0);    

            % Simualte WF in image plane with current DM shape, this is needed for the WF sensing 
%             if k>5
%                 wf2_current = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM10, true, info);
%                 wf2_current = wf2_current * sqrt(normPower);
%             end
            if k==1
                wf2_current = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM10, true, info);
                wf2_current = wf2_current * sqrt(normPower);
            end
    %         if normal_EFC
            if k==1
                immod_flat = abs(wf2_current).^2;
                figure(5)
                immod_flat = immod_flat(N/2-sz_imcam(1)/2:N/2+sz_imcam(1)/2-1,N/2-sz_imcam(2)/2:N/2+sz_imcam(2)/2-1);
%                 immod_flat(Q4G) = max(immod_flat(:));
                imagesc(immod_flat(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
                axis image
                title('Im Cam flat DM')
            end
    %         end


            % Perform WF sensing
            if normal_EFC
                Eab =  EFSensing_RegularEFC_labTest(wf2_current,us_total*poke_amp,info);
            else
                Eab =  EFSensing_EFCwFiber_labTestv2(wf2_current,us_total*poke_amp,info);
%                 if debug
%                     % Check how the sesnsing looks like with flat DM at each
%                     % iteration
%                     wf2_flatDM = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, surf_DM10*0.0, true, info);
%                     wf2_flatDM = wf2_flatDM * sqrt(normPower);
%                     Eab_check =  EFSensing_EFCwFiber_labTest(wf2_flatDM,us_total*poke_amp*0.0,info);
%                     Eab1_fib_check_arr(counttot) = Eab_check(1);
%                     Eab2_fib_check_arr(counttot) = Eab_check(2);
%                     elapsedTime_arr(counttot) = toc;
%                 end
            end
            %

            % Build WF vector
            if normal_EFC
                Eabreg = [Eab(1,:)';Eab(2,:)'; zeros(Nact^2,1)];
            else
                Eabreg = [Eab(1);Eab(2); zeros(Nact^2,1)];                
            end        
            %

            % Check progress
            if ~normal_EFC
                curr_int_est =  sum(sum(abs(Eab).^2))/numel(lam_fracs);
                prev_coupl_SMF = curr_coupl_SMF;
                
%                 us_total2 = vec2mat(us_total,12);
%                 us_total2 = us_total2';

                curr_coupl_SMF = hcstt_GetIntensityFIU(+(us_total)*poke_amp/1e-9,10 );%sum(abs(Eab).^2/totalPower)/numel(lam_fracs);
                coupl_SMF_in_DH = [coupl_SMF_in_DH, curr_coupl_SMF];
                prev_int = prev_coupl_SMF;
                curr_int = curr_coupl_SMF;
                intaux = abs(wf2_current).^2;
                intaux = intaux(N/2-sz_imcam(1)/2:N/2+sz_imcam(1)/2-1,N/2-sz_imcam(2)/2:N/2+sz_imcam(2)/2-1);
                figure(3)
                imagesc(intaux(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
                axis image

            else
                curr_int_est = sum(sum(abs(Eab).^2))/numel(Eab(1,:))/numel(lam_fracs);

                figure(200)
                hcstt_UpdateMultiDM(+us_total)
                im_cam = hcstt_TakeCamImage(true,false,tint)-background;
                imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
                axis image
                colorbar
                drawnow;

                prev_int = curr_int;
                curr_int = mean(im_cam(Q));

                % Plot estimated intensity with measured intensity
                figure(2)
                int = abs(Eab(1,:)).^2+abs(Eab(2,:)).^2;
                plot(1:num_Q,int)
                hold on
                plot(1:num_Q,im_cam(Q))
                hold off
                legend('Estimated','Camera')

                intaux = abs(wf2_current).^2;
                intaux = intaux(N/2-sz_imcam(1)/2:N/2+sz_imcam(1)/2-1,N/2-sz_imcam(2)/2:N/2+sz_imcam(2)/2-1);
                intaux(Q4G) = max(intaux(:));
                figure(3)
                imagesc(intaux(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
                axis image
                im_camaux = im_cam;
                im_camaux(Q) = nan;
                figure(4)
                imagesc(im_camaux(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
                axis image
            end    
            int_in_DH = [int_in_DH, curr_int];
            int_est_in_DH = [int_est_in_DH, curr_int_est];

            fprintf([' Mean MeasIntensity in DH: ', num2str(curr_int)])
            fprintf('\n')
            fprintf([' Mean EstIntensity in DH: ', num2str(curr_int_est)])
            fprintf('\n')

            % Determine whether to continue and/or calculate a new G matrix
            if(prev_int<curr_int)
        %         return;
            elseif(Gcount > Gcountmin)
                recalc_G = true;
            end

            % Sets the new DM poke amplitude based on current dark hole irr. 
%             poke_amp = 1e-9;%sqrt(curr_irr)*lambda0;

            % calculate the G matrix, if needed
            if(recalc_G)
                Gcount = 0;
                disp('Calculating the G matrix for EFC.')
%                 if k==1
%                     load('output\G.mat')
%                 else
                if normal_EFC
                    G = calculateGmatrix_vFiberCouplingOneFibxel_broadband( wfin_noerrors, surf_DM10, Q4G, Nact, ac_spac, poke_amp, infl, lambda0, N , info);
                else
                    G = calculateGmatrix_vFiberCouplingOneFibxel_broadband( wfin_noerrors, surf_DM10, [], Nact, ac_spac, poke_amp, infl, lambda0, N , info);
                end
                save(['output\G.mat'],'G');
%                 end
                recalc_G = false;
            end
            Gcount = Gcount + 1; % Count the times the current G matrix has bee used

            Gsplit = [real(G);imag(G)]; % Splits the G-matrix into real and imaginary parts 

            % Regalarization value, to be updated each iteration?
            numreg = 7;
            numgain = 5;
            regval_arr = logspace(-8,-1,numreg);
            gain_arr = linspace(0.1,4,numgain);
            curr_reg_arr = zeros(1,numreg*numgain) + nan;
            countreg = 1;
            for III=1:numgain
                for II=1:numreg
                %     regval = 0.1; % To be determined
                    regval = regval_arr(II);
                    %

                    Greg = [Gsplit;regval*eye(Nact^2)]; 

                    % Compute the  new actuator heights
                    usII = -1*pinv(Greg)*Eabreg; % 
                    usII = usII*gain_arr(III);
                    %
                    if max(usII)<20
    %                     us_total2 = vec2mat(usII+us_total,12);
    %                     us_total2 = us_total2';

                        if normal_EFC
                            hcstt_UpdateMultiDM(+(us_total+usII)*poke_amp/1e-9)
                            im_cam = hcstt_TakeCamImage(true,false,tint)-background;
                            curr_reg_arr(countreg) = mean(im_cam(Q));
                        else
                            curr_reg_arr(countreg) = hcstt_GetIntensityFIU(+(us_total+usII)*poke_amp/1e-9,10 );
                        end
                    end
                    countreg = countreg + 1;
                end
            end
            [mi,ind_min] = min(curr_reg_arr);
            [ind_minreg,ind_mingain] = ind2sub([numreg,numgain],ind_min);
            regval = regval_arr(ind_minreg);
            gainval = gain_arr(ind_mingain);
            
            regvalfin_arr(k) = regval;
            gainvalfin_arr(k) = gainval;
            
            save([info.outDir,'DMshapes.mat'],'surf_DM10');
            save([info.outDir,'DM1_strokes.mat'],'DM1_strokes');

            Greg = [Gsplit;regval*eye(Nact^2)]; % Final G matrix to be used

            % Compute the  new actuator heights
            us = -1*pinv(Greg)*Eabreg; %  
            us = us*gainval;
            
            if debug
                if k==1
                    Greg = [Gsplit;zeros(Nact^2)];
                    us0 = -1*pinv(Greg)*Eabreg; %
                    save([info.outDir,'us0_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.mat'],'us0');
                end
            end
            %

            %Check if everything is OK
            us_max = max(us*poke_amp*1e9);
            us_max_arr(counttot) = us_max;
            fprintf([' Max us: ', num2str(us_max),'nm'])
            fprintf('\n')
            fprintf([' Regularization value: ', num2str(regval)])
            fprintf([' Gain value: ', num2str(gainval)])
            fprintf('\n')
            fprintf([' Current Suppression ', num2str(int_in_DH(1)/int_in_DH(k))])
            fprintf('\n')

        end
%% update one last time with new solution

        fprintf('Update with last solution');
        k = k+1;
        % Update actuator height with the LMS solution, us
        us_total = us_total + us;

        % Build DM surface from stroke amplitudes 
        count = 0; 
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
        figure(100);
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
%         set(gca,'YDir','normal');
%         set(fig0,'units', 'inches', 'Position', [0 0 5 5])
        if debug
            export_fig([info.outDir,'DM1surf_',num2str(k),'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png']);
        end
        
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
        %

        % Check progress
        if ~normal_EFC
            curr_int_est =  sum(sum(abs(Eab).^2))/numel(lam_fracs);
            prev_coupl_SMF = curr_coupl_SMF;
            curr_coupl_SMF = hcstt_GetIntensityFIU(+(us_total)*poke_amp/1e-9,10 );%sum(abs(Eab).^2/totalPower)/numel(lam_fracs);
            coupl_SMF_in_DH = [coupl_SMF_in_DH, curr_coupl_SMF];
            prev_int = prev_coupl_SMF;
            curr_int = curr_coupl_SMF;
            intaux = abs(wf2_current).^2;
            intaux = intaux(N/2-sz_imcam(1)/2:N/2+sz_imcam(1)/2-1,N/2-sz_imcam(2)/2:N/2+sz_imcam(2)/2-1);
            figure(3)
            imagesc(intaux(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
            axis image
        else
            curr_int_est = sum(sum(abs(Eab).^2))/numel(Eab(1,:))/numel(lam_fracs);

            figure(200)
            hcstt_UpdateMultiDM(+us_total)
            im_cam = hcstt_TakeCamImage(true,false,tint)-background;
            imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
            axis image
            colorbar
            drawnow;

            prev_int = curr_int;
            curr_int = mean(im_cam(Q));

            % Plot estimated intensity with measured intensity
            figure(2)
            int = abs(Eab(1,:)).^2+abs(Eab(2,:)).^2;
            plot(1:num_Q,int)
            hold on
            plot(1:num_Q,im_cam(Q))
            hold off
            legend('Estimated','Camera')

            intaux = abs(wf2_current).^2;
            intaux = intaux(N/2-sz_imcam(1)/2:N/2+sz_imcam(1)/2-1,N/2-sz_imcam(2)/2:N/2+sz_imcam(2)/2-1);
            intaux(Q4G) = max(intaux(:));
            figure(3)
            imagesc(intaux(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
            axis image
            im_camaux = im_cam;
            im_camaux(Q) = nan;
            figure(4)
            imagesc(im_camaux(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
            axis image
        end    
        int_in_DH = [int_in_DH, curr_int];
        int_est_in_DH = [int_est_in_DH, curr_int_est];

        fprintf(['Mean MeasIntensity in DH: ', num2str(curr_int)])
        fprintf([' Mean EstIntensity in DH: ', num2str(curr_int_est)])


    %         close(fig0);    
%% Save data from this EFC run
        
        fig0 = figure(2);
        plot(1:k,int_in_DH/peakInt)
        xlabel('Iteration')
        ylabel('Mean Intensity in DH')
        title(['EFC - Mean Intensity vs it (Suppression of ',num2str(int_in_DH(1)/int_in_DH(k)),')'])
    %     ylim([0.6e-4 1.2e-4])
        % legend('Coupling SMF','Coupling MMF');
        export_fig([outDir,'MeanInt_vs_it',label,'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png'],'-r300');
        close(fig0);

        fig0 = figure(3);
        plot(1:k,int_est_in_DH/peakInt)
        xlabel('Iteration')
        ylabel('Mean EST Intensity in DH')
        title(['EFC - Mean EST Intensity vs it'])
    %     ylim([0 1e-3])
        % legend('Coupling SMF','Coupling MMF');
        export_fig([outDir,'MeanESTInt_vs_it',label,'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png'],'-r300');
        close(fig0);

        if normal_EFC
            fig0 = figure(4);
            imagesc(intaux(Ncam/2-20:Ncam/2+20,Ncam/2-20:Ncam/2+20))
            axis image
            title(['Simulated image to see where Q falls DMconfig',num2str(posII),' apRad',num2str(apRad),])
            export_fig([outDir,'SimulatedImageWhereQFalls',label,'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png'],'-r300');
            close(fig0);
        end

        hcstt_UpdateMultiDM(+us_total)
        figure(6);
        im_cam = hcstt_TakeCamImage(true,false,tint)-background;
        imagesc(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20))
        axis image
        title(['Final Image DMconfig',num2str(posII),'apRad',num2str(apRad),])
        colorbar
        export_fig([outDir,'CamImageFinalImage',label,'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.png'],'-r300');
        im_cam_crop = im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20);

        save([info.outDir,'data_intvsit_dmshapes_',label,'_DMconfig',num2str(posII),'_apRad',num2str(apRad),'.mat'],'us_total','im_cam_crop','int_in_DH','peakInt','int_est_in_DH','regvalfin_arr','gainvalfin_arr');
    end
%     hcstt_test_plotCamImage(im_cam(x_cent_cam-20:x_cent_cam+20,y_cent_cam-20:y_cent_cam+20), [outDir,'CamImage_final','_DMconfig',num2str(posII)], [41,41] );
end
hcstt_DisconnectDevices();

if debug
    if use_fiber
        figure(1)
        plot(1:numel(actxc_fib_check_arr),actxc_fib_check_arr)
        title('Actxcyc sinusoid at each iteration, ie drift')
        xlabel('iteration')
        ylabel('act x cycle')
        export_fig([outDir,'ActxcycDMSinusoidVsIteration.png'],'-r300');
        figure(2)
        plot(1:numel(ang_fib_check_arr),ang_fib_check_arr)
        title('Angle DM sinusoid at each iteration, ie drift')
        xlabel('iteration')
        ylabel('angle')
        export_fig([outDir,'AngleDMSinusoidVsIteration.png'],'-r300');

        % Save us_total as flat:
    %     load('NewFlat_SpeckleOnFiber_May10')
    % %     load('output\NewFlat_testPReviousEFCRun');
        load('ImageSharpening_fminconIt2_Apr1')
        us_total_mat = vec2mat(us_total,Nact);
        us_total_mat = us_total_mat';
        flat_SN = us_total_mat + flat_SN;
        save('output\NewFlat_testPReviousEFCRun.mat','flat_SN')
    end
    
%     figure(3)
%     plot(elapsedTime_arr/60,Eab1_fib_check_arr)
%     hold on
%     plot(elapsedTime_arr/60,Eab2_fib_check_arr)
%     hold off
%     title('Estimated EF at each iteration')
%     xlabel('time (min)')
%     ylabel('Eab')
%     legend('Real{EF}','Imag{EF}')
%     export_fig([outDir,'EabVsIterationv2.png'],'-r300');
end
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


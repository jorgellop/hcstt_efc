%{
DM Writing Function: Write Sinusoid
- Places a Sin Function on Mirror Surface
- Sin calculated using DE_DMMapSin.m function
- Utilizes DE_DMArrayToVect to shape map for write
*** ASSUMES MIRROR CONNECTION ALREADY PRESENT; DOES NOT CLOSE CONNECTION
*** Defaults to DE_DMMapSin output setting 7: only returns heigh
        If desired, hnm can be plotted with imagesc. Refer to DE_DMMapSin
 
******************************************************
- Arguments:
    h0          = Max poke height in nm
    q           = angle of sinusoid
    x0          = actuators per cycle
    alp         = phase delay
    drv_info    = DM info from OPEN_mutliDM
- Returns:
    hnm         = Surface map as matrix in nm
    hV          = Vector of voltage percentages for writing
******************************************************

Compiled By:    Daniel Echeverri
Last Modified:  08/04/2016
%}

function [norm, totalPower,x_fib,y_fib,x_cent_cam,y_cent_cam,x_fib_pix] = ... 
    hcstt_NormalizationFactor_RegularEFC(wfin_noerrors,info,set_pos,totalPowerin,x_planet,y_planet,x_star,y_star)

global cam img Xcr Ycr CExp s drv_inf flat itr himg pm_scale

N = info.N;

[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA_fib,RHO_fib] = cart2pol(X,Y);
fiberDiam_pixII = (info.fiberDiam*info.lambdaOverD);
fibermode0 = sqrt(2/(pi*(fiberDiam_pixII/2)^2))* ...
        exp(-(RHO_fib/(fiberDiam_pixII/2)).^2);

wf2 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wfin_noerrors, zeros(N,N), false, info);
% totalPowerModel = abs(sum(sum(wf2(:,:,1).*fibermode0))).^2;
lambdaOverD_model = info.lambdaOverD;
pix_crop = 10 * lambdaOverD_model;
totalPowerModel = sum(sum(abs(wf2(N/2-pix_crop:N/2+pix_crop,N/2-pix_crop:N/2+pix_crop,1)).^2));

lambdaOverD_bench = 650e-9/(4e-3);
pix = 5.2e-6;
f=100e-3;
x_fib0 = abs(x_star-x_planet)*pix/f/lambdaOverD_bench;
y_fib0 = abs(y_star-y_planet)*pix/f/lambdaOverD_bench;



texp1 = 0.05;
im_cam = hcstt_TakeCamImage(true,true,texp1);
[ma,ind] = max(im_cam(:));
sz = size(im_cam);
[indx,indy] = ind2sub(sz,ind);
PeakPower = sum(sum(im_cam(indx-5:indx+5,indy-5:indy+5)));

texp2 = 4;
im_cam = hcstt_TakeCamImage(true,true,texp2);
sz = size(im_cam);
PeakPower_aux = sum(sum(im_cam(indx-5:indx+5,indy-5:indy+5)));
pix_crop = 10 * lambdaOverD_bench*f/pix;
totalPower = sum(sum(im_cam(sz(1)/2-pix_crop:sz(1)/2+pix_crop,sz(2)/2-pix_crop:sz(2)/2+pix_crop)))-PeakPower_aux...
    +PeakPower* sqrt(texp2/texp1);
%Set minimum exp time  
texp_final = 40;
cam.Timing.Exposure.Set(texp_final);
totalPower = totalPower*sqrt(texp_final/texp2);


disp('Flattening DM');

if(~strcmp('set',set_pos))
    

    %Create StringBuilder to retrieve and hold data; Set cap. to expected size
    sb = System.Text.StringBuilder();
    sb.Capacity = 4100;

    %Create Matrix to hold images: Width(200)xHeight(200)x#Images            
        %Images have been cropped to reduce cube size
    Idat = zeros(400,400,itr);    %Raw cropped images
    %Create column vector to hold Power data
    Pdat = zeros(itr, 1);
    %Create Matrix to hold trial parameters
    ParT = zeros(itr, 4);
    icount = 1;
    REtxt = fopen('temp.txt', 'wt');

    %Subtract 20 from Pdat for easy identification of empty elements later
        %In case fminsearch converges <itr iterations. Do not  want to process
        %or save excess images later
    Pdat        = Pdat - 20;

    % %Reflatten cube for Flat image comparison
    DM_Command = zeros(12);
    FlatCommand = NK_MultiDM_Command(DM_Command, flat);
    JR_UPDATE_MultiDM(drv_inf, FlatCommand);


    % Iterate, Optimize, and Process Data
    disp('Look for maximum output from fiber without coronagraph');

    for i=1:1
        % For aligning to speckle - repeat taking camera images until user input
        global KEY_IS_PRESSED

        KEY_IS_PRESSED = 0;

        gcf

        set(gcf, 'KeyPressFcn', @myKeyPressFcn);
        disp('Take out the coronagraph ');
        disp('When  ');
        while(~KEY_IS_PRESSED)
            % %_______Take Flat Data___________
            %Capture data of flattened DM
            %Acquire image
            cam.Acquisition.Freeze(true);

            %Take Data; Enable state returns to 0 automaticaly when 50 sample done
    %         reading = s.inputSingleScan;
    %         intensity = reading/pm_scale

            % %_______Process and Save Image Data______
            % %Extract image from RAM to array in Matlab
            [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
            % %Reshape array into plottable matrix
            img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
            % %Crop data to (200x200) around user input center; Save image
            Idat(:,:,1)  = img.Data(int16(Ycr-199):int16(Ycr+200),int16(Xcr-199):int16(Xcr+200));
            Iflat        = Idat(:,:,1);
            % %Draw cropped flat image
            set(himg, 'CData', Idat(:,:,1));
            axis tight equal off;
            drawnow;

        end
    end
    %     totalPower = input('Power out fiber? ');
    %     x_star = input('x_star? ');
    %     y_star = input('y_star? ');
    %     x_planet = input('x_planet? ');
    %     y_planet = input('y_planet? ');
end
% totalPower = totalPowerin;
    
% x_fib0in = input('Better guess for x planet in lam/D (press 0 if not)? ');
% if x_fib0in ~= 0
%     x_fib0 = x_fib0in;
% end
% 
% search = input('Search for a new optimal x planet (0 if NO, any other number for YES?) ');
% if search~=0
%     [x_fib,y_fib] = hcstt_FindPosiotionFiberv3(x_fib0,y_fib0);
% else
%     x_fib = x_fib0;
%     y_fib = y_fib0;
% end


norm = totalPower/totalPowerModel;

x_cent_cam = x_star;
y_cent_cam = y_star;
x_fib_pix = abs(x_star-x_planet);
end

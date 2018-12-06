close all
%clear

load('flat.mat');
load('FlatMap.mat');
load('angle_cal');
load('position_cal');

global icount flat REtxt apd Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond
global s pm_scale

%% DM Setup

% % Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera Setup

%Exposure time for camera
CExp    = 5;%ms

% cam = Initialize_Camera(CExp);
    NET.addAssembly(...
        'C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
% end

%Create CCD object handle
cam = uc480.Camera();

%Open and connect to first camera: 0 = 1st cam
char(cam.Init(0));

%Set display mode to bitmap (DiB); Captures to RAM
if ~strcmp(char(cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB)), ...
        'SUCCESS')
    error('Could not set display mode');
end

%Set color mode to 8-bit RAW
if ~strcmp(char(cam.PixelFormat.Set(uc480.Defines.ColorMode.SensorRaw8)), ...
        'SUCCESS')
    error('Could not set pixel format');
end

%Set trigger mode to software (single image acquisition)
if ~strcmp(char(cam.Trigger.Set(uc480.Defines.TriggerMode.Software)), ...
        'SUCCESS')
    error('Could not set trigger format');
end

%Set timing feautres: Pixel clock = frequency; Exposure in ms
cam.Timing.PixelClock.Set(30);
cam.Timing.Exposure.Set(CExp);

%Allocate image memory, define id as img.ID
[ErrChk, img.ID] = cam.Memory.Allocate(true);
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not allocate memory');
end

%Obtain image information
[ErrChk, img.Width, img.Height, img.Bits, img.Pitch] ...
    = cam.Memory.Inquire(img.ID);
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not get image information');
end

%Acquire image
if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
    error('Could not acquire image');
end

%Extract image, store in tmp as vector array
[ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not obtain image data');
end

%Reshape image into WxH matrix, convert tmp to uint8 vector
img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);

%Draw image;
himg = imshow(img.Data, 'Border', 'tight');
title('Sample Image Pre-Iteration');

%Defines Center of cropping Area
Xcr = 454;
Ycr = 344;

%% Power Meter Setup
% 
%Samp = 1;
%PowerMeterSetup(Samp);

%% APD Setup

% % Find a VISA-USB object.
% apd = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0699::0x0368::C014694::0::INSTR', 'Tag', '');
% 
% % Create the VISA-USB object if it does not exist
% % otherwise use the object that was found.
% if isempty(apd)
%     apd = visa('NI', 'USB0::0x0699::0x0368::C014694::0::INSTR');
% else
%     fclose(apd);
%     apd = apd(1);
% end
% 
% % Connect to instrument object, obj1.
% fopen(apd);

%%Setup Femto power meter

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');

s.Rate = 10;
s.DurationInSeconds = 1;

%starting conversion scale (in V/W)
pm_scale = 10^11;

%% Define Initial Parameters
Nactacross = 10;    %Number actuators across any given length

%SET Initial DM SIN PROPERTIES HERE: 
h0  = 10;           %Vary amplitude: 10 to 100, step of 10
x   = 5;         %Fix Distance: Set to meet given speckle
alp = 1.5;         %Vary Phase: 0 to 2pi, step of pi/10
q   = 0;         %Fix Angle: Set to meet given speckle

%SET NUMBER OF ITERATIONS TO PERFORM HERE:
itr = 200;

%Save initial condition as vector for optimization
par0 = [h0, q, x, alp];

%SET PROPERTIES FOR SEARCH
LB = [5.0, -0.3, 3.0, 0.0];              %Lower bounds; same elements as par0
UB = [20.0, 0.3, 10.0, 2*pi];       %Upper bounds; same elements as par0
% sTol = .0000001;                     %Smallest step size
% fTol = 3e-9;                   %Stops iterating when delPower < fTol


%% Flatten DM

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
%% Take background image
take_background = true;
Ncam = 400;
Nact = 12;
tint=0.1;
if(take_background)
    prompt = 'Take out light. Continue? ';
    x = input( prompt );
%     im_cam = zeros(Ncam,Ncam);
%     for II=1:15
%         im_camII = hcstt_TakeCamImage(true,false,tint);
%         im_cam = im_cam + im_camII/15;
%         pause(0.1)
%     end
%     backgroundCam = im_cam;
    backgroundSMF = hcstt_GetIntensityFIU(zeros(Nact,Nact),15,0);
    prompt = 'Put back light on. Continue? ';
    x = input( prompt );
else
    backgroundCam = zeros(Ncam,Ncam);
    backgroundSMF = 0;
end


%% Iterate, Optimize, and Process Data

n_pairs = 4; %initialize with 4 probes (for phase)
%generate equally spaced probes (4 different phis)
probe_phase = [0, 3.14/2, 3.14, 3.14*1.5];
amp_array = linspace(0,25,61)';

%start looping here

amp_cal_debug = zeros(length(amp_array),4);

amp_cal_debug(:,1) = amp_array;

for i=1:1
    % For aligning to speckle - repeat taking camera images until user input
    global KEY_IS_PRESSED

    KEY_IS_PRESSED = 0;

    gcf

    set(gcf, 'KeyPressFcn', @myKeyPressFcn);

    while(~KEY_IS_PRESSED)
        % %_______Take Flat Data___________
        %Capture data of flattened DM
        %Acquire image
        cam.Acquisition.Freeze(true);

        %Take Data; Enable state returns to 0 automaticaly when 50 sample done
        reading = s.inputSingleScan

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


    Pdat(1) = hcstt_GetIntensityFIU(zeros(12,12),10,backgroundSMF);
    Pflat = Pdat(1);

    %Pdat(1) = TakePMData(Samp);
    %Pflat = Pdat(1);

    pause(.5)

    cam.Acquisition.Freeze(true);

    %Take Data; Enable state returns to 0 automaticaly when 50 sample done
    reading = s.inputSingleScan

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

    %PAUSE HERE, GET COORDS
    pause(.5)

%     x_c = input('Center X Coord: ');
%     y_c = input('Center Y Coord: ');
% 
%     x_s = input('Speckle X Coord: ');
%     y_s = input('Speckle Y Coord: ');

%     distance_speckle = sqrt((x_s-x_c)^2+(y_s-y_c)^2);
%     angle_speckle = atan((y_s-y_c)/(x_s-x_c));

    %convert, define location of speckle
%     [c index] = min(abs(distance_speckle-position_cal(:,2)));
%     act = position_cal(index,1);
% 
%     [c index] = min(abs(angle_speckle-angle_cal(:,2)));
%     ang = angle_cal(index,1);
    x_fib_est = 2.58;
    info.backgroundSMF = backgroundSMF;
%     [ang_fib,actxc_fib] = hcstt_FindPosiotionFiberv4(x_fib_est,0,info);
    ang_fib = 0;
    actxc_fib = 2.735;
    for j=1:length(amp_array)
        intensity = zeros(n_pairs,1);

        for n=1:1:n_pairs
%             par_probe = [amp_array(j), ang, act, probe_phase(n)];
            dm_actuators_mat0 = hcstt_DMMapSin(amp_array(j), ang_fib, actxc_fib, probe_phase(n));
            intensity(n) = hcstt_GetIntensityFIU(dm_actuators_mat0,10,backgroundSMF);%DE_supMinimizer_3(par_probe);
        end
        figure(200)
        plot(probe_phase,intensity)
        title(['Probe amp ',num2str(amp_array(j))])
        %this part calculates phase from measurements
        t  = 1:n_pairs;                                                 % Time Vector
        sig  = intensity;                                                 % Signal Vector
        Ts = mean(diff(t));                                         % Sampling Time
        Fs = 1/Ts;                                                  % Sampling Frequency
        Fn = Fs/2;                                                  % Nyquist Frequency
        L  = length(sig);
        fts = fft(sig)/L;                                             % Normalised Fourier Transform
        Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
        Iv = 1:length(Fv);                                          % Index Vector
        amp_fts = abs(fts(Iv))*2;                                   % Spectrum Amplitude
        %phs_fts = angle(fts(Iv));                                   % Spectrum Phase
        
        amp_cal_debug(j,2) = (amp_fts(1)-2*Pflat)/2;
        amp_cal_debug(j,3) = amp_fts(2);
    
    end
    
end

amp_cal_debug(:,4) = amp_cal_debug(:,2)+amp_cal_debug(:,3);

amp_cal_temp = zeros(length(amp_array),2);
amp_cal_temp(:,1) = amp_cal_debug(:,1);
amp_cal_temp(:,2) = amp_cal_debug(:,4);

%% Disconnect all Devices
%Return to MATLAB Directory
% cd('C:\Users\Gary\Documents\MATLAB');

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end

% Close Power Meter
% newp.CloseDevices();


%Close APD
% fclose(apd);
% delete(apd);
% clear apd;

%Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    fprintf('Could not close camera');
else
    fprintf('Successfully Close CCD\n');
end

%%plot
plot(amp_cal_debug(:,1),amp_cal_debug(:,2),amp_cal_debug(:,1),amp_cal_debug(:,3),amp_cal_debug(:,1),amp_cal_debug(:,4))
legend('2','3','4')
amp_cal = amp_cal_debug;
save(['output',filesep,'sn_amp_calDec05.mat'],'amp_cal')
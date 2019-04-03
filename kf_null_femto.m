%Code for testing speckle nulling with a Kalman Filter
%Using the Femto Power Meter
%Note, currently, the software will tell you when you need to change the
%gain on the Femto power meter (and input the new gain so the software can
%convert it). Making this automated is possible with the LUCI-10 driver,
%and will hopefully be in future versions of this code

%INSTRUCTIONS FOR RUNNING THIS CODE
%I. Currently, manually load the amplitude calibration file for the location
%of the speckle that the algorithm is being tested on

%II.Place breakpoints at the following pause points:
%   1. Under "Take Pre-run Zero Data (MAKE SURE LIGHT IS OFF)" (around line
%       215)
%   2. Under "Pause to TURN ON LIGHT and ensure DM is flat" (around line
%       225)
%   3. Under "Take post run off axis data" (around line 489)
%   4. Under "Take post run zero (turn off light)" (around line 500)

%III. Run the code. At the first break point, block the light from going
%into the fiber, and click continue

%IV. At the next breakpoint, unblock the light, make sure the coronagraph
%is not aligned, and maximize the power into the fiber, and click continue

%V. A live camera image should be displayed. Realign the coronagraph and
%inject the test speckle into the fiber. Mark the pixel coordinates of the
%central star and of the speckle. Press a button to continue (you may need
%to first deactivate the cursor if you used it)

%VI. Enter the coordinates you marked down when prompted. Now the speckle
%nulling code should run.

%VII. At the next breakpoint, take out the coronagraph, maximize the
%signal, and click continue.

%VIII. At the next breakpoint, block the light from entering the fiber, and
%click continue. The test should now be complete, with results saved to the
%Results folder.
%% SETUP
close all
%clear

%load calibration files
addpath(genpath('Functions'))
addpath(genpath('Calibrations'))
load('Calibrations/flat.mat');
load('Calibrations/FlatMap.mat');
load('Calibrations/angle_cal');
load('Calibrations/position_cal');

%load amplitude calibration file
load('Calibrations/amp_cal_new_2p4.mat');
amp_cal = amp_cal_new_2p4;

global icount flat REtxt apd Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond
global s pm_scale

%define output dir
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\COO_user\Desktop\HCSTT\Results\';

%Name of TXT file for Raw Results
txtName = ['BandPassNulling' datestr(datetime,30) '.txt'];
Rawnm   = [fldDir,txtName];

%Name of Fits file for Image cube
txtName = ['ImageTrials' datestr(datetime,30) '.fits'];
Idatnm  = [fldDir,txtName];

%Name of Fits file for DM history
txtName = ['DMHistory' datestr(datetime,30) '.fits'];
DMdatnm  = [fldDir,txtName];

%Name of measurement and estimation data file
txtName = ['FilterData' datestr(datetime,30) '.mat'];
filtdatnm  = [fldDir,txtName];

%define number of cycles to run
itr = 30;
num_samples = itr;
num_steps = itr;

%% DM Setup

% % Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera Setup

%Exposure time for camera
CExp    = 10;%ms

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
Xcr = 500;
%Ycr = 780;
Ycr = 550;

%Xcr = 602;
%Ycr = 622;

%% Power Meter Setup

%Setup Femto power meter

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');

s.Rate = 10;
s.DurationInSeconds = 1/10;

%starting conversion scale (in V/W)
pm_scale = 10^11;

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

%Subtract 20 from Pdat for easy identification of empty elements later
    %In case fminsearch converges <itr iterations. Do not  want to process
    %or save excess images later
Pdat        = Pdat - 20;

% %Reflatten cube for Flat image comparison
DM_Command = zeros(12);
FlatCommand = NK_MultiDM_Command(DM_Command, flat);
JR_UPDATE_MultiDM(drv_inf, FlatCommand);

DM_history = zeros(12,12,num_steps);

%% Take Pre-run Zero Data (MAKE SURE LIGHT IS OFF)
pause(.05)

user_cont = 0;
while(~user_cont)

    reading = s.inputSingleScan

    user_cont = input('Continue? 1 or 0: ');
end

pm_scale = 10^input('Exponent? ');

prerun_zero_array = zeros(100,1);
for i=1:1:100
    prerun_zero_array(i) = s.inputSingleScan/pm_scale;
    %prerun_zero_array(i) = TakePMData(Samp);
end
prerun_zero = mean(prerun_zero_array);

%% Pause to TURN ON LIGHT and ensure DM is flat; Note: DM mech response time <100us
pause(.05);

user_cont = 0;
while(~user_cont)
    for i=1:1:50
        pause(.5)
        reading = s.inputSingleScan
    end
    user_cont = input('Continue? 1 or 0: ');
end

pm_scale = 10^input('Exponent? ');
% 
% % Take Prerun Off-Axis Throughput (take out coronagraph and maximize power)

prerun_offaxis_array = zeros(100,1);

%par something with 0 amplitude
for i=1:1:100
   prerun_offaxis_array(i) = s.inputSingleScan/pm_scale;
   %prerun_offaxis_array(i) = TakePMData(Samp);
end

prerun_offaxis = mean(prerun_offaxis_array);

%% Take Flat Data (PUT IN CORONAGRAPH AND ALIGN TO SPECKLE BEFORE PRESSING A KEY)
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


pm_scale = 10^input('Exponent? ');

Pdat(1) = reading/pm_scale;
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

pause(.5)

x_c = input('Center X Coord: ');
y_c = input('Center Y Coord: ');

x_s = input('Speckle X Coord: ');
y_s = input('Speckle Y Coord: ');

distance_speckle = sqrt((x_s-x_c)^2+(y_s-y_c)^2);
angle_speckle = atan((y_s-y_c)/(x_s-x_c));

%convert, define location of speckle
[c index] = min(abs(distance_speckle-position_cal(:,2)));
act = position_cal(index,1);

[c index] = min(abs(angle_speckle-angle_cal(:,2)));
ang = angle_cal(index,1);

%load amplitude calibration file and choose calibration
load('amplitude_calibrations.mat');

%OPTION 1: Use this option generically. Will use the calibration at the
%closest distance that has calibration data. First columns at every
%distance do not need to match up.
act_array = amplitude_calibrations.actuators;
[c index] = min(abs(act-act_array));

dim_amp_cal = size(amplitude_calibrations.calibration);
amp_cal = zeros(dim_amp_cal(2),dim_amp_cal(3));
amp_cal(:,1) = amplitude_calibrations.calibration(index,:,1);
amp_cal(:,2) = amplitude_calibrations.calibration(index,:,2);

% %OPTION 2: Use this option IF AND ONLY IF the calibration at the different
% %distances spans the same range of DM amplitudes. This will linearly
% %interpolate between the calibrations at different distances, assuming that
% %the first columns of every distance calibration match up.
% 
% %find nearest two number of actuators and fraction between
% act_array = amplitude_calibrations.actuators;
% [c index] = min(abs(act-act_array));
% dist_sign = act_array(index)-act;
% if dist_sign<0
%     act_low_idx = index;
%     act_high_idx = index+1;
% else
%     act_low_idx = index-1;
%     act_high_idx = index;    
% end
% 
% act_low = act_array(act_low_idx);
% cal_low = amplitude_calonratopms.calibration(act_low_idx,:,:);
% act_high = act_array(act_high_idx);
% cal_high = amplitude_calibrations.calibration(act_high_idx,:,:);
% 
% fraction = (act-act_low)/(act_high-act_low);
% 
% %create interpolated calibration array
% dim_amp_cal = size(amplitude_calibrations.calibration);
% amp_cal = zeros(dim_amp_cal(2),dim_amp_cal(3));
% amp_cal(:,1) = cal_low(1,:,1);
% amp_cal(:,2) = cal_low(1,:,2)+fraction*(cal_high(1,:,2)-cal_low(1,:,2));

%% Set up savefile
REtxt = fopen(Rawnm, 'wt');

fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
% fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);
fprintf(REtxt, 'Prerun Zero: %05.2e \n', prerun_zero);
fprintf(REtxt, 'Prerun Off Axis Throughput: %05.2e \n', prerun_offaxis);

%Define a counter to allow matching variables to each other
% % %Use fminsearch to iterate and optimize
currentTime = clock;
currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
startingSecond = currentSecond;

icount = 1;     %Start counter at 1

%% Inject fake speckle (OPTIONAL)
%After the realignment, the speckles were dimmer and so there wasn't much
%to null. Add a fake speckle of ~5 nm (in addition to exising speckle)
%and pretend that it's the default flat map

% par_in = [5.0, ang, act, 0.4];
% pow = DE_supMinimizer_3(par_in);
% 
% DM_Map = DE_DMMapSin(par_in(1),par_in(2),par_in(3),par_in(4));
% flat = flat + DM_Map;

%% Run Speckle Nulling

n_pairs = 4; %initialize with 4 probes (for phase)
%generate equally spaced probes (4 different phis)
probe_phase = [0, 3.14/2, 3.14, 3.14*1.5];
probe_amp = 3;

steps = 1:num_steps;

%load interaction matrix
gamma = [1,0;0,1];

%covariance initialization
P_0 = eye(2);

%process noise (on x)
w_k = [.1,0;0,.1]*10^-5;
Q_k = w_k*w_k';

%measurement noise
phase_noise = 1*10^-5;
amplitude_noise = 0.5*10^-5;
R_k = [phase_noise^2,0;0,amplitude_noise^2];

%start looping here

%initialize values
x_k_m = [0;5]; %prior estimates %JLlop: from meas
P_k_m = P_0; %prior covariance
u_k = [0;0]; %control input

%initialize logs
measurements = zeros(num_steps,2);
estimations = zeros(num_steps,2);
inputs = zeros(2,num_steps);

for i=steps
    
    intensity = zeros(n_pairs,1);
    
    H_k = eye(2); %observation matrix
    
    for n=1:1:n_pairs
        par_probe = [probe_amp, ang, act, probe_phase(n)];

        intensity_raw = DE_supMinimizer_3(par_probe);
        
        %simulate gaussian detector noise
        noise_sigma = 0.5*10^-15;
        %draw from gaussian
        noise_amp = normrnd(0,noise_sigma);
        intensity(n) = intensity_raw+noise_amp;
    end
    
    pause(.5)

    %calculate phase and amplitude from measurements
    t  = 1:n_pairs;                                                 % Time Vector
    sig  = intensity;                                                 % Signal Vector
    Ts = mean(diff(t));                                         % Sampling Time
    Fs = 1/Ts;                                                  % Sampling Frequency
    Fn = Fs/2;                                                  % Nyquist Frequency
    L  = length(sig);
    fts = fft(sig)/L;                                             % Normalised Fourier Transform
    Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
    Iv = 1:length(Fv);                                          % Index Vector
    amp_fts = abs(fts(1:3))*2;                                   % Spectrum Amplitude
    phs_fts = angle(fts(1:3));                                   % Spectrum Phase
    
    if phs_fts(2)<0
       phs_fts(2) = phs_fts(2)+2*pi;
    end
    
    amp_meas = amp_fts(1)/2;
    
    %try subtracting a little below lowest known reading
    %offset = 1*10^-10;
    %amp_meas = amp_meas-offset;
    
    %look up closest DM amplitude value
    [c index] = min(abs(amp_meas-amp_cal(:,2)));
    
    %if less than lowest calibrated value, interpolate
    if index <= 2
        amp_solved = amp_meas/amp_cal(2,2)*amp_cal(2,1);
    else
        amp_solved = amp_cal(index,1);
    end
    
    phase_solved = phs_fts(2);
    z_k = [phase_solved; amp_solved]; %measurement

    %Kalman filter
    
    x_k_p = x_k_m+P_k_m*H_k'*inv(H_k*P_k_m*H_k'+R_k)*(z_k-H_k*x_k_m); %update estimate with measurement
    P_k_p = inv(inv(P_k_m)+H_k'*inv(R_k)*H_k); %update covariance
    
    %control input gain (tune)
    K_p = 1;
    
    %control input calculation
    u_k = [-x_k_p(1); -K_p*x_k_p(2)];
    
    %convert control input and actuate DM
    par_in = [u_k(2), ang, act, u_k(1)];
    pow = DE_supMinimizer_3(par_in);
    
    x_k_m = x_k_p+gamma*u_k; %new prior based on control input
    
    %heuristic nonlinear "bottoming out" adjustment
    %amp_offset = 2;
    %x_k_m(2) = x_k_m(2)+amp_offset;
    
    P_k_m = P_k_p+Q_k; %new prior covariance
    
    %set probe amp to expected amplitude, if too low, set to minimum
    %discernible amplitude
    probe_amp = x_k_p(2)/2;
    
    if probe_amp > 3
       probe_amp = 3; 
    end
    
    if probe_amp < 0.05       
        probe_amp = 0.05;  
    end
    
    %update flat map
    DM_Map = DE_DMMapSin(par_in(1),par_in(2),par_in(3),par_in(4));
    flat = flat + DM_Map;
    
    DM_history(:,:,i) = flat;

    %log for analysis
    %real_amp_meas(i) = amp_solved;
    measurements(i,:) = z_k;
    estimations(i,:) = x_k_p;
    inputs(:,i) = u_k;
    
end

figure(1)
plot(steps,measurements(:,1),steps,estimations(:,1))

figure(2)
plot(steps,measurements(:,2),steps,estimations(:,2))

%     [Par, Pmin] = patternsearch(@DE_supMinimizer_3, par0, [], [], [], [], ... 
%          LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));

%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
%Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

% Save all images as single cube
fitswrite(Idat, Idatnm);

%save DM history
fitswrite(DM_history, DMdatnm);


%% Take post run off axis data (take out coronagraph and realign to max)
pause(.05);

user_cont = 0;
while(~user_cont)
    for i=1:1:50
        pause(.5)
        reading = s.inputSingleScan
    end
    user_cont = input('Continue? 1 or 0: ');
end

pm_scale = 10^input('Exponent? ');

postrun_offaxis_array = zeros(100,1);
for i=1:1:100
   postrun_offaxis_array(i) = s.inputSingleScan/pm_scale;
end
postrun_offaxis = mean(postrun_offaxis_array);

fprintf(REtxt, 'Postrun Off Axis Throughput: %05.2e \n', postrun_offaxis);

%% Take post run zero (turn off light)
pause(.05);

user_cont = 0;
while(~user_cont)
    reading = s.inputSingleScan
    user_cont = input('Continue? 1 or 0: ');
end

pm_scale = 10^input('Exponent? ');

postrun_zero_array = zeros(100,1);
for i=1:1:100
   postrun_zero_array(i) = s.inputSingleScan/pm_scale;
end
postrun_zero = mean(postrun_zero_array);
fprintf(REtxt, 'Postrun Zero: %05.2e COMMENTS: 6e-11 offset \n', postrun_zero);

%% SHUTDOWN

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

save(filtdatnm, 'measurements','estimations');
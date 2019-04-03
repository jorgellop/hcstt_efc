%edited script to poke each actuator and take corresponding camera image

close all
%clear

load('flat.mat');
load('FlatMap.mat');

global icount flat REtxt apd Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond

%% DEFINE OUTPUT DIRECTORY AND FILE; DEFINE SETUP PARAMETERS
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\COO_user\Desktop\HCSTT\Results\';

%Name of TXT file for Raw Results
txtName = ['BandPassNulling' datestr(datetime,30) '.txt'];
Rawnm   = [fldDir,txtName];

%Name of Fits file for Image cube
txtName = ['ImageTrials' datestr(datetime,30) '.fits'];
Idatnm  = [fldDir,txtName];

%Number of Samples for PM data set; MUST BE AN INTEGER
% Samp    = 500;

%% DM Setup

% % Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera Setup

%Exposure time for camera
CExp    = 1;%ms

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
Xcr = 520;
Ycr = 780;

%% Power Meter Setup

%Samp = 1;
%PowerMeterSetup(Samp);

% %% APD Setup
% 
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
Idat = zeros(200,200,itr);    %Raw cropped images
%Create column vector to hold Power data
%Pdat = zeros(itr, 1);
%Create Matrix to hold trial parameters
%ParT = zeros(itr, 4);

%Subtract 20 from Pdat for easy identification of empty elements later
    %In case fminsearch converges <itr iterations. Do not  want to process
    %or save excess images later
%Pdat        = Pdat - 20;

% %Reflatten cube for Flat image comparison
DM_Command = zeros(12);
FlatCommand = NK_MultiDM_Command(DM_Command, flat);
JR_UPDATE_MultiDM(drv_inf, FlatCommand);

pause(.1)

%% Take Flat Data
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
   % Pdat(1) = TakePMData(Samp);
   % Pflat = Pdat(1);

    % %_______Process and Save Image Data______
    % %Extract image from RAM to array in Matlab
    [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    % %Reshape array into plottable matrix
    img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
    % %Crop data to (200x200) around user input center; Save image
    temp  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
    Iflat        = temp;
    % %Draw cropped flat image
    set(himg, 'CData', temp);
    axis tight equal off;
    drawnow;
    
end

%% Take poke data

%Define a counter to allow matching variables to each other
icount = 1;     %Start counter at 1 to compensate for flat image
% % %Use fminsearch to iterate and optimize
currentTime = clock;
currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
startingSecond = currentSecond;

poke_array = zeros(12);

poke_command = NK_MultiDM_Command(poke_array, flat);
JR_UPDATE_MultiDM(drv_inf, poke_command);
pause(1);


i_local = 0;
for a=1:1:10
    i_local = i_local+1;    
            %Acquire image
    cam.Acquisition.Freeze(true);
    
    %Measure Time
    currentTime = clock;
    currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
    
    %Take Data; Enable state returns to 0 automaticaly when 50 sample done
    %Samp = 1;
    %Pdat(icount) = TakePMData(Samp);

        %_______Process and Save Image Data_______
    %     %Extract image from RAM to array in Matlab
        [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    %     %Reshape array into plottable matrix
        img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
    %     %Crop data to (200x200) around user input center
        Idat(:,:,icount)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));  

end

poke_amp = 50;


for m=1:1:12
    for n=1:1:12
        
        poke_array = zeros(12);
        poke_array(m,n) = poke_amp;
        
        poke_command = NK_MultiDM_Command(poke_array, flat);
        JR_UPDATE_MultiDM(drv_inf, poke_command);
        pause(1);
        
        for a=1:1:10
            i_local = i_local+1;
            %Acquire image
    cam.Acquisition.Freeze(true);
    
    %Measure Time
    currentTime = clock;
    currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
    
    %Take Data; Enable state returns to 0 automaticaly when 50 sample done
    %Samp = 1;
    %Pdat(icount) = TakePMData(Samp);
        
        %_______Process and Save Image Data_______
    %     %Extract image from RAM to array in Matlab
        [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    %     %Reshape array into plottable matrix
        img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
    %     %Crop data to (200x200) around user input center
        Idat(:,:,i_local)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
        end

    end
end

% %Save all images as single cube
fitswrite(Idat(:,:,1:i_local), Idatnm);



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
%newp.CloseDevices();

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
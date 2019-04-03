%{
Speckling Nulling Proof: Calibrate Speckle Location
- Place Speckles alfkj

******************************************************
- Arguments:
    NONE        User Input on Image       

- Returns:
    NONE        Outputs a .FITS cube with properties:
                    WidthxHeightxImage
                    200  x 200  x (heights -> phase)
- Dependencies:
    OPEN_multiDM()  Create DM connection
    UPDATE_multiDM()Place new DM Map
    CLOSE_multiDM() Disconnect DM
    uc480DotNet.dll Library of CCD Camera Communications
    DE_DMW_Sin      Create Sinusoidal Map for DM
******************************************************


Compiled By:    Daniel Echeverri and Jerry Xuan
Last Modified:  08/09/2016
%}
%close all
clear all

%% DM set-up
% Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera set-up
%Check if camera assembly is present. Import Otherwise
%May need to change search location if library is moved
asm = System.AppDomain.CurrentDomain.GetAssemblies;
if ~any(arrayfun(@(n) strncmpi(char(asm.Get(n-1).FullName), ...
        'uc480DotNet', length('uc480DotNet')), 1:asm.Length))
    NET.addAssembly(...
        'C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
end

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
cam.Timing.Exposure.Set(30);

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
%PROMPTS USER TO SELECT CENTER OF STAR FOR CROPPING (x+-100, y+-100)
%[Xcr, Ycr] = ginput(1);
Xcr = 468;
Ycr = 520;
%% Define DM surface properties
Nactacross = 10;    %Number actuators across any given length
                        %VALUE IS AN ESTIMATE; NEED TO CONFIRM

%SET DM SIN PROPERTIES HERE: Vary 2, fix 2
h0  = 100;                   %Vary amplitude: 10 to 150, step of 10
x   = 4;                  %Fix Distance: Set to meet given speckle
alp = pi/2;                    %Vary Phase: 0 to 2pi, step of pi/10
q   = 4*pi/10;                    %Fix Angle: Set to meet given speckle

%Create Matrix to hold images: Width(200)xHeight(200)x#Images(300)
        %Images have been cropped to optimize cube size
%data_cube = zeros(200,200,numel(alps)*numel(h0s));

%Create new cropped figure: Shows actual saved data LIVE
figure();
%Crop data to (200x200) around user input center
imgCR  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
himgCR = imshow(imgCR);

%Iterate Through Amplitude Values
%i = 1;              %Index for saving in correct cube level
%for h0 = h0s
    %Iterate Through Phases
%    for alp = alps
        %Calculate and place new map on DM surface
        DE_DMW_Sin(h0, q, x, alp, drv_inf);
        %DE_DMW_Flat(0, drv_inf);
        %Acquire image
        cam.Acquisition.Freeze(true);
        %Extract image from RAM to array in Matlab
        [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
        %Reshape array into plottable matrix
        img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
        %Crop data to (200x200) around user chosen center
        imgCR  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
        %Draw Cropped image
        set(himgCR, 'CData', imgCR);
        drawnow;
        %Allow user to see image for an instant
 %       pause(.20);
        %Add data to output cube
%        data_cube(:,:,i) = imgCR;
        %Increment index
%        i = i + 1; 
%    end
%end

%Create name for saving cube as .FITS file
%txtnm = ['h0',num2str(h0s(1)),'_',num2str(h0s(end)),'_', ...
%    num2str(Nactacross/x,2), 'LD_', ...
%    num2str(q*180/pi), 'deg_Phase', ...
%    num2str(alps(1),2), '_', num2str(alps(end),2)];
%fitswrite(data_cube, [txtnm, '.fits']);
%%
%Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    error('Could not close camera');
end
 
% Disable and close Multi-DM driver USB connection
err = CLOSE_multiDM(drv_inf);
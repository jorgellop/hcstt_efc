%% Camera Setup

fldDir  = 'C:\Users\COO_user\Desktop\HCSTT\Results\';

txtName = ['ImageTrials' datestr(datetime,30) '.fits'];
Idatnm  = [fldDir,txtName];

itr = 100;

Idat = zeros(800,800,itr); 

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

%% Loop

for icount=1:itr

    cam.Acquisition.Freeze(true);

    %_______Process and Save Image Data_______
    %     %Extract image from RAM to array in Matlab
        [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
    %     %Reshape array into plottable matrix
        img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
    %     %Crop data to (200x200) around user input center
        Idat(:,:,icount)  = img.Data(200:1000,200:1000);
end

fitswrite(Idat(:,:,1:icount), Idatnm);

%% Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    fprintf('Could not close camera');
else
    fprintf('Successfully Close CCD\n');
end
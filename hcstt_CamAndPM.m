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

global cam img Xcr Ycr CExp s drv_inf flat itr himg pm_scale

addpath(genpath('utils'),genpath('export_scripts'));

hcstt_Initialize(true)
    
disp('Flattening DM');

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

for i=1:1
    % For aligning to speckle - repeat taking camera images until user input
    global KEY_IS_PRESSED

    KEY_IS_PRESSED = 0;

    gcf

    set(gcf, 'KeyPressFcn', @myKeyPressFcn);
        himg = imshow(img.Data, 'Border', 'tight');
        title('Sample Image Pre-Iteration');

    while(~KEY_IS_PRESSED)
        % %_______Take Flat Data___________
        %Capture data of flattened DM
        %Acquire image
        cam.Acquisition.Freeze(true);

        %Take Data; Enable state returns to 0 automaticaly when 50 sample done
        reading = s.inputSingleScan;
        intensity = reading/pm_scale

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

hcstt_DisconnectDevices()
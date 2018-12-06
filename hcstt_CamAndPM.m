%{

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
Idat = zeros(200,200,itr);    %Raw cropped images
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
        

    while(~KEY_IS_PRESSED)
        % %_______Take Flat Data___________
        %Capture data of flattened DM
        %Acquire image
        cam.Acquisition.Freeze(true);

        %Take Data; Enable state returns to 0 automaticaly when 50 sample done
%         reading = s.inputSingleScan;
        intensity = hcstt_GetIntensityFIU(zeros(12,12),1,0);%reading/pm_scale

        % %_______Process and Save Image Data______
        % %Extract image from RAM to array in Matlab
        [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
        % %Reshape array into plottable matrix
        img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
        % %Crop data to (200x200) around user input center; Save image
        Idat(:,:,1)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
        Iflat        = Idat(:,:,1);
        % %Draw cropped flat image
        set(himg, 'CData', Idat(:,:,1));
        axis tight equal off;
        title(['Intensity: ',num2str(intensity)]);

        drawnow;
    end
end

hcstt_DisconnectDevices()
close all

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

%% DM Setup

% % Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera Setup
% 
% %Exposure time for camera
% CExp    = 1;%ms
% 
% % cam = Initialize_Camera(CExp);
%     NET.addAssembly(...
%         'C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
% % end
% 
% %Create CCD object handle
% cam = uc480.Camera();
% 
% %Open and connect to first camera: 0 = 1st cam
% char(cam.Init(0));
% 
% %Set display mode to bitmap (DiB); Captures to RAM
% if ~strcmp(char(cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB)), ...
%         'SUCCESS')
%     error('Could not set display mode');
% end
% 
% %Set color mode to 8-bit RAW
% if ~strcmp(char(cam.PixelFormat.Set(uc480.Defines.ColorMode.SensorRaw8)), ...
%         'SUCCESS')
%     error('Could not set pixel format');
% end
% 
% %Set trigger mode to software (single image acquisition)
% if ~strcmp(char(cam.Trigger.Set(uc480.Defines.TriggerMode.Software)), ...
%         'SUCCESS')
%     error('Could not set trigger format');
% end
% 
% %Set timing feautres: Pixel clock = frequency; Exposure in ms
% cam.Timing.PixelClock.Set(30);
% cam.Timing.Exposure.Set(CExp);
% 
% %Allocate image memory, define id as img.ID
% [ErrChk, img.ID] = cam.Memory.Allocate(true);
% if ~strcmp(char(ErrChk), 'SUCCESS')
%     error('Could not allocate memory');
% end
% 
% %Obtain image information
% [ErrChk, img.Width, img.Height, img.Bits, img.Pitch] ...
%     = cam.Memory.Inquire(img.ID);
% if ~strcmp(char(ErrChk), 'SUCCESS')
%     error('Could not get image information');
% end
% 
% %Acquire image
% if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
%     error('Could not acquire image');
% end
% 
% %Extract image, store in tmp as vector array
% [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
% if ~strcmp(char(ErrChk), 'SUCCESS')
%     error('Could not obtain image data');
% end
% 
% %Reshape image into WxH matrix, convert tmp to uint8 vector
% img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
% 
% %Draw image;
% himg = imshow(img.Data, 'Border', 'tight');
% title('Sample Image Pre-Iteration');
% 
% %Defines Center of cropping Area
% Xcr = 520;
% Ycr = 780;

%% Power Meter Setup



%% Define Initial Parameters

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

%Create Matrix to hold images: Width(200)xHeight(200)x#Images            
    %Images have been cropped to reduce cube size
Idat = zeros(200,200,itr);    %Raw cropped images
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

%% Take Pre-run Zero Data (MAKE SURE LIGHT IS OFF)
pause(.05)
prerun_zero_array = zeros(100,1);
for i=1:1:100
    prerun_zero_array(i) = TakePMData(Samp);
end
prerun_zero = mean(prerun_zero_array);

%% Pause to TURN ON LIGHT and ensure DM is flat; Note: DM mech response time <100us
pause(.05);

% Take Prerun Off-Axis Throughput (take out coronagraph and maximize power)
prerun_offaxis_array = zeros(100,1);

%par something with 0 amplitude
for i=1:1:100
   prerun_offaxis_array(i) = TakePMData(Samp);
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
    Pdat(1) = TakePMData(Samp);
    Pflat = Pdat(1);

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
    drawnow;
    
end

%% Iterate, Optimize, and Process Data

REtxt = fopen(Rawnm, 'wt');

fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
% fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);
fprintf(REtxt, 'Prerun Zero: %05.2e \n', prerun_zero);
fprintf(REtxt, 'Prerun Off Axis Throughput: %05.2e \n', prerun_offaxis);

%Define a counter to allow matching variables to each other
icount = 1;     %Start counter at 1 to compensate for flat image
% % %Use fminsearch to iterate and optimize
currentTime = clock;
currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
startingSecond = currentSecond;

for i=1:10
    power = DE_supMinimizer_3(par0);
end

% for i=1:8
%     minpower = 1000;
%     minpar = [0,5,0,0];
%     
%     par_min = [0, 1.4, 3, 0];
%     par_step = [4, 0.1, 1, 0.4];
%     par_max = [16, 1.6, 8, 6.0];
% 
%     for h0=par_min(1):par_step(1):par_max(1)
%         for q=par_min(2):par_step(2):par_max(2)
%             for x=par_min(3):par_step(3):par_max(3)
%                 for a=par_min(4):par_step(4):par_max(4)
%                     par = [h0,q,x,a];
%                     power = DE_supMinimizer_3(par);   
%                     if(power < minpower)
%                         minpower = power;
%                         minpar = [h0,q,x,a];
%                     end
%                 end
%             end
%         end
%     end
% 
%     for j=1:4
%         if(minpar(j) - par_step(j) > par_min(j))
%             par_min(j) = minpar(j) - par_step(j);
%         end
%         if(minpar(j) + par_step(j) < par_max(j))
%             par_max(j) = minpar(j) + par_step(j);
%         end
%     end
%     
%     par_step = [2, 0.05, 0.5, 0.1];
%     
%     for h0=par_min(1):par_step(1):par_max(1)
%         for x=par_min(3):par_step(3):par_max(3)
%             q = minpar(2);
%             for a=par_min(4):par_step(4):par_max(4)
%                 par = [h0,q,x,a];
%                 power = DE_supMinimizer_3(par);   
%                 if(power < minpower)
%                     minpower = power;
%                     minpar = [h0,q,x,a];
%                 end
%             end
%         end
%     end
% 
%     DM_Map = DE_DMMapSin(minpar(1),minpar(2),minpar(3),minpar(4));
%     flat = flat + DM_Map;    
% 
% end

%minimize power
amp_init = 10;
ang_init = 1.5;
act_init = 5;
pha_init = 1.2;

global fmin_global
fmin_global = Inf;
global xmin_global


for n=1:1:8
    xmin_global = [0,0,5,0];
    opts_coarse = optimoptions(@fmincon,'Algorithm','interior-point','StepTolerance',.05,...
        'DiffMinChange', 0.2);
    problem_coarse = createOptimProblem('fmincon','x0',[20, ang_init,act_init,pha_init],...
        'objective',@DE_supMinimizer_3, 'lb',[10, 1.3, 3.0, 0.0], 'ub',[10, 1.7, 10.0, 6.2], 'options',opts_coarse);

    ms = MultiStart('FunctionTolerance',0,'XTolerance',.05, 'Display','iter',...
        'StartPointsToRun','bounds', 'OutputFcn', @StopIfReached,'PlotFcn',@gsplotbestf);



    [xmin_coarse, fmin_coarse] = run(ms,problem_coarse,18);
    
    opts_fine = optimoptions(@fmincon,'StepTolerance',.05,'DiffMinChange',.05);
    problem_fine = createOptimProblem('fmincon','x0', xmin_global,'objective',@DE_supMinimizer_3,...
        'lb', [10, xmin_global(2)-.05, xmin_global(3)-.5, xmin_global(4)-.1], ...
        'ub', [10, xmin_global(2)+.05, xmin_global(3)+.5, xmin_global(4)+.1], ...
        'options', opts_fine);
    
    [xmin_fine, fmin_fine] = run(ms,problem_fine,18);
    
    %sweep amplitudes
    for amp_sweep=0:5:70
       par_in = [amp_sweep, xmin_global(2), xmin_global(3), xmin_global(4)];
       pow = DE_supMinimizer_3(par_in); 
    end
    
    for amp_step=-2:1:2
       par_in =  [xmin_global(1)+amp_step, xmin_global(2), xmin_global(3), xmin_global(4)];
       pow = DE_supMinimizer_3(par_in);
    end
    
    xmin = xmin_global;
    fmin = fmin_global;
    
    DM_Map = DE_DMMapSin(xmin_global(1),xmin_global(2),xmin_global(3),xmin_global(4));
    flat = flat + DM_Map;    
end

%set to found min
pow = DE_supMinimizer_3(xmin_global);

%     [Par, Pmin] = patternsearch(@DE_supMinimizer_3, par0, [], [], [], [], ... 
%          LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));

%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
%Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

% %Save all images as single cube
fitswrite(Idat(:,:,1:icount), Idatnm);


%% Take post run off axis data (take out coronagraph and realign to max)
pause(.05);
postrun_offaxis_array = zeros(100,1);
for i=1:1:100
   postrun_offaxis_array(i) = TakePMData(Samp); 
end
postrun_offaxis = mean(postrun_offaxis_array);

fprintf(REtxt, 'Postrun Off Axis Throughput: %05.2e \n', postrun_offaxis);

%% Take post run zero (turn off light)
pause(.05);
postrun_zero_array = zeros(100,1);
for i=1:1:100
   postrun_zero_array(i) = TakePMData(Samp);
end
postrun_zero = mean(postrun_zero_array);
fprintf(REtxt, 'Postrun Zero: %05.2e \n', postrun_zero);

%% Disconnect all Devices

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end

%Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    fprintf('Could not close camera');
else
    fprintf('Successfully Close CCD\n');
end
close all
%clear

load('flat.mat');
load('FlatMap.mat');

global icount flat REtxt apd Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond

%% DEFINE OUTPUT DIRECTORY AND FILE; DEFINE SETUP PARAMETERS
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\COO_user\Desktop\HCSTT\Results\';

%Name of TXT file for PM values
PMVnm   = [fldDir,'PM_Values2.txt'];

%Name of TXT file for Raw Results
txtName = ['BandPassNulling' datestr(datetime,30) '.txt'];
Rawnm   = [fldDir,txtName];

%Name of FITS file for Optimum Image
OptImnm = [fldDir,'OptImage.fits'];

%Name of Fits file for Image cube
Idatnm  = [fldDir,'ImageTrials.fits'];

%Number of Samples for PM data set; MUST BE AN INTEGER
% Samp    = 500;

%Exposure time for camera
CExp    = 25;%ms

%% DM Setup

% % Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%% Camera Setup

%cam = Initialize_Camera(CExp);

%Defines Center of cropping Area
Xcr = 468;
Ycr = 520;

%% Power Meter Setup

Samp = 1;
PowerMeterSetup(Samp);

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
itr = 30000;

%Save initial condition as vector for optimization
par0 = [h0, q, x, alp];

%SET PROPERTIES FOR SEARCH
LB = [5.0, -0.3, 3.0, 0.0];              %Lower bounds; same elements as par0
UB = [20.0, 0.3, 10.0, 2*pi];       %Upper bounds; same elements as par0
% sTol = .0000001;                     %Smallest step size
% fTol = 3e-9;                   %Stops iterating when delPower < fTol

%% Take Flat Data

%Create StringBuilder to retrieve and hold data; Set cap. to expected size
sb = System.Text.StringBuilder();
sb.Capacity = 4100;

%Create Matrix to hold images: Width(200)xHeight(200)x#Images            
    %Images have been cropped to reduce cube size
%Idat = zeros(200,200,itr);    %Raw cropped images
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

%Pause to ensure DM is flat; Note: DM mech response time <100us
pause(.05);

% %_______Take Flat Data___________
% %Capture data of flattened DM
% %Acquire image
% cam.Acquisition.Freeze(true);
%Take Data; Enable state returns to 0 automaticaly when 50 sample done
%Pdat(1) = APDMeasurement(apd);
Pdat(1) = TakePMData(Samp);
Pflat = Pdat(1);

% %_______Process and Save Image Data______
% %Extract image from RAM to array in Matlab
% [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
% %Reshape array into plottable matrix
% img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
% %Crop data to (200x200) around user input center; Save image
% Idat(:,:,1)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
% Iflat        = Idat(:,:,1);
% %Draw cropped flat image
% set(himg, 'CData', Idat(:,:,1));
% axis tight equal off;
% drawnow;

%% Iterate, Optimize, and Process Data

REtxt = fopen(Rawnm, 'wt');

fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
% fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);

%Define a counter to allow matching variables to each other
icount = 1;     %Start counter at 1 to compensate for flat image
% % %Use fminsearch to iterate and optimize
currentTime = clock;
currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
startingSecond = currentSecond;


%for i=1:8
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

    %DM_Map = DE_DMMapSin(minpar(1),minpar(2),minpar(3),minpar(4));
    %flat = flat + DM_Map;    

%end

%minimize power
amp_init = 10;
ang_init = 1.5;
act_init = 5;
pha_init = 1.2;

%starting points
%phase
phase_array = (0:.4:6)';
act_array = (3:2:9)';
angle_array = (1.4:.1:1.6)';
amp_array = (10:10:70)';

act_copied = zeros(length(phase_array),length(act_array));
for n=1:1:length(act_array)
    act_copied(:,n) = act_array(n);
end

total_length = length(phase_array)*length(act_array)*length(angle_array);
angle_copied = zeros(length(phase_array)*length(act_array),length(angle_array));
for n=1:1:length(angle_array)
    angle_copied(:,n) = angle_array(n);
end

%stack
angle_stack = reshape(angle_copied,[total_length,1]);

act_stack = reshape(act_copied,[length(phase_array)*length(act_array), 1]);
act_stack = repmat(act_stack,[length(angle_array), 1]);

phase_stack = repmat(phase_array,[length(act_array)*length(angle_array), 1]);

start_array = zeros(total_length,4);

%fill
start_array(:,1) = 20;
start_array(:,2) = angle_stack;
start_array(:,3) = act_stack;
start_array(:,4) = phase_stack;

%shuffle
shuffled_array = start_array(randperm(total_length),:);


%%%%%%%%%%%%%%%%%

tpoints = CustomStartPointSet(shuffled_array);
rpts = RandomStartPointSet('NumStartPoints',200);
allpts = {tpoints,rpts};



global fmin_global
fmin_global = Inf;
global xmin_global


for n=1:1:8
    xmin_global = [0,0,5,0];
    opts_coarse = optimoptions(@fmincon,'Algorithm','interior-point','StepTolerance',.05,...
        'DiffMinChange', 0.1);
    problem_coarse = createOptimProblem('fmincon','x0',[20, ang_init,act_init,pha_init],...
        'objective',@DE_supMinimizer_3, 'lb',[0, 1.3, 3.0, 0.0], 'ub',[20, 1.7, 10.0, 6.2], 'options',opts_coarse);

    ms = MultiStart('FunctionTolerance',0,'XTolerance',.05, 'Display','iter',...
        'StartPointsToRun','bounds', 'OutputFcn', @StopIfReached,'PlotFcn',@gsplotbestf);



    [xmin_coarse, fmin_coarse] = run(ms,problem_coarse,18);
    
    opts_fine = optimoptions(@fmincon,'StepTolerance',.05,'DiffMinChange',.01);
    problem_fine = createOptimProblem('fmincon','x0', xmin_global,'objective',@DE_supMinimizer_3,...
        'lb', [xmin_global(1)-2, xmin_global(2)-.05, xmin_global(3)-.5, xmin_global(4)-.1], ...
        'ub', [xmin_global(1)+2, xmin_global(2)+.05, xmin_global(3)+5, xmin_global(4)+.1], ...
        'options', opts_fine);
    
    [xmin_fine, fmin_fine] = run(ms,problem_fine,18);
    
    xmin = xmin_global;
    fmin = fmin_global;
    
    DM_Map = DE_DMMapSin(xmin_global(1),xmin_global(2),xmin_global(3),xmin_global(4));
    flat = flat + DM_Map;    
end

%for measuring phase effects

% xmin_global = [0,0,5,0];
% opts_coarse = optimoptions(@fmincon,'Algorithm','interior-point','StepTolerance',.05,...
%     'DiffMinChange', 0.1);
% problem_coarse = createOptimProblem('fmincon','x0',[20, ang_init,act_init,pha_init],...
%     'objective',@DE_supMinimizer_3, 'lb',[10, -.3, 3.0, 0.0], 'ub',[70, .3, 10.0, 6.2], 'options',opts_coarse);
% 
% ms = MultiStart('FunctionTolerance',0,'XTolerance',.05, 'Display','iter',...
%     'StartPointsToRun','bounds', 'OutputFcn', @StopIfReached,'PlotFcn',@gsplotbestf);
% 
% 
% 
% [xmin_coarse, fmin_coarse] = run(ms,problem_coarse,72);
% 
% opts_fine = optimoptions(@fmincon,'StepTolerance',.05,'DiffMinChange',.01);
% problem_fine = createOptimProblem('fmincon','x0', xmin_global,'objective',@DE_supMinimizer_3,...
%     'lb', [xmin_global(1)-2, xmin_global(2)-.05, xmin_global(3)-.5, xmin_global(4)-.1], ...
%     'ub', [xmin_global(1)+2, xmin_global(2)+.05, xmin_global(3)+5, xmin_global(4)+.1], ...
%     'options', opts_fine);
% 
% [xmin_fine, fmin_fine] = run(ms,problem_fine,72);
% 
% xmin = xmin_global;
% fmin = fmin_global;
% 
% n = 1;
% for phase = 0:.1:6.2
%     x_input = [xmin_global(1),xmin_global(2),xmin_global(3),phase];
%     
%     phase_sweep_data(n) = DE_supMinimizer_3(x_input);
%     n = n+1;
% end


%now minimize with amplitude

%     [Par, Pmin] = patternsearch(@DE_supMinimizer_3, par0, [], [], [], [], ... 
%          LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));

%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
%Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

%Get index of minimized value
% ind = find(Pdat == Pmin, 1);

%Measure at Minimum
% for i=1:200
%     par0 = [ParT(ind,1), ParT(ind, 2), ParT(ind,3), ParT(ind,4)];
%     Power=DE_supMinimizer_2(par0);
% end

%Calculate Resulting Power Specs
% OptP    = Pmin;
% OptPSup = Pflat/Pmin;

% %Match image to optimum power, show, and save as fits
% OptIm   = Idat(:,:,ind);
% figure()
% imshow(OptIm);
% caxis([0 255]);
% fitswrite(OptIm, OptImnm);
% 
% %Save all images as single cube
% fitswrite(Idat, Idatnm);
% 
% %Save PM values as vector
% PMtxt = fopen(PMVnm, 'wt');
% fprintf(PMtxt, '%d\n', Pdat);
% fclose(PMtxt);
% 
% %_______Save Raw Results___________
% REtxt = fopen(Rawnm, 'wt');
% fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
% % fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
% fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);
% % Save results for each iteration
% for i = 2: icount
%     fprintf(REtxt, 'Iteration: %d     ', i);
%     fprintf(REtxt, 'Parameters: Amp. = %09.5f, Angle = %09.5f, Act. = %09.5f, Phase = %09.5f \n', ...
%         ParT(i,1), ParT(i, 2), ParT(i,3), ParT(i,4));
%     fprintf(REtxt, '       Power: %05.2e      Suppression Ratio: %05.2f \n \n', Pdat(i), Pflat/Pdat(i));
% end
% fprintf(REtxt, '\n \nOptimum Values: Iteration#%d       ', ind);
% fprintf(REtxt, 'Parameters: Amp. = %09.5f, Angle = %09.5f, Act. = %09.5f, Phase = %09.5f \n', ...
%         ParT(ind,1), ParT(ind, 2), ParT(ind,3), ParT(ind,4));
% fprintf(REtxt, 'Power: %05.2e       Suppression Ratio: %05.2f', OptP, OptPSup);
% fclose(REtxt);

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
newp.CloseDevices();

%Close APD
% fclose(apd);
% delete(apd);
% clear apd;

%Close camera
% if ~strcmp(char(cam.Exit), 'SUCCESS')
%     fprintf('Could not close camera');
% else
%     fprintf('Successfully Close CCD\n');
% end
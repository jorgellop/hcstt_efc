%{
Speckling Nulling Optimization: 
- Place Speckles alfkj
***_2 averages 500 PM points in an attempt to decrease noise
***Uses patternsearch() as optimizer instead of fmincon()

***IF CODE FAILS IN SUPMINIMIZER_2() DUE TO INDEXING. POWER CYCLE PM
        Most likely cause by an error in reading from the PM. Can confirm
        via the Log. This error is most easily solved via a power cycle but
        if cycle is impossible. Manually perform several PM reads until PM
        comms are back to normal.
*** Version used to capture images of data for general storage.
    DOES NOT ITERATE. DO NOT TRY TO ITERATE. USE _2

******************************************************
- Arguments:
    NONE        USER MUST SET VARIOUS PARAMETERS:
      fldir     Folder to save data in
      lgnm      Update date at end of string
      h0:q      Initial conditions on DM shape
      UB:LB     Upper and Lower bounds on optimization region
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


Compiled By:    Daniel Echeverri, Jerry Xuan and Michael Randolph
Last Modified:  08/23/2016
%}
% close all
% clear

global icount REtxt Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond

%% DEFINE OUTPUT DIRECTORY AND FILE; DEFINE SETUP PARAMETERS
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\Gary\Desktop\SpeckleNullingTests\GR_111416_2\';
%Set directory to folder with cubes for analysis
% cd(fldDir);

%Name of TXT file for PM values
PMVnm   = [fldDir,'PM_Values.txt'];

%Name of TXT file for Raw Results
Rawnm   = [fldDir,'RawResults.txt'];

%Name of FITS file for Optimum Image
OptImnm = [fldDir,'OptImage.fits'];

%Name of Fits file for Image cube
Idatnm  = [fldDir,'ImageTrials.fits'];

%Number of Samples for PM data set; MUST BE AN INTEGER
Samp    = 500;

%Exposure time for camera
CExp    = 25;%ms

%% DM Setup
% Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);



%% Define Initial Parameters
Nactacross = 10;    %Number actuators across any given length
                        %VALUE IS AN ESTIMATE; NEED TO CONFIRM

%SET Initial DM SIN PROPERTIES HERE: 
h0  = 100;           %Vary amplitude: 10 to 100, step of 10
x   = 4;         %Fix Distance: Set to meet given speckle
alp = 3.1;         %Vary Phase: 0 to 2pi, step of pi/10
q   = 1.6;         %Fix Angle: Set to meet given speckle

%Save initial condition as vector for optimization
par0 = [h0, q, x, alp];

%SET PROPERTIES FOR SEARCH
%LB = [50.0, -0.1, 3.0, 0.0];              %Lower bounds; same elements as par0
%UB = [150.0, 0.1, 11.0, 2*pi];       %Upper bounds; same elements as par0
% sTol = .0000001;                     %Smallest step size
% fTol = 3e-9;                   %Stops iterating when delPower < fTol

DE_DMW_Sin(h0, q, x, alp, drv_inf);

% testing = zeros(12,12);
% testing(5,5) = 100;

%testing = DM_voltage_map;
%DE_DMW_Map(testing, drv_inf);
pause(.05);


%% Disconnect all Devices
%Return to MATLAB Directory
% cd('C:\Users\Gary\Documents\MATLAB');

pause(.05);
pause(.05);

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end

%Close PM
% newp.CloseDevices();
% 
% %Close camera
% if ~strcmp(char(cam.Exit), 'SUCCESS')
%     fprintf('Could not close camera');
% else
%     fprintf('Successfully Close CCD\n');
% end
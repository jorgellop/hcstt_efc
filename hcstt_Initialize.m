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

function hcstt_Initialize(setup_PM)

global cam img Xcr Ycr CExp s drv_inf flat himg pm_scale

load('flat.mat');
load('FlatMap.mat');
load('angle_cal');
load('position_cal');

% Camera Setup
disp('Setting up Camera');
hcstt_InitializeCamera()

% DM Setup
% % Open and initialize Mulit-DM driver USB connection
disp('Setting up DM');
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);

%Flatten DM

%Create StringBuilder to retrieve and hold data; Set cap. to expected size
sb = System.Text.StringBuilder();
sb.Capacity = 4100;

%Create Matrix to hold images: Width(200)xHeight(200)x#Images            
    %Images have been cropped to reduce cube size
% Idat = zeros(400,400,itr);    %Raw cropped images
%Create column vector to hold Power data
% Pdat = zeros(itr, 1);
%Create Matrix to hold trial parameters
% ParT = zeros(itr, 4);
% icount = 1;
% REtxt = fopen('temp.txt', 'wt');

%Subtract 20 from Pdat for easy identification of empty elements later
    %In case fminsearch converges <itr iterations. Do not  want to process
    %or save excess images later
% Pdat        = Pdat - 20;

% %Reflatten cube for Flat image comparison
DM_Command = zeros(12);
FlatCommand = NK_MultiDM_Command(DM_Command, flat);
JR_UPDATE_MultiDM(drv_inf, FlatCommand);


if(setup_PM)
    % PM Setup
    disp('Setting up Powermeter');
    s = daq.createSession('ni');
%     addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');
    addAnalogInputChannel(s, 'SimDev1', 0, 'Voltage');

%     addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');
    addDigitalChannel(s, 'SimDev1', 'Port0/Line0:4', 'OutputOnly');
    outputSingleScan(s, [0, 0, 1,1,0]);
    
    % Setup Femto power meter
    s = daq.createSession('ni');
%     addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
    addAnalogInputChannel(s,'SimDev1', 0, 'Voltage');

    s.Rate = 10;
    s.DurationInSeconds = 1/10;

    %starting conversion scale (in V/W)
    pm_scale = 10^9;
end
end

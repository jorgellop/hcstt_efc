close all
addpath(genpath('utils'));
load('flat.mat');
load('FlatMap.mat');

%% DM Setup

% Open and initialize Mulit-DM driver USB connection
mapping_ID = 6;
[err, drv_inf] = OPEN_multiDM(mapping_ID);



%% Flatten DM

DM_Command = zeros(12);
FlatCommand = NK_MultiDM_Command(DM_Command, flat);
JR_UPDATE_MultiDM(drv_inf, FlatCommand);

%% Pause

w = waitforbuttonpress;

%% Disconnect DM

error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end


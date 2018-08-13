%{
This basic function allow to take an image with a thorlabs detector.
******************************************************
- Arguments: 
  * Exp_time  = Exposure time in miliseconds.

- Returns: 
  * Im   = Image obtained with the detector (Array of double).

- Dependencies:
  * NET.Assembly : structure of function to control the detector.
******************************************************

Compiled By:    JR Delorme
Last Modified:  March 06, 2017
%}

function hcstt_DisconnectDM()

global cam img Xcr Ycr CExp s drv_inf flat

disp('Disconnecting DM')

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end

% Close Power Meter
% newp.CloseDevices();


%Close APD
% fclose(apd);
% delete(apd);
% clear apd;

%Close camera
end
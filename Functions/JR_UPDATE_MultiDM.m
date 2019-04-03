%{
This function apply a command to the DM
- Values capped at 90% of the maximum voltage.

******************************************************
- Arguments:
  drv_inf  = Information about the MultiDM
  command  = Command send to the DM in % of the maximum voltage

- Returns: 
  NONE

- Dependencies:
  NONE
******************************************************

Compiled By:    JR Delorme
Last Modified:  Fabruary 2, 2017
%}

function JR_UPDATE_MultiDM(drv_inf, command)
% Values capped at 90% of the maximum voltage
if sum(command >= 90) > 0
    fprintf('Warning -- JR_MultiDM_Command -- Voltage higher than 90/100 of the maximum voltage.\n')
end
command(command >= 90) = 90;

% Apply command to the DM
[error_code] = UPDATE_multiDM(drv_inf, command);
if error_code == -4
    fprintf('ERROR -- UPDATE_MultiDM -- Data send to driver failed.\n');
elseif error_code == -6
    fprintf('ERROR -- UPDATE_MultiDM -- USB_ID and USB_pointer values unrecognized.\n');
end

end
%{
getIntensityFiberOutput

-
 
******************************************************
- Arguments:
    us          = array with the heights of the DM shape to update
   
- Returns:
    pm_reading         = powermeter reading
******************************************************

Compiled By:    Jorge Llop
Last Modified:  
%}

function pm_reading = hcstt_getIntensityFiberOutput(us, s, pm_exponent, drv_info)

% Write new values to DM
UPDATE_multiDM(drv_info, us);
pause(.05);

pm_reading_raw = s.inputSingleScan;
pm_scale = 10^input(pm_exponent);
pm_reading = pm_reading_raw/pm_scale;

end

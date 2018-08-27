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

function [int] = hcstt_GetIntensityFIUvNoGainUpdate(us,num_avg)

global cam img Xcr Ycr CExp s drv_inf flat itr Idat himg pm_scale


hcstt_UpdateMultiDM(us)
pause(0.2)

for II = 1:num_avg
    int_arr(II) = s.inputSingleScan/pm_scale;
    pause(0.01)
end
int = mean(int_arr);


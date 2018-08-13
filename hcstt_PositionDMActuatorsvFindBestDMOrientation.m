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

function [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,posII)

load('DMAlign_PosActuators_Jul27')
load('DMAlign_RadCenter_Jul27')

scaleR = 0.95; % Correction to the radius

R = R*scaleR;
x = (pos_x-x_c)/R;
y = (pos_y-y_c)/R;

posDM_x0 = x*apRad;
posDM_y0 = y*apRad;

for II=1:12-1
    diff_x(II) = abs(posDM_x0(II+1)-posDM_x0(II));
    diff_y(II) = abs(posDM_y0(II+1)-posDM_y0(II));
end
ac_spac = mean([mean(diff_x),mean(diff_y)]);
%Sort the way it has been observed by comparing the model with the bench
%results
if posII ==1
    posDM_x = sort(posDM_x0,'descend')+N/2;
    posDM_y = sort(posDM_y0,'descend')+N/2;
elseif posII ==2
%     posDM_x = sort(posDM_x0,'ascend')+N/2;
%     posDM_y = sort(posDM_y0,'ascend')+N/2;
    posDM_x = sort(posDM_y0,'descend')+N/2;
    posDM_y = sort(posDM_x0,'descend')+N/2;
elseif posII ==3
    posDM_x = sort(-posDM_x0,'descend')+N/2;
    posDM_y = sort(-posDM_y0,'descend')+N/2;
elseif posII ==4
%     posDM_x = sort(-posDM_x0,'ascend')+N/2;
%     posDM_y = sort(-posDM_y0,'ascend')+N/2;
    posDM_x = sort(-posDM_y0,'descend')+N/2;
    posDM_y = sort(-posDM_x0,'descend')+N/2;
elseif posII ==5
    posDM_x = sort(-posDM_x0,'descend')+N/2;
    posDM_y = sort(posDM_y0,'descend')+N/2;
elseif posII ==6
    posDM_x = sort(posDM_x0,'descend')+N/2;
    posDM_y = sort(-posDM_y0,'descend')+N/2;
elseif posII ==7
    posDM_x = sort(posDM_y0,'descend')+N/2;
    posDM_y = sort(-posDM_x0,'descend')+N/2;
elseif posII ==8
    posDM_x = sort(-posDM_y0,'descend')+N/2;
    posDM_y = sort(posDM_x0,'descend')+N/2;

end
end

%{
hcstt_PositionDMActuatorsvFindBestR

Multiply the radius of the pupil found when imaging the pupil with the
camera to find a better R to match the sinusoid images
%}

function [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestR(N,apRad,scaleR)

load('positionsDMActuators_raw_Feb13v2')
load('alignDMModel_pupilCenterRadius_Feb13')

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
posDM_x = sort(posDM_x0,'descend')+N/2;
posDM_y = sort(posDM_y0,'ascend')+N/2;
end

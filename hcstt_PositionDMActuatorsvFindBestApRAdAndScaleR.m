%{

%}

function [posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestApRAdAndScaleR(N,apRad,scaleR,posII)

load('positionsDMActuators_raw_Feb13v2')
load('alignDMModel_pupilCenterRadius_Feb13')

% scaleR = 0.8509; % Correction to the radius

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
% if posII ==1
%     posDM_x = sort(posDM_x0,'ascend')+N/2;
%     posDM_y = sort(posDM_y0,'ascend')+N/2;
% elseif posII ==2
%     posDM_x = sort(posDM_x0,'ascend')+N/2;
%     posDM_y = sort(-posDM_y0,'ascend')+N/2;
% elseif posII ==3
%     posDM_x = sort(-posDM_x0,'ascend')+N/2;
%     posDM_y = sort(posDM_y0,'ascend')+N/2;
% elseif posII ==4
%     posDM_x = sort(-posDM_x0,'ascend')+N/2;
%     posDM_y = sort(-posDM_y0,'ascend')+N/2;
if posII ==1
    posDM_x = sort(posDM_x0,'descend')+N/2;
    posDM_y = sort(posDM_y0,'descend')+N/2;
elseif posII ==2
    posDM_x = sort(-posDM_x0,'descend')+N/2;
    posDM_y = sort(-posDM_y0,'descend')+N/2;
elseif posII ==3
    posDM_x = sort(posDM_x0,'descend')+N/2;
    posDM_y = sort(-posDM_y0,'descend')+N/2;
elseif posII ==4
    posDM_x = sort(-posDM_x0,'descend')+N/2;
    posDM_y = sort(posDM_y0,'descend')+N/2;
% elseif posII ==9
%     posDM_x = sort(posDM_x0,'descend')+N/2;
%     posDM_y = sort(posDM_y0,'descend')+N/2;
% elseif posII ==10
%     posDM_x = sort(posDM_x0,'descend')+N/2;
%     posDM_y = sort(-posDM_y0,'descend')+N/2;
% elseif posII ==11
%     posDM_x = sort(-posDM_x0,'descend')+N/2;
%     posDM_y = sort(posDM_y0,'descend')+N/2;
% elseif posII ==12
%     posDM_x = sort(-posDM_x0,'descend')+N/2;
%     posDM_y = sort(-posDM_y0,'descend')+N/2;
end
end

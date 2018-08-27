% test_plotDMConfigs.m
%
% plot in a scatter graph the different DM registrations
%
% Jorge Llop - Aug 24, 2018

% close all;
clear all;

load('DMAlign_PosActuators_Jul27')
load('DMAlign_RadCenter_Jul27')

apRad = 142;
N = 1024;
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

% [x,y]=meshgrid(posDM_x0,posDM_y0);
% scatter(x(:),y(:))
% hold on
figure(1)
[x,y]=meshgrid(sort(posDM_y0,'descend'),sort(posDM_x0,'descend'));
scatter(x(:),y(:))
hold on
[x,y]=meshgrid(sort(-posDM_x0,'descend'),sort(-posDM_y0,'descend'));
scatter(x(:),y(:))
[x,y]=meshgrid(sort(-posDM_x0,'descend'),sort(posDM_y0,'descend'));
scatter(x(:),y(:))
[x,y]=meshgrid(sort(posDM_y0,'descend'),sort(-posDM_x0,'descend'));
scatter(x(:),y(:))
hold off
legend('2','3','5','7')
axis image

figure(6)
[x,y]=meshgrid(sort(posDM_y0-1.5,'descend'),sort(posDM_y0,'descend'));
scatter(x(:),y(:))
hold on
% [x,y]=meshgrid(sort(posDM_y0+2,'descend'),sort(posDM_y0,'descend'));
% scatter(x(:),y(:))
% [x,y]=meshgrid(sort(posDM_y0+3,'descend'),sort(posDM_y0,'descend'));
% scatter(x(:),y(:))
[x,y]=meshgrid(sort(posDM_y0+1.5,'descend'),sort(posDM_y0,'descend'));
scatter(x(:),y(:))
% [x,y]=meshgrid(sort(-posDM_x0-1.5,'descend'),sort(-posDM_x0,'descend'));
% scatter(x(:),y(:))
% [x,y]=meshgrid(sort(-posDM_x0+1.5,'descend'),sort(-posDM_x0,'descend'));
% scatter(x(:),y(:))
% [x,y]=meshgrid(sort(posDM_y0,'descend'),sort(-posDM_x0,'descend'));
% scatter(x(:),y(:))
hold off
% legend('2','3','5','7')
axis image

% if posII ==1
%     posDM_x = sort(posDM_x0,'descend')+N/2;
%     posDM_y = sort(posDM_y0,'descend')+N/2;
% elseif posII ==2
% %     posDM_x = sort(posDM_x0,'ascend')+N/2;
% %     posDM_y = sort(posDM_y0,'ascend')+N/2;
%     posDM_x = sort(posDM_y0,'descend')+N/2;
%     posDM_y = sort(posDM_x0,'descend')+N/2;
% elseif posII ==3
%     posDM_x = sort(-posDM_x0,'descend')+N/2;
%     posDM_y = sort(-posDM_y0,'descend')+N/2;
% elseif posII ==4
% %     posDM_x = sort(-posDM_x0,'ascend')+N/2;
% %     posDM_y = sort(-posDM_y0,'ascend')+N/2;
%     posDM_x = sort(-posDM_y0,'descend')+N/2;
%     posDM_y = sort(-posDM_x0,'descend')+N/2;
% elseif posII ==5
%     posDM_x = sort(-posDM_x0,'descend')+N/2;
%     posDM_y = sort(posDM_y0,'descend')+N/2;
% elseif posII ==6
%     posDM_x = sort(posDM_x0,'descend')+N/2;
%     posDM_y = sort(-posDM_y0,'descend')+N/2;
% elseif posII ==7
%     posDM_x = sort(posDM_y0,'descend')+N/2;
%     posDM_y = sort(-posDM_x0,'descend')+N/2;
% elseif posII ==8
%     posDM_x = sort(-posDM_y0,'descend')+N/2;
%     posDM_y = sort(posDM_x0,'descend')+N/2;
% 
% end

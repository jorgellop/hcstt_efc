close all;
clear all;

%{
hcstt_DMAlign_FindPosActuators

Finds position of actuators from the images selected with
hcstt_DMAlign_FindPosActuators.m

Jorge Llop - Dec 2017
%}

load('alignDMModel_ind_ma_arrv2');
load('alignDMModel_ind_arrv2');
load('im_matv2');

Nact = 12;

pos_x = zeros(Nact,Nact);
pos_y = zeros(Nact,Nact);

[I1,I2] = ind2sub(size(im_mat(:,:,1)),ind_ma_arr);

for II = 1:numel(ind_arr)
    pos_x(ind_arr(II)) = I1(II);
    pos_y(ind_arr(II)) = I2(II);
end
    
act_arr = 1:Nact;
mean_x = zeros(Nact,1);
mean_y = zeros(Nact,1);
for II = 1:Nact
    if numel(find(pos_x(II,:))) > 3
        pos_xII = pos_x(II,:)';
        mean_x(II) = mean(pos_xII(find(pos_xII)));
    end
    if numel(find(pos_y(:,II))) > 3
        pos_yII = pos_y(:,II);
        mean_y(II) = mean(pos_yII(find(pos_yII)));
    end
    
end

act_x = find(mean_x);
act_y = find(mean_y);
mean_x = mean_x(find(mean_x));
mean_y = mean_y(find(mean_y));
pos_x_mean = interp1(act_x,mean_x,act_arr,'linear','extrap')
pos_y_mean = interp1(act_y,mean_y,act_arr,'linear','extrap')
pos_x = pos_x_mean;
pos_y = pos_y_mean;
save('positionsDMActuators_raw','pos_x','pos_y')

% act_arr = 1:Nact;
% m_x = zeros(Nact,1);
% m_y = zeros(Nact,1);
% for II = 1:Nact
%     if numel(find(pos_x(:,II))) > 3
%         pos_xII = pos_x(:,II)';
%         pfit_x = polyfit(act_arr(find(pos_xII)),pos_xII(find(pos_xII)),1);  
%         m_x(II) = pfit_x(1);
%         p0_x(II) = pfit_x(2);
%     end
%     if numel(find(pos_y(II,:))) > 3
%         pos_yII = pos_y(II,:);
%         pfit_y = polyfit(act_arr(find(pos_yII)),pos_yII(find(pos_yII)),1);
%         m_y(II) = pfit_y(1);
%         p0_y(II) = pfit_y(2);        
%     end
%     
% end

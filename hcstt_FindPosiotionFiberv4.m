% hcstt_FindPosiotionFiberv4
% 
% Find the position of the fiber wrt to the center of the star. It requires
% a first estimate of the postition, from which this code will scan around the
% neighbourhood and find the maximum FIU output power
%
% Nikita Klimovich, Jorge Llop - Dec 14, 2017

function [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(actxc_est,ang_est,info)
% Set lower and upper bounds of range to look
num_actxc = 7;
num_ang = 7;
actxc_arr = linspace(actxc_est-0.075,actxc_est+0.075,num_actxc);
ang_arr = linspace(ang_est-0.105,ang_est+0.105,num_ang);

% Initialize outputs
fiber_qx = [0 0];
best_amplitude = 0;
ph_arr = linspace(0,3*pi/2,3);

mat = zeros(num_actxc,num_ang);
%Search over the lower and upper bounds
for II=1:num_actxc
    actxc = actxc_arr(II);
    for JJ=1:num_ang
        ang = ang_arr(JJ);
        for KK=1:3
            DM_Command = hcstt_DMMapSin(70,ang,actxc,ph_arr(KK));
    %         amp_max = DE_supMinimizer_3(par0);
            amp_max = hcstt_GetIntensityFIU(DM_Command,3,info.backgroundSMF);
            mat(II,JJ) = mat(II,JJ)+amp_max;
        end
        amp_max = mat(II,JJ);
        % Update fiber location
        if (amp_max > best_amplitude)
            fiber_qx = [II JJ];
            best_amplitude = amp_max;
        end
    end
end

figure(1000)
imagesc(mat)
set(gca, 'XTickLabel', ang_arr)
set(gca, 'YTickLabel', actxc_arr)
axis image
colorbar

% fiber_qx = [0 5];
disp('Fiber Location Search Complete')
actxc_fib = actxc_arr(fiber_qx(1));
ang_fib = ang_arr(fiber_qx(2));
best_amplitude=best_amplitude;
end

% hcstt_Initialize(true);
% [actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4(actxc_est,ang_est)
% hcstt_DisconnectDevices();


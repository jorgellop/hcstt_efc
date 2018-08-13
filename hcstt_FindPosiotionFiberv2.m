% hcstt_FindPosiotionFiberv2
% 
% Find the position of the fiber wrt to the center of the star. It requires
% a first estimate of the postition, from which this code will scan around the
% neighbourhood and find the maximum FIU output power
%
% Nikita Klimovich, Jorge Llop - Dec 14, 2017

function [x_fib,y_fib] = hcstt_FindPosiotionFiberv2(x_fib0,y_fib0)
% Set lower and upper bounds of range to look
LB = [1.0, -0.1, 2.5, 0.0];              %Lower bounds; same elements as par0
UB = [5.0, 0.2, 3.0, 2*pi];       %Upper bounds; same elements as par0

% Initialize outputs
fiber_qx = [0 0];
best_amplitude = 0;
best_phase = 0;

%Search over the lower and upper bounds
for i=LB(2):0.025:UB(2)
    for j=LB(3):0.1:UB(3)
%         par0 = [50 i j 0];
        DM_Command = hcstt_DMMapSin(100,i,j,0);
%         amp_max = DE_supMinimizer_3(par0);
        amp_max = hcstt_GetIntensityFIU(DM_Command,10);
        % Update fiber location
        if (amp_max > best_amplitude)
            fiber_qx = [i j];
            best_amplitude = amp_max;
        end
    end
end

% fiber_qx = [0 5];
disp('Fiber Location Search Complete')
disp(fiber_qx(1))
disp(fiber_qx(2))

x_fib = fiber_qx(1);
y_fib = fiber_qx(2);
best_amplitude
end
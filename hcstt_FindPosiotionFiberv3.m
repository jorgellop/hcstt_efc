% hcstt_FindPosiotionFiberv2
% 
% Find the position of the fiber wrt to the center of the star. It requires
% a first estimate of the postition, from which this code will scan around the
% neighbourhood and find the maximum FIU output power
%
% Nikita Klimovich, Jorge Llop - Dec 14, 2017

function [x_fib,y_fib] = hcstt_FindPosiotionFiberv3(x_fib0,y_fib0)
% Set lower and upper bounds of range to look
disp('Searching for the position of the fiber')
% Initialize outputs
best_amplitude = 0;

apRad2 = 12;
poke_amp = 80e-9;
 
range = 0.5;
incr = 0.01;
x_fib_low = x_fib0 - range/2;
x_fib_up = x_fib0 + range/2;
%Search over the lower and upper bounds
% for i=LB(2):0.025:UB(2)
%     for j=LB(3):0.1:UB(3)
% %         par0 = [50 i j 0];
%         DM_Command = hcstt_DMMapSin(100,i,j,0);
% %         amp_max = DE_supMinimizer_3(par0);
%         amp_max = hcstt_GetIntensityFIU(DM_Command,10);
%         % Update fiber location
%         if (amp_max > best_amplitude)
%             fiber_qx = [i j];
%             best_amplitude = amp_max;
%         end
%     end
% end
for i=x_fib_low:incr:x_fib_up
    cosfct = cos(2*pi*[1:apRad2]/(apRad2) * i ) ;
    a = ones(apRad2,apRad2);
    di = diag(cosfct);
    us = a * di; 
    
    dm_actuators_mat = us' * poke_amp;
    amp_max = hcstt_GetIntensityFIU(+dm_actuators_mat/1e-9,10);
    if (amp_max > best_amplitude)
        best_amplitude = amp_max;
        x_fib = i;
        y_fib = y_fib0;
    end

end
% fiber_qx = [0 5];
disp('Fiber Location Search Complete')
disp(x_fib)
disp(y_fib)

best_amplitude
end
function [actxc_fib,ang_fib] = hcstt_test_getSpatialFreqOfFiber(II)
% Set lower and upper bounds of range to look
num_actxc = 13;
num_ang = 7;
actxc0 = 2.5;
ang0 = 0;
[actxc_arr,ang_arr] = meshgrid(linspace(actxc0-0.11,actxc0+0.11,num_actxc),linspace(ang0-0.115,ang0+0.115,num_ang));
% ang_arr = meshgrid(linspace(ang0-0.115,ang0+0.115,num_ang));
% [aux,ind] = sort(abs(actxc_arr(:)));
% actxc_arr = actxc_arr(ind);
% ang_arr = ang_arr(ind);

actxc_fib = actxc_arr(II);
ang_fib = ang_arr(II);
end
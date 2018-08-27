%%% femto_char.m
%%% determine how different gain levels in the femto read the same signal
%
% this code initializes the digital control for the femto,
% then reads in a stable intensity for 100 scans and saves to an array,
% then switches to the next gain level and repeats with the same intensity,
% and saves both arrays to a file to be analyzed later
%
% Coded in the context of the tests for EFCSMF.
%
%%% Ben Calvin - Aug 23, 2018

s = daq.createSession('ni');

% Define the levels (remember, read as [pin12, pin11, pin10])
lev1 = [0 0 0];
lev2 = [0 0 1];
lev3 = [0 1 0];
lev4 = [0 1 1];
lev5 = [1 0 0];
lev6 = [1 0 1];
lev7 = [1 1 0];
levelList = [lev1; lev2; lev3; lev4; lev5; lev6; lev7];

% Signal in the higher gain
addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');
outputSingleScan(s, [fliplr(lev6), 1, 0])

% Remove the digital channels so we can specify s.Rate and s.Duration
s.removeChannel(1:5)

addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
s.Rate = 10;
s.DurationInSeconds = 1/10;

% Make the array of the lower gain readings first
int1_arr=zeros(size(1:250));
for iii = 1:250
    int1_arr(iii)= s.inputSingleScan;
end

disp(mean(int1_arr))
disp(std(int1_arr))

% When done, remove the analog channel
s.removeChannel(1)

%Signal in the next step higher gain
addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');
outputSingleScan(s, [fliplr(lev7), 1, 0])

s.removeChannel(1:5)

% Same as before
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
s.Rate = 10;
s.DurationInSeconds = 1/10;

% Make the array for the higher gain 
int2_arr=zeros(size(1:250));
for jjj = 1:250
    int2_arr(jjj)= s.inputSingleScan;
end

disp(mean(int2_arr))
disp(std(int2_arr))

% Save them (the L and the 'one' look unfortunately similar)
save('levels6-7.mat', 'int1_arr', 'int2_arr')
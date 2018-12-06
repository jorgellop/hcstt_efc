%{
%}

function [int] = hcstt_GetIntensityFIU(us,num_avg,backgroundCam)

global cam img Xcr Ycr CExp s drv_inf flat itr Idat himg pm_scale

%Load the gain calibration
% load(['output',filesep,'NIPhotodiodeGainCalibration'])
gain_arr = ones(1,6)*9.97;
%Update DM with us, ie the shape
hcstt_UpdateMultiDM(us)

level0 = 2; % The gain is calibrated wrt level 2
pause(0.05)

for II = 1:num_avg
    if II==1
        reading = s.inputSingleScan;
        pm_scale=hcstt_getPMNewLevel(reading,pm_scale);
        level = log10(pm_scale) - 4;
        gain = prod(gain_arr(level0:level-1))/10^(numel(gain_arr(level0:level-1)));
        int_arr(II) = s.inputSingleScan/pm_scale*gain-backgroundCam;
    else
        int_arr(II) = s.inputSingleScan/pm_scale*gain-backgroundCam;
    end
    pause(0.1)
end
int = mean(int_arr);
% if pm_scale<0.01


%%% I don't know how MATLAB networking works, so I'm just writing the
%%% function down here -Ben

function scale=hcstt_getPMNewLevel(reading, old_scale)
%%% Assumes that there are only the first three digital channels set up and
%%% that switching the H/L or AC/DC is not needed.

    global cam img Xcr Ycr CExp s drv_inf flat itr Idat himg pm_scale
    highVolt=7.5;
    if pm_scale==1e11
        lowVolt=0.03;
    else
        lowVolt=.075;
    end
    scale = old_scale;
    while reading>highVolt || reading<lowVolt
        
        
        level = log10(scale) - 4;
        disp(['Changing Level, current level: ', num2str(level),', reading: ',num2str(reading)])

        s = daq.createSession('ni');
        addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');

        addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');

        % Read pins as [pin12, pin11, pin10], use fliplr when outputScaning
        if level==1
            lev = [0 0 0];
        elseif level==2
            lev = [0 0 1];
        elseif level==3
            lev = [0 1 0];
        elseif level==4
            lev = [0 1 1];
        elseif level==5
            lev = [1 0 0];
        elseif level==6
            lev = [1 0 1];
        elseif level==7
            lev = [1 1 0];
        end
        outputSingleScan(s, [fliplr(lev),1,0]);
        
        lev1 = [0 0 0];
        lev2 = [0 0 1];
        lev3 = [0 1 0];
        lev4 = [0 1 1];
        lev5 = [1 0 0];
        lev6 = [1 0 1];
        lev7 = [1 1 0];
      
        levelList = [lev1; lev2; lev3; lev4; lev5; lev6; lev7];

        %%% I'm defining the voltage reading to be accurate between
        %%% 0.075 and 7.5 V. Outside these, the level wants to switch
        if reading>highVolt
            if not(level==1)
                level = level-1;
                scale = power(10, (level+4));
                outputSingleScan(s, [fliplr(levelList(level, :)),1,0]);
            else
                disp('Power too high.')
            end

        elseif reading<lowVolt
            if not(level==7)
                level = level+1;
                scale = power(10, (level+4));
                outputSingleScan(s, [fliplr(levelList(level, :)),1,0]);
            else
                disp('Gain is at maximum, scale=10^11, and reading is not high enough')
                break
            end
        end
        s = daq.createSession('ni');
        addAnalogInputChannel(s,'Dev1', 0, 'Voltage');

        s.Rate = 10;
        s.DurationInSeconds = 1/10;

        reading = s.inputSingleScan;
        disp(['Level changed to: ', num2str(level),', reading: ',num2str(reading)])
    end

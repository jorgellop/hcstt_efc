%%% femto_auto.m
%%% automatic switching of the femto gain switch
% this code initializes the digital control for the femto,
% reads power from the femto power meter,
% then decides if the gain setting should be changed.
%
% in the case it doesn't need to be changed, we read the power again.
% in the case it does need to be changed, increase
%
% Coded in the context of the tests for EFCSMF.
%
%%% Ben Calvin - Aug 17, 2018

%%                     INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% defining the gain level in the newLevel function
level=1; %starts at the least sensitive and will min out first
pm_scale = 10^5;

%%% ACDC 
% this shouldn't need to change. It's the supply voltage regulator
% true means DC and false means AC
global ACDC s
ACDC = 1;

%%% Gain amplifier
% modes are 'high speed' (0) and 'low noise' (1)
% high speed is for when the signal needs to be boosted further
noise = 0;


%%                        SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = daq.createSession('ni');
addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');

addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');

outputSingleScan(s, [0,0,0,1,0]);


%%                        LIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
user_cont=1;
highVolt=7.5;
lowVolt=0.075;

while(user_cont)
    for iii=1:120
        pause(.4)
        reading = s.inputSingleScan;
        
        %%% I'm defining the voltage reading to be accurate between
        %   0.075 and 7.5 V. Outside these, the level gets switched
        if reading>highVolt
            if not(level==1)
                level = level-1;
                pm_scale = newLevel(level, noise);
            elseif level==1 && noise==0
                noise=1;
                level = level+1;
                pm_scale = newLevel(level, noise);
            else
                disp('Power too high. Pausing')
                break
            end
        
        elseif reading<lowVolt
            if level==7 && noise==0
                
            elseif level==7 && noise==1
                noise=0;
                level = level-1;
                pm_scale = newLevel(level, noise);
            else
                level = level+1;
                pm_scale = newLevel(level, noise);
            end
        end
        
        power = reading/pm_scale;
        disp([num2str(power), ' Watts  ', num2str(reading), ' Volts'])
        
    end
    user_cont = user_cont-1;
    if user_cont==0
        user_cont=input('Enter how many more minutes you want: ');
    end
end


%%                 SWITCH THE GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pm_scale=newLevel(level, noise)

% This function determines the new power meter gain scaling and
% alters the signal sent to the femto accordingly

    global ACDC s;
    
    disp('Changing Level')
    
    %%% defining the gain levels
    lev1 = [0 0 0];
    lev2 = [0 0 1];
    lev3 = [0 1 0];
    lev4 = [0 1 1];
    lev5 = [1 0 0];
    lev6 = [1 0 1];
    lev7 = [1 1 0];

    levelList = [lev1; lev2; lev3; lev4; lev5; lev6; lev7];
    currentlevel = levelList(level, :);

%%  Determine new gain scaling
    
    pm_scale = power(10, (level+2));
    
    if ~noise
        pm_scale = pm_scale*100;
    end
    
%%  CHANGE THE CABLE INPUT
    pin10 = currentlevel(3);
    pin11 = currentlevel(2);
    pin12 = currentlevel(1);
    pin13 = ACDC;
    pin14 = noise;
    
    %%% send the signal through the cable
    outputSingleScan(s,[pin10, pin11, pin12, pin13, pin14])

end


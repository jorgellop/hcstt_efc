close all

global icount flat REtxt apd Rawnm Pdat Idat ParT Pflat sb newp cam drv_inf img Xcr Ycr Samp currentTime currentSecond startingSecond
global s pm_scale

%% DEFINE OUTPUT DIRECTORY AND FILE; DEFINE SETUP PARAMETERS
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\COO_user\Desktop\HCSTT\Results\NK_broadband\';

%Name of TXT file for Raw Results
txtName = ['BandPassNulling' datestr(datetime,30) '.txt'];
Rawnm   = [fldDir,txtName];

%Name of Fits file for Image cube
txtName = ['ImageTrials' datestr(datetime,30) '.fits'];
Idatnm  = [fldDir,txtName];

%% Initialize System

hcstt_Initialize(1);

%% Define Initial Parameters
Nactacross = 10;    %Number actuators across any given length

%SET Initial DM SIN PROPERTIES HERE: 
h0  = 3;           %Vary amplitude: 10 to 100, step of 10
q   = 1.6;         %Fix Angle: Set to meet given speckle
x   = 5;         %Fix Distance: Set to meet given speckle
alp = 1.5;         %Vary Phase: 0 to 2pi, step of pi/10

%SET NUMBER OF ITERATIONS TO PERFORM HERE:
itr = 200;

%Save initial condition as vector for optimization
par0 = [h0, q, x, alp];

%SET PROPERTIES FOR SEARCH
LB = [1.0, 1.2, 3.0, 0.0];              %Lower bounds; same elements as par0
UB = [20.0, 2.0, 10.0, 2*pi];       %Upper bounds; same elements as par0

%% Flatten DM

%Create StringBuilder to retrieve and hold data; Set cap. to expected size
sb = System.Text.StringBuilder();
sb.Capacity = 4100;

%Create Matrix to hold images: Width(200)xHeight(200)x#Images            
    %Images have been cropped to reduce cube size
Idat = zeros(200,200,itr);    %Raw cropped images
%Create column vector to hold Power data
Pdat = zeros(itr, 1);
%Create Matrix to hold trial parameters
ParT = zeros(itr, 4);

%Subtract 20 from Pdat for easy identification of empty elements later
    %In case fminsearch converges <itr iterations. Do not  want to process
    %or save excess images later
Pdat        = Pdat - 20;

% %Reflatten cube for Flat image comparison
DM_Command = zeros(12);
FlatCommand = NK_MultiDM_Command(DM_Command, flat);
JR_UPDATE_MultiDM(drv_inf, FlatCommand);

%% Take Pre-run Zero Data (MAKE SURE LIGHT IS OFF)
pause(.05)

user_cont = 0;
while(~user_cont)

    reading = s.inputSingleScan
    user_cont = input('Light is OFF? Continue? ');
end

pm_scale = 10^input('Exponent?');

prerun_zero_array = zeros(100,1);
for i=1:1:100
    prerun_zero_array(i) = s.inputSingleScan/pm_scale;
end
prerun_zero = mean(prerun_zero_array);

%% Pause to TURN ON LIGHT and ensure DM is flat; Note: DM mech response time <100us
pause(.05);

user_cont = 0;
while(~user_cont)

    reading = s.inputSingleScan
    user_cont = input('Light is ON? Continue? ');
end

prerun_offaxis = mean(prerun_offaxis_array);

%% Take Flat Data (PUT IN CORONAGRAPH AND ALIGN TO SPECKLE BEFORE PRESSING A KEY)
% For aligning to speckle - repeat taking camera images until user input
global KEY_IS_PRESSED

Pdat(1) = reading/pm_scale;
Pflat = Pdat(1);
%% Iterate, Optimize, and Process Data

REtxt = fopen(Rawnm, 'wt');

fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
% fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);
fprintf(REtxt, 'Prerun Zero: %05.2e \n', prerun_zero);
fprintf(REtxt, 'Prerun Off Axis Throughput: %05.2e \n', prerun_offaxis);

%Define a counter to allow matching variables to each other
icount = 1;     %Start counter at 1 to compensate for flat image
% % %Use fminsearch to iterate and optimize
currentTime = clock;
currentSecond = currentTime(6) + 60*currentTime(5) + 3600*currentTime(4);
startingSecond = currentSecond;

par0(1) = 0;
for i=1:5
    power = DE_supMinimizer_3(par0);
end

%% Find Location

% Set lower and upper bounds of range to look
LB = [1.0, -0.3, 3.0, 0.0];              %Lower bounds; same elements as par0
UB = [5.0, 0.3, 10.0, 2*pi];       %Upper bounds; same elements as par0

% Initialize outputs
fiber_qx = [0 0];
best_amplitude = 0;
best_phase = 0;

Search over the lower and upper bounds
for i=LB(2):0.1:UB(2)
    for j=LB(3):1:UB(3)
        par0 = [50 i j 0];
        amp_max = DE_supMinimizer_3(par0);
        
        % Update fiber location
        if (amp_max > best_amplitude)
            fiber_qx = [i j];
            best_amplitude = amp_max;
        end
    end
end

fiber_qx = [-0.2 4];
disp('Fiber Location Search Complete')
disp(fiber_qx(1))
disp(fiber_qx(2))
      
%% Exhaustive search around fiber location

% Put 8 layers of antispeckles
for i=1:8
    % Measure power with flat DM for comparison
    par0 = [0 fiber_qx(1) fiber_qx(2) 0];
    starting_amplitude = DE_supMinimizer_3(par0);
    
    % Search
    
    % Initialize outputs
    fiber_qx_fine = [fiber_qx(1) fiber_qx(2)];
    best_amplitude_fine = starting_amplitude;
    best_phase_fine = 0;
    amp = 10/i;
    if (amp < 1)
        amp = 1;
    end
    
    % Fine search
    for q_index = (fiber_qx(1) - 0.15):0.075:(fiber_qx(1) + 0.15)
        for x_index = (fiber_qx(2) - 1.5):0.5:(fiber_qx(2) + 1.5)
            % Measure at 5 phases
            phase_measurements = [0 0 0 0 0 0 0];
            phase_index = 1;
            for k=0:2*pi/7:12*pi/7
                par0 = [amp q_index x_index k];
                phase_measurements(phase_index) = DE_supMinimizer_3(par0);
                phase_index = phase_index + 1;
            end

            % Fit to sinusoid to find minimum
            x = [0 2*pi/7 4*pi/7 6*pi/7 8*pi/7 10*pi/7 12*pi/7];
            y_max = max(phase_measurements);
            y_min = min(phase_measurements);
            y_range = y_max - y_min;
            y_z = phase_measurements - y_max + (y_range/2);
            per = 2*pi;
            y_mean = mean(phase_measurements);

            fit = @(b,x) b(1).*(sin(x + b(2))) + b(3);
            fcn = @(b) sum((fit(b,x) - phase_measurements).^2);
            sinusoid_fit = fminsearch(fcn, [y_range; -1; y_mean]);

            % Calculate Residuals
            fitted_values = sinusoid_fit(1).*sin(x + sinusoid_fit(2)) + sinusoid_fit(3);
            ss_tot = sum((phase_measurements - mean(phase_measurements)).^2);
            ss_res = sum((fitted_values - phase_measurements).^2);
            coefficient = 1 - ss_res/ss_tot
            
            xp = linspace(min(x), max(x));
            figure(1)
            plot(x,phase_measurements,xp,fit(sinusoid_fit,xp))
            xlim([0 2*pi])
            grid
            
            % Find minimum amplitude and phase
            amp_min = sinusoid_fit(3) - abs(sinusoid_fit(1));

            if(sinusoid_fit(1) < 0)
                sinusoid_fit(2) = sinusoid_fit(2) + pi;
            end
            
            phase_min = 3*pi/2 - sinusoid_fit(2);
            if (phase_min < 0)
                phase_min = phase_min + 2*pi;
            end        

            % Update fiber location
            if (amp_min < best_amplitude_fine)
                fiber_qx_fine = [q_index x_index];
                best_amplitude_fine = amp_min;
                best_phase_fine = phase_min;
            end
        end            
    end
    
    % Raise the amplitude to maximum
    par0 = [0 fiber_qx_fine(1) fiber_qx_fine(2) best_phase_fine];
    
    best_amp = 0;
    best_power = starting_amplitude;
    for amp=0:0.25:(4*amp)
        par0 = [amp fiber_qx_fine(1) fiber_qx_fine(2) best_phase_fine];
        current_power = DE_supMinimizer_3(par0);
        if (current_power < best_power)
            best_power = current_power;
            best_amp = amp;
        end
    end    
    
    % Create new flat
    DM_Map = DE_DMMapSin((best_amp),par0(2),par0(3),par0(4));
    flat = flat + DM_Map;     
   
    % Sit on new flat for a bit
    for j=1:10
        par0(1) = 0;
        disp('New Flat')
        Power = DE_supMinimizer_3(par0)
    end
    
end

%%


% for i=1:8
%     minpower = 1000;
%     minpar = [0,5,0,0];
%     
%     par_min = [0, 1.4, 3, 0];
%     par_step = [4, 0.1, 1, 0.4];
%     par_max = [16, 1.6, 8, 6.0];
% 
%     for h0=par_min(1):par_step(1):par_max(1)
%         for q=par_min(2):par_step(2):par_max(2)
%             for x=par_min(3):par_step(3):par_max(3)
%                 for a=par_min(4):par_step(4):par_max(4)
%                     par = [h0,q,x,a];
%                     power = DE_supMinimizer_3(par);   
%                     if(power < minpower)
%                         minpower = power;
%                         minpar = [h0,q,x,a];
%                     end
%                 end
%             end
%         end
%     end
% 
%     for j=1:4
%         if(minpar(j) - par_step(j) > par_min(j))
%             par_min(j) = minpar(j) - par_step(j);
%         end
%         if(minpar(j) + par_step(j) < par_max(j))
%             par_max(j) = minpar(j) + par_step(j);
%         end
%     end
%     
%     par_step = [2, 0.05, 0.5, 0.1];
%     
%     for h0=par_min(1):par_step(1):par_max(1)
%         for x=par_min(3):par_step(3):par_max(3)
%             q = minpar(2);
%             for a=par_min(4):par_step(4):par_max(4)
%                 par = [h0,q,x,a];
%                 power = DE_supMinimizer_3(par);   
%                 if(power < minpower)
%                     minpower = power;
%                     minpar = [h0,q,x,a];
%                 end
%             end
%         end
%     end
% 
%     DM_Map = DE_DMMapSin(minpar(1),minpar(2),minpar(3),minpar(4));
%     flat = flat + DM_Map;    
% 
% end

%minimize power
% amp_init = 10;
% ang_init = 1.5;
% act_init = 5;
% pha_init = 1.2;
% 
% global fmin_global
% fmin_global = Inf;
% global xmin_global
% 
% 
% for n=1:1:8
%     xmin_global = [0,0,5,0];
%     opts_coarse = optimoptions(@fmincon,'Algorithm','interior-point','StepTolerance',.05,...
%         'DiffMinChange', 0.2);
%     problem_coarse = createOptimProblem('fmincon','x0',[20, ang_init,act_init,pha_init],...
%         'objective',@DE_supMinimizer_3, 'lb',[10, 1.4, 3.0, 0.0], 'ub',[10, 1.7, 10.0, 6.2], 'options',opts_coarse);
% 
%     ms = MultiStart('FunctionTolerance',0,'XTolerance',.05, 'Display','iter',...
%         'StartPointsToRun','bounds','PlotFcn',@gsplotbestf);
% 
% 
% 
%     [xmin_coarse, fmin_coarse] = run(ms,problem_coarse,18);
%     
%     opts_fine = optimoptions(@fmincon,'StepTolerance',.05,'DiffMinChange',.05);
%     problem_fine = createOptimProblem('fmincon','x0', xmin_global,'objective',@DE_supMinimizer_3,...
%         'lb', [10, xmin_global(2)-.05, xmin_global(3)-.5, xmin_global(4)-.1], ...
%         'ub', [10, xmin_global(2)+.05, xmin_global(3)+.5, xmin_global(4)+.1], ...
%         'options', opts_fine);
%     
%     [xmin_fine, fmin_fine] = run(ms,problem_fine,18);
%     
%     %sweep amplitudes
%     for amp_sweep=0:5:90
%        par_in = [amp_sweep, xmin_global(2), xmin_global(3), xmin_global(4)];
%        pow = DE_supMinimizer_3(par_in); 
%     end
%     
%     for amp_step=-2:1:2
%        par_in =  [xmin_global(1)+amp_step, xmin_global(2), xmin_global(3), xmin_global(4)];
%        pow = DE_supMinimizer_3(par_in);
%     end
%     
%     xmin = xmin_global;
%     fmin = fmin_global;
%     
%     DM_Map = DE_DMMapSin(xmin_global(1),xmin_global(2),xmin_global(3),xmin_global(4));
%     flat = flat + DM_Map;    
% end
% 
% %set to found min
% pow = DE_supMinimizer_3(xmin_global);

% [Par, Pmin] = patternsearch(@DE_supMinimizer_3, par0, [], [], [], [], ... 
%      LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));
% 
%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
%Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

% %Save all images as single cube
%fitswrite(Idat(:,:,1:icount), Idatnm);


%% Take post run off axis data (take out coronagraph and realign to max)
% pause(.05);
% 
% user_cont = 0;
% while(~user_cont)
%     for i=1:1:50
%         pause(.5)
%         reading = s.inputSingleScan
%     end
%     user_cont = input('Continue? 1 or 0');
% end
% 
% pm_scale = 10^input('Exponent?');
% 
% postrun_offaxis_array = zeros(100,1);
% for i=1:1:100
%    postrun_offaxis_array(i) = s.inputSingleScan/pm_scale;
% end
% postrun_offaxis = mean(postrun_offaxis_array);
% 
% fprintf(REtxt, 'Postrun Off Axis Throughput: %05.2e \n', postrun_offaxis);

%% Take post run zero (turn off light)
pause(.05);

postrun_zero = mean(postrun_zero_array);
fprintf(REtxt, 'Postrun Zero: %05.2e \n', postrun_zero);

%% Disconnect all Devices

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end


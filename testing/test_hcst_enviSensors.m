% test_hcst_enviSensors.m
%
%
%
% Grady Morrissey - July 10, 2019


numtry = 100;
for II=1:numtry
    %get sensor data
    [humi, temp, etcetera] = hcst_getEnviData();
    disp(['Humidity: ',num2str(humi)])
    disp(['Temp: ',num2str(temp)])
    disp(['Etc...: ',num2str(etce)])
    
%     % code a simple graphical display
%     figure(201); 
    
    %pause
    pause(10)
end
    
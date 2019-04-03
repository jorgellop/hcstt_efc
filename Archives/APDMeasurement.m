function frequency = APDMeasurement(apd)

%% Instrument Configuration and Control

% Communicating with instrument object, obj1.
output = query(apd, 'COUNTERFreq:CH1Value?');
frequency = str2double(output);
pause(1.5);
% 
% counts = zeros(100,1);
% for i=1:30
%     output = query(apd, 'MEASUrement:IMMed:VALue?');
%     frequency = str2double(output);
%     if (frequency < 100)
%         counts(i) = frequency;
%     end
% end
% frequency = sum(counts);

end     


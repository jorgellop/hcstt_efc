function power = TakePMData(Samp)

global newp sb

% Clear PM data store of any old data
newp.Write(4, 'PM:DS:CLear');

% Take Data
newp.Write(4, 'PM:DS:ENable 1');
pause(.1+Samp/1000);        %Collection takes Samp*.0001 seconds

% Process and Save PM Data
sb.Clear();         %Ensure buffer is empty
% Place data in PM internal readable storage
newp.Write(4, ['PM:DS:Get? +' num2str(Samp)]);
pause(.05);         %Wait for data to finish storing in PM

% Loop to read arbitrary length data; Reads until no more data in PM
fl = 1;             %Flag to end reading           
str = '';           %variable to store data
while fl            %Loops until flag is tripped
    sb.Clear();
    newp.ReadBinary(4, sb);             %Read 4096 characters from PM
    pause(.05);
    str = [str char(sb.ToString())];   %Append new read values
    % Change flag when the 'End of Data' marker is found
    fl = ~logical(strcmpi(str(end-12:end), ['End of Data' 13 10]));
end

% Remove Header and footer from Data
h = strfind(str,'End of Header');
str = str((h+13):end-13); 
% disp(str);

% Save first PM value as flat; Take mean of Samp points since signal
        %is noisy without the Analog and Digital Filter
power = mean(str2num(str));

end

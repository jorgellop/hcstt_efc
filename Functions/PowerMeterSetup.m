function PowerMeterSetup(Samp)

global newp 

% %Check if Log File exists; rename. Seems to prevent errors
% daN = clock();
% lgnm = sprintf('C:\\Users\\Gary\\Documents\\Log\\Log%04.0f-%02.0f-%02.0f.txt',daN(1), daN(2), daN(3));
% if exist(lgnm, 'file')
%     fprintf('Log File already existed. Renamed to append time\n');
%     da = num2str(daN, '%0.0f');
%     dt = strcat(lgnm(1:end-4),'-',da(15:16),da(19:20),da(23:24));
%     movefile(lgnm, [dt '.txt']);
% end

% Check if PM assembly is present. Import Otherwise
% May need to change search location if library is moved
asm = System.AppDomain.CurrentDomain.GetAssemblies;
if ~any(arrayfun(@(n) strncmpi(char(asm.Get(n-1).FullName), ...
        'UsbDllWrap.dll', length('UsbDllWrap.dll')), 1:asm.Length))
    anet = NET.addAssembly( ...
        'C:\Program Files\Newport\Newport USB Driver\Bin\UsbDllWrap.dll');
end

% Create instance of PM; TRUE = start logging
newp = Newport.USBComm.USB(true);

% Begin Connection to PM; Make sure it is on and connected
if ~newp.OpenDevices()
    error('Failed to connect to Power Meter')
end

% Set Buffer Behavior to 0 (Fixed Size)
newp.Write(4, 'PM:DS:BUFfer 0');

% Clear data store of any old data
newp.Write(4, 'PM:DS:CLear');

% Set Sample rate to 10kHz (.1 ms Interval)
newp.Write(4, 'PM:DS:INTerval 1');

% Set Size to 50; will take 50 samples
newp.Write(4, ['PM:DS:SIZE ' num2str(Samp)]);

end

%% To Close Power Meter
% newp.CloseDevices();
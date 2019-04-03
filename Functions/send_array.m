function send_array(driver_info,input_array)
%input a 144 command array for dm, in voltages (0-300)
max_V=250;


 amplitudes = reshape(input_array,144,1);
% amplitudes = reshape(input_array,12,12);

amplitudes = max(amplitudes,0);
amplitudes = min(amplitudes,max_V);
UPDATE_multiDM(driver_info, amplitudes/3);
% error_code = send_voltages_to_channels(channels, amplitudes/3,length(channels));
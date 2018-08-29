% hcstt_calibrateNIPhotodiodeGainChange.m
%
% Calibrate the readings from the NI photodiode wrt to a specific gain
% setting
%
% Jorge Llop - Aug 23, 2018

for II=1:6
    load(['/Users/jllopsay/Documents/GitHub/hcstt_efc/output/NIPhotoDiodeCalibration/levels',num2str(II),'-',num2str(II+1)]);
    gain_arr(II) = mean(int2_arr)/mean(int1_arr);
    level_arr(II) = II;
end
save('output/NIPhotodiodeGainCalibration.mat','level_arr','gain_arr')

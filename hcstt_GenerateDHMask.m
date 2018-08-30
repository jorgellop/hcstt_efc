% hcstt_GenerateDHMask
% 
% Generate DH mask to match the throughput of a SMF fiber
%
% Jorge Llop - Aug 29, 2018

function [Q,totalPowerCam] = hcstt_GenerateDHMask(im_mod,X,Y,info)
lambdaOverD = info.lambdaOverD;
N = info.N;
RHO = info.RHO;

sidepix = 30;
im_mod_crop = im_mod(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
totalPower = sum(abs(im_mod_crop(:)).^2);
throughput = 0.8119;
throughput_arr = [];
numtry  = 35;
diam_arr = linspace(0.5,3,numtry);
for JJ = 1:numtry
     Q_pixJJ = exp( -(RHO/(diam_arr(JJ)/2*lambdaOverD)).^100 );
     throughput_arr = [throughput_arr; sum(sum(abs(im_mod.*Q_pixJJ).^2))/totalPower];
end
[THETA_fib,RHO_fib] = cart2pol(X - info.x_fib_pix ,Y - info.y_fib_pix);
diamII = interp1(throughput_arr,diam_arr,throughput);
Q = exp(-(RHO_fib/(diamII/2*lambdaOverD)).^100);
Qcent = exp(-(RHO/(diamII/2*lambdaOverD)).^100);
totalPowerCam = sum(sum(abs(im_mod.*Qcent).^2))*info.normPowerRegEFC;
end
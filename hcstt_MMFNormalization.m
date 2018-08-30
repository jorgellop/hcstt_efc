% hcstt_MMFNormalization
% 
% 
%
% Jorge Llop - Aug 30, 2018
N = 1024;
lambdaOverD = N/apRad/2;
[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X,Y);

info.useGPU = false;
info.RHO = RHO;
info.THETA = THETA;
info.N = N;
info.lam_arr =  635e-9;
info.useApodizer = false;
info.useGPU = false; 
info.FPM = exp(1i*4*THETA);
info.apRad = 142;
info.lambdaOverD = lambdaOverD;
info.LPM = exp(-(RHO/(0.85*apRad)).^1000);

wf1 = exp( -(RHO/(apRad)).^100 );
wf2 = prescription_DM1toImage_compact_vFiberCoupling_broadband( wf1, zeros(N,N), false, info);

sidepix = 30;
wf2_crop = wf2(N/2-sidepix+1:N/2+sidepix+1,N/2-sidepix+1:N/2+sidepix+1);
totalPower = sum(abs(wf2_crop(:)).^2);
throughput = 0.8119;
throughput_arr = [];
numtry  = 35;
diam_arr = linspace(0.5,3,numtry);
for JJ = 1:numtry
     Q_pixJJ = exp( -(RHO/(diam_arr(JJ)/2*lambdaOverD)).^100 );
     throughput_arr = [throughput_arr; sum(sum(abs(wf2.*Q_pixJJ).^2))/totalPower];
end
% [THETA_fib,RHO_fib] = cart2pol(X - info.x_fib_pix ,Y - info.y_fib_pix);
diamII = interp1(throughput_arr,diam_arr,throughput);
% Q = exp(-(RHO_fib/(diamII/2*lambdaOverD)).^100);
Qcent = exp(-(RHO/(diamII/2*lambdaOverD)).^100);
% totalPowerCam = sum(sum(abs(wf2.*Qcent).^2))*info.normPowerRegEFC;
%%
hcstt_Initialize(false)
take_background = true;
Ncam = 400;

tint = 0.05;
keeptrying = true;

while keeptrying
    if(take_background)
        prompt = 'Take out light. Continue? ';
        x = input( prompt );
        im_cam = zeros(400,400);
        for II=1:15
            im_camII = hcstt_TakeCamImage(true,false,tint);
            im_cam = im_cam + im_camII/15;
            pause(0.1)
        end
        backgroundCam = im_cam;
        prompt = 'Put back light on. Continue? ';
        x = input( prompt );
    else
        backgroundCam = zeros(Ncam,Ncam);
    end

    for II=1:15
        im_camII = hcstt_TakeCamImage(true,false,tint);
        im_cam = im_cam + im_camII/15;
        pause(0.1)
    end
    DH = im_cam(find(Qcent));
    if min(DH)<10
        tint = tint*2;
    elseif max(DH)>200
        tint = tint/2;
    else
        kepptrying = false;
    end
end

totalPowerMMF_calibration = sum(DH);
tint_MMF_calibration = tint;
diamMMF = diamII;
save('calibrationMMF_Aug30','totalPowerMMF_calibration','tint_MMF_calibration','diamMMF')
hcstt_DisconnectDevices()

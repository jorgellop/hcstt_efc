% test_LOWFEvsSMFCoupling
%
% Code to plot the couling into a fiber for different LOWFEs. Done for
% Jacques Delorme work on KPIC
%
% Jorge Llop - Feb 19, 2019
close all;clear all;
addpath(genpath('utils'),genpath('export_scripts'));

N =1024;
fiberDiam = 2*0.71; % Fiber diam. (lambda_0/D)

apRad = N/8; % Aperture radius (samples)
lambdaOverD = N/apRad/2; % lambda/D (samples) 

fiberDiam_pix = (fiberDiam*lambdaOverD);
[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X, Y);
[THETA45,~] = cart2pol(X+Y, -X+Y);
fibermode0 = sqrt(2/(pi*(fiberDiam_pix/2)^2))* ...
        exp(-(RHO/(fiberDiam_pix/2)).^2);

figure(1)
imagesc(fibermode0)

wf0 = zeros(N);
wf0(RHO<apRad) = 1+0i;
power_tot = sum(wf0(:));

% apply DM1
% wfP = wfP.*fftshift(exp(1i*4*pi*Z/633e-9));
wfF = myfft2(wf0);

figure(2)
imagesc(abs(wfF).^2)
power_tot = sum(sum(abs(wfF).^2));
coupl0 = abs(sum(sum(wfF.*fibermode0))).^2/power_tot;

numZernikes=11;
numAmp = 100;
minAmp = 0;
maxAmp = 0.5; %waves
wfe_arr = linspace(minAmp,maxAmp,numAmp);
fig0 = figure('visible','off','color','none');
plot(wfe_arr,ones(1,numAmp)*coupl0,'LineStyle','--','Color','k')
hold on
% good_ind = [2,4,5,7,9,11];
good_ind = [2,3];

%Initial Defocus
Z00 = generateZernike( 6, apRad, RHO, THETA  );
Z00 = Z00*100/1600;
for II = 2:numZernikes
    if ismember(II,good_ind)
        Z0 =generateZernike( II, apRad, RHO, THETA  );
        coupl = [];
        for JJ=1:numAmp
            wfe = wfe_arr(JJ);
            Z = Z0*wfe + Z00;
            wfP = (wf0).*(exp(1i*2*pi*Z));
            wfF = myfft2(wfP);
            power_tot = sum(sum(abs(wfF).^2));
            coupl = [coupl, abs(sum(sum(wfF.*fibermode0))).^2/power_tot];
        end
        plot(wfe_arr,coupl,'LineWidth',4)
        hold on
    end
end
for II = 2:numZernikes
    if ismember(II,good_ind)
        Z0 =generateZernike( II, apRad, RHO, THETA45  );
        coupl = [];
        for JJ=1:numAmp
            wfe = wfe_arr(JJ);
            Z = Z0*wfe + Z00;
            wfP = (wf0).*(exp(1i*2*pi*Z));
            wfF = myfft2(wfP);
            power_tot = sum(sum(abs(wfF).^2));
            coupl = [coupl, abs(sum(sum(wfF.*fibermode0))).^2/power_tot];
        end
        plot(wfe_arr,coupl,'LineWidth',4)
        hold on
    end
end
title(['SMF Coupling w/ LOWFE - Initial AstigY 100 nmRMS'],'Interpreter','latex')
xlabel('WFE [wavesRMS]','Interpreter','latex');
ylabel('SMF coupling','Interpreter','latex');
%         xlim([0 2])
ylim([0 1])
set(gca,'FontSize',20)
% lgd = legend('No WFE','Tip/Tilt','Defocus','Astigmatism','Coma','Trefoil','Spherical');
lgd = legend('No WFE','Tip/Tilt X','Tip/Tilt Y','Tip/Tilt X+45deg','Tip/Tilt Y+45deg');
title(lgd,'Zernike Modes')
hold off
figure(fig0)
export_fig(['output/SMFCoupling_LOWFE_AstigY100nmrms.png']);

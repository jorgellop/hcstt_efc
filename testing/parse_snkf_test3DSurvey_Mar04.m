% parse_snkf_test3DSurvey_Mar04.m
%
% Measure which filter is best from 3 D survey on weeekend Mar 1, 2019
%
% Rebecca Zhang - Mar 04, 2019


clear all; close all;

path2files = '';
numfiles;

numw = 5;
numphn=numw;
numamn=numw;
w_arr = linspace(-5,-3,numw);
phn_arr = linspace(-5,-3,numphn);
amn_arr = linspace(-5,-3,numamn);

foldernm = '';

w = w_arr(1);
for IIphn=1:numphn
    for IIamn=1:numamn
        phn = phn_arr(IIphn);
        amn = phn_arr(IIphn);
        nmfl = [foldernm,'snkf_w',num2str(w),'_phn',num2str(phn),'_amn',num2str(amn)];
        load(nmfl)
        me_NI_unfilt = median(intensity_arr/normI);
        me_NI_kf = median(intensityKF_arr/normI);
        improv_mat(ind_j, ind_i) = me_NI_unfilt/me_NI_kf;
        stddev_NI_unfilt = std(intensity_arr/normI);
    end
end
fig0 = figure('visible','off','color','w','pos',[10 10 500 500]);
imagesc(improv_mat)
title('')
xlabel()
ylabel()
set(gca,'FontSize',15)
figure(fig0)
export_fig(['.png');

export_fig
fileSMF = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/EFC_wFiber_LabDemonstration_EFCSMF_laserSource_1config_Aug31v4/data_intvsit_dmshapes__EFCSMF_laserSource_1config_Aug31v4_DMconfig1_apRad142';
fileReg = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/EFC_wFiber_LabDemonstration_RegEFC_laserSource_1config_Aug31/data_intvsit_dmshapes__RegEFC_laserSource_1config_Aug31_DMconfig1_apRad142';
fileBroad = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/EFC_wFiber_LabDemonstration_EFCSMF_broadband_1config_Aug31/data_intvsit_dmshapes__EFCSMF_broadband_1config_Aug31_DMconfig1_apRad142';

load(fileSMF)
numit = numel(coupl_SMF_in_DH);
figure(301)
plot(1:numit,log10(coupl_SMF_in_DH/peakIntEFCSMF),'linewidth',3,'color','r');%,'color',[1 0.1 0.4])
hold on
plot(1:numit,log10(coupl_MMF_in_DH./totalPowerMMF_arr),'linewidth',3,'color','b');%,'color',[0.1 0.4 0.4])
plot(1:numit,log10(coupl_pixGauss_in_DH./totalPowerPixGauss_arr),'linewidth',3,'color','k');%,'color',[0.4 0.1 0.4])
hold off
ylabel('Raw Contrast (Log Scale)')
xlabel('Iteration')
% title(['EFCSMF Log Scale - 2.6 actxcycl - Suppression: ',num2str(int_in_DH(1)/int_in_DH(numel(int_in_DH)))])
title(['HCST-T EFCSMF - 635nm Monochromatic Light'])
legend('SMF Raw Contrast','Pixel Aperture Raw Contrast','Mode Filtered Intensity Raw Contrast');
set(gca,'FontSize',14);
ax = gca;
ax.YLim = [-5.5 -2.7];
export_fig('RC_it_EFCSMF_Aug31.png','-r300');

load(fileReg)
figure(302)
numit = numel(coupl_SMF_in_DH);
plot(1:numit,log10(coupl_SMF_in_DH/peakIntEFCSMF),'linewidth',3,'color','r');%,'color',[1 0.1 0.4])
hold on
plot(1:numit,log10(int_in_DH/peakIntRegEFC),'linewidth',3,'color',[255/255, 165/255, 0])
plot(1:numit,log10(coupl_MMF_in_DH./totalPowerMMF_arr),'linewidth',3,'color','b');%,'color',[0.1 0.4 0.4])
plot(1:numit,log10(coupl_pixGauss_in_DH./totalPowerPixGauss_arr),'linewidth',3,'color','k');%,'color',[0.4 0.1 0.4])
hold off
ylabel('Raw Contrast (Log Scale)')
xlabel('Iteration')
% title(['EFCSMF Log Scale - 2.6 actxcycl - Suppression: ',num2str(int_in_DH(1)/int_in_DH(numel(int_in_DH)))])
title(['HCST-T Regular EFC - 635nm Monochromatic Light'])
legend('SMF Raw Contrast','Mean DH Raw Contrast','Pixel Aperture Raw Contrast','Mode Filtered Intensity Raw Contrast');
set(gca,'FontSize',14);
ax = gca;
ax.YLim = [-5.5 -2.7];
export_fig('RC_it_regEFC_Aug31.png','-r300');

load(fileBroad)
numit = numel(coupl_SMF_in_DH);
figure(303)
plot(1:numit,log10(coupl_SMF_in_DH/peakIntEFCSMF),'linewidth',3,'color','r');%,'color',[1 0.1 0.4])
hold on
plot(1:numit,log10(coupl_MMF_in_DH./totalPowerMMF_arr),'linewidth',3,'color','b');%,'color',[0.1 0.4 0.4])
plot(1:numit,log10(coupl_pixGauss_in_DH./totalPowerPixGauss_arr),'linewidth',3,'color','k');%,'color',[0.4 0.1 0.4])
hold off
ylabel('Raw Contrast (Log Scale)')
xlabel('Iteration')
% title(['EFCSMF Log Scale - 2.6 actxcycl - Suppression: ',num2str(int_in_DH(1)/int_in_DH(numel(int_in_DH)))])
title(['HCST-T EFCSMF - 625nm 8% Broadband Light'])
legend('SMF Raw Contrast','Pixel Aperture Raw Contrast','Mode Filtered Intensity Raw Contrast');
set(gca,'FontSize',14);
ax = gca;
ax.YLim = [-5.5 -2.7];
export_fig('RC_it_broadband_Aug31.png','-r300');

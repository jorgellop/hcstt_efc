pathSimuSMF = '/Users/jllopsay/Documents/MATLAB/HCSTR_fiber_proper_model/output/HCST_FIU_broadband_12actrng16_EFCSMF_Nov13_EFCwDM1_VC6_N1024_circ_LSout0.9_lam6.5e-07_BW10_surf_error_PSDfit/EFC_fiberConsistentWLab_IWA3.5_OWA4.5/';
pathSimuRegEFC = '/Users/jllopsay/Documents/MATLAB/HCSTR_fiber_proper_model/output/HCST_FIU_broadband_12actrng16_RegEFC_Nov13_EFCwDM1_VC6_N1024_circ_LSout0.9_lam6.5e-07_BW10_surf_error_PSDfit/EFC_fiberConsistentWLab_IWA3.5_OWA4.5/';

N = 1024;
apRad = N/8; % Aperture radius (samples)
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);

load([pathSimuSMF,'DMshapes'])
    fig0 = figure('visible','off','color','w');
imagesc(xvals/apRad,yvals/apRad,surf_DM1*1e9);
colormap(gray(256));hcb=colorbar;
ylabel(hcb, 'Surface height (nm)')
%caxis([0 10e-9]);
axis image;%axis off;% 
axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
%text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
% title('DM1 surface height (nm)');
set(gca,'FontSize', 35,...
                'TickDir','out',...
                'TickLength',[.02 .02]);
set(gca,'YDir','normal');
set(fig0,'units', 'inches', 'Position', [0 0 10 10])
figure(fig0)
export_fig(['DM1surf_SimuSMF.png']);
close(fig0)

load([pathSimuRegEFC,'DMshapes'])
fig0 = figure('visible','off','color','none');
imagesc(xvals/apRad,yvals/apRad,surf_DM1*1e9);
colormap(gray(256));hcb=colorbar;
ylabel(hcb, 'Surface height (nm)')
%caxis([0 10e-9]);
axis image;%axis off;% 
axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
%text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
% title('DM1 surface height (nm)');
set(gca,'FontSize', 35,...
                'TickDir','out',...
                'TickLength',[.02 .02]);
set(gca,'YDir','normal');
set(fig0,'units', 'inches', 'Position', [0 0 10 10])
export_fig(['DM1surf_SimuRegEFC.png']);
figure(fig0)
close(fig0)
%%
load('normI4Paper')
% load([pathSimuSMF,'bm'])
load([pathSimuRegEFC,'bm'])
% load(['/Users/jllopsay/Documents/MATLAB/HCSTR_fiber_proper_model/output/HCST_FIU_broadband_12actrng6_RegularEFC_Oct01_EFCwDM1_VC4_N1024_circ_LSout0.9_lam6.35e-07_BW10_surf_error_PSDfit/EFC_fiberConsistentWLab_IWA3.5_OWA4.5/','bm'])
info.xvals = xvals;
info.yvals = yvals; 
info.lambda0 = 650e-9; 
info.numOfWavelengths = 3;
info.normI = normI;
info.N = N;

% plotBBpsf4Paper(bm,'BBPSF_SimuEFCSMF_Jan29',info)
plotBBpsf4Paper(bm,'BBPSF_SimuRegEFC_Jan29',info)

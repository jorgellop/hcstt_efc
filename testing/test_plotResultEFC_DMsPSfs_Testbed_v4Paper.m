addpath(genpath('utils'),genpath('export_scripts'));
pathOutput = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/';
fileSMF = 'EFC_wFiber_LabDemonstration_EFCSMF_laserSource_1config_Aug31v4/data_intvsit_dmshapes__EFCSMF_laserSource_1config_Aug31v4_DMconfig1_apRad142';
fileReg = 'EFC_wFiber_LabDemonstration_RegEFC_laserSource_1config_Aug31/data_intvsit_dmshapes__RegEFC_laserSource_1config_Aug31_DMconfig1_apRad142';
fileBroad = 'EFC_wFiber_LabDemonstration_EFCSMF_broadband_1config_Aug31/data_intvsit_dmshapes__EFCSMF_broadband_1config_Aug31_DMconfig1_apRad142';

disp('SMF troughput')
load([pathOutput,fileSMF])
N = 1024;
Nact = 12;
apRad = 142;
lambda0 = 635e-9;
[posDM_x,posDM_y,ac_spac] = hcstt_PositionDMActuatorsvFindBestDMOrientation(N,apRad,1);
infl = loadInfluenceFunction( 'influence_dm5v2.fits', ac_spac ); % Influence function. Need of a model for actual DM
DM1_strokes = zeros(N,N);
count = 0; 
for ix = 1:Nact
    xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
    for iy = 1:Nact
        count = count + 1;
        ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        DM1_strokes(xpos,ypos) = us_total(count)*1e-9 + DM1_strokes(xpos,ypos);
    end
end
surf_DM10 = conv2(DM1_strokes,infl,'same');
[X,Y] = meshgrid(-N/2:N/2-1); 
xvals = X(1,:);yvals = Y(:,1);

%%
fig0 = figure('visible','off','color','none');
imagesc(xvals/apRad,yvals/apRad,surf_DM10*1e9);
colormap(gray(256));hcb=colorbar;
ylabel(hcb, 'Surface height (nm)')
%caxis([0 10e-9]);
axis image;%axis off;% 
axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
%text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
% title('DM1 surface height (nm)');
set(gca,'FontSize', 23,...
                'TickDir','out',...
                'TickLength',[.02 .02]);
set(gca,'YDir','normal');
set(fig0,'units', 'inches', 'Position', [0 0 10 10])
figure(fig0)
export_fig(['DM1surf_LabSMFv2.png']);
close(fig0)

load([pathOutput,fileBroad])
DM1_strokes = zeros(N,N);
count = 0; 
for ix = 1:Nact
    xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
    for iy = 1:Nact
        count = count + 1;
        ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        DM1_strokes(xpos,ypos) = us_total(count)*1e-9 + DM1_strokes(xpos,ypos);
    end
end
surf_DM10 = conv2(DM1_strokes,infl,'same');
fig0 = figure('visible','off','color','none');
imagesc(xvals/apRad,yvals/apRad,surf_DM10*1e9);
colormap(gray(256));hcb=colorbar;
ylabel(hcb, 'Surface height (nm)')
%caxis([0 10e-9]);
axis image;%axis off;% 
axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
%text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
% title('DM1 surface height (nm)');
set(gca,'FontSize', 23,...
                'TickDir','out',...
                'TickLength',[.02 .02]);
set(gca,'YDir','normal');
set(fig0,'units', 'inches', 'Position', [0 0 10 10])
export_fig(['DM1surf_LabBroadbandv2.png']);
figure(fig0)
close(fig0)

%%
load([pathOutput,fileReg])
DM1_strokes = zeros(N,N);
count = 0; 
for ix = 1:Nact
    xpos = round(posDM_x(ix));%round(N/2+1+(ix-Nact/2-0.5)*ac_spac);
    for iy = 1:Nact
        count = count + 1;
        ypos = round(posDM_y(iy));%round(N/2+1+(iy-Nact/2-0.5)*ac_spac);
        DM1_strokes(xpos,ypos) = us_total(count)*1e-9 + DM1_strokes(xpos,ypos);
    end
end
surf_DM10 = conv2(DM1_strokes,infl,'same');
fig0 = figure('visible','off','color','none');
imagesc(xvals/apRad,yvals/apRad,surf_DM10*1e9);
colormap(gray(256));hcb=colorbar;
ylabel(hcb, 'Surface height (nm)')
%caxis([0 10e-9]);
axis image;%axis off;% 
axis([-1.1 1.1 -1.1 1.1]);set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
%text(-1.05,0.98,'{\bf(c)}','FontSize',12,'Color','w')
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
% title('DM1 surface height (nm)');
set(gca,'FontSize', 23,...
                'TickDir','out',...
                'TickLength',[.02 .02]);
set(gca,'YDir','normal');
set(fig0,'units', 'inches', 'Position', [0 0 10 10])
figure(fig0)
export_fig(['DM1surf_LabRegEFCv2.png']);
close(fig0)
%%
% axismax = 255/peakIntRegEFC;
% axismax = 255/peakIntRegEFC;
    load([pathOutput,fileSMF])
hcstt_test_plotCamImagev4Paper(im_cam_crop, 'Camera Image - Flat DM',['testing/CamImage_finalEFCSMF'], [41,41] );
    load([pathOutput,fileReg])
hcstt_test_plotCamImagev4Paper(im_cam_crop, 'Camera Image - Flat DM',['testing/CamImage_finalRegEFC'], [41,41] );
    load([pathOutput,fileBroad])
hcstt_test_plotCamImagev4Paper(im_cam_crop, 'Camera Image - Flat DM',['testing/CamImage_finalBroadband'], [41,41] );

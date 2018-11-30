% test_MarechalApproxDMrms
%
%
%
% Jorge Llop - Sept 26, 2018
addpath(genpath('utils'),genpath('export_scripts'));
pathOutput = '/Users/jllopsay/Documents/GitHub/hcstt_efc/output/';
fileSMF = 'EFC_wFiber_LabDemonstration_EFCSMF_laserSource_1config_Aug31v4/data_intvsit_dmshapes__EFCSMF_laserSource_1config_Aug31v4_DMconfig1_apRad142';
fileReg = 'EFC_wFiber_LabDemonstration_RegEFC_laserSource_1config_Aug31/data_intvsit_dmshapes__RegEFC_laserSource_1config_Aug31_DMconfig1_apRad142';
fileBroad = 'EFC_wFiber_LabDemonstration_EFCSMF_broadband_1config_Aug31/data_intvsit_dmshapes__EFCSMF_broadband_1config_Aug31_DMconfig1_apRad142';

disp('SMF troughput')
load([pathOutput,fileSMF])
N = 1024;
[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X,Y);
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
surf_DM10(RHO>apRad) = nan;
figure(1)
imagesc(surf_DM10)
axis image
dm = surf_DM10(RHO<apRad);
1-exp(-2*pi*std(dm/lambda0))
dmrms = std(dm)

disp('RegEFC troughput')
load([pathOutput,fileReg])
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
surf_DM10(RHO>apRad) = nan;
figure(1)
imagesc(surf_DM10)
axis image
[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X,Y);
dm = surf_DM10(RHO<apRad);
1-exp(-2*pi*std(dm/lambda0))
dmrms = std(dm)

disp('BB troughput')
load([pathOutput,fileBroad])
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
surf_DM10(RHO>apRad) = nan;
figure(1)
imagesc(surf_DM10)
axis image
[X,Y] = meshgrid(-N/2:N/2-1); 
[THETA,RHO] = cart2pol(X,Y);
dm = surf_DM10(RHO<apRad);
1-exp(-2*pi*std(dm/625e-9))
dmrms = std(dm)

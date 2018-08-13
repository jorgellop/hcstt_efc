% hcstt_BenchModelPixelPitchMatching
% 
% Align model and bench to have the same pixel pitch
%
% Jorge Llop - Jan 31, 2018

clear all;
close all;

addpath(genpath('utils'),genpath('export_scripts'));

load('alignDMModel_ind_ma_arrvPlot0207');
load('alignDMModel_ind_arrvPlot0207');
load('im_matv2');

sidepix = 40;
im_mat2 = zeros(830,830);

for II=1:numel(ind_arr) 
    imaux = im_mat(:,:,ind_arr(II))/max(max(im_mat(:,:,ind_arr(II))));
    [ma,ind_ma] = max(imaux(:));
    [ind_ma_I,ind_ma_J] = ind2sub(size(imaux),ind_ma);
    imcrop = zeros(830,830);
    imcrop(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix) = imaux(ind_ma_I-sidepix:ind_ma_I+sidepix,ind_ma_J-sidepix:ind_ma_J+sidepix);
    im_mat2 = im_mat2+imcrop;
end
imagesc(im_mat2)
axis image
export_fig(['utils/AlignmentDMModel/DMAllActsPupil'],'-r300');
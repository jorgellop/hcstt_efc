clear all;
close all;
%{
hcstt_DMAlign_FindPosActuators

Selects from im_mat(830,830,144) the best images from the DM pupil image

Jorge Llop - Dec 2017
%}
load('im_matv2');

sz = size(im_mat);
im_mat2 = zeros(sz(1),sz(2));
max_arr = zeros(sz(3),1);
% for II = 1:sz(3)
%     imagesc(im_mat(:,:,II));
%     axis image
%     prompt = 'Accept image to find maximum? \n';
%     accept = input(prompt);
% end
for II = 1:sz(3)
%     im_mat2 = im_mat2 + im_mat(:,:,II);
    max_arr(II) = max(max(im_mat(:,:,II)));
end
plot(1:sz(3),max_arr)
prompt = 'What is the cutoff value? \n';
cutoff = input(prompt);
% cutoff = 50;
ind = find(max_arr>50);
ind_arr = zeros(numel(ind),1);
ind_ma_arr = zeros(numel(ind),1);
for II = 1:numel(ind)
%     im_mat2 = im_mat2 + im_mat(:,:,ind(II))/max(max(im_mat(:,:,ind(II))));
%     im_mat2 = im_mat2 + im_mat(:,:,ind(II));
    imII = im_mat(:,:,ind(II));
    imagesc(imII)
    prompt = 'Accept image to find maximum? \n';
    accept = input(prompt);
    if accept ~= 0
        [ma,ind_ma] = max(imII(:));
        ind_arr(II) = ind(II);
        ind_ma_arr(II) = ind_ma;
    end
end
ind_arr = ind_arr(find(ind_arr));
ind_ma_arr = ind_ma_arr(find(ind_ma_arr));
save('utils/AlignmentDMModel/alignDMModel_ind_arrvPlot0207','ind_arr')% imagesc(im_mat2)
save('utils/AlignmentDMModel/alignDMModel_ind_ma_arrvPlot0207','ind_ma_arr')% imagesc(im_mat2)
% axis image
% colorbar;

%% Fit circle to points from pupil image

% load('measuredPupilPoints')
% x = measuredPupilPoints(:,1);
% y = measuredPupilPoints(:,2);
% 
% [xc yx R] = circfit(x,y)
%% Find center of pupil and radius

% load('alignDMMdeol_pupil')
% im_pupil = im_cam;
% imagesc(im_pupil)
% axis image
% prompt = 'Point 1.1 \n';
% p11 = input(prompt);
% prompt = 'Point 1.2 \n';
% p12 = input(prompt);
% prompt = 'Point 2.1 \n';
% p21 = input(prompt);
% prompt = 'Point 2.2 \n';
% p22 = input(prompt);
% prompt = 'Point 3.1 \n';
% p31 = input(prompt);
% prompt = 'Point 3.2 \n';
% p32 = input(prompt);

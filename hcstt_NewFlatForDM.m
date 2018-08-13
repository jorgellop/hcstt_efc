% hcstt_NewFlatForDM
%
% Load new flat array for DM
%
% Jorge Llop - Feb 22, 2018

function [] = hcstt_NewFlatForDM(flat_file)
global cam img Xcr Ycr CExp s drv_inf flat himg pm_scale

load(flat_file)
flat = flat_SN;
end
% hcstt_AddToFlatForDM
%
% Add new flat to old flat for DM
%
% Jorge Llop - Feb 25, 2019

function [] = hcstt_AddToFlatForDM(flat_new)
global cam img Xcr Ycr CExp s drv_inf flat himg pm_scale

flat = flat + flat_new;
end
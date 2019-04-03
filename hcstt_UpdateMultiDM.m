%{
%}

function hcstt_UpdateMultiDM(h)

global cam img Xcr Ycr CExp s drv_inf flat

    FlatCommand = NK_MultiDM_Command(h, flat);
    JR_UPDATE_MultiDM(drv_inf, FlatCommand);
    pause(0.1);
end

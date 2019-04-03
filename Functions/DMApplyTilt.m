%{
DM Surface Mapping and Applying Function: Tilt
- h values are deformation in nanometers
*** Assumes N=12 for mirror grid size

******************************************************
- Arguments:
    ho          = Max poke height
    axis           = 0 for tilt around x axis, 1 for tilt around y axis
- Returns:
    h           = Surface Map as Matrix of Deformations in nm
******************************************************
%}
function h = DMApplyTilt(ho,axis)

    global flat drv_inf
    
    global fmin_global xmin_global
    
    %Map tilt
    h = zeros(12);
    if axis==0
        %xtilt
        for i=1:12
           h(:,i) = ho/12*i-ho/2;
        end

    elseif axis==1
        %ytilt
        for i=1:12
           h(i,:) = ho/12*i-ho/2;
        end

    end

    %Place tilt on mirror
    DM_Command = h;
    FlatCommand = NK_MultiDM_Command(DM_Command, flat);
    JR_UPDATE_MultiDM(drv_inf, FlatCommand);
    
    pause(1);
end
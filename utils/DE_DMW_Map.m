function [hnm, hV] = DE_DMW_Map(map, drv_info)

% Calculate Sinusoidal map in nm
hnm = map;
% Convert Map to vector for writing
hV  = DE_DMArrayToVect(hnm);

% Write new values to DM
UPDATE_multiDM(drv_info, hV);
end
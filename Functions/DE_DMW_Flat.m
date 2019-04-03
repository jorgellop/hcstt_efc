%{
DM Writing Function: Write Flat
- Places Flat Front on Mirror Surface
- Can piston whole surface by given distance
- Utilizes DE_DMArrayToVect to shape map for write
*** ASSUMES MIRROR CONNECTION ALREADY PRESENT; DOES NOT CLOSE CONNECTION
*** Assumes 144 actuators

******************************************************
- Arguments
    h0          = Piston Distance in nm
    drv_info    = DM info from OPEN_mutliDM
- Returns:
    hnm         = Surface map as matrix in nm
    hV          = Vector of voltage percentages for writing
******************************************************


Compiled By:    Daniel Echeverri
Last Modified:  08/04/2016
%}

function [hnm, hV] = DE_DMW_Flat(h0, drv_info)

hnm = zeros(12,12) + h0;        %Create map with piston distance
hV  = DE_DMArrayToVect(hnm);    %Convert map to vector for writing

UPDATE_multiDM(drv_info, hV);   %Write values to DM
end

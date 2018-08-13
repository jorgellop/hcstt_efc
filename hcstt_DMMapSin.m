%{
DM Surface Mapping Function: Sinusoidal
- Translated and adapted from Jerry Xuan Python Code
    -maps.py
- h values are deformation in nanometers
*** Assumes N=12 for mirror grid size

******************************************************
* Outputs various items depending on "out" argument
    *** All return h matrix
    out     Outputs
    0       display h, plot h, create .txt file
    1       display h
    2       plot h
    3       displayh, plot h
    4       create .txt file
    #       ONLY RETURNS h
- Arguments:
    ho          = Max poke height
    q           = angle of sinusoid
    xo          = actuators per cycle
    alp         = phase delay
    out         = output selector
    filename    = name for save file (.txt)
- Returns:
    h           = Surface Map as Matrix of Deformations in nm
******************************************************

Compiled By:    Daniel Echeverri
Last Modified:  08/04/2016
%}
function h = hcstt_DMMapSin(ho, q, xo, alp)

N   = 12;                           %DM axis base value (12x12)
Y   = repmat([-N/2:(N/2-1)],12,1);  %Matrix with ascending rows
X   = Y.';                          %Matrix with ascending columns
RHO = sqrt(X.^2 + Y.^2);            %Matrix of element dist. from cent.
TTA = atan2(Y,X);                   %Matrix of element angle from cent.

% Equation for corresponding actuator height
h   = 0.5*ho*(1+cos(2*pi.*RHO.*cos(TTA-q)/xo+alp));

end





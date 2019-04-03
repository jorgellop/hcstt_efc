%{
DM Writing Function: Convert a map in nm to a vector in % of the maxV
- Reshape the map in nm as a Vector.
- Convert nm in V.
- Convert to percentage of max V.
- Values capped at 80% of the maximum voltage.

******************************************************
- Arguments:
  Map             = Input map as 2D matrix of nm deformations

- Returns: 
  MultiDM_Command = Vector in % of the maxV directly usable by UPDATE_mutliDM

- Dependencies:
  NONE
******************************************************

Compiled By:    JR Delorme
Last Modified:  Fabruary 2, 2017
%}

function MultiDM_Command = JR_MultiDM_Command(map)

% Reshape the map in nm as a Vector
MultiDM_Command = reshape(flipud(map),[],1);  

% Convert nm in V || Equation from DM characteristics sheet. 
MultiDM_Command = - 30.16 + 5.704e-5 * sqrt(7.01264e9 * MultiDM_Command + 2.79527e11);

% Convert to percentage of max V (212V for this MultiDM)
MultiDM_Command = 100*(MultiDM_Command/212);

% Values capped at 80% of the maximum voltage
if sum(MultiDM_Command >= 80) > 0
    fprintf('Warning -- JR_MultiDM_Command -- Voltage higher than 80/100 of the maximum voltage.\n')
end
MultiDM_Command( MultiDM_Command >= 80) = 80;

end
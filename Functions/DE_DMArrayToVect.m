%{
DM Writing Function: Convert Map to Array
- Converts a full matrix map in nm to array in V
- Output will be a %V of max V
- Return can be fed directly into UPDATE_multiDM function
- Values capped at 80%

******************************************************
- Arguments:
    map         = Input map as 2D matrix of nm deformations
- Returns: 
    hV          = Vector of V percents ready for writing via UPDATE_mutliDM
******************************************************

Compiled By:    Daniel Echeverri
Last Modified:  08/04/2017
%}

function hV = DE_DMArrayToVect(map)

h = reshape(flipud(map),[],1);  %Value in nm and as vector

%Calculate Value in V. Eq. from DM characteristics sheet
%hV = .0000140619*(sqrt(142228e6 *h + 94738e9) - 9733357); %Value in V
hV = -30.1571 + 0.0000570399*sqrt(7.01264e9 * h + 2.79527e11);
hV = 100*(hV/212);              %Convert to percentage of max V
for i=1:length(hV)
    if (hV(i) > 80)
        hV(i) = 80;
    end
end

end
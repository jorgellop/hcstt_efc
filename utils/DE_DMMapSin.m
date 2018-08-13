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
function h = DE_DMMapSin(ho, q, xo, alp)

N   = 12;                           %DM axis base value (12x12)
Y   = repmat([-N/2:(N/2-1)],12,1);  %Matrix with ascending rows
X   = Y.';                          %Matrix with ascending columns
RHO = sqrt(X.^2 + Y.^2);            %Matrix of element dist. from cent.
TTA = atan2(Y,X);                   %Matrix of element angle from cent.

% fnm = strcat(filename, '.txt'); %function name for saving as .txt file

% Equation for corresponding actuator height
h   = 0.5*ho*(1+cos(2*pi.*RHO.*cos(TTA-q)/xo+alp));

% if out == 0 %display, plot, .txt file
%     disp(h);
%     imagesc(h);
%     colormap(gray);
%     colorbar();
%     
%     hVE = reshape(h.', [], 1);            %Reshape h to vector for printing
%     % Create formatted text file of actuator values; traverse across rows
%     fid = fopen( fnm, 'wt');
%     for i = 1:length(hVE)
%         if (i ~= length(hVE))
%             fprintf(fid, '%5.3f, ', hVE(i));
%         else
%             fprintf(fid, '%5.3f', hVE(i));
%         end
%     end
%     fclose(fid);
% elseif out == 1 %Just Display
%     disp(h);
% elseif out == 2 %Just Plot
%     % Plot Surface Map as 2d array with color for deformation level
%     imagesc(h);
%     colormap(gray);
%     colorbar();
% elseif out ==3 %Display and Plot
%     disp(h);
%     imagesc(h);
%     colormap(gray);
%     colorbar();
% elseif out ==4 %Just .txt file
%     hVE = reshape(h.', [], 1);            %Reshape h to vector for printing
%     % Create formatted text file of actuator values; traverse across rows
%     fid = fopen(fnm, 'wt');
%     for i = 1:length(hVE)
%         if (i ~= length(hVE))
%             fprintf(fid, '%5.3f, ', hVE(i));
%         else
%             fprintf(fid, '%5.3f', hVE(i));
%         end
%     end
%     fclose(fid);
% end

end





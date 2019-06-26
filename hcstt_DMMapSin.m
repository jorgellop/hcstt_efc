%{
%}
function h = hcstt_DMMapSin(ho, q, xo, alp)

N   = 12;                           %DM axis base value (12x12)
Y   = repmat([-N/2:(N/2-1)],N,1);  %Matrix with ascending rows
X   = Y.';                          %Matrix with ascending columns
RHO = sqrt(X.^2 + Y.^2);            %Matrix of element dist. from cent.
TTA = atan2(Y,X);                   %Matrix of element angle from cent.

% Equation for corresponding actuator height
h   = 0.5*ho*(1+cos(2*pi.*RHO.*cos(TTA-q)/xo+alp));
% imagesc(h)
end



% h = hcstt_DMMapSin(1, 0, 20, 0);
% imagesc(h)
% axis image
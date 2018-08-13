function [img, pixx] = multi(wl, nx, prm)
% function [img, pixx] = multi(wl, nx, prm)
% Example of a Proper function to run in parallel
% img  = image frame
% pixx = distance between pixels x (m)
% wl   = wavelength (nm)
% nx   = number of pixels x
% prm  = parameter structure

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2015may22 gmg New routine

% Create coordinates
  cz   = prm.cz;                    % Zernike coefficient vector
  iz   = prm.iz - 1;                % Zernike index vector (Fricker)
  ny   = nx;
  dx   = 0.0001;                    % distance between points x (m)
  dy   = dx    ;                    % distance between points y (m)
  znr  = (nx - 1) * dx / 2.0;       % Zernike normalization radius (m)
  zny  = (ny - 1) * dy / 2.0;
  if zny > znr
    znr  = zny;
  end
  zx   = [0 : nx - 1] * dx - znr;   % Zernike coordinates x
  zy   = [0 : ny - 1] * dy - znr;   % Zernike coordinates y
  [zgx, zgy] = meshgrid(zx, zy);

% Convert the (x, y) coordinates of each array point to polar coordinates
% cc   = angle coordinate about Z axis (radians)
% cr   = radius coordinate (m)
  [cc, cr] = cart2pol(zgx, zgy);
  cr   = cr / znr;                  % Normalize radius coordinates
  ccv  = cc(:);                     % zernfun2 wants vector inputs
  crv  = cr(:);                     % zernfun2 wants vector inputs
  crv(crv > 1.0) = 0.0;             % zernfun2 doesn't like radius > 1.0

% Calculate the Zernike functions using Paul Fricker, MathWorks, routine
  zf   = zernfun2(iz, crv, ccv) * cz;
  img  = reshape(zf, ny, nx);
  img(cr > 1.0) = NaN;
  pixx = dx;

  pause on;
  pause(10);                    % pause for 10 seconds
  pause off;
end                             % function multi

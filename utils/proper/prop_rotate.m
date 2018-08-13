function arro = prop_rotate(arri, angd, meth, extr, sx, sy, cx, cy)
% function arro = prop_rotate(arri, angd, meth, extr, sx, sy, cx, cy)
% Rotate and shift an array via interpolation (bilinear by default).
% Returns a rotated and shifted array with the same dimensions as the
% input array.
% Outputs:
% arro = 2D output array
% Required inputs:
% arri = array to be rotated
% angd = angle to rotate array counter-clockwise (degrees)
% Optional inputs:
% meth = interpolation method
%      = 'nearest' = nearest neighbor interpolation
%      = 'linear'  = bilinear interpolation (default)
%      = 'spline'  = spline interpolation
%      = 'cubic'   = bicubic interpolation as long as the data is
%                    uniformly spaced, otherwise the same as 'spline'
%                  = same as IDL interpolate routine with CUBIC = -0.5
% extr = scalar value for arro outside of the domain created by arri
%        (default = 0.0)
% sx   = amount to shift arro in x direction (pixels) (default = 0.0)
% sy   = amount to shift arro in y direction (pixels) (default = 0.0)
% cx   = pixel coordinates of arri center (default fix(nx / 2) + 1)
% cy   = pixel coordinates of arri center (default fix(ny / 2) + 1)

% 2005 Feb     jek  created idl routine
% 2014 Aug 18  gmg  Matlab translation; works correctly for non-square arrays

  if nargin < 2
    error('Proper:PROP_ROTATE', ...
          'Must specify input array and rotation angle.\n');
  elseif nargin < 3
    meth = 'linear';
  end

  if nargin < 4
    extr = 0.0;
  end

  if nargin < 6
    sx   = 0.0;
    sy   = 0.0;
  end

  [ny, nx] = size(arri);
  if nargin < 8
    cx = fix(nx / 2) + 1;
    cy = fix(ny / 2) + 1;
  end

  [aox, aoy] = meshgrid([1 : nx], [1 : ny]);    % output array coordinates
% Minus sign for angle is needed for compatibility with IDL prop_rotate.pro
  cosa = cosd(-angd);
  sina = sind(-angd);
% Transform output array coordinates to input array coordinates
  aix  = (aox - cx - sx) * cosa - (aoy - cy - sy) * sina + cx;
  aiy  = (aox - cx - sx) * sina + (aoy - cy - sy) * cosa + cy;
  arro = interp2(arri, aix, aiy, meth, extr);
end                     % function prop_rotate

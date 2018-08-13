function [zc, fit] = prop_fit_zernikes(wf, msk, apr, nz, bc, cx, cy, shr)
%        [zc, fit] = prop_fit_zernikes(wf, msk, apr, nz, bc, cx, cy, shr)
% Fit circular Zernike polynomials to a 2D error map.  The user provides
% the error map, a 2D mask (zero or one) that defines valid data values
% in the map, and the radius of the aperture.
% zc  = Fitted Zernike coefficients, ordered such that zc(1) is the
%       first Zernike polynomial (Z1, piston).
%       The units are the same as those in wf.
% fit = 2D map of the Zernike polynomials fitted to wf.
%       This map can be directly subtracted from wf, for instance.
% wf  = A 2D array containing the aberration map.  The returned Zernike
%       coefficients will have the same data units as this map.
% msk = 2D mask array indicating which corresponding values in wf are
%       valid.  A value of 1 indicates a valid point, 0 is a bad point
%       or is obscured.
% apr = Aperture radius of the beam in the map in pixels.  The Zernike
%       polynomials are normalized for a circle of this radius.
% nz  = Maximum number of Zernike polynomials to fit (1 to nz using the
%       Noll ordering scheme).  This is arbitrary for unobscured
%       polynomials.  For obscured ones, the max allowed is 22.
% bc  = Obscuration ratio of the central obscuration radius to the
%       aperture radius.  Specifying this value will cause Zernike
%       polynomials normalized for an obscured aperture to be fit rather
%       than unobscured Zernikes, which is the default. (optional)
% cx  = Specifies the center of the wavefront in the wavefront array in
% cy    pixels, with the center of the lower-left pixel being (0.0, 0.0).
%       By default, the wavefront center is at the center of the array
%       (nx / 2, ny / 2). (optional)
% shr = Amount to "shrink" the map before fitting it.  This can
%       significantly increase the speed and lower the memory
%       requirements when fitting large maps.  In effect, this indicates
%       the number of pixels between samples in wf.
%       This must be an integer value (2 or greater). (optional)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 2005 Feb     jek  created idl routine
% 2016 Mar 28  gmg  Matlab translation

  itmx = 10;            % maximum number of iterations

  if nargin < 5
    bc   = 0.0;
  end

  if bc > 0.0 & nz > 22
    error('Proper:PROP_FIT_ZERNIKES', ...
          'Limited to first 22 Obscured Zernikes.\n');
  end

  [ny, nx] = size(wf);
  if nargin < 6
    cx = fix(nx / 2) + 1;               % center pixel x (pixels)
    cy = fix(ny / 2) + 1;               % center pixel y (pixels)
  end

  if nargin == 8        % apply shrink
    nx0  = nx;
    ny0  = ny;
    nx   = nx0 / shr;
    ny   = ny0 / shr;
    [xx, yy] = meshgrid(nx, ny);
    xx   = xx - fix(nx / 2) - 1;
    yy   = yy - fix(ny / 2) - 1;
    xx   = xx * shr + fix(nx0 / 2) + 1;
    yy   = yy * shr + fix(ny0 / 2) + 1;
% bandpass filter the data prior to interpolation
% since Matlab does not have anything which replicates
% the IDL smooth function, we have to do it manually!
    wf   = prop_smooth(wf, shr);
    wf   = interp2(wf, xx, yy);
% bandpass filter the data prior to interpolation
    msk  = prop_smooth(msk, shr);
    msk  = interp2(msk, xx, yy);
    apr  = apr / shr;
    cx   = cx / shr;
    cy   = cy / shr;
  end

% Create coordinate arrays (pixels)
  [px, py] = meshgrid((1 : nx) - cx, (1 : ny) - cy);
  r    = sqrt(px.^2 + py.^2) / apr;     % normalized radius
  t    = atan2(py, px);                 % azimuth angle (radians)

  ab   = zeros(ny, nx, nz);     % 2D array values for each Zernike
  ab(:, :, 1) = 1.0;
  sr02 = sqrt( 2);
  sr03 = sqrt( 3);
  sr05 = sqrt( 5);
  sr06 = sqrt( 6);
  sr07 = sqrt( 7);
  sr10 = sqrt(10);

  if bc > 0.0
    r2   = r .* r;
    r3   = r .* r2;
    r4   = r .* r3;
    r5   = r .* r4;
    r6   = r .* r5;
    if nz >  1
      ab(:,:, 2) = 2        * cos(  t) .* r  ...
                                / sqrt(bc^2 + 1);
    end
    if nz >  2
      ab(:,:, 3) = 2        * sin(  t) .* r  ...
                                / sqrt(bc^2 + 1);
    end
    if nz >  3
      ab(:,:, 4) =     sr03 * (1 + bc^2 - 2*r2) ...
                 / (bc^2 - 1);
    end
    if nz >  4
      ab(:,:, 5) =     sr06 * sin(2*t) .* r2 ...
                                / sqrt(bc^4 + bc^2 + 1);
    end
    if nz >  5
      ab(:,:, 6) =     sr06 * cos(2*t) .* r2 ...
                                / sqrt(bc^4 + bc^2 + 1);
    end
    if nz >  6
      ab(:,:, 7) = 2 * sr02 * sin(  t) .* r  ...
                 .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
                 / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
    end
    if nz >  7
      ab(:,:, 8) = 2 * sr02 * cos(  t) .* r  ...
                 .* (2 + 2*bc^4 - 3*r2 + (2 - 3*r2) * bc^2) ...
                 / (bc^2 - 1)   / sqrt(bc^6 + 5*bc^4 + 5*bc^2 + 1);
    end
    if nz >  8
      ab(:,:, 9) = 2 * sr02 * sin(3*t) .* r3 ...
                                / sqrt(bc^6 + bc^4 + bc^2 + 1);
    end
    if nz >  9
      ab(:,:,10) = 2 * sr02 * cos(3*t) .* r3 ...
                                / sqrt(bc^6 + bc^4 + bc^2 + 1);
    end
    if nz > 10
      ab(:,:,11) =     sr05 * (1 + bc^4 - 6*r2 + 6*r4 + (4 - 6*r2)*bc^2) ...
                 / (bc^2 - 1)^2;
    end
    if nz > 11
      ab(:,:,12) =     sr10 * cos(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
                 + (3 - 4*r2) * (bc^4 + bc^2)) ...
                 / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
                 / (bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
    end
    if nz > 12
      ab(:,:,13) =     sr10 * sin(2*t) .* r2 .* (3 + 3*bc^6 - 4*r2 ...
                 + (3 - 4*r2) * (bc^4 + bc^2)) ...
                 / (bc^2 - 1)   / sqrt(bc^4 + bc^2 + 1) ...
                 / (bc^8 + 4*bc^6 + 10*bc^4 + 4*bc^2 + 1);
    end
    if nz > 13
      ab(:,:,14) =     sr10 * cos(4*t) .* r4 ...
                                / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
    end
    if nz > 14
      ab(:,:,15) =     sr10 * sin(4*t) .* r4 ...
                                / sqrt(bc^8 + bc^6 + bc^4 + bc^2 + 1);
    end
    if nz > 15
      ab(:,:,16) = 2 * sr03 * cos(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
                 + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
                 + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
                 / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
                 / (bc^6 + 9*bc^4 + 9*bc^2 + 1);
    end
    if nz > 16
      ab(:,:,17) = 2 * sr03 * sin(  t) .* r  .* (3 + 3*bc^8 - 12*r2 ...
                 + 10*r4 - 12*(r2 - 1)*bc^6 + 2*(15 - 24*r2 + 5*r4)*bc^4 ...
                 + 4*(3 - 12*r2 + 10*r4)*bc^2) ...
                 / (bc^2 - 1)^2 / sqrt(bc^4 + 4*bc^2 + 1) ...
                 / (bc^6 + 9*bc^4 + 9*bc^2 + 1);
    end
    if nz > 17
      ab(:,:,18) = 2 * sr03 * cos(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
                 + 12*(4 - 5*r2) * (bc^6 + bc^4  + bc^2))...
                 / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
                 / (bc^12 + 4*bc^10 + 10*bc^8 + 20*bc^6 + 10*bc^4 + 4*bc^2 + 1);
    end
    if nz > 18
      ab(:,:,19) = 2 * sr03 * sin(3*t) .* r3 .* (4 + 4*bc^8 - 5*r2 ...
                 + 12*(4 - 5*r2) * (bc^6 + bc^4  + bc^2))...
                 / (bc^2 - 1)   / sqrt(bc^6 + bc^4 + bc^2 + 1) ...
                 / (bc^12 + 4*bc^10 + 10*bc^8 + 20*bc^6 + 10*bc^4 + 4*bc^2 + 1);
    end
    if nz > 19
      ab(:,:,20) = 2 * sr03 * cos(5*t) .* r5 ...
                                / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
    end
    if nz > 20
      ab(:,:,21) = 2 * sr03 * sin(5*t) .* r5 ...
                                / sqrt(bc^10 + bc^8 + bc^6 + bc^4 + bc^2 + 1);
    end
    if nz > 21
      ab(:,:,22) =     sr07 * (1 + bc^6 - 12*r2 + 30*r4 - 20*r6 ...
                 + (9 - 36*r2 + 30*r4) * bc^2 + (9 - 12*r2) * bc^4) ...
                 / (bc^2 - 1)^3;
    end
  else
    zca = prop_noll_zernikes(nz);       % cell array of Zernike equations
    for iz = 1 : nz
      ab(:,:,iz) = eval(char(zca(iz)));
    end
  end

  fiti = wf;            % 2D wavefront array fit at current iteration
  err  = sqrt(sum(sum(fiti.^2 .* msk)));
  nmsk = sum(sum(msk));
  zc   = zeros(nz, 1);

  for it = 1 : itmx
    fiti = fiti .* msk;
    err0 = err;
    for iz = 1 : nz
      dzc  = sum(sum(ab(:, :, iz) .* fiti)) / nmsk;
      zc(iz) = zc(iz) + dzc;
    end
    fiti = zeros(ny, nx);
    for iz = 1 : nz
      fiti = fiti + zc(iz) * ab(:, :, iz);
    end
    fit  = fiti;
    fiti = fiti .* msk;
    fiti = wf - fiti;
    err  = sqrt(sum(sum(fiti.^2 .* msk)));
    fprintf(1, '*** Iteration  %2d  err: %24.15e\n', it, err);
    if err == 0.0 | (err0 - err) / err0 < 0.01 | err0 < err
      break;
    end
  end
end                     % function prop_fit_zernikes

function [bm,map]=prop_errormap(bm,flnm,type,sx,sy,unts,mul,magn,rot,cx,cy,smpl)
% [bm,map]=prop_errormap(bm,flnm,type,sx,sy,unts,mul,magn,rot,cx,cy,smpl)
% Read in an amplitude, surface, or wavefront error map from a
% Flexible Image Transport System (FITS) file.
% One of the types (amplitude, mirror surface, or wavefront) should be
% specified in order to properly apply the map to the wavefront.
% The amplitude map must range from 0 to 1.
% For mirror surface or wavefront error maps, the map values are assumed
% to be in meters, unless the unts switch is used to specify the units.
% The map will be interpolated to match the current wavefront sampling
% if necessary.
% Outputs:
% bm   = wavefront structure with map applied
% map  = 2D output array
% Required inputs:
% bm   = the current wavefront structure
% flnm = file name of FITS file, including extension, containing map
% Optional inputs:
% type = type of map
%        1 = file contains an amplitude error map
%        2 = file contains a mirror surface height error map (m)
%            A positive value indicates a surface point higher than the
%            mean surface.  The map will be multiplied by -2 to convert
%            it to a wavefront map to account for reflection and
%            wavefront delay (a low region on the surface causes a
%            positive increase in the phase relative to the mean).
%        3 = file contains a wavefront error map (m) (default)
% sx   = amount to shift map in x direction (m) (default = 0.0)
% sy   = amount to shift map in y direction (m) (default = 0.0)
% unts = 0 = amplitude
%        1 = meters (default)
%        2 = millimeters
%        3 = micrometers
%        4 = nanometers
% mul  = multiplies the map amplitude by the specified factor (default = 1.0)
% magn = spatially magnify the map by this factor from its default size
% rot  = Counter-clockwise rotation of map, after any resampling and
%        shifting (degrees) (default = 0.0)
% cx   = pixel coordinates of map center
% cy   = pixel coordinates of map center
% smpl = sampling of map (m)

% 2005 Feb     jek  created idl routine
% 2014 Aug 18  gmg  Matlab translation
% 2014 Sep 02  gmg  Changed prop_shift_center call to allow odd size arrays

  if exist(flnm) ~= 2           % file does not exist or is not FITS
    error('Proper:PROP_ERRORMAP', ...
          'File %s does not exist or is not a FITS file.\n', flnm);
  end

  if nargin > 5 & type == 1 & unts ~= 0
    error('Proper:PROP_ERRORMAP', ...
          'Cannot specify units for an amplitude map.\n');
  end

  if nargin < 5
    sx   = 0.0;
    sy   = 0.0;
  end

  if nargin < 11
    info = fitsinfo(flnm);      % FITS file header info
    [lna1, ina1] = ismember('NAXIS1', info.PrimaryData.Keywords(:, 1));
    [lna2, ina2] = ismember('NAXIS2', info.PrimaryData.Keywords(:, 1));
    cx   = fix(info.PrimaryData.Keywords{ina1, 2} / 2) + 1;
    cy   = fix(info.PrimaryData.Keywords{ina2, 2} / 2) + 1;
  end

  if nargin > 11
    map  = prop_readmap(bm, flnm, sx, sy, cx, cy, smpl);
  else
    map  = prop_readmap(bm, flnm, sx, sy, cx, cy);
  end

  if nargin > 7                 % apply magnification and/or rotation
    map  = prop_shift_center(map);
    if nargin > 8               % apply rotation
      map  = prop_rotate(map, rot, 'cubic', 0.0);
    end
    if magn ~= 1.0              % apply magnification
      map  = prop_magnify(map, magn, 0, size(map, 2));
    end
    map  = prop_shift_center(map, 1);
  end

  if nargin > 5                 % apply units
    if     unts == 2
      map  = map * 1.0e-3;
    elseif unts == 3
      map  = map * 1.0e-6;
    elseif unts == 4
      map  = map * 1.0e-9;
    end
  end

  if nargin > 6                 % apply multiplication factor
    map  = map * mul;
  end

  if nargin > 2 & type == 1     % amplitude error map
    bm.wf = bm.wf .* map;
  elseif nargin > 2 & type == 2 % mirror surface height error map
    bm.wf = bm.wf .* exp(complex(0.0, -4.0 * pi * map / bm.wl));
  else                          % wavefront error map
    bm.wf = bm.wf .* exp(complex(0.0,  2.0 * pi * map / bm.wl));
  end

  map  = prop_shift_center(map);
end                     % function prop_errormap

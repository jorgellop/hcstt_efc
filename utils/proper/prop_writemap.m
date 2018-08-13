function prop_writemap(map, flnm, type, pixr, smpl)
% function prop_writemap(map, flnm, type, pixr, smpl)
% Write a wavefront, surface, or amplitude aberration map to a
% Flexible Image Transport System (FITS) file.
% map  = 2D input array (required)
% flnm = file name of FITS file, including extension (required)
% type = type of map (optional)
%        1 = amplitude
%        2 = mirror surface (not wavefront) (m)
%        3 = wavefront (m) (default)
% pixr = Specifies the beam radius in the map in pixels. (optional)
%        If specified, the value of smpl (if provided) is ignored.
%        When this file is read by prop_errormap, the map will be
%        resampled as necessary to match the sampling of the beam.
% smpl = Map sampling (m) (optional)
%        Ignored if pixr is specified.

% 2005 Feb     jek  created idl routine
% 2014 Aug 05  gmg Matlab translation

  if exist(flnm) == 2           % file already exists
    error('Proper:PROP_WRITEMAP', ...
          'File %s already exists.\n', flnm);
  end
  import matlab.io.*;
  fptr = fits.createFile(flnm);

  [ny, nx] = size(map);         % number of rows and columns in map
  fits.createImg(fptr, 'double_img', [ny nx]);
  fits.deleteKey(fptr, 'EXTEND');
  fits.writeDate(fptr);
  fits.writeImg(fptr, map);

  if nargin < 3 | type == 3     % wavefront map
    mptp = 'wavefront';
  elseif type == 1              % amplitude map
    mptp = 'amplitude';
  elseif type == 2              % mirror map
    mptp = 'mirror';
  end
  fits.writeKey(fptr, 'MAPTYPE', mptp, ' error map type');

  fits.writeKey(fptr, 'X_UNIT', 'meters', ' X & Y units');
  if nargin < 3 | type ~= 1
    fits.writeKey(fptr, 'Z_UNIT', 'meters', ' Error units');
  end

  if nargin == 4 & pixr > 0
    fits.writeKey(fptr, 'RADPIX', pixr, ' beam radius in pixels');
  elseif nargin > 4
    fits.writeKey(fptr, 'PIXSIZE', smpl, ' spacing in meters', 11);
  end

  icpx = int16(floor(nx / 2) + 1);      % index of center pixel x
  icpy = int16(floor(ny / 2) + 1);      % index of center pixel y
  fits.writeKey(fptr, 'XC_PIX', icpx, ' Center X pixel coordinate');
  fits.writeKey(fptr, 'YC_PIX', icpy, ' Center Y pixel coordinate');

  fits.closeFile(fptr);
end                     % function prop_writemap

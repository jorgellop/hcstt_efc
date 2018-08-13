function [map] = prop_readmap(bm, flnm, stx, sty, mpcx, mpcy, smpl)
% function [map] = prop_readmap(bm, flnm, stx, sty, mpcx, mpcy, smpl)
% Read in a surface, wavefront, or amplitude error map from a
% Flexible Image Transport System (FITS) file, scaling if necessary.
% map  = output array
% bm   = beam structure
% flnm = file name of FITS file
% stx  = amount to shift map in x (m)
% sty  = amount to shift map in y (m)
% mpcx = map center pixel x (default is fix(mnx / 2) + 1)
% mpcy = map center pixel y (default is fix(mny / 2) + 1)
% smpl = sampling of map (m) (will override any sampling specified in
%        the file header; must be specified if header does not specify
%        sampling using the PIXSIZE value).
%        Note: if the header value RADPIX is specified (the radius of
%        the beam in the map in pixels), then this will override any
%        other sampling specifiers, either in the header or using smpl.
% Intended for internal use by Proper routines.
% Users should call either prop_errormap or prop_psd_errormap.

% 2005 Feb     jek  created idl routine
% 2014 Jul 17  gmg Matlab translation
% 2014 Sep 02  gmg  Changed fftshift to ifftshift to allow for odd size arrays

  if nargin < 4                 % Set default map shifts to 0.0
    stx  = 0.0;
    sty  = 0.0;
  end

  info = fitsinfo(flnm);        % FITS file header info
  map  = fitsread(flnm);        % FITS file array

  if nargin < 6                 % Set default map center indices
    [mny, mnx] = size(map);
    mpcx = fix(mnx / 2.0) + 1;
    mpcy = fix(mny / 2.0) + 1;
  end

% If the radius of the beam (in pixels) in the map is specified in
% the header (RADPIX value), this will override any other sampling
% specifications, either from the header (PIXSIZE value) or procedure
% call (smpl value).  If smpl is set, it overrides PIXSIZE.

  [lrdp, irdp] = ismember('RADPIX' , info.PrimaryData.Keywords(:, 1));
  [lpxs, ipxs] = ismember('PIXSIZE', info.PrimaryData.Keywords(:, 1));
  if lrdp > 0                   % radpix is defined
    pxsz = prop_get_beamradius(bm) / info.PrimaryData.Keywords{irdp, 2};
  elseif nargin > 6             % smpl is defined
    pxsz = smpl;
  elseif lpxs > 0               % pixsize is defined
    pxsz = info.PrimaryData.Keywords{ipxs, 2};
  else
    error('Proper:PROP_READMAP', ...
          'No pixel scale specified in %s header', flnm);
  end

% Resample map to current wavefront grid spacing
  map  = prop_resamplemap(bm, map, pxsz, mpcx, mpcy, stx, sty);

% Shift center of map to (1, 1)
  map  = ifftshift(map);
end                     % function prop_readmap

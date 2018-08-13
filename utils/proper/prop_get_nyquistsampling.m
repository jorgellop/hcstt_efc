function [dx] = prop_get_nyquistsampling(bm, wl)
% [dx] = prop_get_nyquistsampling(bm, wl)
% Return the Nyquist sampling interval for the current wavefront.
% This routine determines the Nyquist sampling interval for
% the current beam, which is the focal_ratio * wavelength / 2.
% dx   = Nyquist sampling interval (m)
% bm   = beam structure
% wl   = wavelength (m) (optional)
%        By default, the current wavefront wavelength is used.
%        This parameter can be used when you want to know the Nyquist
%        sampling for a wavelength other than for the current wavefront.

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  if nargin > 1
    dx   = bm.fr *    wl / 2.0;
  else
    dx   = bm.fr * bm.wl / 2.0;
  end
end                     % function prop_get_nyquistsampling

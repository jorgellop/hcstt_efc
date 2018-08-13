function [wl] = prop_get_wavelength(bm)
% [wl] = prop_get_wavelength(bm)
% Return the wavelength of the current beam.
% bm   = beam structure
% wl   = current beam's wavelength (m)

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  wl   = bm.wl;
end                     % function prop_get_wavelength

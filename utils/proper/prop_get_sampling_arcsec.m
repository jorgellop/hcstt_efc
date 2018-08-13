function da = prop_get_sampling_arcsec(bm)
% [da] = prop_get_sampling_arcsec(bm)
% Return the current wavefront sampling interval in arcseconds / pixel
% This routine is only valid when the current wavefront is at focus.
% bm   = beam structure
% da   = sampling interval (arcseconds)

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation

  fl   = prop_get_fratio(bm) * bm.diam;
  da   = prop_get_sampling(bm) * 360. * 3600. / (2. * pi * fl);
end                     % function prop_get_sampling_arcsec

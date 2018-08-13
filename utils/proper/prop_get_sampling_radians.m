function dr = prop_get_sampling_radians(bm)
% [dr] = prop_get_sampling_radians(bm)
% Return the current wavefront sampling interval in radians / pixel
% This routine is only valid when the current wavefront is at focus.
% bm   = beam structure
% dr   = sampling interval (radians)

% 2005 Feb     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  fl   = prop_get_fratio(bm) * bm.diam;
  dr   = prop_get_sampling(bm) / fl;
end                     % function prop_get_sampling_radians

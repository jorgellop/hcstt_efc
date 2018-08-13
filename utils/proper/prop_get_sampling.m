function dx = prop_get_sampling(bm)
% [dx] = prop_get_sampling(bm)
% Return the current wavefront sampling interval in meters / pixel
% bm   = beam structure
% dx   = sampling interval (m)

% 2005 Feb     jek  created idl routine
% 2014 May 07  gmg  Matlab translation

  dx   = bm.dx;
end                     % function prop_get_sampling

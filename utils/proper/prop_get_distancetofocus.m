function [dist] = prop_get_distancetofocus(bm)
% [dist] = prop_get_distancetofocus(bm)
% Determine the distance to focus from the current location (m)
% bm   = beam structure
% dist = distance to focus (m)

% 2005 Feb     jek  created idl routine
% 2014 Jul 29  gmg  Matlab translation

  dist = bm.w0_pz - bm.pz;
end                     % function prop_get_distancetofocus

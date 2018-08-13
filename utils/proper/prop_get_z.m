function [pz] = prop_get_z(bm)
% [pz] = prop_get_z(bm)
% bm   = beam structure
% pz   = distance from the initialization of the wavefront to the
%        current surface (m)

% 2013 Nov     jek  created idl routine
% 2016 Feb 24  gmg  Matlab translation

  pz   = bm.pz;
end                     % function prop_get_z
